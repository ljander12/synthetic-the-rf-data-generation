function [T, info] = diagnose_rmi_wavefield_consistency(phantom_cases)
%DIAGNOSE_RMI_WAVEFIELD_CONSISTENCY Compare Eq.6 residual across wavefields.
%
% USAGE (from repository/Analysis/):
%   [T,info] = diagnose_rmi_wavefield_consistency();
%   [T,info] = diagnose_rmi_wavefield_consistency('large_inclusions');
%   [T,info] = diagnose_rmi_wavefield_consistency( ...
%       {'large_inclusions','low_contrast','high_contrast','chaoran_baseline'});
%
% QUESTION: Does using USTB's actual filtered 200-Hz wavefield change the
% relationship between the Choi Eq.(6) 1-D residual and SWS error?
%
% Three inputs for SAME fixed PG/AC speeds (no new SWS estimates):
%  1. Original broadband time-domain residual (already saved).
%  2. Harmonic residual from USTB's UNFILTERED frequency component.
%  3. Harmonic residual from ALL USTB directional filter channels, with
%     sum of squared model errors / sum of squared x-curvatures across
%     channels. This uses all channels, not one hand-selected direction.
%
% Harmonic-equivalent PDE: u_tt = -omega^2 U, u_xx = U_xx. For a complex
% phasor U(x,z), the Eq.(6) time-sum residual is proportional to
%   r = sum_angles |(-omega^2 U) - c_est^2*U_xx|^2 / sum_angles |U_xx|^2
% for steady-state harmonic motion. The factor 1/2 cancels. 1D x-only
% dynamics remain a MODEL ASSUMPTION, not proven to hold for every angle.
%
% IMPORTANT: The old residual uses the full displacement time series; the
% second and third share a harmonic input, allowing a controlled before/
% after directional-filter comparison. No GT is used to make residuals;
% GT is loaded only at the evaluation stage (AUC, |SWS-GT| > 0.5 m/s).
% The 0.8/0.2 RMI weights and RMI normalization are NOT involved here.
% This is a diagnostic, not a new proposed metric or threshold.
%
% PREREQUISITES (READ ONLY):
%   Data/<case>_displacement_with_quality.mat (or existing alternative)
%   Results/method_comparison/batch/<case>/method_comparison.mat
%   Results/sws_measurement_quality/<case>/phase_consistency_diagnostic/
%       phase_quality_results.mat
%   Results/sws_measurement_quality/<case>/rmi_1d_diagnostic/
%       rmi_1d_results.mat
%
% OUTPUTS (NEW DIRECTORY ONLY):
%   Results/sws_measurement_quality/rmi_wavefield_consistency/
%       wavefield_auc_all_regions.csv
%       wavefield_auc_global.csv
%       wavefield_auc_global.png
%       <case>/wavefield_residual_maps.png
%       <case>/wavefield_residual_results.mat
%
% MATLAB R2020b+ with USTB on disk. No k-Wave/Field II simulation.
% V2 (2026-11-12): Reconstructed AC has isolated singular-value outliers
% on some datasets. Reconstructed PG is strictly checked; saved AC maps
% are checked strictly against each other. A separate reconstructed-AC
% comparison is diagnostic only because reconstructed AC is NOT used in
% residuals or AUC. This avoids rejecting valid analysis due to a max.
% The MATLAB/USTB run has NOT been executed in the authoring environment.

    if nargin < 1 || isempty(phantom_cases)
        phantom_cases = {'large_inclusions','low_contrast', ...
                         'high_contrast','chaoran_baseline'};
    end
    phantom_cases = cellstr(string(phantom_cases));
    phantom_cases = phantom_cases(:);
    if isempty(phantom_cases), error('Supply at least one phantom.'); end

    here = fileparts(mfilename('fullpath'));
    if isfolder(fullfile(here,'Results')), repo = here;
    else, repo = fileparts(here); end
    assert(isfolder(fullfile(repo,'Results')), ...
        'Place this file inside <repo>/Analysis/.');

    % Add USTB BEFORE loading any USTB MATLAB objects from data files.
    if isempty(which('uff.elastography'))
        tb = fullfile(fileparts(repo),'elastography','ElastographyToolbox');
        if isfolder(tb), addpath(genpath(tb)); rehash; end
    end
    assert(~isempty(which('uff.elastography')), ...
        ['USTB uff.elastography not found. Add the repository''s ' ...
         'ElastographyToolbox/USTB to the MATLAB path.']);
    assert(~isempty(which('postprocess.elastography.extract_freq_component')) && ...
           ~isempty(which('postprocess.elastography.apply_filter')) && ...
           ~isempty(which('postprocess.elastography.phase_gradient_method')) && ...
           ~isempty(which('postprocess.elastography.autocorrelation_SWS_method')), ...
           'The required USTB postprocessing classes are not on the path.');

    root = fullfile(repo,'Results','sws_measurement_quality');
    outdir = fullfile(root,'rmi_wavefield_consistency');
    if ~isfolder(outdir), mkdir(outdir); end
    err_threshold = 0.5; % Predefined SWS accuracy threshold, m/s
    modes = {'BroadbandTime','HarmonicUnfiltered','HarmonicDirectional12'};
    regions = {'Global','InclusionCore','BackgroundInterior'};
    methods = {'PG','AC'};
    T = table();
    info = struct();
    info.phantom_cases = phantom_cases;
    info.error_threshold_mps = err_threshold;
    info.normalization_used = false;
    info.residual_uses_GT = false;
    info.filtered_mode = 'sum_sq_errors_over_sum_sq_x_curvatures_all_angles';
    info.recomputed_PG_strictly_verified = true;
    info.saved_AC_maps_strictly_verified = true;
    info.recomputed_AC_diagnostic_only = true;

    selftest_wave_operator();
    fprintf('\n=== RMI Eq.(6): UNFILTERED vs USTB FILTERED WAVEFIELD ===\n');
    fprintf('No GT used in residuals. Same saved PG/AC maps throughout.\n');
    fprintf('AUC: higher is better, 0.5 = chance, <0.5 = reversed ranking.\n');

    for ci=1:numel(phantom_cases)
        cas = phantom_cases{ci};
        fprintf('\n--- Case %d/%d: %s ---\n',ci,numel(phantom_cases),cas);
        compare_path = fullfile(repo,'Results','method_comparison','batch', ...
            cas,'method_comparison.mat');
        phase_path = fullfile(root,cas,'phase_consistency_diagnostic', ...
            'phase_quality_results.mat');
        rmi_path = fullfile(root,cas,'rmi_1d_diagnostic', ...
            'rmi_1d_results.mat');
        assert(isfile(compare_path), 'Missing method comparison: %s',compare_path);
        assert(isfile(phase_path), 'Missing phase data: %s',phase_path);
        assert(isfile(rmi_path), 'Missing original RMI: %s',rmi_path);
        M=load(compare_path,'PG','AC','x_axis','z_axis', ...
            'PRF','target_frequency','filter_speed');
        P=load(phase_path,'PG','AC','common','x','z');
        R=load(rmi_path,'resPG','resAC','cPG','cAC','x','z','prf');
        needed={'PG','AC','x_axis','z_axis', ...
            'PRF','target_frequency','filter_speed'};
        assert(all(isfield(M,needed)), ...
            'Missing fields in saved method comparison for %s.',cas);
        assert(all(isfield(P,{'PG','AC','common','x','z'})) && ...
            all(isfield(R,{'resPG','resAC','cPG','cAC','x','z','prf'})), ...
            'Saved RMI/phase maps missing expected fields for %s.',cas);

        x=double(M.x_axis(:)); z=double(M.z_axis(:));
        [nz,nx]=size(M.PG);
        assert(nx>=5 && nz>=5 && all(diff(x)>0) && all(diff(z)>0), ...
            'Invalid scan dimensions/axes.');
        expected=[nz nx];
        assert(isequal(size(M.PG),expected,size(M.AC), ...
            size(P.PG),size(P.AC),size(R.resPG),size(R.resAC), ...
            size(R.cPG),size(R.cAC),size(P.common)), ...
            'Saved image dimensions mismatch.');
        assert(max(abs(double(P.x(:))-x))<1e-10 && ...
               max(abs(double(P.z(:))-z))<1e-10 && ...
               max(abs(double(R.x(:))-x))<1e-10 && ...
               max(abs(double(R.z(:))-z))<1e-10, ...
               'Saved data use different scan axes.');

        check_agreement(double(M.PG),double(P.PG), ...
            'method comparison PG vs phase PG',cas);
        check_agreement(double(M.AC),double(P.AC), ...
            'method comparison AC vs phase AC',cas);
        check_agreement(double(R.cPG),double(P.PG), ...
            'saved RMI PG vs phase PG',cas);
        check_agreement(double(R.cAC),double(P.AC), ...
            'saved RMI AC vs phase AC',cas);

        datafile = locate_data(repo,here,cas);
        D=load(datafile,'b_data','displacement_estimation_autocorr', ...
            'shear_params','prf');
        assert(all(isfield(D, ...
            {'b_data','displacement_estimation_autocorr','shear_params'})), ...
            'Original dataset lacks required displacement/scan metadata.');
        scan=D.b_data.scan;
        assert(max(abs(double(scan.x_axis(:))-x))<1e-10 && ...
               max(abs(double(scan.z_axis(:))-z))<1e-10, ...
               'Ultrasound scan axes differ from saved SWS maps.');
        PRF=double(M.PRF); f0=double(M.target_frequency);
        filter_speed=double(M.filter_speed);
        assert(isscalar(PRF) && PRF>2*f0 && isfinite(PRF) && ...
            isscalar(f0) && f0>0 && isfinite(f0) && ...
            isscalar(filter_speed) && filter_speed>0 && isfinite(filter_speed), ...
            'Saved filter/PRF/frequency metadata invalid.');
        assert(abs(double(D.shear_params.source_freq)-f0)<1e-6, ...
            'Original excitation differs from saved target frequency.');
        assert(isscalar(R.prf) && abs(double(R.prf)-PRF)<1e-6, ...
            'Original RMI PRF differs from method comparison PRF.');
        if isfield(D,'prf') && ~isempty(D.prf)
            assert(abs(double(D.prf)-PRF)<1e-6, ...
                'Original PRF differs from saved method comparison.');
        end
        v=D.displacement_estimation_autocorr.data;
        assert(isnumeric(v) && isreal(v) && mod(numel(v),nz*nx)==0, ...
            'Invalid displacement size or type.');
        nt=numel(v)/(nz*nx);
        assert(nt>=5,'Too few displacement time frames.');
        u=-reshape(single(v),nz,nx,nt); % EXACT sign convention in batch SWS.
        clear D v;

        fprintf('Input: %dx%dx%d | PRF %.1f | f0 %.1f | filter %.3f m/s\n', ...
            nz,nx,nt,PRF,f0,filter_speed);
        elasto=uff.elastography();
        elasto.PRF=PRF;
        elasto.target_frequency=f0;
        elasto.scan=scan;
        elasto.autocorrelation_field=uff.data_field();
        elasto.autocorrelation_field.set_component('u', ...
            reshape(u,nz,nx,1,nt));
        clear u;
        nang=12;
        angles=2*pi*linspace(0,1-1/nang,nang);
        elasto.directional_filtering_angles=angles;
        extractor=postprocess.elastography.extract_freq_component(elasto);
        extractor.target_frequency=f0;
        extractor.component='u';
        extractor.run();
        raw_phasor=read_complex_field(elasto.unfiltered_phase_image_xz, ...
            nz,nx,'unfiltered_phase_image_xz');
        assert(size(raw_phasor,3)==1, ...
            'Expected one unfiltered harmonic phasor.');

        omega=2*pi*f0;
        elasto.radial_filter=uff.spatial_filter('gaussian_charite','radial');
        elasto.radial_filter.sigma=filter_speed/omega;
        elasto.angular_filter=uff.spatial_filter('gaussian_charite','angular');
        elasto.angular_filter.deltaTheta=median(diff(angles));
        filt=postprocess.elastography.apply_filter(elasto);
        filt.run();
        filtered_phasors=read_complex_field( ...
            elasto.directional_filtered_field_xz,nz,nx, ...
            'directional_filtered_field_xz');
        assert(size(filtered_phasors,3)==nang, ...
            'Expected %d directional fields, got %d.', ...
            nang,size(filtered_phasors,3));

        % RECOMPUTE ONLY TO CHECK that the SWS pipeline is identical.
        % The residuals below STILL use saved PG/AC speed maps.
        pg=postprocess.elastography.phase_gradient_method(elasto);
        pg.run();
        rePG=squeeze(double(elasto.sws_map.u));
        ac=postprocess.elastography.autocorrelation_SWS_method(elasto);
        ac.run();
        reAC=squeeze(double(elasto.sws_map.u));
        check_agreement(rePG,double(P.PG), ...
            'reconstructed PG vs saved PG',cas);
        % Reconstructed AC is NOT an input to any residual/AUC. AC can
        % produce isolated large values near singular spatial estimates,
        % even when nearly all valid pixels agree to 1e-5 m/s. Keep this
        % as an explicit diagnostic, never silently ignore disagreement.
        % A valid reconstruction of PG provides the strict check of
        % USTB extraction, directional filtering and signal provenance.
        acMask=logical(P.common) & isfinite(double(P.AC)) & ...
            double(P.AC)>0 & double(P.AC)<10;
        check_agreement(reAC,double(P.AC), ...
            'reconstructed AC vs saved AC (valid common support)', ...
            cas,false,acMask);
        clear elasto extractor filt pg ac;

        cs={double(P.PG),double(P.AC)};
        rsTime={double(R.resPG),double(R.resAC)};
        [rRaw,denRaw]=harmonic_wave_residual(raw_phasor,x,omega,cs);
        [rFiltered,denFilt]=harmonic_wave_residual(filtered_phasors,x,omega,cs);
        clear raw_phasor filtered_phasors;

        % EVALUATION STAGE: GT and ROI labels enter only here.
        V=load(compare_path,'GT','roi');
        assert(all(isfield(V,{'GT','roi'})) && ...
            all(isfield(V.roi, ...
            {'image_mask','inclusion_core','background_core'})), ...
            'Saved GT/ROI missing for %s.',cas);
        GT=double(V.GT);
        assert(isequal(size(GT),expected, ...
            size(V.roi.image_mask),size(V.roi.inclusion_core), ...
            size(V.roi.background_core)), ...
            'GT or ROI dimensions differ from SWS maps.');
        image=logical(V.roi.image_mask);
        incl=logical(V.roi.inclusion_core);
        back=logical(V.roi.background_core);
        common=image & logical(P.common) & isfinite(GT) & ...
            isfinite(cs{1}) & isfinite(cs{2}) & ...
            cs{1}>0 & cs{1}<10 & cs{2}>0 & cs{2}<10;
        % One fixed common support across all modes AND both estimators.
        for mi=1:2
            common=common & isfinite(rsTime{mi}) & ...
                isfinite(rRaw{mi}) & isfinite(rFiltered{mi});
        end
        assert(nnz(common)>100, ...
            'Too few common finite pixels. Check filtered-field derivative support.');
        masks={common,common & incl,common & back};
        fprintf('SAME-PIXEL support: Global %d | Inclusion %d | Background %d\n', ...
            nnz(masks{1}),nnz(masks{2}),nnz(masks{3}));
        if nnz(common)<0.1*nnz(image)
            warning('Only %.1f%% of image ROI has common support.', ...
                100*nnz(common)/nnz(image));
        end

        out_case=fullfile(outdir,cas);
        if ~isfolder(out_case),mkdir(out_case);end
        mapPairs=cell(2,3);
        for mi=1:2
            mapPairs(mi,:)={rsTime{mi},rRaw{mi},rFiltered{mi}};
            er=abs(cs{mi}-GT);
            bad=er>err_threshold;
            for ri=1:numel(regions)
                mask=masks{ri};
                nBad=nnz(mask&bad); nGood=nnz(mask&~bad);
                for si=1:numel(modes)
                    res=mapPairs{mi,si};
                    rv=res(mask);
                    bm=bad(mask);
                    auc=auc_good_has_smaller_residual(rv,bm);
                    medGood=finite_median(rv(~bm));
                    medBad=finite_median(rv(bm));
                    row=table(string(cas),string(regions{ri}), ...
                        string(methods{mi}),string(modes{si}), ...
                        nnz(mask),nGood,nBad,auc,medGood,medBad, ...
                        'VariableNames',{'Phantom','Region','Estimator', ...
                        'Wavefield','N','NGood','NBad','AUC', ...
                        'MedianResidualGood','MedianResidualBad'});
                    T=[T;row]; %#ok<AGROW>
                end
            end
        end
        plot_case(out_case,x,z,GT,cs{1},common, ...
            rsTime{1},rRaw{1},rFiltered{1});
        metadata=struct('PRF',PRF,'f0',f0, ...
            'filter_speed',filter_speed,'nang',nang, ...
            'nx',nx,'nz',nz,'nt',nt, ...
            'old_time_domain_from_saved_RMI',true, ...
            'residual_model','Choi Eq6 1D x-only, harmonic complex phasor', ...
            'normalization','none', ...
            'common_valid_count',nnz(common));
        caseRows=T(T.Phantom==string(cas),:);
        save(fullfile(out_case,'wavefield_residual_results.mat'), ...
            'caseRows','metadata','common','GT','x','z', ...
            'denRaw','denFilt','mapPairs','-v7.3');
        disp(caseRows(caseRows.Region=="Global", ...
            {'Estimator','Wavefield','NGood','NBad','AUC', ...
             'MedianResidualGood','MedianResidualBad'}));
    end

    G=T(T.Region=="Global",:);
    writetable(T,fullfile(outdir,'wavefield_auc_all_regions.csv'));
    writetable(G,fullfile(outdir,'wavefield_auc_global.csv'));
    fig=figure('Color','w','Visible','off','Position',[100 100 1250 510]);
    labels=strings(numel(phantom_cases)*2,1);
    Y=nan(numel(labels),numel(modes));
    ir=0;
    for ci=1:numel(phantom_cases)
        for mi=1:2
            ir=ir+1;
            labels(ir)=string(phantom_cases{ci})+" / "+string(methods{mi});
            for si=1:numel(modes)
                ix=G.Phantom==string(phantom_cases{ci}) & ...
                   G.Estimator==string(methods{mi}) & ...
                   G.Wavefield==string(modes{si});
                if nnz(ix)==1, Y(ir,si)=G.AUC(ix); end
            end
        end
    end
    bar(Y,'grouped'); grid on;hold on;
    yline(0.5,'k--','AUC 0.5 (chance)');
    xticks(1:numel(labels));xticklabels(cellstr(labels));xtickangle(30);
    ylabel('AUC: correct SWS has smaller Eq.(6) residual');ylim([0 1]);
    legend({'Broadband time (old)','Harmonic, unfiltered', ...
            'Harmonic, 12 directions'},'Location','northoutside', ...
            'Orientation','horizontal');
    title('Same SWS, same valid pixels | error threshold 0.5 m/s');
    exportgraphics(fig,fullfile(outdir,'wavefield_auc_global.png'), ...
        'Resolution',160);close(fig);
    save(fullfile(outdir,'wavefield_summary.mat'),'T','info');
    fprintf('\nSaved: %s\n',outdir);
    fprintf(['INTERPRETATION: Harmonic unfiltered vs directional is the ' ...
       'controlled filtering contrast. Broadband vs harmonic also ' ...
       'changes temporal frequency content.\n']);
    fprintf(['WARNING: Eq.(6) assumes x-only propagation; 12 angular ' ...
       'channels may violate that assumption. No GT-based calibration ' ...
       'or threshold was fitted.\n']);
end

function f=locate_data(repo,analysis,cas)
    fn=[cas '_displacement_with_quality.mat'];
    c={fullfile(repo,'Data',cas,fn),fullfile(repo,'Data',fn), ...
        fullfile(repo,fn),fullfile(analysis,fn)};
    f='';
    for ii=1:numel(c)
        if isfile(c{ii}),f=c{ii};break;end
    end
    assert(~isempty(f),'Could not locate displacement data %s.',fn);
end

function check_agreement(a,b,what,cas,strict,validMask)
    % strict=true for saved-map provenance and reconstructed PG.
    % strict=false ONLY for recomputed AC, which is not used downstream.
    if nargin<5,strict=true;end
    if nargin<6 || isempty(validMask),validMask=true(size(a));end
    assert(isequal(size(a),size(b),size(validMask)), ...
        'Shape mismatch for %s (%s).',what,cas);
    keep=isfinite(a)&isfinite(b)&logical(validMask);
    assert(nnz(keep)>100,'Insufficient overlap: %s (%s).',what,cas);
    diffabs=abs(a(keep)-b(keep));
    p95=sample_percentile(diffabs,95);
    p99=sample_percentile(diffabs,99);
    dmax=max(diffabs);
    nLarge=nnz(diffabs>=0.2);
    fracLarge=nLarge/numel(diffabs);
    fprintf(['Consistency %s: N %d, P95 %.6g, P99 %.6g, max %.6g m/s, ' ...
        '|diff|>=0.2: %d (%.4f%%)\n'], ...
        what,numel(diffabs),p95,p99,dmax,nLarge,100*fracLarge);
    if strict
        assert(p95<0.02 && dmax<0.2, ...
            ['SWS maps disagree too strongly (%s, %s). ' ...
             'Do not interpret residual test until USTB settings/versions match.'], ...
             cas,what);
    else
        if p95>=0.02 || fracLarge>0.01
            warning(['Recomputed AC disagrees on a material fraction of ' ...
                'valid pixels for %s; check AC settings. ' ...
                'The residual diagnostic uses saved AC exclusively, ' ...
                'and strict reconstructed PG validation is still required.'],cas);
        elseif dmax>=0.2
            fprintf(['NOTE: isolated reconstructed AC outliers observed; ' ...
                'these values are NOT used for residuals or AUC.\n']);
        end
    end
end

function [res,denOut]=harmonic_wave_residual(U,x,omega,cMaps)
    % Wave equation Eq.(6), complex harmonic field, all angular channels.
    [nz,nx,nd]=size(U);
    assert(numel(x)==nx && nx>=5 && nd>=1, ...
        'Unexpected harmonic wavefield dimensions.');
    den=zeros(nz,nx);
    nums={zeros(nz,nx),zeros(nz,nx)};
    dxr=x(3:end)-x(2:end-1);
    dxl=x(2:end-1)-x(1:end-2);
    span=x(3:end)-x(1:end-2);
    c2={cMaps{1}.^2,cMaps{2}.^2};
    for d=1:nd
        field=double(U(:,:,d));
        fwd=bsxfun(@rdivide,field(:,3:end)-field(:,2:end-1),dxr');
        bwd=bsxfun(@rdivide,field(:,2:end-1)-field(:,1:end-2),dxl');
        xx=bsxfun(@rdivide,2*(fwd-bwd),span');
        utt=-(omega^2)*field(:,2:end-1);
        den(:,2:end-1)=den(:,2:end-1)+abs(xx).^2;
        for mi=1:2
            err=utt-c2{mi}(:,2:end-1).*xx;
            nums{mi}(:,2:end-1)=nums{mi}(:,2:end-1)+abs(err).^2;
        end
    end
    positives=den(isfinite(den)&den>0);
    assert(~isempty(positives),'No positive harmonic spatial curvature.');
    threshold=median(positives)*1e-12;
    support=isfinite(den)&den>threshold;
    support(:,[1 end])=false;
    res=cell(1,2);
    for mi=1:2
        v=nums{mi}./den;
        v(~support|~isfinite(v))=NaN;
        res{mi}=v;
    end
    denOut=den;
    denOut(~support)=NaN;
end

function U=read_complex_field(v,nz,nx,label)
    a=[];
    if isnumeric(v)
        a=v;
    elseif isstruct(v)
        if isfield(v,'u'),a=v.u;
        elseif isfield(v,'data'),a=v.data;
        end
    elseif isobject(v)
        if ismethod(v,'get_component')
            try, a=v.get_component('u');catch, end
        end
        if isempty(a) && isprop(v,'u')
            try, a=v.u;catch, end
        end
        if isempty(a) && isprop(v,'data')
            try, a=v.data;catch, end
        end
    end
    if isstruct(a) && isfield(a,'u'),a=a.u;end
    if iscell(a) && numel(a)==1,a=a{1};end
    assert(isnumeric(a) && ~isempty(a), ...
        'Cannot read USTB field %s (%s).',label,class(v));
    if isvector(a)
        assert(mod(numel(a),nz*nx)==0, ...
            'Vectorized %s has incorrect length.',label);
        a=reshape(a,nz,nx,[]);
    else
        dims=size(a);
        assert(dims(1)==nz && dims(2)==nx, ...
            'USTB %s shape %s != [%d %d].', ...
            label,mat2str(dims),nz,nx);
        a=reshape(a,nz,nx,[]);
    end
    U=double(a);
    % Complex amplitude is necessary. Wrapped phase alone is not adequate.
    assert(~isreal(U) && any(abs(imag(U(:)))>0), ...
        ['Field %s is real-valued. Cannot compute PDE residual ' ...
         'from a phase-only field.'],label);
end

function a=auc_good_has_smaller_residual(res,bad)
    if nnz(bad)==0 || nnz(~bad)==0 || numel(res)<2
        a=NaN;return;
    end
    q=-double(res(:)); good=~logical(bad(:));
    [s,ix]=sort(q);
    good=double(good(ix));
    group=cumsum([true;diff(s)~=0]);
    ng=accumarray(group,good);
    nb=accumarray(group,1-good);
    below=[0;cumsum(nb(1:end-1))];
    a=sum(ng.*(below+0.5*nb))/(sum(ng)*sum(nb));
end

function m=finite_median(v)
    v=v(isfinite(v));
    if isempty(v),m=NaN;else,m=median(v);end
end

function p=sample_percentile(v,pct)
    v=sort(double(v(isfinite(v))));
    assert(~isempty(v),'No finite samples in percentile computation.');
    k=1+(numel(v)-1)*pct/100;
    lo=floor(k);hi=ceil(k);
    p=v(lo)+(k-lo)*(v(hi)-v(lo));
end

function plot_case(out_case,x,z,GT,PG,mask,rtime,rraw,rfilt)
    maps={abs(PG-GT),rtime,rraw,rfilt,rfilt./rraw,double(mask)};
    titles={'Absolute PG SWS error [m/s]', ...
        'Old broadband residual (log10)', ...
        'Unfiltered 200 Hz residual (log10)', ...
        'Directional 200 Hz residual (log10)', ...
        'Filtered / unfiltered residual (log10)', ...
        'COMMON evaluation mask'};
    fig=figure('Color','w','Visible','off','Position',[100 100 1450 780]);
    tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
    for j=1:6
        nexttile;
        a=double(maps{j});
        if j>=2 && j<=5
            a(a<=0)=NaN;
            a=log10(a);
        end
        a(~mask)=NaN;
        imagesc(x*1e3,z*1e3,a);
        axis image;set(gca,'YDir','normal');colorbar;
        xlabel('x [mm]');ylabel('z [mm]');
        title(titles{j});
    end
    exportgraphics(fig,fullfile(out_case,'wavefield_residual_maps.png'), ...
        'Resolution',160);
    close(fig);
end

function selftest_wave_operator()
    % Analytic harmonic wave: correct speed must give a lower PDE residual.
    x=(-150:150)'*0.00015;
    c=2.7;f=200;omega=2*pi*f;
    U=reshape(exp(1i*(omega/c)*x.'),1,[],1);
    [rr,~]=harmonic_wave_residual(U,x,omega, ...
        {c*ones(1,numel(x)),(1.6*c)*ones(1,numel(x))});
    i=floor(numel(x)/2);
    assert(rr{1}(1,i) < 0.01*rr{2}(1,i), ...
        'Harmonic PDE synthetic self-test failed.');
    a=auc_good_has_smaller_residual([1;2;3;4],[false;false;true;true]);
    assert(abs(a-1)<1e-12,'AUC self-test failed.');
end
