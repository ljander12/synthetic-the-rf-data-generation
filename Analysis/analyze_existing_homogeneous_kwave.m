function T = analyze_existing_homogeneous_kwave(source_name)
%ANALYZE_EXISTING_HOMOGENEOUS_KWAVE Reuse saved k-Wave GT displacement.
% No k-Wave simulation, no Field II, no RF/beamforming, no new quality metric.
%
% INSTALL: synthetic-the-rf-data-generation/Analysis/
% RUN:
%   T = analyze_existing_homogeneous_kwave();
%   T = analyze_existing_homogeneous_kwave('homogeneous_ground_truth_theta0.mat');
%
% INPUTS (read only): homogeneous_ground_truth*.mat,
%   homogeneous_kwave_checkpoint.mat, and optional saved PG/AC metadata.
% OUTPUT: Results/sws_estimator_diagnostic/existing_homogeneous_kwave/<source>/
%
% IMPORTANT: disp_gt is a COARSE k-Wave-derived 2-D displacement image,
% not 1024x512 RF-estimated displacement. Work on native axes; no
% interpolation or artificial claim of higher spatial resolution.

    if nargin < 1 || isempty(source_name)
        source_name = 'homogeneous_ground_truth.mat';
    end
    source_name = char(string(source_name));
    [~,source_stem,ext] = fileparts(source_name);
    if isempty(ext),source_name = [source_name '.mat'];end
    root = fileparts(fileparts(mfilename('fullpath')));
    assert(isfolder(fullfile(root,'Analysis')) && isfolder(fullfile(root,'Results')), ...
        'Put this script in the repository Analysis folder.');
    addpath(genpath(root));
    toolbox_dir = fullfile(fileparts(root),'elastography','ElastographyToolbox');
    assert(isfolder(toolbox_dir), ...
        'USTB ElastographyToolbox missing: %s',toolbox_dir);
    addpath(genpath(toolbox_dir));

    for c = {'uff.elastography','uff.linear_scan', ...
             'postprocess.elastography.phase_gradient_method', ...
             'postprocess.elastography.autocorrelation_SWS_method'}
        assert(~isempty(which(c{1})), 'Class missing on MATLAB path: %s',c{1});
    end

    data_file = locate_mat(root,source_name);
    assert(~isempty(data_file), 'Cannot find %s below %s.',source_name,root);
    checkpoint = fullfile(fileparts(data_file),'homogeneous_kwave_checkpoint.mat');
    if ~isfile(checkpoint)
        checkpoint = locate_mat(root,'homogeneous_kwave_checkpoint.mat');
    end
    assert(~isempty(checkpoint) && isfile(checkpoint), ...
        'Cannot find homogeneous_kwave_checkpoint.mat for the spatial grid.');

    D = load(data_file,'disp_gt','shear_params','vec_T','prf', ...
        'array_length','scene_depth','prb_theta');
    K = load(checkpoint,'kgrid');
    req = {'disp_gt','shear_params','vec_T','prf', ...
           'array_length','scene_depth'};
    assert(all(isfield(D,req)) && isfield(K,'kgrid'), ...
        'Saved ground truth/checkpoint missing required fields.');
    p = D.shear_params;
    assert(isfield(p,'c_shear_bkg') && isfield(p,'source_freq'), ...
        'shear_params missing c_shear_bkg/source_freq.');
    c_true = double(p.c_shear_bkg);
    f0 = double(p.source_freq);
    PRF = double(D.prf);
    assert(all(isfinite([c_true,f0,PRF])) && c_true>0 && f0>0 && PRF>2*f0, ...
        'Invalid material speed/source frequency/PRF.');
    if isfield(p,'c_shear_incl') && ...
            abs(double(p.c_shear_incl)-c_true)>1e-7
        error(['This MAT file has c_shear_incl != c_shear_bkg; ' ...
               'NOT safe to assume homogeneous material.']);
    end

    A = double(D.disp_gt);
    [Nz,Nx,Nraw] = size(A);
    assert(Nraw>=25 && Nz>=20 && Nx>=20, ...
        'Expected a spatial-by-spatial-by-time displacement cube.');
    assert(all(isfinite(A(:))), 'disp_gt contains nonfinite values.');
    kgrid = K.kgrid;
    dx = double(kgrid.dx);
    dz = double(kgrid.dz);
    array_length = double(D.array_length);
    scene_depth = double(D.scene_depth);
    x_expected = (-array_length/2:dx:array_length/2).';
    z_expected = (5e-3:dz:scene_depth).';
    assert(numel(x_expected)==Nx && numel(z_expected)==Nz, ...
        ['Ground-truth axes inferred from the original get_GT_kwave ' ...
         'sampling do not match disp_gt dimensions. ' ...
         'Expected %dx%d, got %dx%d; inspect projection settings.'], ...
         numel(z_expected),numel(x_expected),Nz,Nx);
    x = x_expected;
    z = z_expected;
    clear K kgrid;
    tt = double(D.vec_T(:));
    assert(numel(tt)==Nraw && all(diff(tt)>0), ...
        'vec_T does not match displacement temporal dimension.');
    dt = median(diff(tt));
    if abs(PRF*dt-1)>0.02
        warning('vec_T interval %.6g s differs from 1/PRF %.6g s.', ...
            dt, 1/PRF);
    end

    % Follow previous GT-to-SWS comparison: packet-average 4 increments.
    % This is an analysis input assumption, not a fitted parameter.
    packet_size = 4;
    Nt = Nraw-packet_size+1;
    u = zeros(Nz,Nx,Nt,'single');
    for i = 1:Nt
        u(:,:,i) = single(mean(A(:,:,i:i+packet_size-1),3));
    end
    clear A D;

    % Use the original method-comparison filter-speed metadata when saved.
    infofile = fullfile(root,'Results','method_comparison','batch', ...
        'large_inclusions','method_comparison.mat');
    if ~isfile(infofile)
        infofile = fullfile(root,'Results','method_comparison','batch', ...
            'high_contrast','method_comparison.mat');
    end
    if isfile(infofile)
        M = load(infofile,'filter_speed','PRF','target_frequency');
        assert(isfield(M,'filter_speed'), ...
            'Method-comparison metadata missing filter_speed.');
        filter_speed = double(M.filter_speed);
        if isfield(M,'PRF') && abs(double(M.PRF)-PRF)>1e-6
            warning('Reference PRF differs: using homogeneous PRF %.2f Hz.',PRF);
        end
        if isfield(M,'target_frequency') && ...
                abs(double(M.target_frequency)-f0)>1e-6
            warning('Reference f0 differs: using homogeneous source f0 %.2f Hz.',f0);
        end
    else
        filter_speed = 3.2; % Legacy original compare_SWS_methods_batch fallback
        warning(['Could not find saved method-comparison settings. ' ...
                 'Using documented filter_speed=3.2 m/s; check comparability.']);
    end
    assert(isfinite(filter_speed) && filter_speed>0, 'Invalid filter speed.');

    out_dir = fullfile(root,'Results','sws_estimator_diagnostic', ...
        'existing_homogeneous_kwave',source_stem);
    if ~isfolder(out_dir),mkdir(out_dir);end

    fprintf('\n=== SAVED HOMOGENEOUS k-WAVE GT -> PG/AC ===\n');
    fprintf('Data: %s\n',data_file);
    fprintf('Checkpoint: %s\n',checkpoint);
    fprintf('GT grid: %d axial x %d lateral x %d frames\n',Nz,Nx,Nraw);
    fprintf('Native spacing: dx=%.4f mm, dz=%.4f mm\n',dx*1e3,dz*1e3);
    fprintf('Material c=%.4f m/s; f0=%.2f Hz; PRF=%.2f Hz\n', ...
        c_true,f0,PRF);
    fprintf('Lambda=%.3f mm; dx/lambda=%.3f; frames after packet=%d\n', ...
        1e3*c_true/f0,dx*f0/c_true,Nt);
    fprintf('No interpolation or new simulation.\n');
    if dx*f0/c_true > 0.3
        warning('Spatial sampling is relatively coarse: interpret PG/AC carefully.');
    end

    scan = uff.linear_scan('x_axis',x,'z_axis',z);
    [PG,AC] = same_pg_ac_pipeline(u,scan,PRF,f0,filter_speed);
    amp = sqrt(mean(double(u).^2,3));
    clear u;
    assert(isequal(size(PG),[Nz,Nx]) && isequal(size(AC),[Nz,Nx]), ...
        'SWS dimensions do not match native spatial grid.');

    GT = c_true*ones(Nz,Nx);
    margin_x = max(2,round(0.15*Nx));
    margin_z = max(2,round(0.15*Nz));
    roi = false(Nz,Nx);
    roi(1+margin_z:end-margin_z,1+margin_x:end-margin_x)=true;
    names = ["PG";"AC"];
    maps = {PG,AC};
    med=zeros(2,1); bias=zeros(2,1); rmse=zeros(2,1);
    mae=zeros(2,1); valid=zeros(2,1);
    for i=1:2
        ok = roi & isfinite(maps{i}) & maps{i}>0 & maps{i}<10;
        valid(i)=nnz(ok)/nnz(roi);
        if any(ok(:))
            delta=maps{i}(ok)-c_true;
            med(i)=median(maps{i}(ok));
            bias(i)=mean(delta);
            rmse(i)=sqrt(mean(delta.^2));
            mae(i)=mean(abs(delta));
        else
            med(i)=NaN;bias(i)=NaN;rmse(i)=NaN;mae(i)=NaN;
        end
    end
    T=table(names,repmat(c_true,2,1),med,bias,rmse,mae, ...
        100*rmse/c_true,valid,'VariableNames', ...
        {'Estimator','TrueSWS_mps','MedianSWS_mps','Bias_mps', ...
         'RMSE_mps','MAE_mps','RelativeRMSE_percent','ValidFraction'});
    disp(T);
    writetable(T,fullfile(out_dir,'existing_homogeneous_summary.csv'));

    fig=figure('Color','w','Position',[80 80 1340 820]);
    tl=tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
    plotmaps={GT,PG,AC,amp,PG-GT,AC-GT};
    labels={'Material GT SWS','PG reconstruction','AC reconstruction', ...
            'k-Wave displacement RMS','PG error','AC error'};
    for k=1:6
        ax=nexttile(tl,k);
        imagesc(ax,x*1e3,z*1e3,plotmaps{k});
        set(ax,'YDir','normal');axis(ax,'image');colorbar(ax);
        title(ax,labels{k});xlabel(ax,'Lateral [mm]');
        ylabel(ax,'Depth [mm]');
        if k<=3,clim(ax,[0,8]);end
        if k>=5,clim(ax,[-3,3]);end
        if any(k==[2 3 5 6])
            hold(ax,'on');
            xx=x([1+margin_x,end-margin_x])*1e3;
            zz=z([1+margin_z,end-margin_z])*1e3;
            plot(ax,xx([1 2 2 1 1]),zz([1 1 2 2 1]), ...
                'w--','LineWidth',1);
        end
    end
    title(tl,sprintf('Existing homogeneous k-Wave | c_{GT}=%.3f m/s | %s', ...
        c_true,source_stem),'Interpreter','none');
    exportgraphics(fig,fullfile(out_dir,'existing_homogeneous_maps.png'), ...
        'Resolution',160);
    close(fig);

    meta=struct('source',data_file,'checkpoint',checkpoint, ...
        'packet_size',packet_size,'PRF_Hz',PRF,'f0_Hz',f0, ...
        'filter_speed_mps',filter_speed, ...
        'warning_note', ['Native GT spatial grid; ' ...
            'no spatial interpolation, no RF estimation. ' ...
            'A poor result alone does not prove a code bug.']);
    save(fullfile(out_dir,'existing_homogeneous_results.mat'), ...
        'PG','AC','GT','amp','roi','x','z','T','meta','-v7.3');
    fprintf('Saved results in: %s\n',out_dir);
end

function path = locate_mat(root,filename)
    path='';
    prioritized={fullfile(root,'Data',filename),fullfile(root,filename), ...
        fullfile(root,'Results',filename)};
    for k=1:numel(prioritized)
        if isfile(prioritized{k}),path=prioritized{k};return;end
    end
    found=dir(fullfile(root,'**',filename));
    if ~isempty(found)
        path=fullfile(found(1).folder,found(1).name);
        if numel(found)>1
            warning('Found %d matches for %s; using %s', ...
                numel(found),filename,path);
        end
    end
end

function [PG,AC] = same_pg_ac_pipeline(u,scan,PRF,f0,filter_speed)
    elasto=uff.elastography();
    elasto.PRF=PRF;
    elasto.target_frequency=f0;
    elasto.scan=scan;
    elasto.autocorrelation_field=uff.data_field();
    elasto.autocorrelation_field.set_component('u', ...
        reshape(u,size(u,1),size(u,2),1,size(u,3)));
    angles=12;
    elasto.directional_filtering_angles= ...
        2*pi*linspace(0,1-1/angles,angles);
    extractor=postprocess.elastography.extract_freq_component(elasto);
    extractor.target_frequency=f0;
    extractor.component='u';
    extractor.run();
    omega=2*pi*f0;
    elasto.radial_filter=uff.spatial_filter('gaussian_charite','radial');
    elasto.radial_filter.sigma=filter_speed/omega;
    elasto.angular_filter=uff.spatial_filter('gaussian_charite','angular');
    elasto.angular_filter.deltaTheta= ...
        median(diff(elasto.directional_filtering_angles));
    filter_operator=postprocess.elastography.apply_filter(elasto);
    filter_operator.run();
    pg=postprocess.elastography.phase_gradient_method(elasto);
    pg.run();
    PG=squeeze(double(elasto.sws_map.u));
    ac=postprocess.elastography.autocorrelation_SWS_method(elasto);
    ac.run();
    AC=squeeze(double(elasto.sws_map.u));
end
