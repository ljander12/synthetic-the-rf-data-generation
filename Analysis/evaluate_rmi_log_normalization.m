function [T, summary, calibration] = evaluate_rmi_log_normalization()
%EVALUATE_RMI_LOG_NORMALIZATION  Controlled normalization improvement of RMI-1D.
%
% PURPOSE
%   Compare the ORIGINAL Choi-2015 Eq.(6) 1D-adapted RMI against a
%   candidate with smooth logarithmic residual normalization. Neither
%   k-Wave, RF, displacement, SWS estimates, nor ORIGINAL RMI files change.
%
% ORIGINAL:
%   Qres_old = clip((res_max-res)/(res_max-res_min), 0, 1)
%   RMI_old = 0.8*Qres_old + 0.2*Qamp
%
% CANDIDATE (same residual, same amplitude, same 0.8/0.2 weights):
%   l = log10(res)
%   Qres_log = 1/(1+exp((l-mu)/scale))
%   RMI_log = 0.8*Qres_log + 0.2*Qamp
%
%   mu and scale use UNLABELED residuals from:
%       large_inclusions, low_contrast, chaoran_baseline
%   high_contrast is NOT used for the candidate's scaling.
%   We use the pooled median as mu; P10 and P90 define scale such that
%   their log-span roughly spans logistic scores 0.9 to 0.1.
%   GT and per-case error labels NEVER enter this normalization.
%
% KEY SCIENTIFIC POINT
%   Log normalization is strictly monotonic in the residual. Therefore
%   Qres_log has exactly the same ranking/AUC as -raw residual, except
%   machine-precision ties. It CANNOT fix reversed ranking observed in
%   large_inclusions (or in certain inclusion-only regions).
%   Adding Qamp to form RMI_log CAN change the ranking/AUC.
%
% AUC: probability that correctly estimated SWS is ranked as more reliable
% than incorrectly estimated SWS; 0.5 chance; 1 best; <0.5 reversed.
% Correct/incorrect is |estimated SWS - material GT| <=/> 0.5 m/s.
% GT is accessed ONLY in the separate EVALUATION stage.
% No quality metric is optimized or thresholded against GT here.
%
% This is an EXPLORATORY comparison. All listed phantom cases have been
% investigated previously; high_contrast is normalization-held-out, NOT
% fully independent/untouched validation.
%
% INSTALL: <repo>/Analysis/evaluate_rmi_log_normalization.m
% RUN:
%   [T, summary, calibration] = evaluate_rmi_log_normalization();
%
% OUTPUTS: Results/sws_measurement_quality/rmi_log_normalization_candidate/
%   rmi_log_all_regions.csv
%   rmi_log_global_summary.csv
%   rmi_log_inclusion_summary.csv
%   rmi_log_global_auc.png
%   rmi_log_saturation.png
%   rmi_log_maps_high_contrast.png
%   rmi_log_maps_large_inclusions.png
%   <phantom>/rmi_log_candidate_maps.mat
%   rmi_log_summary.mat
%
% MATLAB R2020b+; no Statistics Toolbox required.

    cases = {'large_inclusions','low_contrast', ...
             'chaoran_baseline','high_contrast'};
    calibCases = {'large_inclusions','low_contrast','chaoran_baseline'};
    methods = {'PG','AC'};
    scoreNames = {'RawResidual','OriginalQres','LogQres', ...
                  'Qamp','OriginalRMI','LogRMI'};
    regions = {'Global','InclusionCore','BackgroundInterior'};
    epsilon = 0.5;  % FIXED SWS-error definition [m/s]

    here = fileparts(mfilename('fullpath'));
    if isfolder(fullfile(here,'Results'))
        repo = here;
    else
        repo = fileparts(here);
    end
    root = fullfile(repo,'Results','sws_measurement_quality');
    assert(isfolder(root),'Missing directory: %s',root);
    outdir = fullfile(root,'rmi_log_normalization_candidate');
    if ~isfolder(outdir), mkdir(outdir); end

    A = load(fullfile(root,'rmi_1d_calibration.mat'),'limits');
    assert(isfield(A,'limits'),'Original RMI calibration missing limits.');
    lim = A.limits;
    assert(isfield(lim,'algorithm') && ...
        strcmp(lim.algorithm,'Choi2015_Eq6_RMI_1D_adaptation'), ...
        'The saved RMI is not the expected Eq.(6) 1D adaptation.');
    assert(all(isfield(lim,{'res_min','res_max','u_min','u_max'})), ...
        'Incomplete original RMI normalization limits.');
    assert(lim.res_max>lim.res_min && lim.u_max>lim.u_min, ...
        'Invalid original calibration limits.');

    fprintf('\n=== RMI-1D: LOG NORMALIZATION CONTROL ===\n');
    fprintf('Original RMI formula and PG/AC SWS are NOT modified.\n');
    fprintf('Candidate residual normalization calibrated without GT.\n');
    fprintf('Calibration phantoms: %s, %s, %s\n',calibCases{:});
    fprintf('No high_contrast data used for new normalization parameters.\n');

    % STAGE 1: UNSUPERVISED CALIBRATION. DO NOT LOAD GT OR ROI MASKS.
    % Equal sampling cap per phantom/estimator avoids overweighting one.
    trainLog = [];
    for pi=1:numel(calibCases)
        phantom=calibCases{pi};
        f=component_file(root,phantom);
        D=load(f,'resPG','resAC','cPG','cAC','urep');
        assert(all(isfield(D,{'resPG','resAC','cPG','cAC','urep'})), ...
            'Missing raw residual or speed for %s',phantom);
        rlist={D.resPG,D.resAC};
        clist={D.cPG,D.cAC};
        for mi=1:2
            r=double(rlist{mi});
            c=double(clist{mi});
            u=double(D.urep);
            valid=isfinite(r) & (r>0) & isfinite(c) & ...
                (c>0) & (c<10) & isfinite(u) & (u>=0);
            assert(nnz(valid)>1000,'Insufficient residual samples: %s',phantom);
            v=log10(r(valid));
            maxSamples=100000;
            stride=max(1,ceil(numel(v)/maxSamples));
            v=v(1:stride:end);
            trainLog=[trainLog;v(:)]; %#ok<AGROW>
            fprintf('  calibration %s / %s: %d samples\n', ...
                phantom,methods{mi},numel(v));
        end
    end
    assert(numel(trainLog)>10000,'Insufficient calibration data.');
    mu=percentile_simple(trainLog,50);
    p10=percentile_simple(trainLog,10);
    p90=percentile_simple(trainLog,90);
    scale=(p90-p10)/(2*log(9));
    assert(isfinite(mu) && isfinite(scale) && scale>1e-5, ...
        'Residual log-distribution has no usable spread.');

    calibration = struct();
    calibration.algorithm='RMI_1D_smooth_log_residual_candidate';
    calibration.mu_log10_residual=mu;
    calibration.scale_log10_residual=scale;
    calibration.p10_log10_residual=p10;
    calibration.p90_log10_residual=p90;
    calibration.calibration_cases=calibCases;
    calibration.calibration_sample_count=numel(trainLog);
    calibration.error_threshold_for_evaluation_mps=epsilon;
    calibration.original_RMI_limits=lim;
    calibration.uses_ground_truth_for_normalization=false;
    calibration.exploratory=true;
    fprintf('LOG parameters (no GT): mu=%.6g, scale=%.6g, P10=%.6g, P90=%.6g\n', ...
        mu,scale,p10,p90);

    % STAGE 2: EVALUATION. ONLY HERE do GT, ground-truth ROI and error
    % labels enter; they affect statistics/plots, NOT quality scores.
    selftest_auc();
    T=table();
    for pi=1:numel(cases)
        phantom=cases{pi};
        R=load(component_file(root,phantom), ...
            'resPG','resAC','urep','RmiPG','RmiAC', ...
            'cPG','cAC','GT','basic','masks','regions','limits','x','z');
        needed={'resPG','resAC','urep','RmiPG','RmiAC', ...
            'cPG','cAC','GT','basic','masks','regions','limits'};
        assert(all(isfield(R,needed)),'Incomplete RMI file for %s.',phantom);
        for k={'res_min','res_max','u_min','u_max'}
            key=k{1};
            assert(abs(double(R.limits.(key))-double(lim.(key))) < ...
                1e-8*max(1,abs(double(lim.(key)))), ...
                'Original RMI calibration differs for %s: %s',phantom,key);
        end
        gt=double(R.GT);
        sz=size(gt);
        rlist={double(R.resPG),double(R.resAC)};
        clist={double(R.cPG),double(R.cAC)};
        storedList={double(R.RmiPG),double(R.RmiAC)};
        amp=double(R.urep);
        assert(isequal(sz,size(amp),size(R.basic), ...
            size(rlist{1}),size(rlist{2}), ...
            size(clist{1}),size(clist{2}), ...
            size(storedList{1}),size(storedList{2})), ...
            'Inconsistent map dimensions: %s',phantom);

        qamp=clip01((amp-double(lim.u_min))/ ...
            (double(lim.u_max)-double(lim.u_min)));
        qold=cell(1,2);
        qlog=cell(1,2);
        original=cell(1,2);
        candidate=cell(1,2);
        shared=logical(R.basic) & isfinite(gt) & isfinite(qamp);
        for mi=1:2
            r=rlist{mi};
            qold{mi}=clip01((double(lim.res_max)-r)/ ...
                (double(lim.res_max)-double(lim.res_min)));
            qlog{mi}=nan(sz);
            positive=isfinite(r) & r>0;
            qlog{mi}(positive)=stable_reliability( ...
                (log10(r(positive))-mu)/scale);
            original{mi}=0.8*qold{mi}+0.2*qamp;
            candidate{mi}=0.8*qlog{mi}+0.2*qamp;
            valid=logical(R.basic) & isfinite(gt) & ...
                isfinite(clist{mi}) & (clist{mi}>0) & (clist{mi}<10) & ...
                isfinite(r) & positive & isfinite(qamp) & ...
                isfinite(storedList{mi}) & isfinite(candidate{mi});
            assert(nnz(valid)>1000,'Too few pixels in %s/%s', ...
                phantom,methods{mi});
            discrepancy=max(abs(original{mi}(valid)-storedList{mi}(valid)));
            fprintf('%s / %s: original RMI reconstruction max diff %.3g\n', ...
                phantom,methods{mi},discrepancy);
            assert(discrepancy<1e-6, ...
                'Cannot reproduce original RMI for %s/%s',phantom,methods{mi});
            shared=shared & valid;
        end
        assert(nnz(shared)>1000,'Too few SAME pixels in %s',phantom);

        roiMasks=cell(1,numel(regions));
        roiMasks{1}=true(sz);
        savedRegionNames=string(R.regions(:));
        for ri=2:numel(regions)
            index=find(savedRegionNames==string(regions{ri}),1);
            assert(~isempty(index),'Region %s missing for %s',regions{ri},phantom);
            roiMasks{ri}=logical(R.masks{index});
            assert(isequal(size(roiMasks{ri}),sz),'Region mask size mismatch.');
        end

        fprintf('\n%s: %d same pixels.\n',phantom,nnz(shared));
        for mi=1:2
            err=abs(clist{mi}-gt);
            scores={-rlist{mi},qold{mi},qlog{mi},qamp, ...
                    storedList{mi},candidate{mi}};
            for ri=1:numel(regions)
                valid=shared & roiMasks{ri};
                if nnz(valid)<30,continue;end
                bad=err(valid)>epsilon;
                for si=1:numel(scoreNames)
                    score=double(scores{si}(valid));
                    auc=auc_ties(score,~bad);
                    if si==1
                        sat=NaN;
                    else
                        sat=mean(score>=1-1e-10);
                    end
                    row=table(string(phantom),string(methods{mi}), ...
                        string(regions{ri}),string(scoreNames{si}), ...
                        numel(score),nnz(bad),nnz(~bad),auc,sat, ...
                        'VariableNames',{'Phantom','Estimator','Region', ...
                        'Score','N','NBad','NGood','AUC','FractionAtOne'});
                    T=[T;row]; %#ok<AGROW>
                end
            end
            globalSel=shared;
            b=err(globalSel)>epsilon;
            aucRaw=auc_ties(-rlist{mi}(globalSel),~b);
            aucLog=auc_ties(qlog{mi}(globalSel),~b);
            if isfinite(aucRaw) && isfinite(aucLog)
                gap=abs(aucRaw-aucLog);
                fprintf('  %s: raw AUC %.5f, log Qres AUC %.5f (gap %.3g)\n', ...
                    methods{mi},aucRaw,aucLog,gap);
                assert(gap<1e-4, ...
                    'Log normalization unexpectedly changes raw residual ranking.');
            end
        end

        caseOutput=fullfile(outdir,phantom);
        if ~isfolder(caseOutput),mkdir(caseOutput);end
        QresLogPG=qlog{1}; %#ok<NASGU>
        QresLogAC=qlog{2}; %#ok<NASGU>
        RmiLogPG=candidate{1}; %#ok<NASGU>
        RmiLogAC=candidate{2}; %#ok<NASGU>
        save(fullfile(caseOutput,'rmi_log_candidate_maps.mat'), ...
            'QresLogPG','QresLogAC','RmiLogPG','RmiLogAC', ...
            'calibration','-v7.3');
        if ismember(phantom,{'high_contrast','large_inclusions'})
            plot_case_maps(R,shared,qold,qlog,storedList,candidate, ...
                outdir,phantom);
        end
    end

    summary=table();
    fields={'RawResidual','OriginalQres','LogQres', ...
            'Qamp','OriginalRMI','LogRMI'};
    for pi=1:numel(cases)
        for mi=1:2
            subset=T(T.Phantom==string(cases{pi}) & ...
                T.Estimator==string(methods{mi}) & T.Region=="Global",:);
            assert(height(subset)==6,'Expected six scores for %s/%s', ...
                cases{pi},methods{mi});
            vals=nan(1,6);
            sats=nan(1,6);
            for k=1:6
                ix=find(subset.Score==string(fields{k}),1);
                assert(~isempty(ix),'Missing %s',fields{k});
                vals(k)=subset.AUC(ix);
                sats(k)=subset.FractionAtOne(ix);
            end
            row=table(string(cases{pi}),string(methods{mi}), ...
                vals(1),vals(2),vals(3),vals(4),vals(5),vals(6), ...
                sats(2),sats(3),vals(6)-vals(5),vals(6)-vals(3), ...
                'VariableNames',{'Phantom','Estimator','AUCRaw', ...
                'AUCOldQres','AUCLogQres','AUCQamp','AUCOldRMI', ...
                'AUCLogRMI','OldQresSatAtOne','LogQresSatAtOne', ...
                'DeltaRMIAUC','AmplitudeMixingEffect'});
            summary=[summary;row]; %#ok<AGROW>
        end
    end

    writetable(T,fullfile(outdir,'rmi_log_all_regions.csv'));
    writetable(summary,fullfile(outdir,'rmi_log_global_summary.csv'));
    inclusion=T(T.Region=="InclusionCore",:);
    writetable(inclusion,fullfile(outdir,'rmi_log_inclusion_summary.csv'));
    save(fullfile(outdir,'rmi_log_summary.mat'), ...
        'T','summary','calibration','epsilon','cases','-v7');
    plot_global_auc(summary,outdir);
    plot_saturation(summary,outdir);

    fprintf('\n=== GLOBAL AUC: SAME PIXELS / HIGHER = BETTER ===\n');
    disp(summary(:,{'Phantom','Estimator','AUCRaw','AUCOldQres', ...
        'AUCLogQres','AUCOldRMI','AUCLogRMI','DeltaRMIAUC', ...
        'OldQresSatAtOne','LogQresSatAtOne'}));
    fprintf('\n=== INCLUSION CORE AUC (NaN = only good or only bad pixels) ===\n');
    disp(inclusion(:,{'Phantom','Estimator','Score','NBad','NGood','AUC'}));
    fprintf('\nOutputs: %s\n',outdir);
    fprintf(['INTERPRETATION: LogQres must have raw-residual AUC; ' ...
             'RMI_log may differ because Qamp is mixed in.\n']);
    fprintf(['WARNING: large_inclusions has reversed raw-residual ranking. ' ...
             'A monotonic normalization CANNOT fix this.\n']);
    fprintf(['NO new threshold chosen; NO GT used for calibrating Q; ' ...
             'NO old files overwritten.\n']);
    fprintf(['EXPLORATORY: cases previously inspected; NOT independent ' ...
             'out-of-sample proof.\n']);
end

function name=component_file(root,phantom)
    name=fullfile(root,phantom,'rmi_1d_diagnostic','rmi_1d_results.mat');
    assert(isfile(name),'Missing stored RMI components: %s',name);
end

function q=clip01(x)
    q=max(0,min(1,double(x)));
end

function q=stable_reliability(z)
% Stable evaluation of 1/(1+exp(z)); strictly decreasing in z.
    z=double(z);
    q=zeros(size(z));
    hi=z>=0;
    t=exp(-z(hi));
    q(hi)=t./(1+t);
    q(~hi)=1./(1+exp(z(~hi)));
end

function p=percentile_simple(v,pct)
    v=sort(double(v(isfinite(v(:)))));
    assert(~isempty(v),'Empty percentile sample.');
    u=1+(numel(v)-1)*pct/100;
    i=floor(u);
    j=ceil(u);
    p=v(i)+(u-i)*(v(j)-v(i));
end

function selftest_auc()
    assert(abs(auc_ties([1 2 3 4],[false false true true])-1)<1e-12, ...
        'AUC direction self-test failed.');
    assert(abs(auc_ties([1 1 2 2],[false true false true])-0.5)<1e-12, ...
        'AUC tie self-test failed.');
    assert(isnan(auc_ties([1 2],[true true])), ...
        'AUC all-good case must be undefined.');
end

function auc=auc_ties(s,good)
% Exact rank-based AUC, midrank ties; no Statistics Toolbox.
    s=double(s(:));
    good=logical(good(:));
    ok=isfinite(s);
    s=s(ok); good=good(ok);
    ng=nnz(good);
    nb=nnz(~good);
    if ng==0 || nb==0
        auc=NaN;
        return;
    end
    [ss,ix]=sort(s);
    begins=find([true;diff(ss)~=0]);
    ends=[begins(2:end)-1;numel(ss)];
    gid=cumsum([true;diff(ss)~=0]);
    ranks=(begins(gid)+ends(gid))/2;
    r=zeros(numel(s),1);
    r(ix)=ranks;
    auc=(sum(r(good))-ng*(ng+1)/2)/(ng*nb);
    auc=min(1,max(0,auc));
end

function plot_global_auc(S,outdir)
    f=figure('Visible','off','Color','w','Position',[65 65 1550 600]);
    y=[S.AUCRaw,S.AUCOldQres,S.AUCLogQres, ...
       S.AUCOldRMI,S.AUCLogRMI];
    bar(y,'grouped');
    ylim([0,1]);grid on;hold on;
    yline(0.5,'k--','Random AUC 0.5');hold off;
    xticks(1:height(S));
    names=cell(height(S),1);
    for i=1:height(S)
        names{i}=sprintf('%s / %s', ...
            strrep(char(S.Phantom(i)),'_',' '),char(S.Estimator(i)));
    end
    xticklabels(names);xtickangle(25);
    ylabel('AUC: correct SWS ranked above wrong SWS');
    title('RMI residual normalization: same pixels, error > 0.5 m/s');
    legend({'Raw residual','Old Qres (clipped)','Log Qres', ...
        'Original RMI','Log-normalized RMI'}, ...
        'Location','northoutside','Orientation','horizontal');
    exportgraphics(f,fullfile(outdir,'rmi_log_global_auc.png'), ...
        'Resolution',155);
    close(f);
end

function plot_saturation(S,outdir)
    f=figure('Visible','off','Color','w','Position',[65 65 1300 530]);
    bar(100*[S.OldQresSatAtOne,S.LogQresSatAtOne], 'grouped');
    ylim([0 100]);grid on;
    xticks(1:height(S));
    names=cell(height(S),1);
    for i=1:height(S)
        names{i}=sprintf('%s / %s', ...
            strrep(char(S.Phantom(i)),'_',' '),char(S.Estimator(i)));
    end
    xticklabels(names);xtickangle(25);
    ylabel('Qres close to 1 [%]');
    legend({'Original clipped Qres','Smooth log Qres'}, ...
        'Location','northoutside','Orientation','horizontal');
    title('Residual-score saturation: same-pixel comparison');
    exportgraphics(f,fullfile(outdir,'rmi_log_saturation.png'), ...
        'Resolution',155);
    close(f);
end

function plot_case_maps(R,shared,qold,qlog,original,candidate,outdir,phantom)
    gt=double(R.GT);
    c=double(R.cPG);
    inc=gt>min(gt(isfinite(gt)))+0.2;
    x=double(R.x(:))*1000;
    z=double(R.z(:))*1000;
    err=abs(c-gt);
    f=figure('Visible','off','Color','w','Position',[65 50 1550 800]);
    tl=tiledlayout(f,2,3,'TileSpacing','compact','Padding','compact');
    vals={err,qold{1},qlog{1},original{1},candidate{1}, ...
          candidate{1}-original{1}};
    labels={'PG | absolute SWS error [m/s]', ...
        'Original clipped Qres','New smooth log Qres', ...
        'Original RMI-1D','Candidate log RMI-1D', ...
        'Candidate RMI - original RMI'};
    limits={[0 4],[0 1],[0 1],[0 1],[0 1],[-1 1]};
    for i=1:6
        ax=nexttile(tl);
        v=double(vals{i});
        v(~shared)=NaN;
        imagesc(ax,x,z,v);
        axis(ax,'image');set(ax,'YDir','normal');
        caxis(ax,limits{i});
        colorbar(ax);
        hold(ax,'on');
        contour(ax,x,z,double(inc),[0.5 0.5], ...
            'w','LineWidth',1.2);
        hold(ax,'off');
        title(ax,labels{i},'Interpreter','none');
        xlabel(ax,'x [mm]');ylabel(ax,'z [mm]');
    end
    title(tl,sprintf('%s | normalized residual diagnostic (PG) | white = GT inclusion', ...
        strrep(phantom,'_',' ')),'Interpreter','none');
    exportgraphics(f,fullfile(outdir, ...
        ['rmi_log_maps_' phantom '.png']),'Resolution',145);
    close(f);
end
