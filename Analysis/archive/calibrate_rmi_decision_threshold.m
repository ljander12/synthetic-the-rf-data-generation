function [T, selected, calibration] = calibrate_rmi_decision_threshold()
%CALIBRATE_RMI_DECISION_THRESHOLD  Fixed RMI-1D reliability cutoff.
%
% Calibrates ONE absolute decision threshold for existing RMI-1D Q:
%     accept SWS if Q_RMI >= qstar; reject SWS if Q_RMI < qstar.
%
% This is THRESHOLD CALIBRATION, not a fourth quality metric, and does
% not alter the existing RMI residual/amplitude weighting or SWS maps.
% Ground truth is used ONLY for training selection and offline evaluation.
%
% Data split set BEFORE execution:
%  TRAIN: low_contrast, chaoran_baseline (PG and AC, equally weighted).
%  EVALUATE: high_contrast (same fixed cutoff; previously inspected).
%  DESCRIPTIVE: large_inclusions (used for RMI score normalization).
%
% Primary error: |c_est-c_GT| > 0.5 m/s. Also report 1.0 m/s.
% Cutoff optimization: maximize mean (bad-pixel recall - good-pixel
% rejection rate) across the four TRAIN phantom/estimator combinations,
% subject to mean acceptance fraction between 20 and 90 percent.
% No test-case data are used when choosing qstar. This cutoff maximizes
% balanced discrimination; it does NOT by itself guarantee low error risk.
%
% Input: previously saved MAT files ONLY:
%   Results/sws_measurement_quality/<phantom>/
%     rmi_1d_diagnostic/rmi_1d_results.mat
%     phase_consistency_diagnostic/phase_quality_results.mat
%
% Output (does NOT overwrite original results):
%   Results/sws_measurement_quality/rmi_threshold_calibration/
%     fixed_rmi_threshold.mat
%     calibration_curve.csv
%     threshold_evaluation.csv
%     rmi_threshold_calibration.png
%     rmi_threshold_outcomes.png
%     high_contrast_threshold_maps.png
%
% Run (put in Analysis/):
%   [T, selected, calibration] = calibrate_rmi_decision_threshold();
%
% IMPORTANT:
% - This is EXPLORATORY internal validation, not truly unseen data:
%   all these cases have already been examined during development.
% - Pixels are spatially correlated; phantom-level sample size is tiny.
% - RMI residual/amplitude normalization was derived on large_inclusions.
% - No fixed cutoff is guaranteed reliable; inspect coverage, bad recall,
%   rejection of good pixels AND false confidence among accepted pixels.
%
% MATLAB R2020b+.

    repo = find_repo(fileparts(mfilename('fullpath')));
    root = fullfile(repo, 'Results', 'sws_measurement_quality');
    out = fullfile(root,'rmi_threshold_calibration');
    assert(isfolder(root),'Missing Results/sws_measurement_quality folder.');
    if ~isfolder(out), mkdir(out); end

    TRAIN = {'low_contrast','chaoran_baseline'};
    TEST = {'high_contrast'};
    NORMALIZATION_SAMPLE = {'large_inclusions'};
    estimators = {'PG','AC'};
    regions = {'Global','InclusionCore','BackgroundInterior'};
    eps_primary = 0.5;
    eps_eval = [0.5,1.0];
    q_grid = (0:0.005:1).';  % thresholds FIXED in advance
    min_coverage = 0.20;
    max_coverage = 0.90;

    fprintf('\n=== RMI-1D: FIXED DECISION THRESHOLD ===\n');
    fprintf('Train: %s, %s | Test: %s\n', TRAIN{1},TRAIN{2},TEST{1});
    fprintf('Descriptive only: %s (RMI normalization source).\n',NORMALIZATION_SAMPLE{1});
    fprintf('Event: |estimated SWS - material GT| > %.2f m/s.\n',eps_primary);
    fprintf('ONE threshold shared across PG, AC and all phantoms.\n');
    fprintf('This is a cutoff on existing RMI, NOT a new quality metric.\n');

    % TRAIN STAGE. Only the two TRAIN phantoms can influence qstar.
    nGroups = numel(TRAIN)*numel(estimators);
    recall = nan(numel(q_grid),nGroups);
    falseRejection = nan(numel(q_grid),nGroups);
    coverage = nan(numel(q_grid),nGroups);
    g = 0;
    for pi=1:numel(TRAIN)
        D = load_one_case(root,TRAIN{pi});
        for mi=1:2
            g=g+1;
            q = D.rmi{mi}(D.common);
            e = abs(D.speed{mi}(D.common)-D.gt(D.common));
            bad = e>eps_primary;
            assert(any(bad) && any(~bad), ...
                'Train %s / %s has no bad OR no good pixels at eps=%.2f.',...
                TRAIN{pi},estimators{mi},eps_primary);
            [nAcceptBad,nAcceptGood] = accepted_counts(q,bad,q_grid);
            nBad=nnz(bad); nGood=nnz(~bad); nAll=numel(q);
            recall(:,g)=1-nAcceptBad/nBad;
            falseRejection(:,g)=1-nAcceptGood/nGood;
            coverage(:,g)=(nAcceptBad+nAcceptGood)/nAll;
            fprintf('Train %-18s %-2s  N=%d  event=%.1f %%\n', ...
                TRAIN{pi},estimators{mi},nAll,100*nBad/nAll);
        end
    end

    meanRecall = mean(recall,2);
    meanFalseReject = mean(falseRejection,2);
    meanCoverage = mean(coverage,2);
    youden = meanRecall-meanFalseReject;
    eligible = meanCoverage>=min_coverage & meanCoverage<=max_coverage & ...
        isfinite(youden);
    if ~any(eligible)
        error('No q grid values retain 20-90%% mean coverage in TRAIN.');
    end
    ranked=youden;
    ranked(~eligible)=-Inf;
    best=max(ranked);
    % Tie-break: smallest threshold => more pixels accepted.
    ix=find(abs(ranked-best)<1e-12,1,'first');
    qstar=q_grid(ix);
    if best<=0
        warning(['No useful positive balanced-discrimination score on TRAIN. ' ...
                 'Any selected cutoff is exploratory; do NOT deploy.']);
    end

    calibration=table(q_grid,meanRecall,meanFalseReject,meanCoverage,...
        youden,eligible,'VariableNames',{'Qthreshold',...
        'MeanBadDetection','MeanGoodFalseRejection','MeanCoverage',...
        'MeanYoudenJ','EligibleCoverage'});
    writetable(calibration,fullfile(out,'calibration_curve.csv'));

    selected=struct();
    selected.qstar=qstar;
    selected.error_threshold_mps=eps_primary;
    selected.calibration_phantoms=TRAIN;
    selected.test_phantoms=TEST;
    selected.rmi_normalization_sample=NORMALIZATION_SAMPLE;
    selected.estimators=estimators;
    selected.min_mean_coverage=min_coverage;
    selected.max_mean_coverage=max_coverage;
    selected.objective='Max macro Youden J (bad recall - good false rejection)';
    selected.qgrid_step=q_grid(2)-q_grid(1);
    selected.train_mean_bad_detection=meanRecall(ix);
    selected.train_mean_good_false_rejection=meanFalseReject(ix);
    selected.train_mean_coverage=meanCoverage(ix);
    selected.train_mean_youden=youden(ix);
    selected.certified=false;
    selected.validation_note=['Exploratory development split only. High contrast ' ...
        'was previously inspected; large_inclusions normalized the RMI.'];
    fprintf('\nLOCKED qstar=%.3f on TRAIN only.\n',qstar);
    fprintf('TRAIN bad detection=%.1f%%, good false rejection=%.1f%%, coverage=%.1f%%\n',...
        100*meanRecall(ix),100*meanFalseReject(ix),100*meanCoverage(ix));
    fprintf('TRAIN macro Youden J=%.3f\n',youden(ix));

    % EVALUATION STAGE: fixed qstar (no refitting, no per-case cutoff).
    T=table();
    cases = [TRAIN TEST NORMALIZATION_SAMPLE];
    labels = [repmat({'TRAIN'},1,numel(TRAIN)),...
              repmat({'TEST_PREVIOUSLY_SEEN'},1,numel(TEST)),...
              repmat({'RMI_NORMALIZATION_IN_SAMPLE'},1,numel(NORMALIZATION_SAMPLE))];
    highContrastView=struct();
    for pi=1:numel(cases)
        D=load_one_case(root,cases{pi});
        for mi=1:2
            qmap=D.rmi{mi};
            errmap=abs(D.speed{mi}-D.gt);
            for ri=1:numel(regions)
                switch regions{ri}
                    case 'Global', mask=D.common;
                    case 'InclusionCore', mask=D.inclusion;
                    case 'BackgroundInterior', mask=D.background;
                end
                q=qmap(mask);
                e=errmap(mask);
                if isempty(q),continue;end
                keep=q>=qstar;
                reject=~keep;
                for ei=1:numel(eps_eval)
                    epsilon=eps_eval(ei);
                    bad=e>epsilon;
                    good=~bad;
                    nb=nnz(bad); ng=nnz(good); na=nnz(keep);
                    nr=nnz(reject);
                    if nb>0, badDetected=nnz(reject & bad)/nb;
                    else, badDetected=NaN;end
                    if ng>0, goodRejected=nnz(reject & good)/ng;
                    else, goodRejected=NaN;end
                    if na>0, risk=nnz(keep & bad)/na;
                    else, risk=NaN;end
                    if nr>0, rejectPrecision=nnz(reject & bad)/nr;
                    else, rejectPrecision=NaN;end
                    if na>0, accMAE=mean(e(keep));
                    else, accMAE=NaN;end
                    row=table(string(cases{pi}),string(labels{pi}),...
                        string(estimators{mi}),string(regions{ri}),...
                        epsilon,qstar,numel(q),nb,ng,nnz(keep),nnz(reject),...
                        mean(bad),mean(keep),badDetected,goodRejected,...
                        risk,rejectPrecision,accMAE,mean(e),...
                        'VariableNames',{'Phantom','Split','Estimator','Region',...
                        'ErrorLimit_mps','Qthreshold','N_Pixels','N_Bad','N_Good',...
                        'N_Accepted','N_Rejected','BadPrevalence',...
                        'Coverage','BadDetectionRecall','GoodFalseRejection',...
                        'FalseConfidenceAmongAccepted','BadPrecisionAmongRejected',...
                        'MAE_Accepted_mps','MAE_All_mps'});
                    T=[T;row]; %#ok<AGROW>
                end
            end
            if strcmp(cases{pi},'high_contrast') && mi==1
                highContrastView.gt=D.gt;
                highContrastView.error=errmap;
                highContrastView.q=qmap;
                highContrastView.common=D.common;
                highContrastView.inclusion=D.inclusion;
                highContrastView.x=D.x;
                highContrastView.z=D.z;
            end
        end
        fprintf('Evaluated %-20s (%s)\n',cases{pi},labels{pi});
    end
    writetable(T,fullfile(out,'threshold_evaluation.csv'));

    calibration_file=fullfile(out,'fixed_rmi_threshold.mat');
    save(calibration_file,'selected','TRAIN','TEST','NORMALIZATION_SAMPLE',...
        'eps_primary','eps_eval','q_grid','calibration','-v7');
    plot_calibration(calibration,qstar,...
        fullfile(out,'rmi_threshold_calibration.png'));
    plot_outcomes(T,fullfile(out,'rmi_threshold_outcomes.png'));
    if isfield(highContrastView,'gt')
        plot_high_contrast(highContrastView,qstar,eps_primary,...
            fullfile(out,'high_contrast_threshold_maps.png'));
    end
    save(fullfile(out,'threshold_diagnostic_results.mat'),...
        'T','selected','calibration','-v7.3');

    fprintf('\n=== LOCKED THRESHOLD RESULTS: GLOBAL, ERROR > 0.5 m/s ===\n');
    rows=T(T.Region=="Global" & T.ErrorLimit_mps==eps_primary,:);
    disp(rows(:,{'Phantom','Split','Estimator','Coverage',...
        'BadDetectionRecall','GoodFalseRejection',...
        'FalseConfidenceAmongAccepted'}));
    fprintf('\n=== HIGH CONTRAST INCLUSION (all or most may be bad) ===\n');
    rows=T(T.Phantom=="high_contrast" & T.Region=="InclusionCore" & ...
        T.ErrorLimit_mps==eps_primary,:);
    disp(rows(:,{'Estimator','N_Bad','N_Good','Coverage',...
        'BadDetectionRecall','FalseConfidenceAmongAccepted'}));
    fprintf('\nCalibration and fixed-cutoff test saved to: %s\n',out);
    fprintf('NOTE: no independent untouched phantom; do not claim clinical calibration.\n');
    fprintf('NOTE: do not update production RMI or SWS code from this run alone.\n');
    fprintf('NOTE: low coverage or high false confidence means cutoff is inadequate.\n');
end

function repo=find_repo(here)
    if isfolder(fullfile(here,'Results'))
        repo=here;
    else
        repo=fileparts(here);
    end
    assert(isfolder(fullfile(repo,'Results')),...
        'Script must be placed in repo/Analysis/ or repo root.');
end

function D=load_one_case(root,phantom)
    p=fullfile(root,phantom);
    rmiFile=fullfile(p,'rmi_1d_diagnostic','rmi_1d_results.mat');
    phaseFile=fullfile(p,'phase_consistency_diagnostic','phase_quality_results.mat');
    assert(isfile(rmiFile),'Missing RMI source: %s',rmiFile);
    assert(isfile(phaseFile),'Missing Q/ROI source: %s',phaseFile);
    R=load(rmiFile,'RmiPG','RmiAC','cPG','cAC','GT',...
        'basic','regions','masks','x','z');
    P=load(phaseFile,'common','GT');
    reqR={'RmiPG','RmiAC','cPG','cAC','GT',...
        'basic','regions','masks','x','z'};
    assert(all(isfield(R,reqR)) && all(isfield(P,{'common','GT'})),...
        'Missing expected fields in MAT sources for %s',phantom);
    D.gt=double(R.GT);
    D.speed={double(R.cPG),double(R.cAC)};
    D.rmi={double(R.RmiPG),double(R.RmiAC)};
    D.x=double(R.x(:));
    D.z=double(R.z(:));
    sz=size(D.gt);
    assert(isequal(size(P.GT),sz,size(R.basic),sz,size(P.common),sz,...
        size(D.rmi{1}),sz,size(D.rmi{2}),sz,...
        size(D.speed{1}),sz,size(D.speed{2}),sz),...
        'Map dimensions mismatch for %s.',phantom);
    assert(numel(D.x)==sz(2) && numel(D.z)==sz(1),...
        'Grid dimensions mismatch for %s.',phantom);
    finiteGT=isfinite(D.gt) & isfinite(double(P.GT));
    assert(nnz(finiteGT)>0 && ...
        max(abs(D.gt(finiteGT)-double(P.GT(finiteGT))))<1e-6,...
        'GT differs between saved source files for %s.',phantom);
    D.common=logical(R.basic) & logical(P.common) & finiteGT;
    for mi=1:2
        D.common=D.common & isfinite(D.speed{mi}) & ...
            isfinite(D.rmi{mi}) & (D.rmi{mi}>=-1e-7) & ...
            (D.rmi{mi}<=1+1e-7);
        D.rmi{mi}=max(0,min(1,D.rmi{mi}));
    end
    assert(nnz(D.common)>1000,'Too few common pixels for %s.',phantom);
    assert(numel(R.regions)==numel(R.masks),...
        'Region names/masks mismatch for %s.',phantom);
    names=string(R.regions(:));
    i=find(names=="InclusionCore",1);
    b=find(names=="BackgroundInterior",1);
    assert(~isempty(i) && ~isempty(b),...
        'InclusionCore / BackgroundInterior missing in %s.',phantom);
    D.inclusion=D.common & logical(R.masks{i});
    D.background=D.common & logical(R.masks{b});
    assert(isequal(size(D.inclusion),sz,size(D.background),sz),...
        'ROI size mismatch for %s.',phantom);
    assert(~any(D.inclusion(:)&D.background(:)),...
        'ROI overlap for %s.',phantom);
end

function [acceptBad,acceptGood]=accepted_counts(q,bad,qgrid)
% Fast fixed-threshold counts without ground-truth-dependent selection.
    q=double(q(:));bad=logical(bad(:));
    assert(numel(q)==numel(bad) && all(isfinite(q)) && ...
        all(q>=0 & q<=1),'Invalid scores/labels.');
    edges=[qgrid(:); Inf];
    cb=histcounts(q(bad),edges).';
    cg=histcounts(q(~bad),edges).';
    acceptBad=flipud(cumsum(flipud(cb)));
    acceptGood=flipud(cumsum(flipud(cg)));
    assert(numel(acceptBad)==numel(qgrid));
    assert(acceptBad(1)==nnz(bad) && acceptGood(1)==nnz(~bad));
end

function plot_calibration(C,qstar,outfile)
    f=figure('Color','w','Position',[80 80 1100 480],'Visible','off');
    subplot(1,2,1);
    plot(C.Qthreshold,100*C.MeanBadDetection,'LineWidth',1.8);hold on;
    plot(C.Qthreshold,100*C.MeanGoodFalseRejection,'LineWidth',1.8);
    xline(qstar,'k--','LineWidth',1.2);
    xlabel('RMI threshold q');ylabel('Percent [%]');
    title('TRAIN: bad detection vs false rejection of good pixels');
    legend({'Detect wrong SWS','Reject correct SWS','Selected q'},...
        'Location','northwest');grid on;ylim([0 100]);
    subplot(1,2,2);
    plot(C.Qthreshold,C.MeanYoudenJ,'LineWidth',1.7);hold on;
    plot(C.Qthreshold,100*C.MeanCoverage/100,'LineWidth',1.7);
    xline(qstar,'k--','LineWidth',1.2);
    xlabel('RMI threshold q');ylabel('Score / fraction');
    title('TRAIN: balanced discrimination and accepted coverage');
    legend({'Macro Youden J','Mean coverage','Selected q'},...
        'Location','best');grid on;
    sgtitle(sprintf('RMI-1D cutoff selected ONLY on TRAIN | q*=%.3f',qstar));
    exportgraphics(f,outfile,'Resolution',160);close(f);
end

function plot_outcomes(T,outfile)
    rows=T(T.Region=="Global" & T.ErrorLimit_mps==0.5,:);
    n=height(rows);
    labels=cell(n,1);
    for k=1:n
        labels{k}=sprintf('%s / %s',char(rows.Phantom(k)),char(rows.Estimator(k)));
    end
    f=figure('Color','w','Position',[80 80 1200 780],'Visible','off');
    subplot(2,1,1);
    bar(100*[rows.BadDetectionRecall,rows.GoodFalseRejection]);
    ylabel('Percent [%]');title('Detection tradeoff (higher bad recall, lower good rejection)');
    legend({'Wrong SWS rejected','Correct SWS rejected'},'Location','best');
    grid on;xticks(1:n);xticklabels(labels);xtickangle(20);
    subplot(2,1,2);
    bar(100*[rows.Coverage,rows.FalseConfidenceAmongAccepted]);
    ylabel('Percent [%]');title('Coverage and error rate among accepted pixels');
    legend({'Accepted fraction','Error > 0.5 m/s among accepted'},...
        'Location','best');grid on;xticks(1:n);xticklabels(labels);xtickangle(20);
    sgtitle('ONE fixed RMI-1D decision cutoff on all cases (not per-case tuning)');
    exportgraphics(f,outfile,'Resolution',160);close(f);
end

function plot_high_contrast(V,qstar,eps_mps,outfile)
    f=figure('Color','w','Position',[70 80 1180 520],'Visible','off');
    x=V.x*1e3;z=V.z*1e3;
    good=V.common & (V.error<=eps_mps);
    bad=V.common & (V.error>eps_mps);
    accept=V.common & (V.q>=qstar);
    wrongAccepted=bad & accept;
    panels={V.error,V.q,double(accept),double(wrongAccepted)};
    titles={'Absolute SWS error [m/s]','Existing RMI-1D Q',...
        sprintf('Accepted: Q >= %.3f',qstar),'Wrong SWS still accepted'};
    for k=1:4
        subplot(1,4,k);
        img=panels{k};img(~V.common)=NaN;
        imagesc(x,z,img);axis image;axis xy;
        colormap(gca,'parula');colorbar;
        if k==1,caxis([0 4]);elseif k==2,caxis([0 1]);else,caxis([0 1]);end
        hold on;contour(x,z,double(V.inclusion),[0.5 0.5],'w','LineWidth',1.2);
        hold off;xlabel('x [mm]');ylabel('z [mm]');title(titles{k});
    end
    sgtitle(sprintf('HIGH CONTRAST / PG | fixed q*=%.3f | error > %.1f m/s',qstar,eps_mps));
    exportgraphics(f,outfile,'Resolution',160);close(f);
end
