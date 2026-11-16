function [T, globalSummary] = diagnose_rmi_component_discrimination(phantom_cases)
%DIAGNOSE_RMI_COMPONENT_DISCRIMINATION  Why does calibrated RMI-1D fail?
%
% Diagnostic of the EXISTING RMI-1D, not a fourth quality metric.
% Uses ONLY saved maps from validate_sws_rmi.m; no new reconstruction,
% k-Wave simulation, filter changes, or overwriting existing results.
%
% Comparison: correctly and incorrectly reconstructed SWS pixels, using
% ground truth solely for OFFLINE validation:
%     bad = abs(c_est - c_GT) > 0.5 m/s
%
% Four scores are compared:
%   Raw Eq.(6) 1D wave residual (LOWER = better)
%   Qres = clip((res_max - residual)/(res_max - res_min),0,1)
%   Qamp = clip((amp - u_min)/(u_max - u_min),0,1)
%   RMI = 0.8*Qres + 0.2*Qamp
%
% AUC = probability that a randomly chosen CORRECT pixel has better score
% than a randomly chosen WRONG pixel (ties count as 0.5). For raw residual,
% lower is better.  AUC = 0.5: chance; 1: perfect; <0.5: reversed ranking.
% AUC is a VALIDATION STATISTIC, not a new quality metric.
%
% Interpretation:
%   Raw-residual AUC > clipped-Qres AUC => clipping loses ranking ability.
%   Qres AUC high but RMI AUC low => amplitude weighting may be harmful.
%   Both raw residual and Qamp AUC low => retuning normalization alone
%   cannot fix the underlying lack of error-discriminating information.
%
% IMPORTANT: (1) GT-based region masks are for validation only.
% (2) In an all-bad or all-good region AUC is undefined: NaN, NOT zero.
% (3) large_inclusions is IN-SAMPLE for existing RMI normalization.
% (4) All phantoms have been examined in development: NOT independent test.
% (5) Spatially neighboring pixels are dependent; per-pixel AUC is NOT
%     a statistically independent estimate of out-of-sample performance.
%
% Put this file in <repo>/Analysis/ and run:
%   [T, summary] = diagnose_rmi_component_discrimination();
%
% Outputs under Results/sws_measurement_quality/rmi_component_discrimination/
%   component_discrimination_all.csv
%   component_auc_global.csv
%   component_auc_global.png
%   component_score_distributions_<phantom>.png (four figures)
%   component_discrimination_results.mat
%
% MATLAB R2020b+ recommended; no Statistics Toolbox required.

    if nargin < 1 || isempty(phantom_cases)
        phantom_cases = {'large_inclusions','high_contrast', ...
                         'low_contrast','chaoran_baseline'};
    end
    phantom_cases = cellstr(string(phantom_cases(:)));
    here = fileparts(mfilename('fullpath'));
    if isfolder(fullfile(here,'Results'))
        repo = here;
    else
        repo = fileparts(here);
    end
    root = fullfile(repo,'Results','sws_measurement_quality');
    assert(isfolder(root),'No RMI results folder: %s',root);
    output = fullfile(root,'rmi_component_discrimination');
    if ~isfolder(output), mkdir(output); end

    C = load(fullfile(root,'rmi_1d_calibration.mat'),'limits');
    assert(isfield(C,'limits'),'Missing RMI calibration limits.');
    lim = C.limits;
    neededLimits = {'res_min','res_max','u_min','u_max', ...
                    'reference_phantom','algorithm'};
    assert(all(isfield(lim,neededLimits)) && ...
        strcmp(lim.algorithm,'Choi2015_Eq6_RMI_1D_adaptation') && ...
        lim.res_max>lim.res_min && lim.u_max>lim.u_min, ...
        'Unsupported or invalid RMI-1D calibration file.');

    epsilon = 0.5;  % LOCKED GT error definition for this analysis [m/s]
    methods = {'PG','AC'};
    scoreNames = {'Raw residual','Qres','Qamp','RMI-1D'};
    regionNames = {'Global','InclusionCore','BackgroundInterior'};
    T = table();
    selftest_auc();
    plotCache = cell(numel(phantom_cases),2);

    fprintf('\n=== RMI-1D COMPONENT DISCRIMINATION ===\n');
    fprintf('Error = |estimated SWS - GT| > %.2f m/s. GT used ONLY for validation.\n',epsilon);
    fprintf('Original RMI normalization: %s; fixed Eq.(6) adaptation.\n',lim.reference_phantom);
    fprintf('Residual bounds [%.6g, %.6g], amplitude bounds [%.6g, %.6g].\n', ...
        lim.res_min,lim.res_max,lim.u_min,lim.u_max);
    fprintf('NO refitting, NO changes to existing SWS or RMI.\n');

    for pi=1:numel(phantom_cases)
        phantom = phantom_cases{pi};
        file = fullfile(root,phantom,'rmi_1d_diagnostic','rmi_1d_results.mat');
        assert(isfile(file),'Missing saved RMI components: %s',file);
        R = load(file,'resPG','resAC','urep','RmiPG','RmiAC', ...
            'cPG','cAC','GT','basic','masks','regions','limits','x','z');
        required = {'resPG','resAC','urep','RmiPG','RmiAC', ...
            'cPG','cAC','GT','basic','masks','regions','limits'};
        assert(all(isfield(R,required)),'Incomplete RMI MAT file for %s.',phantom);
        assert(all(isfield(R.limits,neededLimits)), ...
            'RMI file for %s has incomplete limits.',phantom);
        for f={'res_min','res_max','u_min','u_max'}
            name=f{1};
            assert(abs(double(R.limits.(name))-double(lim.(name))) < ...
                1e-8*max(1,abs(double(lim.(name)))), ...
                'Inconsistent calibration in %s: %s',phantom,name);
        end

        gt=double(R.GT);
        sz=size(gt);
        res={double(R.resPG),double(R.resAC)};
        speed={double(R.cPG),double(R.cAC)};
        rmiStored={double(R.RmiPG),double(R.RmiAC)};
        amp=double(R.urep);
        assert(isequal(size(amp),sz,size(R.basic),sz, ...
             size(res{1}),sz,size(res{2}),sz, ...
             size(speed{1}),sz,size(speed{2}),sz, ...
             size(rmiStored{1}),sz,size(rmiStored{2}),sz), ...
             'Map dimensions do not match in %s.',phantom);
        assert(numel(R.regions)==numel(R.masks), ...
            'Region mask/name count mismatch for %s.',phantom);
        rnames=string(R.regions(:));
        regionMasks=cell(1,numel(regionNames));
        regionMasks{1}=true(sz);
        for ri=2:numel(regionNames)
            idx=find(rnames==string(regionNames{ri}),1);
            assert(~isempty(idx),'Missing %s ROI in %s',regionNames{ri},phantom);
            assert(isequal(size(R.masks{idx}),sz), ...
                'Invalid %s ROI map size.',regionNames{ri});
            regionMasks{ri}=logical(R.masks{idx});
        end

        qamp=clip01((amp-double(lim.u_min))/ ...
            (double(lim.u_max)-double(lim.u_min)));
        qres=cell(1,2);
        allValid=logical(R.basic) & isfinite(gt) & isfinite(qamp);
        for mi=1:2
            qres{mi}=clip01((double(lim.res_max)-res{mi})/ ...
                (double(lim.res_max)-double(lim.res_min)));
            reconstructed=0.8*qres{mi}+0.2*qamp;
            valid=allValid & isfinite(speed{mi}) & isfinite(res{mi}) & ...
                isfinite(qres{mi}) & isfinite(rmiStored{mi}) & ...
                speed{mi}>0 & speed{mi}<10;
            assert(nnz(valid)>100,'Too few RMI pixels for %s / %s', ...
                phantom,methods{mi});
            delta=max(abs(reconstructed(valid)-rmiStored{mi}(valid)));
            fprintf('%s/%s RMI reconstruction max diff: %.3g\n',phantom,methods{mi},delta);
            assert(delta<1e-6, ...
                'Reconstructed RMI does not match saved RMI in %s/%s.', ...
                phantom,methods{mi});
            allValid=valid;
        end
        assert(nnz(allValid)>1000,'Too few same-pixel comparisons in %s.',phantom);

        fprintf('\n%s: %d shared valid pixels\n',phantom,nnz(allValid));
        for mi=1:2
            errorMap=abs(speed{mi}-gt);
            scoreMaps={res{mi},qres{mi},qamp,rmiStored{mi}};
            for ri=1:numel(regionNames)
                valid=allValid & regionMasks{ri};
                if nnz(valid)<30
                    fprintf('Skipping %s / %s: too few pixels.\n', ...
                        methods{mi},regionNames{ri});
                    continue;
                end
                err=errorMap(valid);
                bad=err>epsilon;
                for qi=1:numel(scoreNames)
                    s=scoreMaps{qi}(valid);
                    if qi==1
                        auc=auc_ties(-s,~bad); % lower raw residual is better
                        direction="LowerIsBetter";
                        atZero=NaN; atOne=NaN;
                    else
                        auc=auc_ties(s,~bad);
                        direction="HigherIsBetter";
                        atZero=mean(s<=1e-10);
                        atOne=mean(s>=1-1e-10);
                    end
                    goodMedian=nanmedian_local(s(~bad));
                    badMedian=nanmedian_local(s(bad));
                    row=table(string(phantom),string(methods{mi}), ...
                        string(regionNames{ri}),string(scoreNames{qi}), ...
                        direction,numel(s),nnz(bad),nnz(~bad), ...
                        mean(bad),auc,goodMedian,badMedian,atZero,atOne, ...
                        'VariableNames',{'Phantom','Estimator','Region', ...
                        'Component','ScoreDirection','N','NBad','NGood', ...
                        'BadPrevalence','AUC','GoodMedian','BadMedian', ...
                        'FractionAtZero','FractionAtOne'});
                    T=[T;row]; %#ok<AGROW>
                end
            end
            plotCache{pi,mi}=struct('res',res{mi}(allValid), ...
                'qres',qres{mi}(allValid),'qamp',qamp(allValid), ...
                'rmi',rmiStored{mi}(allValid), ...
                'bad',errorMap(allValid)>epsilon);
        end
        plot_case_distributions(output,phantom,plotCache(pi,:),methods);
    end

    globalSummary=table();
    for pi=1:numel(phantom_cases)
        for mi=1:2
            sel=T.Phantom==string(phantom_cases{pi}) & ...
                T.Estimator==string(methods{mi}) & T.Region=="Global";
            V=T(sel,:);
            assert(height(V)==4, ...
                'Expected 4 components for global %s/%s', ...
                phantom_cases{pi},methods{mi});
            aucs=nan(1,4); saturations=nan(1,4);
            for qi=1:4
                k=find(V.Component==string(scoreNames{qi}),1);
                assert(~isempty(k),'Missing component %s',scoreNames{qi});
                aucs(qi)=V.AUC(k);
                saturations(qi)=V.FractionAtOne(k);
            end
            row=table(string(phantom_cases{pi}),string(methods{mi}), ...
                V.NBad(1),V.NGood(1),V.BadPrevalence(1), ...
                aucs(1),aucs(2),aucs(3),aucs(4), ...
                aucs(1)-aucs(2),saturations(2),saturations(3), ...
                'VariableNames',{'Phantom','Estimator','NBad','NGood', ...
                'BadPrevalence','AUCRawResidual','AUCQres','AUCQamp', ...
                'AUCRMI','AUCClippingLoss','QresSaturationAtOne', ...
                'QampSaturationAtOne'});
            globalSummary=[globalSummary;row]; %#ok<AGROW>
        end
    end

    writetable(T,fullfile(output,'component_discrimination_all.csv'));
    writetable(globalSummary,fullfile(output,'component_auc_global.csv'));
    plot_global_auc(globalSummary,output);
    save(fullfile(output,'component_discrimination_results.mat'), ...
        'T','globalSummary','lim','epsilon','phantom_cases','-v7');

    fprintf('\n=== GLOBAL COMPONENT AUC (higher = better) ===\n');
    disp(globalSummary(:,{'Phantom','Estimator','AUCRawResidual', ...
        'AUCQres','AUCQamp','AUCRMI','AUCClippingLoss', ...
        'QresSaturationAtOne'}));
    fprintf('\n=== INCLUSION-CORE COMPONENT AUC ===\n');
    inc=T(T.Region=="InclusionCore",:);
    disp(inc(:,{'Phantom','Estimator','Component','NBad', ...
        'NGood','AUC','GoodMedian','BadMedian'}));
    fprintf('\nSaved diagnostics in %s\n',output);
    fprintf('AUC=NaN when a region has no correct OR no incorrect pixels.\n');
    fprintf('Raw residual and Qres use identical underlying data: AUC differences reflect clipping/ties.\n');
    fprintf('This is EXPLORATORY, not untouched phantom-level validation.\n');
end

function selftest_auc()
    % Minimal checks for AUC direction and exact ties.
    assert(abs(auc_ties([1 2 3 4],[false false true true])-1)<1e-12, ...
        'AUC orientation self-test failed.');
    assert(abs(auc_ties([1 1 2 2],[false true false true])-0.5)<1e-12, ...
        'AUC ties self-test failed.');
    assert(isnan(auc_ties([1 2 3],[true true true])), ...
        'All-good region must give undefined AUC.');
end

function q=clip01(v)
    q=max(0,min(1,double(v)));
end

function m=nanmedian_local(v)
    v=double(v(isfinite(v)));
    if isempty(v),m=NaN; else,m=median(v);end
end

function auc=auc_ties(s,good)
% Efficient rank-based ROC AUC with exact tie handling; no toolbox needed.
    s=double(s(:)); good=logical(good(:));
    ok=isfinite(s);
    s=s(ok);good=good(ok);
    ng=nnz(good);nb=nnz(~good);
    if ng==0 || nb==0
        auc=NaN;return;
    end
    [sorted,ix]=sort(s);
    starts=find([true;diff(sorted)~=0]);
    ends=[starts(2:end)-1;numel(s)];
    groups=cumsum([true;diff(sorted)~=0]);
    rankSorted=(starts(groups)+ends(groups))/2;
    rankOriginal=zeros(numel(s),1);
    rankOriginal(ix)=rankSorted;
    auc=(sum(rankOriginal(good))-ng*(ng+1)/2)/(ng*nb);
    assert(auc>=-1e-9 && auc<=1+1e-9,'AUC numeric error.');
    auc=min(1,max(0,auc));
end

function plot_case_distributions(out,phantom,record,methods)
    f=figure('Visible','off','Color','w','Position',[80 60 1560 750]);
    tl=tiledlayout(f,2,4,'TileSpacing','compact','Padding','compact');
    labels={'Raw residual (log10)','Qres','Qamp','RMI-1D'};
    dataField={'res','qres','qamp','rmi'};
    for mi=1:2
        d=record{mi};
        for qi=1:4
            ax=nexttile(tl);
            values=double(d.(dataField{qi})(:));
            bad=logical(d.bad(:));
            if qi==1
                values=log10(max(values,1e-16));
                edgeLow=quantile_simple(values,2);
                edgeHigh=quantile_simple(values,98);
                if ~isfinite(edgeLow) || ~isfinite(edgeHigh) || ...
                        edgeHigh<=edgeLow
                    edgeLow=min(values)-0.5;
                    edgeHigh=max(values)+0.5;
                end
                edges=linspace(edgeLow,edgeHigh,35);
                xlabel(ax,'log10(raw residual)');
            else
                edges=linspace(0,1,41);
                xlabel(ax,'Existing quality score');
            end
            centers=(edges(1:end-1)+edges(2:end))/2;
            a=histcounts(values(~bad),edges,'Normalization','probability');
            b=histcounts(values(bad),edges,'Normalization','probability');
            % Normalize within each class, so unequal prevalence is not
            % confused with component separability.
            plot(ax,centers,a,'LineWidth',1.4); hold(ax,'on');
            plot(ax,centers,b,'LineWidth',1.4); hold(ax,'off');
            grid(ax,'on');
            ylabel(ax,'Fraction within class / bin');
            title(ax,sprintf('%s / %s',methods{mi},labels{qi}), ...
                'Interpreter','none');
            if nnz(~bad)==0
                text(ax,.05,.87,'NO correct SWS pixels', ...
                    'Units','normalized','FontSize',9);
            elseif nnz(bad)==0
                text(ax,.05,.87,'NO incorrect SWS pixels', ...
                    'Units','normalized','FontSize',9);
            end
            if mi==1 && qi==1
                legend(ax,{'Correct SWS','Wrong SWS'},'Location','best');
            end
        end
    end
    title(tl,sprintf('%s | residual vs normalized scores | error > 0.5 m/s', ...
        strrep(phantom,'_',' ')));
    exportgraphics(f,fullfile(out, ...
        ['component_score_distributions_' phantom '.png']), ...
        'Resolution',150);
    close(f);
end

function p=quantile_simple(v,pct)
    v=sort(double(v(isfinite(v))));
    if isempty(v),p=NaN;return;end
    k=1+(numel(v)-1)*pct/100;
    i=floor(k);j=ceil(k);
    p=v(i)+(k-i)*(v(j)-v(i));
end

function plot_global_auc(S,out)
    f=figure('Visible','off','Color','w','Position',[70 70 1450 620]);
    y=[S.AUCRawResidual,S.AUCQres,S.AUCQamp,S.AUCRMI];
    bar(y,'grouped');
    ylim([0 1]); grid on;hold on;
    yline(0.5,'k--','Chance (0.5)');
    hold off;
    xticks(1:height(S));
    tick=cell(height(S),1);
    for i=1:height(S)
        tick{i}=sprintf('%s / %s',char(strrep(S.Phantom(i),'_',' ')), ...
            char(S.Estimator(i)));
    end
    xticklabels(tick);xtickangle(22);
    legend({'Raw residual','Qres (clipped)','Qamp','RMI-1D'}, ...
        'Location','northoutside','Orientation','horizontal');
    ylabel('AUC: correct SWS ranked above incorrect SWS');
    title('GLOBAL | error > 0.5 m/s | 0.5 = chance; 1 = perfect');
    exportgraphics(f,fullfile(out,'component_auc_global.png'), ...
        'Resolution',160);
    close(f);
end
