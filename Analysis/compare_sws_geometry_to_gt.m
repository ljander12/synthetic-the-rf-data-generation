function T = compare_sws_geometry_to_gt(phantom_cases)
%COMPARE_SWS_GEOMETRY_TO_GT  Inspect actual PG/AC SWS versus material ground truth.
%
% PURPOSE: Determine whether the apparent inclusion shape in reconstructed
% SWS maps agrees with the physical inclusion geometry. This is a
% VALIDATION/DIAGNOSTIC, not a fourth quality metric or a new SWS estimator.
%
% QUICK START (place in <repo>/Analysis):
%   T = compare_sws_geometry_to_gt();
%   T = compare_sws_geometry_to_gt({'large_inclusions','high_contrast'});
%
% READ ONLY INPUTS (already generated; no k-Wave/Field II/USTB required):
%   Preferred: Results/method_comparison/batch/<case>/method_comparison.mat
%     Fields: PG, AC, GT, x_axis, z_axis, roi, common_valid.
%   Fallback: Results/sws_measurement_quality/<case>/rmi_1d_diagnostic/
%             rmi_1d_results.mat
%     Fields: cPG, cAC, GT, x, z, basic.
%   Optional cross-check: Data/<case>_displacement_with_quality.mat
%     Load ONLY sws_gt, shear_params, prb_theta (no large RF matrices).
%
% OUTPUT (does not replace or modify any original files):
%   Results/sws_estimator_diagnostic/sws_vs_gt_geometry/<case>/
%       sws_gt_overlay.png, sws_geometry_classification.png,
%       sws_geometry_summary.csv, sws_geometry_results.mat
%   Combined: Results/sws_estimator_diagnostic/sws_vs_gt_geometry/
%       sws_geometry_summary_all.csv
%
% IMPORTANT: The thresholded PG/AC areas are APPARENT contrast regions,
% not physical inclusion boundaries. A low Dice score does not by itself
% prove that the phantom or estimator is incorrectly programmed.
%
% No new simulation, RF generation, beamforming, or SWS reconstruction.

if nargin < 1 || isempty(phantom_cases)
    phantom_cases = {'large_inclusions','high_contrast'};
end
if ischar(phantom_cases) || isstring(phantom_cases)
    phantom_cases = cellstr(phantom_cases);
end
assert(iscell(phantom_cases) && ~isempty(phantom_cases), ...
    'Provide a phantom name or cell array of names.');

here = fileparts(mfilename('fullpath'));
if isfolder(fullfile(here,'Data')) && isfolder(fullfile(here,'Results'))
    repo = here;
else
    repo = fileparts(here);
end
assert(isfolder(fullfile(repo,'Results')), ...
    'Cannot find Results/. Place this script in <repo>/Analysis/.');

rootout = fullfile(repo,'Results','sws_estimator_diagnostic','sws_vs_gt_geometry');
if ~isfolder(rootout), mkdir(rootout); end
T = table();

for ci = 1:numel(phantom_cases)
    p = char(phantom_cases{ci});
    fprintf('\n====== ACTUAL SWS versus GT: %s ======\n',p);
    comparefile = fullfile(repo,'Results','method_comparison', ...
        'batch',p,'method_comparison.mat');
    rmifile = fullfile(repo,'Results','sws_measurement_quality', ...
        p,'rmi_1d_diagnostic','rmi_1d_results.mat');

    if isfile(comparefile)
        S = load(comparefile,'PG','AC','GT','x_axis','z_axis', ...
            'roi','common_valid');
        req = {'PG','AC','GT','x_axis','z_axis'};
        assert(all(isfield(S,req)), ...
            'Method comparison file lacks required SWS maps/axes: %s',comparefile);
        GT = double(S.GT); PG = double(S.PG); AC = double(S.AC);
        x = double(S.x_axis(:)); z = double(S.z_axis(:));
        source = comparefile;
        if isfield(S,'common_valid')
            valid_saved = logical(S.common_valid);
        else
            valid_saved = true(size(GT));
        end
        if isfield(S,'roi') && isfield(S.roi,'image_mask')
            roi_mask = logical(S.roi.image_mask);
        else
            roi_mask = true(size(GT));
        end
        clear S;
    elseif isfile(rmifile)
        S = load(rmifile,'cPG','cAC','GT','x','z','basic');
        req = {'cPG','cAC','GT','x','z'};
        assert(all(isfield(S,req)), ...
            'RMI result file lacks required SWS maps/axes: %s',rmifile);
        GT = double(S.GT); PG = double(S.cPG); AC = double(S.cAC);
        x = double(S.x(:)); z = double(S.z(:));
        source = rmifile;
        if isfield(S,'basic')
            valid_saved = logical(S.basic);
            roi_mask = logical(S.basic);
        else
            valid_saved = true(size(GT));
            roi_mask = true(size(GT));
        end
        clear S;
        warning('Method comparison file not found; using saved RMI input SWS maps for %s.',p);
    else
        error(['No saved PG/AC maps found for %s. Expected:\n%s\nor\n%s\n' ...
            'Run compare_SWS_methods_batch for this case if needed.'], ...
            p,comparefile,rmifile);
    end

    sz = size(GT);
    assert(isequal(size(PG),sz,size(AC),sz,size(valid_saved),sz, ...
        size(roi_mask),sz), ...
        'PG, AC, GT or mask sizes do not match for %s.',p);
    assert(numel(x)==sz(2) && numel(z)==sz(1), ...
        'Coordinate vector length mismatch for %s.',p);
    assert(all(diff(x)>0) && all(diff(z)>0), ...
        'Coordinate axes must be strictly increasing.');
    assert(max(abs([x(:);z(:)]))<1, ...
        'Expected ultrasound axes in METRES, but values look too large.');

    % Choose the original displacement file for optional small metadata
    % variables. Do NOT load uff objects or the large displacement cube.
    basename = [p '_displacement_with_quality.mat'];
    candidates = {fullfile(repo,'Data',p,basename), ...
        fullfile(repo,'Data',basename),fullfile(repo,basename), ...
        fullfile(repo,'Analysis',basename)};
    hit = find(cellfun(@isfile,candidates),1,'first');
    originalGTmatch = NaN; angledeg = NaN;
    c_bkg = NaN; c_incl = NaN;
    if ~isempty(hit)
        datafile = candidates{hit};
        W = whos('-file',datafile);
        names = {W.name};
        wanted = {'shear_params','sws_gt','prb_theta'};
        wanted = wanted(ismember(wanted,names));
        meta = load(datafile,wanted{:});
        if isfield(meta,'shear_params')
            sp = meta.shear_params;
            if isfield(sp,'c_shear_bkg'),c_bkg=double(sp.c_shear_bkg);end
            if isfield(sp,'c_shear_incl'),c_incl=double(sp.c_shear_incl);end
        end
        if isfield(meta,'prb_theta')
            th = double(meta.prb_theta);
            angledeg = rad2deg(th(1));
        end
        if isfield(meta,'sws_gt')
            srcgt = double(meta.sws_gt);
            if isequal(size(srcgt),sz)
                a = isfinite(GT) & isfinite(srcgt);
                if any(a(:))
                    originalGTmatch = max(abs(GT(a)-srcgt(a)));
                    if originalGTmatch>1e-6
                        warning(['GT in SWS results differs from GT in original ' ...
                            'dataset (max %.4g m/s). Inspect provenance!'], ...
                            originalGTmatch);
                    end
                end
            else
                warning('Original sws_gt has different grid size; cross-check skipped.');
            end
        end
        clear meta W;
    end

    % If original metadata is unavailable, derive two GT material levels
    % from the actual GT map. No estimator is involved in this threshold.
    gtFinite = GT(isfinite(GT));
    assert(~isempty(gtFinite), 'Ground-truth SWS contains no finite pixels.');
    if ~(isscalar(c_bkg) && isfinite(c_bkg) && ...
            isscalar(c_incl) && isfinite(c_incl) && c_bkg~=c_incl)
        levels = unique(gtFinite);
        assert(numel(levels)>1, ...
            'GT does not contain distinct material speeds; cannot outline inclusion.');
        c_bkg = double(levels(1));
        c_incl = double(levels(end));
        warning(['No distinct background/inclusion speeds in metadata for %s. ' ...
            'Using GT minimum/maximum to form contour threshold; ' ...
            'review material assignment.'],p);
    end
    boundary_cut = (c_bkg+c_incl)/2;
    % Explicitly select the class nearest the inclusion speed.
    gtInc = isfinite(GT) & (abs(GT-c_incl)<abs(GT-c_bkg));
    assert(any(gtInc(:)) && any(~gtInc(:)), ...
        'GT inclusion class was not found: check material parameters.');

    valid = isfinite(GT) & isfinite(PG) & isfinite(AC) & ...
        valid_saved & roi_mask & PG>0 & PG<10 & AC>0 & AC<10;
    assert(nnz(valid)>100, 'Fewer than 100 common valid ROI pixels for %s.',p);
    fprintf('SWS source: %s\n',source);
    fprintf('Material: background=%.4f, inclusion=%.4f m/s; threshold=%.4f\n', ...
        c_bkg,c_incl,boundary_cut);
    if isfinite(angledeg)
        fprintf('Probe angle in original dataset: %.1f degrees\n',angledeg);
    else
        fprintf('Probe angle: not stored in original displacement MAT\n');
    end
    if isfinite(originalGTmatch)
        fprintf('Saved GT vs original displacement GT max difference: %.4g m/s\n', ...
            originalGTmatch);
    end
    fprintf('GT inclusion pixels (all/ROI): %d / %d\n', ...
        nnz(gtInc),nnz(gtInc & valid));
    fprintf('COMMON VALID ROI: %d pixels\n',nnz(valid));

    outdir = fullfile(rootout,p);
    if ~isfolder(outdir),mkdir(outdir);end

    % Actual SWS maps with identical colour axis; WHITE = true inclusion.
    cmax = max(8,ceil(max([c_bkg,c_incl])));
    swsLimits = [0 cmax];
    ePG=PG-GT; eAC=AC-GT;
    fig = figure('Color','w','Visible','off','Position',[100 80 1500 850]);
    tl = tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
    draw_map(nexttile(tl,1),x,z,GT,swsLimits,gtInc,'Material ground truth');
    draw_map(nexttile(tl,2),x,z,PG,swsLimits,gtInc,'PG SWS | white = GT contour');
    draw_map(nexttile(tl,3),x,z,AC,swsLimits,gtInc,'AC SWS | white = GT contour');
    draw_map(nexttile(tl,4),x,z,double(valid),[0 1],[], ...
        'Fixed common evaluation mask');
    errlim = max(1, ceil(max(abs([ePG(valid);eAC(valid)]))));
    errlim = min(errlim,10);
    draw_map(nexttile(tl,5),x,z,ePG,[-errlim errlim],gtInc, ...
        'PG - GT [m/s]');
    draw_map(nexttile(tl,6),x,z,eAC,[-errlim errlim],gtInc, ...
        'AC - GT [m/s]');
    title(tl,sprintf('%s: ORIGINAL SWS maps vs TRUE geometry', ...
        strrep(p,'_',' ')),'Interpreter','none');
    exportgraphics(fig,fullfile(outdir,'sws_gt_overlay.png'),'Resolution',160);
    close(fig);

    % Thresholded maps show only where each estimator appears inclusion-like.
    % These are NOT ground-truth boundaries and are not used in SWS RMSE.
    if c_incl>c_bkg
        pgInc = PG>boundary_cut;
        acInc = AC>boundary_cut;
    else
        pgInc = PG<boundary_cut;
        acInc = AC<boundary_cut;
    end
    pgInc = pgInc & valid;
    acInc = acInc & valid;
    fig2=figure('Color','w','Visible','off','Position',[150 130 1400 460]);
    tl2=tiledlayout(fig2,1,3,'Padding','compact','TileSpacing','compact');
    draw_map(nexttile(tl2,1),x,z,double(gtInc & valid),[0 1], ...
        gtInc,'GT inclusion (in evaluation ROI)');
    draw_map(nexttile(tl2,2),x,z,double(pgInc),[0 1],gtInc, ...
        'Apparent PG inclusion (threshold)');
    draw_map(nexttile(tl2,3),x,z,double(acInc),[0 1],gtInc, ...
        'Apparent AC inclusion (threshold)');
    title(tl2,sprintf('Contrast classification only (threshold %.3f m/s) | %s', ...
        boundary_cut,strrep(p,'_',' ')),'Interpreter','none');
    exportgraphics(fig2,fullfile(outdir,'sws_geometry_classification.png'), ...
        'Resolution',160);
    close(fig2);

    % Numeric diagnostic: common mask for both estimators; no quality score.
    names = ["PG";"AC"];
    maps = {PG,AC};
    rows = table();
    for mi=1:2
        M=maps{mi};
        E=M-GT;
        g = valid;
        inc = valid & gtInc;
        bg  = valid & ~gtInc;
        assert(nnz(inc)>10 && nnz(bg)>10, ...
            'Not enough inclusion/background pixels for %s.',p);
        predicted = (M>boundary_cut);
        if c_incl<c_bkg, predicted = (M<boundary_cut); end
        predicted = predicted & valid;
        overlap=nnz(predicted & gtInc);
        dice = 2*overlap / max(1,(nnz(predicted)+nnz(gtInc & valid)));
        J = nnz(predicted & gtInc)/max(1,nnz((predicted | gtInc) & valid));
        meanInc=mean(M(inc));meanBg=mean(M(bg));
        contrastRecovery = (meanInc-meanBg)/(c_incl-c_bkg);
        r = table(string(p),names(mi),nnz(g),nnz(inc),nnz(bg), ...
            sqrt(mean(E(g).^2)),mean(E(g)),mean(abs(E(g))), ...
            sqrt(mean(E(inc).^2)),sqrt(mean(E(bg).^2)), ...
            mean(E(inc)),mean(E(bg)),meanInc,meanBg, ...
            contrastRecovery,dice,J,nnz(predicted), ...
            c_bkg,c_incl,boundary_cut,angledeg,originalGTmatch, ...
            string(source), ...
            'VariableNames',{'Phantom','Estimator','N_Valid','N_GT_Inclusion', ...
            'N_GT_Background','RMSE_mps','Bias_mps','MAE_mps', ...
            'InclusionRMSE_mps','BackgroundRMSE_mps', ...
            'InclusionBias_mps','BackgroundBias_mps', ...
            'MeanEstimatedInclusion_mps','MeanEstimatedBackground_mps', ...
            'ContrastRecovery','ThresholdDice','ThresholdIoU', ...
            'PredictedInclusionPixels','GT_Background_mps', ...
            'GT_Inclusion_mps','ClassificationThreshold_mps', ...
            'ProbeAngle_deg','OriginalGT_MaxDiff_mps','SourceFile'});
        rows = [rows;r]; %#ok<AGROW>
    end
    disp(rows(:,{'Estimator','RMSE_mps','InclusionRMSE_mps', ...
        'BackgroundRMSE_mps','ContrastRecovery','ThresholdDice'}));
    fprintf(['Note: ThresholdDice measures agreement of a simple speed-based ' ...
        'classification, NOT anatomical edge precision.\n']);
    writetable(rows,fullfile(outdir,'sws_geometry_summary.csv'));
    save(fullfile(outdir,'sws_geometry_results.mat'), ...
        'GT','PG','AC','x','z','valid','gtInc','pgInc','acInc', ...
        'rows','c_bkg','c_incl','boundary_cut','source','-v7.3');
    fprintf('Saved diagnostic to: %s\n',outdir);
    T = [T;rows]; %#ok<AGROW>
end
writetable(T,fullfile(rootout,'sws_geometry_summary_all.csv'));
fprintf('\nDONE. Original simulations and reconstruction results were NOT modified.\n');
end

function draw_map(ax,x,z,M,lims,outline,ttl)
imagesc(ax,x*1e3,z*1e3,M);
axis(ax,'image');
set(ax,'YDir','normal');
clim(ax,lims);
colorbar(ax);
title(ax,ttl,'Interpreter','none');
xlabel(ax,'Lateral [mm]');
ylabel(ax,'Depth [mm]');
if ~isempty(outline) && any(outline(:)) && any(~outline(:))
    hold(ax,'on');
    contour(ax,x*1e3,z*1e3,double(outline),[0.5 0.5], ...
        'LineColor','w','LineWidth',1.5);
    hold(ax,'off');
end
end
