function T = validate_sws_quality_regions(phantom_case)
%VALIDATE_SWS_QUALITY_REGIONS Regional evaluation of two quality-map pilots.
%
% READ ONLY:
%  Results/sws_measurement_quality/<case>/phase_consistency_diagnostic/phase_quality_results.mat
%  Results/method_comparison/batch/<case>/method_comparison.mat
%
% OUTPUT (new folder only):
%  Results/sws_measurement_quality/<case>/regional_quality_diagnostic/
%       regional_summary.csv, regional_summary.mat, roi_regions.png,
%       regional_spearman.png
%
% Q metrics use displacement only. GT/ROI are used ONLY for offline validation.
% Pixel-level associations are descriptive because nearby pixels correlate.
% Block-level correlations are computed on spatially aggregated 1-mm cells.
%
% Usage:
%   validate_sws_quality_regions('large_inclusions')
%   validate_sws_quality_regions('high_contrast')

    if nargin < 1 || isempty(phantom_case)
        phantom_case = 'large_inclusions';
    end
    phantom_case = char(string(phantom_case));
    here = fileparts(mfilename('fullpath'));
    if isfolder(fullfile(here,'Results'))
        repo = here;
    else
        repo = fileparts(here);
    end

    phase_file = fullfile(repo,'Results','sws_measurement_quality', ...
        phantom_case,'phase_consistency_diagnostic','phase_quality_results.mat');
    comparison_file = fullfile(repo,'Results','method_comparison', ...
        'batch',phantom_case,'method_comparison.mat');
    assert(isfile(phase_file), 'Missing phase results: %s',phase_file);
    assert(isfile(comparison_file), 'Missing method comparison: %s',comparison_file);

    A = load(phase_file,'Qold','Qphase','common','GT','PG','AC','x','z');
    B = load(comparison_file,'roi','GT');
    needed = {'Qold','Qphase','common','GT','PG','AC','x','z'};
    assert(all(isfield(A,needed)), 'Phase results MAT lacks expected variables.');
    assert(isfield(B,'roi') && isfield(B,'GT'), ...
        'Method comparison MAT lacks roi or GT.');

    GT = double(A.GT);
    PG = double(A.PG);
    AC = double(A.AC);
    Qold = double(A.Qold);
    Qphase = double(A.Qphase);
    common = logical(A.common);
    roi = B.roi;
    x = double(A.x(:)); z = double(A.z(:));
    [Nz,Nx] = size(GT);
    assert(isequal(size(GT),size(PG),size(AC),size(Qold), ...
        size(Qphase),size(common),size(B.GT)), ...
        'SWS/Q images do not have matching sizes.');
    assert(max(abs(GT(:)-double(B.GT(:))),[],'omitnan') < 1e-6, ...
        'GT from phase diagnostic differs from batch comparison GT.');
    assert(all(isfield(roi,{'inclusion_core','background_core','image_mask'})), ...
        'ROI structure lacks required masks.');
    assert(isequal(size(roi.image_mask),[Nz Nx]) && ...
        isequal(size(roi.inclusion_core),[Nz Nx]) && ...
        isequal(size(roi.background_core),[Nz Nx]), ...
        'ROI masks have wrong dimensions.');
    assert(numel(x)==Nx && numel(z)==Nz,'Axis dimensions disagree.');

    % Fixed before viewing regional validation results.
    edge_width_mm = 3.5; % lateral strip measured from scan outer boundary
    block_width_mm = 1.0; % aggregate into 1x1 mm spatial blocks
    quality_threshold = 0.8;
    min_region_pixels = 30;

    [X,Z] = meshgrid(x,z);
    valid = common & isfinite(GT) & isfinite(PG) & isfinite(AC) & ...
        isfinite(Qold) & isfinite(Qphase);
    assert(nnz(valid)>100,'No sufficient common pixels available.');

    % Existing masks for inclusion/background; separate lateral-boundary
    % region is explicitly limited to GT-background material.
    finiteGT = GT(isfinite(GT));
    bkg_speed = min(finiteGT);
    gt_bkg = isfinite(GT) & abs(GT - bkg_speed) < 1e-5;
    near_side = X <= (min(x)+edge_width_mm*1e-3) | ...
                X >= (max(x)-edge_width_mm*1e-3);

    masks = {valid & logical(roi.inclusion_core), ...
             valid & logical(roi.background_core) & ~near_side, ...
             valid & logical(roi.image_mask) & gt_bkg & near_side};
    region_names = {'InclusionCore','BackgroundInterior','LateralBoundary'};
    for r=1:numel(masks)
        fprintf('%-21s: %d valid pixels\n',region_names{r},nnz(masks{r}));
        if nnz(masks{r}) < min_region_pixels
            warning('Region %s has too few pixels; statistics will be NaN.', ...
                region_names{r});
        end
    end
    assert(nnz(masks{1})>=min_region_pixels && ...
        nnz(masks{2})>=min_region_pixels, ...
        'Too few valid inclusion/background pixels.');

    errors = {abs(PG-GT),abs(AC-GT)};
    methods = {'PG','AC'};
    Qs = {Qold,Qphase};
    q_names = {'CorrelationQ','PhasePlanarityQ'};
    T = table();

    % Group pixels into fixed spatial cells to reduce (not eliminate)
    % pseudo-replication from neighboring ultrasound pixels.
    bx = floor((X-min(x))/(block_width_mm*1e-3));
    bz = floor((Z-min(z))/(block_width_mm*1e-3));
    nxblock = max(bx(:))+1;
    block_id = bz*nxblock+bx+1;

    for r=1:numel(masks)
        mask = masks{r};
        for qi=1:2
            for m=1:2
                q = double(Qs{qi}(mask));
                e = double(errors{m}(mask));
                ids = double(block_id(mask));
                if numel(q)<min_region_pixels
                    values = nan(1,11);
                else
                    hi = q>=quality_threshold;
                    lo = q<quality_threshold;
                    [blkq,blke] = grouped_means(ids,q,e);
                    values = [numel(q), mean(q), std(q), ...
                        mean(q>=1-1e-6), ...
                        spearman_safe(q,e), ...
                        numel(blkq), spearman_safe(blkq,blke), ...
                        nnz(hi), nnz(lo), ...
                        safe_mean(e(hi)),safe_mean(e(lo))];
                end
                row = table(string(phantom_case),string(region_names{r}), ...
                    string(methods{m}),string(q_names{qi}), ...
                    values(1),values(2),values(3),values(4), ...
                    values(5),values(6),values(7),values(8),values(9), ...
                    values(10),values(11), ...
                    'VariableNames',{'Phantom','Region','Estimator','Quality', ...
                    'ValidPixels','MeanQ','StdQ','FracQSaturated', ...
                    'PixelSpearman','SpatialBlocks','BlockSpearman', ...
                    'N_Qge08','N_Qlt08','MAE_Qge08_mps','MAE_Qlt08_mps'});
                T = [T;row]; %#ok<AGROW>
            end
        end
    end
    disp(T);

    outdir = fullfile(repo,'Results','sws_measurement_quality', ...
        phantom_case,'regional_quality_diagnostic');
    if ~isfolder(outdir), mkdir(outdir); end
    writetable(T,fullfile(outdir,'regional_summary.csv'));
    save(fullfile(outdir,'regional_summary.mat'),'T','masks', ...
        'region_names','edge_width_mm','block_width_mm', ...
        'quality_threshold','-v7.3');

    % Plot region definitions. White pixels are not included in any region.
    region_map = nan(Nz,Nx);
    for r=1:3
        region_map(masks{r})=r;
    end
    f1=figure('Color','w','Name',['Regional masks: ' phantom_case]);
    imagesc(x*1e3,z*1e3,region_map);
    set(gca,'YDir','normal'); axis image;
    colormap(gca,[0.28 0.64 0.43;0.27 0.49 0.82;0.86 0.38 0.23]);
    clim([0.5 3.5]);
    c=colorbar; c.Ticks=1:3; c.TickLabels=region_names;
    xlabel('Lateral [mm]');ylabel('Depth [mm]');
    title(['Validation regions: ' strrep(phantom_case,'_',' ')]);
    exportgraphics(f1,fullfile(outdir,'roi_regions.png'),'Resolution',160);

    % Summarize regional associations. Negative Spearman is desirable.
    f2=figure('Color','w','Name',['Regional Q evaluation: ' phantom_case]);
    tl=tiledlayout(f2,1,2,'TileSpacing','compact','Padding','compact');
    metric_labels={'Correlation Q / PG','Correlation Q / AC', ...
                   'Phase Q / PG','Phase Q / AC'};
    for k=1:2
        ax=nexttile(tl);
        V=nan(3,4);
        for r=1:3
            for qi=1:2
                for m=1:2
                    idx=strcmp(T.Region,region_names{r}) & ...
                        strcmp(T.Quality,q_names{qi}) & ...
                        strcmp(T.Estimator,methods{m});
                    col=(qi-1)*2+m;
                    if nnz(idx)==1
                        if k==1, V(r,col)=T.PixelSpearman(idx);
                        else, V(r,col)=T.BlockSpearman(idx); end
                    end
                end
            end
        end
        imagesc(ax,V,[-1 1]);
        colormap(ax,parula);
        set(ax,'XTick',1:4,'XTickLabel',metric_labels, ...
            'XTickLabelRotation',25,'YTick',1:3, ...
            'YTickLabel',region_names);
        colorbar(ax);
        if k==1, title(ax,'Pixel-level Spearman');
        else, title(ax,'1-mm block Spearman'); end
        for r=1:3
            for j=1:4
                if isfinite(V(r,j))
                    text(ax,j,r,sprintf('%+.2f',V(r,j)), ...
                        'HorizontalAlignment','center','Color','w', ...
                        'FontWeight','bold');
                end
            end
        end
    end
    title(tl,['Regional association: ' strrep(phantom_case,'_',' ')]);
    exportgraphics(f2,fullfile(outdir,'regional_spearman.png'), ...
        'Resolution',160);

    fprintf('\nSaved regional outputs to: %s\n',outdir);
    fprintf(['Negative Spearman(Q, |SWS-GT|) is desirable, but neither ' ...
        'pixel Spearman nor block Spearman proves reliable confidence.\n']);
    fprintf('Q >= %.1f groups: inspect N_Qge08 and N_Qlt08 before comparing MAE.\n', ...
        quality_threshold);
end

function v = safe_mean(a)
    if isempty(a), v=NaN; else, v=mean(a,'omitnan'); end
end

function rho = spearman_safe(a,b)
    a=a(:);b=b(:);
    good=isfinite(a)&isfinite(b);
    a=a(good);b=b(good);
    if numel(a)<10 || numel(unique(a))<2 || numel(unique(b))<2
        rho=NaN;return;
    end
    rho=corr(a,b,'Type','Spearman','Rows','complete');
end

function [q_b,e_b]=grouped_means(ids,q,e)
    [~,~,g]=unique(ids);
    sums_q=accumarray(g,q,[],@mean);
    sums_e=accumarray(g,e,[],@mean);
    n=accumarray(g,1);
    keep=n>=5;
    q_b=sums_q(keep);
    e_b=sums_e(keep);
end
