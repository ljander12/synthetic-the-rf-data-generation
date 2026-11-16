function summary = test_sws_measurement_quality(phantom_cases)
%TEST_SWS_MEASUREMENT_QUALITY Pilot of Schoen et al. (2025) SWS quality Q.
%
% Research source:
%   Schoen et al., Ultrasound in Medicine & Biology 51(7), 1070-1077 (2025).
%   DOI: 10.1016/j.ultrasmedbio.2025.03.003
%
% IMPORTANT: This is an EXPLORATORY ADAPTATION, not an exact reproduction.
%   - Spatially adjacent lateral displacement waveforms are compared using
%     ZERO-LAG, temporally mean-centered normalized cross-correlation.
%   - The correlation magnitude is transformed with the published piecewise
%     Q(C12) mapping. Lag selection and original proprietary processing may
%     differ from this prototype.
%   - Input Q requires NO ground truth. Ground truth is used only for
%     independent validation against SWS reconstruction error.
%   - Harmonic 200 Hz motion can produce near-perfect correlations even when
%     SWS is biased. Check the printed Q saturation fraction!
%
% Place file in PROJECT_ROOT/Analysis/ and call in MATLAB:
%   test_sws_measurement_quality('large_inclusions')
%   test_sws_measurement_quality(["large_inclusions","high_contrast", ...
%       "low_contrast","chaoran_baseline"])
%
% Required existing files per phantom:
%   Data/[case]/[case]_displacement_with_quality.mat
%   Results/method_comparison/batch/[case]/method_comparison.mat
%
% Writes NEW outputs only to:
%   Results/sws_measurement_quality/[case]/
% Does not change the existing MATLAB data or analysis scripts.

    if nargin < 1 || isempty(phantom_cases)
        phantom_cases = "large_inclusions";
    end
    phantom_cases = string(phantom_cases(:));
    script_dir = fileparts(mfilename('fullpath'));
    if isfolder(fullfile(script_dir, 'Results'))
        repo_dir = script_dir;
    else
        repo_dir = fileparts(script_dir);
    end
    results_dir = fullfile(repo_dir, 'Results');

    % IMPORTANT: displacement_estimation_autocorr was saved as an
    % uff.beamformed_data object. Its class must be available BEFORE load().
    % If USTB is installed elsewhere, enter its root directory below.
    ustb_custom_dir = '';  % e.g., 'C:\path\to\USTB' (containing +uff)
    ensure_ustb_beamformed_class(repo_dir, ustb_custom_dir);

    out_root = fullfile(results_dir, 'sws_measurement_quality');
    if ~isfolder(out_root), mkdir(out_root); end

    summary = table();
    for icase = 1:numel(phantom_cases)
        name = char(phantom_cases(icase));
        fprintf('\n===== QUALITY PILOT: %s =====\n', name);

        comparison_file = fullfile(results_dir, 'method_comparison', ...
            'batch', name, 'method_comparison.mat');
        assert(isfile(comparison_file), ...
            'Missing saved SWS maps: %s', comparison_file);
        S = load(comparison_file, 'GT', 'PG', 'AC', 'roi', ...
            'x_axis', 'z_axis');
        assert(all(isfield(S, {'GT','PG','AC','roi','x_axis','z_axis'})), ...
            'Comparison MAT is missing required variables.');
        GT = double(S.GT);
        PG = double(S.PG);
        AC = double(S.AC);
        x_axis = double(S.x_axis(:));
        z_axis = double(S.z_axis(:));
        [Nz, Nx] = size(GT);
        assert(isequal(size(PG),size(GT),size(AC)), ...
            'SWS map dimensions differ.');
        assert(Nz == numel(z_axis) && Nx == numel(x_axis), ...
            'Scan axes do not match SWS maps.');

        filename = [name '_displacement_with_quality.mat'];
        candidate_files = {fullfile(repo_dir,'Data',name,filename), ...
            fullfile(repo_dir,'Data',filename), ...
            fullfile(repo_dir,'Analysis',filename), ...
            fullfile(repo_dir,filename)};
        displacement_file = '';
        for j = 1:numel(candidate_files)
            if isfile(candidate_files{j})
                displacement_file = candidate_files{j};
                break;
            end
        end
        assert(~isempty(displacement_file), ...
            'Cannot locate displacement dataset %s.', filename);
        fprintf('Reading displacement from %s\n', displacement_file);

        % Only load the estimated displacement object; do not load b_data
        % or k-Wave wavefields into memory unnecessarily.
        D = load(displacement_file, 'displacement_estimation_autocorr');
        assert(isfield(D,'displacement_estimation_autocorr'), ...
            'Dataset lacks displacement_estimation_autocorr.');
        if ~isobject(D.displacement_estimation_autocorr) && ...
                ~isstruct(D.displacement_estimation_autocorr)
            error(['The displacement was not restored as a USTB object ' ...
                '(loaded as %s). Ensure the same compatible USTB version ' ...
                'is on the MATLAB path, then rerun.'], ...
                class(D.displacement_estimation_autocorr));
        end
        u_raw = D.displacement_estimation_autocorr.data;
        assert(isnumeric(u_raw) && isreal(u_raw), ...
            'Expected real numeric displacement data.');
        assert(mod(numel(u_raw), Nz*Nx) == 0, ...
            'Displacement dimensions do not match SWS grid.');
        Nt = numel(u_raw)/(Nz*Nx);
        assert(Nt >= 5, 'Need at least five displacement frames.');
        u = reshape(single(u_raw), Nz, Nx, Nt);
        clear D u_raw;
        fprintf('Displacement size: %d x %d x %d\n', Nz, Nx, Nt);

        [Q, C12dB, rho] = quality_map_from_displacement(u);
        clear u;

        % Common evaluation ROI. GT is not used in quality_map_from_displacement.
        valid = logical(S.roi.image_mask) & isfinite(GT) & ...
            isfinite(PG) & PG > 0 & PG < 10 & ...
            isfinite(AC) & AC > 0 & AC < 10 & isfinite(Q);
        assert(nnz(valid) > 10, ...
            'Too few common valid Q / SWS pixels for %s.', name);
        qv = Q(valid);
        fprintf('Valid pixels: %d; Q >= 0.63: %.2f%%; Q=1: %.2f%%\n', ...
            nnz(valid),100*mean(qv>=0.63),100*mean(qv>=1-1e-10));
        if mean(qv>=1-1e-10) > 0.90
            warning(['More than 90%% of valid pixels have Q=1. ' ...
                'The published transformation is saturated for this ' ...
                'harmonic dataset and may not discriminate SWS errors.']);
        end

        out_dir = fullfile(out_root, name);
        if ~isfolder(out_dir), mkdir(out_dir); end
        sws_maps = {PG, AC};
        methods = {"PG", "AC"};
        case_rows = table();
        for k = 1:2
            err = abs(sws_maps{k} - GT);
            [r_spear, mae_hi, mae_lo] = evaluate_quality(Q, err, valid);
            row = table(string(name), methods{k}, nnz(valid), ...
                mean(qv,'omitnan'), mean(qv>=0.45), ...
                mean(qv>=0.63), mean(qv>=1-1e-10), ...
                r_spear, mae_hi, mae_lo, ...
                'VariableNames', {'Phantom','Method','ValidPixels', ...
                'MeanQ','FractionQAbove045','FractionQAbove063', ...
                'FractionQSaturated','Spearman_Q_vs_AbsError', ...
                'MAE_QHigh_mps','MAE_QLow_mps'});
            if isempty(case_rows)
                case_rows = row;
            else
                case_rows = [case_rows;row]; %#ok<AGROW>
            end
        end
        disp(case_rows);
        writetable(case_rows, fullfile(out_dir,'quality_validation.csv'));
        save(fullfile(out_dir,'quality_validation.mat'), ...
            'Q','C12dB','rho','GT','PG','AC','valid', ...
            'x_axis','z_axis','case_rows','-v7.3');

        fig = figure('Name',['SWS Q pilot: ' name], ...
            'Color','w','NumberTitle','off');
        layout = tiledlayout(fig,1,3,'TileSpacing','compact', ...
            'Padding','compact');
        ax1 = nexttile(layout);
        imagesc(ax1,x_axis*1e3,z_axis*1e3,Q);
        axis(ax1,'image'); set(ax1,'YDir','normal');
        clim(ax1,[0 1]); colorbar(ax1);
        title(ax1,'Displacement-based Q');
        xlabel(ax1,'Lateral [mm]'); ylabel(ax1,'Depth [mm]');

        ax2 = nexttile(layout);
        imagesc(ax2,x_axis*1e3,z_axis*1e3,abs(PG-GT));
        axis(ax2,'image'); set(ax2,'YDir','normal');
        clim(ax2,[0 4]); colorbar(ax2);
        title(ax2,'|PG - GT| [m/s]');
        xlabel(ax2,'Lateral [mm]'); ylabel(ax2,'Depth [mm]');

        ax3 = nexttile(layout);
        imagesc(ax3,x_axis*1e3,z_axis*1e3,abs(AC-GT));
        axis(ax3,'image'); set(ax3,'YDir','normal');
        clim(ax3,[0 4]); colorbar(ax3);
        title(ax3,'|AC - GT| [m/s]');
        xlabel(ax3,'Lateral [mm]'); ylabel(ax3,'Depth [mm]');
        title(layout,['SWS measurement quality pilot: ' strrep(name,'_',' ')]);
        exportgraphics(fig,fullfile(out_dir,'quality_vs_error.png'), ...
            'Resolution',160);

        fig2 = figure('Name',['Q versus SWS error: ' name], ...
            'Color','w','NumberTitle','off');
        err_pg = abs(PG(valid)-GT(valid));
        err_ac = abs(AC(valid)-GT(valid));
        nplot = min(numel(qv),5000);
        idx = unique(round(linspace(1,numel(qv),nplot)));
        plot(qv(idx),err_pg(idx),'.','MarkerSize',4);
        hold on;
        plot(qv(idx),err_ac(idx),'.','MarkerSize',4);
        hold off;
        xlim([0 1]); xlabel('Q (no GT used)');
        ylabel('Absolute SWS error [m/s] (GT used for validation)');
        legend('Phase gradient','Spatial autocorrelation', ...
            'Location','best'); grid on;
        title(['Q versus SWS error: ' strrep(name,'_',' ')]);
        exportgraphics(fig2,fullfile(out_dir,'quality_vs_error_scatter.png'), ...
            'Resolution',160);

        if isempty(summary), summary=case_rows;
        else, summary=[summary;case_rows]; %#ok<AGROW>
        end
        fprintf('Saved pilot outputs to %s\n',out_dir);
    end
    writetable(summary,fullfile(out_root,'batch_summary.csv'));
    fprintf('\n===== COMPLETED QUALITY PILOT =====\n');
    disp(summary);
end

function ensure_ustb_beamformed_class(repo_dir, custom_dir)
% Find the USTB root folder containing the +uff package before loading data.
    if ~isempty(which('uff.beamformed_data'))
        fprintf('USTB class found: %s\n', which('uff.beamformed_data'));
        return;
    end

    home_dir = getenv('USERPROFILE');
    candidates = {custom_dir, ...
        fullfile(home_dir,'USTB'), ...
        fullfile(home_dir,'Documents','MATLAB','USTB'), ...
        fullfile(home_dir,'Documents','GitHub','USTB'), ...
        fullfile(fileparts(repo_dir),'USTB'), ...
        fullfile(fileparts(fileparts(repo_dir)),'USTB'), ...
        fullfile(repo_dir,'USTB')};
    for k = 1:numel(candidates)
        root = candidates{k};
        if isempty(root) || ~isfolder(root), continue; end
        addpath(genpath(root));
        if ~isempty(which('uff.beamformed_data'))
            fprintf('USTB loaded from %s\n', root);
            return;
        end
    end

    error(['Cannot find class uff.beamformed_data. USTB must be on the ' ...
        'MATLAB path BEFORE loading the displacement MAT file. ' ...
        'Set ustb_custom_dir at the beginning of this script to the ' ...
        'USTB folder containing +uff, or run addpath(genpath(...)) ' ...
        'manually. Check with: which(''uff.beamformed_data'').']);
end

function [Q, C12dB, rho_map] = quality_map_from_displacement(u)
% First-pass zero-lag, mean-centered correlation between adjacent lateral
% displacement waveforms; no ground truth enters this function.
    [Nz,Nx,Nt] = size(u);
    assert(Nx >= 3, 'Need at least three lateral pixels.');
    sa=zeros(Nz,Nx-1); sb=sa; saa=sa; sbb=sa; sab=sa;
    for t=1:Nt
        a=double(u(:,1:end-1,t));
        b=double(u(:,2:end,t));
        sa=sa+a; sb=sb+b;
        saa=saa+a.^2; sbb=sbb+b.^2; sab=sab+a.*b;
    end
    numerator = sab - sa.*sb/Nt;
    va = max(0,saa - sa.^2/Nt);
    vb = max(0,sbb - sb.^2/Nt);
    denominator = sqrt(va.*vb);
    rho = numerator ./ denominator;
    rho(denominator <= eps(max(denominator(:)))) = NaN;
    rho = max(-1,min(1,rho));

    % Tiny-amplitude waveforms can have misleading correlations. Use a
    % deliberately conservative numerical RMS floor (0.01% of median).
    pair_rms=sqrt((saa+sbb)/(2*Nt));
    good_rms=pair_rms(isfinite(pair_rms) & pair_rms>0);
    if ~isempty(good_rms)
        rho(pair_rms < median(good_rms)*1e-4)=NaN;
    end

    % Published Q mapping; C12 here is based on *our* zero-lag correlation
    % implementation, not guaranteed to match the original signal processing.
    C_pair=20*log10(max(abs(rho),realmin));
    Q_pair=nan(Nz,Nx-1);
    Q_pair(C_pair >= -10)=1;
    mid=C_pair > -50 & C_pair < -10;
    Q_pair(mid)=cos((pi/2)*(C_pair(mid)+10)/40);
    Q_pair(C_pair <= -50)=0;

    % Assign pairwise quality to pixels conservatively (worst neighbor).
    Q=nan(Nz,Nx); C12dB=nan(Nz,Nx); rho_map=nan(Nz,Nx);
    Q(:,1)=Q_pair(:,1); Q(:,end)=Q_pair(:,end);
    Q(:,2:end-1)=min(Q_pair(:,1:end-1),Q_pair(:,2:end));
    C12dB(:,1)=C_pair(:,1); C12dB(:,end)=C_pair(:,end);
    C12dB(:,2:end-1)=min(C_pair(:,1:end-1),C_pair(:,2:end));
    rho_map(:,1)=rho(:,1); rho_map(:,end)=rho(:,end);
    rho_map(:,2:end-1)=0.5*(rho(:,1:end-1)+rho(:,2:end));
end

function [r, mae_hi, mae_lo] = evaluate_quality(Q, err, valid)
% Spearman correlation should be NEGATIVE when Q predicts low error.
    qv=double(Q(valid));
    ev=double(err(valid));
    isvalid=isfinite(qv) & isfinite(ev);
    qv=qv(isvalid); ev=ev(isvalid);
    if numel(qv)<10 || std(qv) < 1e-10 || std(ev) < 1e-10
        r=NaN;
    else
        r=corr(qv,ev,'Type','Spearman');
    end
    hi=qv>=0.63;
    if any(hi), mae_hi=mean(ev(hi)); else, mae_hi=NaN; end
    if any(~hi), mae_lo=mean(ev(~hi)); else, mae_lo=NaN; end
end
