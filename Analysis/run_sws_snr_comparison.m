function [measurements, summary] = run_sws_snr_comparison(mode, phantom_case)
%RUN_SWS_SNR_COMPARISON  Beamformed-data SNR test of two SWS estimators.
%
% Put this file in synthetic-the-rf-data-generation/Analysis/.
%
% QUICK PILOT (one 20 dB noise realization, plus clean reference):
%   run_sws_snr_comparison
%   run_sws_snr_comparison('pilot', 'large_inclusions')
%
% FULL STUDY (40, 30, 20, 10 dB, ten realizations per level):
%   run_sws_snr_comparison('full', 'large_inclusions')
%
% Notes:
% - Adds Gaussian noise to BEAMFORMED RF/IQ data, before displacement estimation.
% - Re-runs USTB autocorrelation displacement estimation for each realization.
% - Uses the SAME frequency extraction and directional filtering for both
%   phase-gradient (PG) and spatial-autocorrelation (AC) SWS estimators.
% - Uses the same stored ground truth and ROI masks for both methods.
% - No k-Wave or Field II re-simulation; original data/results are READ ONLY.
% - Direct-KDE gCNR is APPROXIMATE, based on deterministic subsampling
%   for feasible runtime; unlike RMSE/CNR it is not based on every ROI pixel.
% - Existing compare_SWS_methods_batch.m outputs are NOT overwritten.
%
% Requires: original dataset; Results/<case>/<case>_SWS_quality_results.mat;
% USTB/ElastographyToolbox; Statistics Toolbox (ksdensity).

if nargin < 1 || isempty(mode), mode = 'pilot'; end
if nargin < 2 || isempty(phantom_case), phantom_case = 'large_inclusions'; end
mode = lower(char(string(mode)));
phantom_case = char(string(phantom_case));
assert(any(strcmp(mode, {'pilot', 'full'})), ...
    'mode must be ''pilot'' or ''full''.');

if strcmp(mode, 'pilot')
    snr_dB_levels = 20;
    n_realizations = 1;
else
    snr_dB_levels = [40 30 20 10];
    n_realizations = 10;
end

% Deterministic KDE subsampling. Increase for greater gCNR precision.
kde_max_pixels_per_roi = 5000;
kde_eval_points = 500;

script_dir = fileparts(mfilename('fullpath'));
repo_dir = find_repository(script_dir);
master_dir = fileparts(repo_dir);
results_root = fullfile(repo_dir, 'Results');

% Local repo and the sibling USTB toolkit (same convention as working batch).
addpath(genpath(repo_dir));
toolbox_dir = fullfile(master_dir, 'elastography', 'ElastographyToolbox');
if isfolder(toolbox_dir), addpath(genpath(toolbox_dir)); end

assert(~isempty(which('uff.elastography')), ...
    'Cannot find USTB uff.elastography. Add USTB/ElastographyToolbox to MATLAB path.');
assert(~isempty(which('postprocess.autocorrelation_displacement_estimation')), ...
    'Cannot find postprocess.autocorrelation_displacement_estimation.');
assert(~isempty(which('postprocess.elastography.phase_gradient_method')), ...
    'Cannot find phase_gradient_method.');
assert(~isempty(which('postprocess.elastography.autocorrelation_SWS_method')), ...
    'Cannot find autocorrelation_SWS_method.');
assert(~isempty(which('ksdensity')), ...
    'Statistics and Machine Learning Toolbox (ksdensity) is required.');

baseline_file = fullfile(results_root, phantom_case, ...
    [phantom_case '_SWS_quality_results.mat']);
assert(isfile(baseline_file), 'Baseline SWS result missing: %s', baseline_file);
B = load(baseline_file, 'results', 'sws_gt');
assert(isfield(B, 'sws_gt') && isfield(B, 'results') && ...
    isfield(B.results, 'roi'), ...
    'Baseline must contain sws_gt and results.roi: %s', baseline_file);
GT = double(B.sws_gt);
roi = B.results.roi;
assert(all(isfield(roi, {'image_mask', 'inclusion_core', 'background_core'})), ...
    'Baseline results.roi is missing a required mask.');

source_file = locate_data(repo_dir, script_dir, phantom_case);
% Object-valued MAT variables require the USTB path, added above.
D = load(source_file, 'b_data', 'channel_data', 'displacement_estimation_autocorr', ...
    'shear_params', 'autocorr_params', 'prf', 'sws_gt');
assert(all(isfield(D, {'b_data','channel_data', ...
    'displacement_estimation_autocorr','shear_params'})), ...
    'Dataset must contain b_data, channel_data, displacement estimate and shear_params.');
if isfield(D, 'sws_gt')
    assert(isequaln(double(D.sws_gt),GT), ...
        'Ground truth in baseline results differs from source dataset.');
end

scan = D.b_data.scan;
x_axis = double(scan.x_axis(:));
z_axis = double(scan.z_axis(:));
Nz = numel(z_axis);
Nx = numel(x_axis);
assert(isequal(size(GT), [Nz Nx]), 'Ground truth does not match scan axes.');
assert(isequal(size(roi.image_mask),size(GT)) && ...
       isequal(size(roi.inclusion_core),size(GT)) && ...
       isequal(size(roi.background_core),size(GT)), ...
       'ROI masks do not match ground-truth dimensions.');

if isfield(D, 'prf') && ~isempty(D.prf)
    PRF = double(D.prf);
elseif isprop(D.channel_data, 'PRF')
    PRF = double(D.channel_data.PRF);
else
    PRF = 3000;
    warning('Legacy dataset: PRF absent. Assuming PRF = 3000 Hz.');
end
assert(isscalar(PRF) && isfinite(PRF) && PRF>0, 'Invalid PRF.');
frequency = double(D.shear_params.source_freq);
assert(isfinite(frequency) && isscalar(frequency) && frequency>0, ...
    'Invalid source frequency.');
filter_speed = 3.2;
if isfield(B.results,'filter_speed') && ~isempty(B.results.filter_speed)
    filter_speed = double(B.results.filter_speed);
end

if isfield(D, 'autocorr_params') && ~isempty(D.autocorr_params)
    params = D.autocorr_params;
else
    params = struct('z_gate',4, 'x_gate',2, 'packet_size',4);
    warning('autocorr_params absent; using legacy values [4,2,4].');
end
assert(all(isfield(params,{'z_gate','x_gate','packet_size'})), ...
    'autocorr_params lacks gate/packet_size fields.');

output_dir = fullfile(results_root,'snr_method_comparison',phantom_case,mode);
if ~isfolder(output_dir), mkdir(output_dir); end

fprintf('\n=== SWS INPUT-SNR COMPARISON ===\n');
fprintf('Phantom: %s | mode: %s | f: %.1f Hz | PRF: %.1f Hz\n', ...
    phantom_case,mode,frequency,PRF);
fprintf('Noise added to beamformed data, NOT to SWS maps.\n');
fprintf('Dataset: %s\nOutput: %s\n',source_file,output_dir);

clean_b = D.b_data.data;
assert(isnumeric(clean_b) && ~isempty(clean_b), ...
    'b_data.data must be a nonempty numeric matrix.');
signal_power = mean(abs(double(clean_b(:))).^2);
assert(isfinite(signal_power) && signal_power>0, ...
    'Clean beamformed data have invalid signal power.');

% Use saved CLEAN maps when possible to avoid recomputation. Fall back to
% reconstruction from the stored clean displacement (not a new simulation).
clean_file = fullfile(results_root,'method_comparison','batch',phantom_case, ...
    'method_comparison.mat');
if isfile(clean_file)
    C = load(clean_file,'PG','AC','GT','PRF','target_frequency','filter_speed');
    assert(all(isfield(C,{'PG','AC'})), 'Saved method maps are incomplete.');
    if isfield(C,'GT')
        assert(isequaln(double(C.GT),GT), ...
            'Clean method comparison uses different ground truth.');
    end
    if isfield(C,'PRF') && abs(double(C.PRF)-PRF)>1e-8
        error('Saved clean map uses different PRF; rerun method comparison.');
    end
    if isfield(C,'target_frequency') && ...
            abs(double(C.target_frequency)-frequency)>1e-8
        error('Saved clean map uses a different excitation frequency.');
    end
    if isfield(C,'filter_speed') && abs(double(C.filter_speed)-filter_speed)>1e-8
        error('Saved clean map uses a different spatial filter parameter.');
    end
    PG_clean = double(C.PG);
    AC_clean = double(C.AC);
    fprintf('Using saved clean PG/AC maps.\n');
else
    u_raw = double(D.displacement_estimation_autocorr.data);
    assert(mod(numel(u_raw),Nz*Nx)==0, 'Clean displacement dimensions invalid.');
    u = -reshape(u_raw,Nz,Nx,[]);
    [PG_clean,AC_clean] = reconstruct_two_methods(u,scan,PRF,frequency,filter_speed);
    clear u u_raw;
    fprintf('Reconstructed clean PG/AC maps from stored clean displacement.\n');
end
assert(isequal(size(PG_clean),size(GT)) && ...
       isequal(size(AC_clean),size(GT)), ...
       'Clean SWS map dimensions do not match ground truth.');

baseline_valid = valid_sws(PG_clean,AC_clean,GT);
assert(nnz(baseline_valid & logical(roi.inclusion_core))>1 && ...
       nnz(baseline_valid & logical(roi.background_core))>1, ...
       'Too few valid clean ROI samples.');

measurements = table();
methods = ["Phase gradient", "Spatial autocorrelation"];

% The clean rows provide reference points. Their SNR is Inf, not an input
% SNR measurement, and they are intentionally excluded from SNR plots.
measurements = add_two_rows(measurements,phantom_case,Inf,0, ...
    methods,PG_clean,AC_clean,GT,roi,baseline_valid, ...
    kde_max_pixels_per_roi,kde_eval_points);
write_outputs(measurements,output_dir);

for ii = 1:numel(snr_dB_levels)
    snr_dB = snr_dB_levels(ii);
    noise_power = signal_power / 10^(snr_dB/10);
    for rr = 1:n_realizations
        fprintf('SNR = %g dB | realization %d/%d\n', ...
            snr_dB,rr,n_realizations);
        rng(20000 + 100*ii + rr,'twister');
        if isreal(clean_b)
            noise = sqrt(noise_power).*randn(size(clean_b));
        else
            noise = sqrt(noise_power/2).* ...
                (randn(size(clean_b))+1i*randn(size(clean_b)));
        end
        D.b_data.data = clean_b + noise;
        clear noise;

        estimator = postprocess.autocorrelation_displacement_estimation();
        estimator.channel_data = D.channel_data;
        estimator.input = D.b_data;
        estimator.z_gate = params.z_gate;
        estimator.x_gate = params.x_gate;
        estimator.packet_size = params.packet_size;
        displacement_noisy = estimator.go();

        u_raw = double(displacement_noisy.data);
        assert(mod(numel(u_raw),Nz*Nx)==0, ...
            'Noisy displacement has unexpected dimensions.');
        u_noisy = -reshape(u_raw,Nz,Nx,[]);
        clear u_raw displacement_noisy estimator;
        [PG,AC] = reconstruct_two_methods(u_noisy,scan,PRF,frequency,filter_speed);
        clear u_noisy;
        assert(isequal(size(PG),size(GT)) && isequal(size(AC),size(GT)), ...
            'SWS map dimension mismatch.');

        % Same valid pixel intersection for BOTH methods within each trial.
        common_valid = baseline_valid & valid_sws(PG,AC,GT);
        assert(nnz(common_valid & logical(roi.inclusion_core))>1 && ...
               nnz(common_valid & logical(roi.background_core))>1, ...
            'Too few valid common ROI pixels at SNR %g dB, trial %d.',snr_dB,rr);
        measurements = add_two_rows(measurements,phantom_case,snr_dB,rr, ...
            methods,PG,AC,GT,roi,common_valid, ...
            kde_max_pixels_per_roi,kde_eval_points);

        % Save continuously so interrupted full sweeps retain completed trials.
        write_outputs(measurements,output_dir);
        fprintf('  PG RMSE %.4f | AC RMSE %.4f m/s\n', ...
            measurements.GlobalRMSE_mps(end-1), ...
            measurements.GlobalRMSE_mps(end));
    end
end
D.b_data.data = clean_b;
summary = make_summary(measurements);
writetable(summary,fullfile(output_dir,'summary.csv'));
save(fullfile(output_dir,'snr_results.mat'), ...
    'measurements','summary','phantom_case','mode','snr_dB_levels', ...
    'n_realizations','frequency','PRF','filter_speed', ...
    'kde_max_pixels_per_roi','kde_eval_points');
plot_summary(summary,measurements,output_dir,phantom_case);

fprintf('\nSNR COMPARISON COMPLETE.\nCSV: %s\n', ...
    fullfile(output_dir,'measurements.csv'));
disp(summary);
end

function repo_dir = find_repository(folder)
repo_dir = folder;
for i=1:6
    if isfolder(fullfile(repo_dir,'Results')) && ...
            isfolder(fullfile(repo_dir,'Analysis'))
        return;
    end
    parent = fileparts(repo_dir);
    if strcmp(parent,repo_dir), break; end
    repo_dir = parent;
end
error(['Cannot locate repository. Place this file inside the ' ...
    'synthetic-the-rf-data-generation project.']);
end

function data_file = locate_data(repo_dir,script_dir,case_name)
name = [case_name '_displacement_with_quality.mat'];
options = { ...
    fullfile(repo_dir,'Data',case_name,name), ...
    fullfile(repo_dir,'Data',name), ...
    fullfile(repo_dir,'Analysis',name), ...
    fullfile(repo_dir,name), ...
    fullfile(script_dir,name)};
data_file = '';
for k=1:numel(options)
    if isfile(options{k}), data_file=options{k}; return; end
end
error('Dataset %s not found in expected repo locations.',name);
end

function [PG,AC] = reconstruct_two_methods(u,scan,PRF,f,filter_speed)
elasto=uff.elastography();
elasto.PRF=PRF;
elasto.target_frequency=f;
elasto.scan=scan;
elasto.autocorrelation_field=uff.data_field();
elasto.autocorrelation_field.set_component('u', ...
    reshape(u,size(u,1),size(u,2),1,size(u,3)));
num_angles=12;
elasto.directional_filtering_angles = ...
    2*pi*linspace(0,1-1/num_angles,num_angles);
extractor=postprocess.elastography.extract_freq_component(elasto);
extractor.target_frequency=f;
extractor.component='u';
extractor.run();
omega=2*pi*f;
elasto.radial_filter=uff.spatial_filter('gaussian_charite','radial');
elasto.radial_filter.sigma=filter_speed/omega;
elasto.angular_filter=uff.spatial_filter('gaussian_charite','angular');
elasto.angular_filter.deltaTheta=median(diff(elasto.directional_filtering_angles));
filter_operator=postprocess.elastography.apply_filter(elasto);
filter_operator.run();
pg=postprocess.elastography.phase_gradient_method(elasto);
pg.run();
PG=squeeze(double(elasto.sws_map.u));
ac=postprocess.elastography.autocorrelation_SWS_method(elasto);
ac.packetSize=4; % Explicit documented default used in the existing comparison.
ac.run();
AC=squeeze(double(elasto.sws_map.u));
end

function valid=valid_sws(PG,AC,GT)
valid=isfinite(GT) & isfinite(PG) & PG>0 & PG<10 & ...
    isfinite(AC) & AC>0 & AC<10;
end

function T=add_two_rows(T,case_name,snr,trial,methods,PG,AC,GT,roi,valid,kde_max,ke)
maps={PG,AC};
mask_global=valid & logical(roi.image_mask);
mask_incl=valid & logical(roi.inclusion_core);
mask_bkg=valid & logical(roi.background_core);
assert(nnz(mask_global)>1 && nnz(mask_incl)>1 && nnz(mask_bkg)>1, ...
    'Too few common valid pixels for metrics.');
for k=1:2
    m=compute_metrics(maps{k},GT,mask_global,mask_incl,mask_bkg,kde_max,ke);
    new=table(string(case_name),double(snr),double(trial),methods(k), ...
        m.global_rmse,m.inclusion_rmse,m.background_rmse, ...
        m.global_bias,m.inclusion_bias,m.background_bias, ...
        m.CNR,m.gCNR,m.contrast_recovery, ...
        double(nnz(mask_global)),double(nnz(mask_incl)),double(nnz(mask_bkg)), ...
        'VariableNames',{'Phantom','SNR_dB','Realization','Method', ...
        'GlobalRMSE_mps','InclusionRMSE_mps','BackgroundRMSE_mps', ...
        'GlobalBias_mps','InclusionBias_mps','BackgroundBias_mps', ...
        'CNR','gCNR_KDE_subsampled','ContrastRecovery', ...
        'GlobalPixels','InclusionPixels','BackgroundPixels'});
    if isempty(T), T=new; else, T=[T;new]; end %#ok<AGROW>
end
end

function m=compute_metrics(map,GT,gmask,imask,bmask,kde_max,ke)
err_g=map(gmask)-GT(gmask);
err_i=map(imask)-GT(imask);
err_b=map(bmask)-GT(bmask);
incl=double(map(imask));
bkg=double(map(bmask));
m.global_rmse=sqrt(mean(err_g.^2));
m.inclusion_rmse=sqrt(mean(err_i.^2));
m.background_rmse=sqrt(mean(err_b.^2));
m.global_bias=mean(err_g);
m.inclusion_bias=mean(err_i);
m.background_bias=mean(err_b);
denom=hypot(std(incl),std(bkg));
m.CNR=NaN;
if isfinite(denom) && denom>0
    m.CNR=abs(mean(incl)-mean(bkg))/denom;
end
true_contrast=mean(GT(imask))-mean(GT(bmask));
m.contrast_recovery=NaN;
if isfinite(true_contrast) && abs(true_contrast)>eps
    m.contrast_recovery=(mean(incl)-mean(bkg))/true_contrast;
end

% Approximate DIRECT KDE gCNR on systematic pixel subsamples.
% Full-resolution RMSE/bias/CNR/contrast recovery do NOT use this subset.
incl=incl(round(linspace(1,numel(incl),min(kde_max,numel(incl)))));
bkg=bkg(round(linspace(1,numel(bkg),min(kde_max,numel(bkg)))));
xmin=min([incl;bkg]);
xmax=max([incl;bkg]);
m.gCNR=NaN;
if isfinite(xmin) && isfinite(xmax)
    if xmax==xmin
        m.gCNR=0;
    elseif xmax>xmin
        x=linspace(xmin,xmax,ke);
        fi=ksdensity(incl,x);
        fb=ksdensity(bkg,x);
        ni=trapz(x,fi); nb=trapz(x,fb);
        if ni>0 && nb>0
            m.gCNR=1-trapz(x,min(fi/ni,fb/nb));
        end
    end
end
end

function write_outputs(T,out_dir)
writetable(T,fullfile(out_dir,'measurements.csv'));
end

function summary=make_summary(T)
summary=groupsummary(T,{'SNR_dB','Method'}, {'mean','std'}, ...
    {'GlobalRMSE_mps','InclusionRMSE_mps','BackgroundRMSE_mps', ...
    'GlobalBias_mps','InclusionBias_mps','BackgroundBias_mps', ...
    'CNR','gCNR_KDE_subsampled','ContrastRecovery'});
end

function plot_summary(S,T,out_dir,case_name)
levels=sort(unique(T.SNR_dB(isfinite(T.SNR_dB))));
methods=["Phase gradient","Spatial autocorrelation"];
fig=figure('Name',['SWS SNR study | ' case_name],'Color','w');
l=tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
fields={'GlobalRMSE_mps','CNR'};
labels={'Global SWS RMSE [m/s]','CNR'};
for p=1:2
    ax=nexttile(l); hold(ax,'on');
    for j=1:2
        y=nan(size(levels)); ysd=y;
        for k=1:numel(levels)
            mask=S.SNR_dB==levels(k) & S.Method==methods(j);
            if any(mask)
                mean_values = S.(['mean_' fields{p}]);
                std_values = S.(['std_' fields{p}]);
                y(k)=mean_values(mask);
                ysd(k)=std_values(mask);
            end
        end
        errorbar(ax,levels,y,ysd,'-o','LineWidth',1.5, ...
            'DisplayName',char(methods(j)));
    end
    xlabel(ax,'Beamformed-data SNR [dB]');
    ylabel(ax,labels{p});
    grid(ax,'on');
    legend(ax,'Location','best');
end
title(l,['SWS reconstruction sensitivity | ' strrep(case_name,'_',' ')]);
exportgraphics(fig,fullfile(out_dir,'SNR_method_comparison.png'),'Resolution',200);
savefig(fig,fullfile(out_dir,'SNR_method_comparison.fig'));
end
