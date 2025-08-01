function [amp_all, nmrFit, nmrFit_ppm, dyn, detrend_peaks, fitted_peaks] = ...
    dynamicSummary(raw_path, dyn_path, BHs, save_fig_flag, rbc_axis_lim)
% dynamicSummary: produce detrended oscillations and oscillation fits. 
%       detects and determines if dyn is Voigt or Lorentzian
% Inputs:
%   raw_path: raw data file either .dat or .7
%   dyn_path: dynamic spectroscopy matlab structure or structure location.
%   BHs: either a vector containing the time of inhale/exhale or the filepath
%       for a saved matlab vector with the BH information. Default
%       values are [2 10].
%   save_fig_flag: logical input for saving figures. Figures are saved in the
%       same folder as dyn_path.
%   rbc_axis_lim: Sets axis limits on rbc oscillation plots. If 0,
%       automatic limits are used.
% Outputs:
%   amp_all: structure with fit amplitudes, heartrate, and dynamic SNR
%   nmrFit: spectroscopy fit object (Hz)
%   nmrFit_ppm: spectroscopy fit object (ppm)
%   dyn: dynamic spectroscopy matlab structure
%   detrend_peaks: structure containing detrended oscillations
%   fitted_peaks: structure containing fits to detrended oscillations

% Example inputs
% raw_path = '~/Desktop/subject/subject_id/meas_MID00026_FID69955_Xe_fid_DynamicSpec_high_flip.dat';
% 
% dyn_name = 'dynV5_2';
% 
% dyn_path = [path(1:end-7),'4 - Spectra/',dyn_name,'.mat'];
% 
% BHs = [3 10];
% save_fig_flag = 0;
% 
% figname = '';
% fig_path = fileparts(dyn_path);

%%
if isempty(raw_path)
    disp('Please select the raw dat file')
    raw_path = filepath();
    dyn_folder = fileparts(raw_path);
else 
    dyn_folder = fileparts(raw_path);
end 

if isempty(dyn_path)
    disp('Please select the dyn variable')
    dyn_path = filepath();
    load(dyn_path)
elseif ischar(dyn_path)
    load(dyn_path)
    [dyn_folder] = fileparts(dyn_path);
else
    dyn = dyn_path;
end

[BHstart, BHend] = findBHs(dyn.t(:,1), BHs);

if ~isfield(dyn,'fwhmG') % if the data is not Voigt, set all FWHM_G to 0
    dyn.fwhmG = zeros(size(dyn.t));
    dyn.fwhmL = dyn.fwhm;
    fitType = 'L';
else
    fitType = 'V';
end

%% Calculate Oscillations Amplitudes and Static Parameters

% amplitude for the peaks fit
b = highpassfilter(length(dyn.area(:,1)));
[amp_peaks, detrend_peaks, fitted_peaks] = calculateOscillationAmps(dyn,BHs,b, 'peaks', 'oscType', 'peaks_lu');
close all
if ~isfield(dyn,'snrs')
    dyn.snrs = zeros(1,length(dyn.t));
end 
amp_peaks.snr = mean(dyn.snrs(BHstart:BHend));

amp_all.peaks = amp_peaks;

% % Display oscillation results
% disp([10 'file: ',raw_path])
% disp([10 'RBC Dynamic Oscillation Amplitudes:'])
% disp('   Area      Freq      FWHM      Phase      HR      SNR');
% disp([sprintf('%8.1f', amp.area*100), ' ' ...
%     sprintf('%8.2f',amp.freq),  '  ' ...
%     sprintf('%8.2f',amp.fwhm), '  ' ...
%     sprintf('%8.1f',amp.phase), '   '...
%     sprintf('%6.0f',amp.hr),'  ',...
%     sprintf('%7.1f',amp.snr)]);

%% Calculate Oscillations Amplitudes and Static Parameters
[nmrFit, ~, nmrFit_ppm, ~] = calculateStaticSpectroscopy(raw_path, BHs, fitType);

%% Plots

%%%%%%%%%%%%%%%%%%%%%% Peaks Fit figures %%%%%%%%%%%%%%%%%%%%%%%%
% Plotting dynamics, oscillations, dynamic SNR, static spectroscopy 
% for the peaks fit, peaks is only for oscillations but keeping all the
% figures for peaks fit to prepare the slides seamlessly

figName = 'peaks';
dynSave_path = [dyn_folder,'/dyn',fitType,figName,'.tif'];
plotDynamics(dyn, BHs, save_fig_flag, dynSave_path);

oscSave_path = [dyn_folder,'/dyn',fitType,'_oscs_new',figName,'.tif'];
plotOscillations(dyn, BHs, detrend_peaks, fitted_peaks, rbc_axis_lim, save_fig_flag, oscSave_path)

snrSave_path = [dyn_folder,'/dyn',fitType,'_error',figName,'.tif'];
plotDynSNR(dyn, BHs, save_fig_flag, snrSave_path)

staticSave_path = [dyn_folder,'/static',fitType,figName,'.tif'];
plotStaticSpectroscopy(raw_path, nmrFit, save_fig_flag, staticSave_path)
close all;
