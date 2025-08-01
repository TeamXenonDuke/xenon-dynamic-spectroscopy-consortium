function dyn = fitDynamicSpec(raw_path, dyn_save_name)
% Fit dynamic spectroscopy.
% raw_path: full filepath of the raw scanner data
% save_name: name of .mat file for saving. The output will be saved in
%            the same folder ans the raw data

peaks = 3;
nToAvg = 5; % 1 for high flip
skipSize = 1;
gasFrames = 20;

[raw_folder, ~, ~] = fileparts(raw_path);
% Read in twix or P file and define associated variables
[raw_fids, dwell_time, npts, tr, xeFreqMHz, rf_excitation] = readRawDyn(raw_path);
% For high flip, also comment out SIFT
% raw_fids = raw_fids(10:end,5:end);

% SIFT raw fids
disFrames = size(raw_fids,2) - gasFrames;
raw_fids = raw_fids(:, 2:disFrames);
raw_fids = SIFT(raw_fids, dwell_time, tr);
% raw_fids(:,901:end-1) = []; % shorten long aquisitions

% Separate fids from gas frames
gas_fid = mean(raw_fids(:, end-1:end), 2);

% Separate dissolved frames
fids = raw_fids(:, 1:(end -2));

npts = size(raw_fids, 1); % Number of samples
t = dwell_time * (0:(npts - 1))'; % FID time vector
nFrames = size(fids, 2); % Number of disolved frames
t_tr = tr * (1:nFrames); % tr time vector
% average data for high SNR fit
% avgdata = mean(fids(:,150:250),2);
avgdata = mean(fids(:, 100:200), 2);
switch peaks
    case 3
        %         Three peak fit
        % perform findpeaks to guide initial area and frequency guess.
        % finds 3 peaks (barrier, RBC, gas) in descending order
        avgdata_fft = abs(fftshift(fft(avgdata)));
        % Finds frequency of peak
        freq_ppm = linspace(-1/(2 * dwell_time), 1/(2 * dwell_time), length(avgdata)) / xeFreqMHz; % in ppm
        freq_ppm_dis = freq_ppm(-50 < freq_ppm & freq_ppm < 50);
        freq_ppm_gas = freq_ppm(-175 > freq_ppm);
        % find peaks in dissolved region
        [pks_dis, locs_dis, ~] = findpeaks(avgdata_fft(-50 < freq_ppm & freq_ppm < 50), ...
            'SortStr', 'descend', 'NPeaks', 2);
        % find the gas peak
        [pk_gas, loc_gas, ~] = findpeaks(avgdata_fft(freq_ppm < -175), ...
            'SortStr', 'descend', 'NPeaks', 1);

        % Flag of Using Junlan's guessing
        Using_Junlan_guessing = false;

        % checks if there are 2 dissolved peaks that are approximately
        % spaced apart. Otherwise, use hardcoded guess
        if Using_Junlan_guessing && (length(pks_dis) == 2 && (abs(freq_ppm_dis(locs_dis(2)) ...
                -freq_ppm_dis(locs_dis(1))) < 25 && abs(freq_ppm_dis(locs_dis(2)) ...
                -freq_ppm_dis(locs_dis(1))) > 15))
            pk_membrane = pks_dis(1);
            % Set the barrier freq
            membrane_freq = freq_ppm_dis(locs_dis(1)); % in ppm
            % find the smallest peak near [-50, 50] ppm range (rbc peak)
            pk_rbc = pks_dis(2);
            % Set the rbc freq
            rbc_freq = freq_ppm_dis(locs_dis(2)); % in ppm
            % Set the rbc freq
            gas_freq = freq_ppm_gas(loc_gas); % in ppm
        else
            % Calculate the initial guess based on the read-in ppm 
            % Using Bas equation
            rbc_freq = 217.2-rf_excitation;
            membrane_freq = 197.7-rf_excitation;
            gas_freq = 0-rf_excitation;
  
            pk_rbc = 1; pk_membrane = 1; pk_gas=1;
         
        end
        % Three peak fit initial guesses
        % Set starting guess of frequency.
        freq_orig = [rbc_freq, membrane_freq, gas_freq] * xeFreqMHz; % in Hz
        % freq_orig = [0 -12 -218]*xeFreqMHz+7500; For 90 degree spoiled
        fwhmL_orig = [8.8, 5.0, 1.2] * xeFreqMHz;
        fwhmG_orig = [0, 6.1, 0] * xeFreqMHz;
        % Set the starting guess of area.
        area_orig = [pk_rbc, pk_membrane, pk_gas];
        phase_orig = [0, 0, 0];
        titles = {'RBC', 'Barrier', 'Gas'};
    case 4
        % Four peak fit
        area_orig = [1, 1, 1, 1];
        freq_orig = [-20, -285, -393, -3600] * 2;
        fwhm_orig = [215, 200, 150, 70];
        phase_orig = [0, 0, 0, 0];
        titles = {'RBC', 'Barrier 1', 'Barrier 2', 'Gas'};
    case 5
        % Five peak fit
        area_orig = [1, 1, 1, 1, 1];
        freq_orig = [-20, -285, -393, -3465, -3842] * 2;
        fwhm_orig = [215, 200, 150, 70, 30];
        phase_orig = [0, 0, 0, 0, 0];
        titles = {'RBC', 'Barrier 1', 'Barrier 2', 'Gas 1', 'Gas 2'};
end

%% High SNR and fit for starting guesses
nmrFit = NMR_TimeFit_v(avgdata, t, area_orig, freq_orig, fwhmL_orig, fwhmG_orig, phase_orig, [], []);
nmrFit = nmrFit.fitTimeDomainSignal();

% Create guesses from fitted data
area = nmrFit.area;
freq = nmrFit.freq;
fwhmL = nmrFit.fwhm;
fwhmG = nmrFit.fwhmG;
phase = nmrFit.phase;

% Define reference frequency for ppm conversion
ref_freq = nmrFit.freq(end);

% Fit spectra
nComp = length(area);
startingTimePoints = 1:skipSize:(nFrames - nToAvg);
nTimePoints = length(startingTimePoints);

area_dyn = zeros(nTimePoints, nComp);
freq_dyn = zeros(nTimePoints, nComp);
fwhmL_dyn = zeros(nTimePoints, nComp);
fwhmG_dyn = zeros(nTimePoints, nComp);
phase_dyn = zeros(nTimePoints, nComp);
snrs = zeros(nTimePoints, 1);
residualSpectrum = zeros(npts, nTimePoints);
tzero = min(t(:));

t_dyn = repmat(t_tr(startingTimePoints)', [1, nComp]);

disp('    Area     Freq      FWHM      FWHMg     Phase');
for k = 1:peaks
    disp([sprintf('%8.3f', nmrFit.area(k)/sum(nmrFit.area(2))), ' ', ...
        sprintf('%8.1f', (nmrFit.freq(k) - ref_freq)/xeFreqMHz), '  ', ...
        sprintf('%8.1f', nmrFit.fwhm(k)/xeFreqMHz), '  ', ...
        sprintf('%8.1f', nmrFit.fwhmG(k)/xeFreqMHz), '  ', ...
        sprintf('%9.1f', nmrFit.phase(k))]);

end

avgdata = movmean(fids, nToAvg, 2, 'Endpoints', 'discard');

tic
parfor iTimePoint = 1:nTimePoints
    nmrFit = NMR_TimeFit_v(avgdata(:, iTimePoint), t, area, freq, fwhmL, fwhmG, phase, [], []);
    nmrFit = nmrFit.fitTimeDomainSignal();

    area_dyn(iTimePoint, :) = nmrFit.area(:);
    freq_dyn(iTimePoint, :) = (nmrFit.freq(:) - ref_freq) / xeFreqMHz; %ppm units
    fwhmL_dyn(iTimePoint, :) = nmrFit.fwhm(:) / xeFreqMHz; %ppm units
    fwhmG_dyn(iTimePoint, :) = nmrFit.fwhmG(:) / xeFreqMHz; %ppm units
    phase_dyn(iTimePoint, :) = nmrFit.phase(:);

    fittedSignal = nmrFit.calcTimeDomainSignal(t);
    fittedSpectrum = dwell_time * fftshift(fft(fittedSignal));
    residualSpectrum(:, iTimePoint) = (nmrFit.spectralDomainSignal - fittedSpectrum);
    snrs(iTimePoint) = snr(fittedSignal, avgdata(:, iTimePoint)-fittedSignal);
end
toc

% Create Structure for Easy Saving of Data
fnames = {'area', 'freq', 'fwhmL', 'fwhmG', 'phase', 't', 'resSpec', 'ref_freq', 'snrs'};

dyn.(fnames{1}) = area_dyn;
dyn.(fnames{2}) = freq_dyn;
dyn.(fnames{3}) = fwhmL_dyn;
dyn.(fnames{4}) = fwhmG_dyn;
dyn.(fnames{5}) = phase_dyn;
dyn.(fnames{6}) = t_dyn;
dyn.(fnames{7}) = residualSpectrum;
dyn.(fnames{8}) = ref_freq;
dyn.(fnames{9}) = snrs;

dyn.processed_date = datestr(now);
dyn.raw_file = raw_path;
% save .mat file of the dyn struct
save(fullfile(raw_folder, 'Spectroscopy', dyn_save_name), 'dyn');
