function [nmrFit, nmrFit_ppm, fids, tr] = calculateStaticSpectroscopy(raw_path, BHs, hr, model)

[fids, dwell_time, npts, tr, xeFreqMHz, rf_excitation] = readRawDyn(raw_path);
nFrames = size(fids,2); % Consortium standard: 500 dis., 20 gas (520 total)
t_tr = tr*(1:nFrames);
nAvg = ceil(60/hr/tr);
[BHstart, ~] = findBHs(t_tr, BHs);

%% Calculate Static Spectral Parameters
t = dwell_time*(0:(npts-1))';
% Analyze only a fraction of the total data
avgdata = mean(fids(:,BHstart:BHstart+nAvg),2); % average 1s of data
% perform findpeaks to guide initial area and frequency guess. 
% finds 3 peaks (barrier, RBC, gas) in descending order
avgdata_fft = abs(fftshift(fft(avgdata)));
% Finds frequency of peak
freq_ppm = linspace(-1/(2*dwell_time), 1/(2*dwell_time), length(avgdata))/xeFreqMHz; % in ppm
freq_ppm_dis = freq_ppm(-50 < freq_ppm & freq_ppm < 50);
freq_ppm_gas = freq_ppm(-175 > freq_ppm);
% find peaks in dissolved region
[pks_dis, locs_dis, ~] = findpeaks(avgdata_fft(-50 < freq_ppm & freq_ppm < 50), ...
                'SortStr','descend','NPeaks', 2);
% find the gas peak
[pk_gas, loc_gas, ~] = findpeaks(avgdata_fft(freq_ppm < -175), ...
                'SortStr','descend','NPeaks', 1);
% Flag of Using Junlan's guessing
Using_Junlan_guessing = false;

% checks if there are 2 dissolved peaks that are approximately 
% spaced apart. Otherwise, use hardcoded guess
if Using_Junlan_guessing && (length(pks_dis) == 2 && (abs(freq_ppm_dis(locs_dis(2))...
        - freq_ppm_dis(locs_dis(1))) < 25 && abs(freq_ppm_dis(locs_dis(2))...
        - freq_ppm_dis(locs_dis(1))) > 15) )
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
freq_orig = [rbc_freq membrane_freq gas_freq]*xeFreqMHz; % in Hz
% freq_orig = [0 -12 -218]*xeFreqMHz+7500; For 90 degree spoiled
fwhmL_orig = [8.8 5.0 1.2]*xeFreqMHz;
fwhmG_orig = [0 6.1 0]*xeFreqMHz;
% Set the starting guess of area.
area_orig = [pk_rbc, pk_membrane, pk_gas];
phase_orig = [0 0 0];
peaks = 3; ngas = 3; nbar = 2;

switch model
    case 'V' 
        nmrFit = NMR_TimeFit_v(avgdata, t, area_orig,freq_orig,fwhmL_orig,fwhmG_orig,phase_orig,[],[]);
    case 'L'
        nmrFit = NMR_TimeFit(avgdata, t, area_orig,freq_orig,fwhm_orig,phase_orig,[],[]);
end 
nmrFit = nmrFit.fitTimeDomainSignal();

%% Calculate Peak SNRs
% historical "tail of fit residuals" SNR method
timeFit=nmrFit.calcTimeDomainSignal(nmrFit.t);  % calculate fitted data in time domain
timeRes=timeFit - nmrFit.timeDomainSignal;
n25pct = round(length(timeRes)/4);
std25 = std(timeRes(end-n25pct:end)); % take standard deviation of tail of residuals
SNRresid = nmrFit.area/std25; % will contain a value for each peak (RBC, bar, gas)

% Doreen's "simulated noise frame (snf)" SNR method:
BH_fids = fids(:,BHstart:BHstart+nAvg);
fid_tail = length(BH_fids)/2; % focus calculation 2nd half of fids
cropped_fids = BH_fids(fid_tail+1:end,:); 
mean_fid = mean(cropped_fids,2);
diff = mean_fid-cropped_fids; % array of difference fids
diff = diff - mean (diff); % subtract off any DC bias
noiseDis = std(real(diff));
SNR_frames = nmrFit.area'./noiseDis; % SNR of each frame
SNRsnf = mean(SNR_frames,2)*sqrt(nAvg); % simulated noise frame SNR

%% Change units from Hz to PPM
ref_freq = nmrFit.freq(ngas);

positive_phase = nmrFit.phase - nmrFit.phase(2);
positive_phase(positive_phase < 0) = positive_phase(positive_phase < 0) + 360;

nmrFit_ppm.area = nmrFit.area/sum(nmrFit.area(nbar));
nmrFit_ppm.freq = (nmrFit.freq-ref_freq)/xeFreqMHz;
nmrFit_ppm.fwhm = nmrFit.fwhm/xeFreqMHz;
switch model
case 'V'
    nmrFit_ppm.fwhmG = nmrFit.fwhmG/xeFreqMHz;
end 
nmrFit_ppm.phase = positive_phase;
nmrFit_ppm.SNRresid = SNRresid;
nmrFit_ppm.SNRsnf = SNRsnf;

%% Display Results
switch model
case 'V'
    disp([10 'Static Parameters:'])
    disp('    Area     Freq      FWHM      FWHMg     Phase     SNRresid   SNRsnf');
    for k = 1:peaks
        disp([sprintf('%8.2f', nmrFit.area(k)/sum(nmrFit.area(nbar))), ' ' ...
            sprintf('%8.1f',(nmrFit.freq(k)-ref_freq)/xeFreqMHz),  '  ' ...
            sprintf('%8.1f',nmrFit.fwhm(k)/xeFreqMHz), '  ' ...
            sprintf('%8.1f',nmrFit.fwhmG(k)/xeFreqMHz), '  ' ...
            sprintf('%9.1f',positive_phase(k)), '  '...
            sprintf('%9.1f',SNRresid(k)) '  '...
            sprintf('%9.1f',SNRsnf(k))]);
    end
case 'L' 
    disp([10 'Static Parameters:'])
    disp('    Area     Freq      FWHM       Phase     SNRresid   SNRsnf');
    nmrFit.phase(nmrFit.phase < 0) = nmrFit.phase(nmrFit.phase < 0)+360;
    for k = 1:peaks
        disp([sprintf('%8.2f', nmrFit.area(k)/sum(nmrFit.area(nbar))), ' ' ...
            sprintf('%8.1f',(nmrFit.freq(k)-ref_freq)/xeFreqMHz),  '  ' ...
            sprintf('%8.1f',nmrFit.fwhm(k)/xeFreqMHz), '  ' ...
            sprintf('%9.1f',positive_phase(k)), '  '...
            sprintf('%9.1f',SNRresid(k)) '  '...
            sprintf('%9.1f',SNRsnf(k))]);
    end
end 

end 
