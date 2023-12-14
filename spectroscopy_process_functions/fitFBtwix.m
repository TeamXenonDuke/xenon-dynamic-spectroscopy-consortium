function dyn = fitFBtwix(twix_path,peaks,nToAvg)

% raw_path: full filepath of the raw scanner data
% save_name: name of .mat file for saving. The output will be saved in
%            the same folder ans the raw data

skipSize = 1;
linebroadening = 0;
calFrames = 21;

[raw_folder, ~, scanner] = fileparts(twix_path);
% Read in twix or P file and define associated variables
[raw_fids, dwell_time, npts, tr, xeFreqMHz] = readRawDyn(twix_path);

% SIFT raw fids
raw_fids = raw_fids(:,2:end-calFrames);
raw_fids = SIFT(raw_fids, dwell_time, tr);
% raw_fids(:,901:end-1) = []; % shorten long aquisitions

switch peaks
    case 3 
        % Three peak fit
        area_orig = [1 1 1];
        freq_orig = [30 -650 -7430];
        fwhm_orig = [215 200 40];
        phase_orig = [0 0 0];
        titles = {'RBC' 'Barrier' 'Gas 1' 'Gas 2'};
    case 4
        % Four peak fit
        area_orig = [1 1 1 1];
        freq_orig = [-20 -285 -3600 -3650]*2;
        fwhm_orig = [215 200 70 70];
        phase_orig = [0 0 0 0];
        titles = {'RBC' 'Barrier 1' 'Barrier 2' 'Gas 1' 'Gas 2'};
    case 5
        % Five peak fit
        area_orig = [1 1 1 1 1];
        freq_orig = [-20 -285 -393 -3465 -3842]*2;
        fwhm_orig = [215 200 150 70 30];
        phase_orig = [0 0 0 0 0 ];
end 

fids = raw_fids;

% Separate fids from gas frames
gas_fids = fids(:,end-1:end);
gas_fid = mean(gas_fids,2);


% Separate dissolved frames
fids = fids(:,1:(end-21));

% Create array of sample times (sec)
npts = size(fids,1);                   % Number of samples                  % Receiver bandwidth (kHz)
t = dwell_time*(0:(npts-1))';
nFrames = size(fids,2);
t_tr = tr*((1:nFrames)-1);

% Fit gas sample

%% Fit gas-phase data to 2 peaks (airway, alveolar) for a proper reference frequency

% Constants
ox_shift = (4.5/5.5)*(0.13)*(273/298)*0.917;
xe_shift = (0.3/5.5)*(273/298)*(0.548);
susc_shift = -9.06/3;
tot_shift = -xe_shift;
tot_shift_hz = tot_shift*xeFreqMHz;

gas_fit_guess = [
    1           0           36.0          0; % Component #1
    1           -10           49.7          0; % Component #2
    ];

% Fit dedicated gas
gasFit= NMR_TimeFit(gas_fid,t,...
    gas_fit_guess(:,1),gas_fit_guess(:,2),...
    gas_fit_guess(:,3),gas_fit_guess(:,4),...
    linebroadening,0);
gasFit= gasFit.fitTimeDomainSignal();

% Make two peaks reasonable in frequency
badGasFreq_low = gasFit.freq < -200;
badGasFreq_high =  gasFit.freq > 100;
freqChange = 70;
while(any(badGasFreq_low) || any(badGasFreq_high))
    if(any(badGasFreq_low))
        if(sum(badGasFreq_low) ==2)
            gasFit= gasFit.fitTimeDomainSignal();
            disp('Coulndt find freq! All too low...')
            break            
%             error('Coulndt find freq! All too low...');
        else
            if(badGasFreq_low(2))
                gas_fit_guess(1,2) = gasFit.freq(1);
                gas_fit_guess(2,2) = gasFit.freq(1)-freqChange;
            else
                gasFit= gasFit.fitTimeDomainSignal();
                disp('Impossible! to have higher freq too low')
                break
%                 error('Impossible! to have higher freq too low')
            end
        end
    else
        if(sum(badGasFreq_high)==2)
            gasFit= gasFit.fitTimeDomainSignal();
            disp('Couldnt find freqs! all too high...')
            break
%             error('Couldnt find freqs! all too high...')
        else
            if(badGasFreq_high(2))
                gasFit= gasFit.fitTimeDomainSignal();
                disp('Impossible! to have lower freq too high')
                break
%                 error('Impossible! to have lower freq too high')
            else
                gas_fit_guess(1,2) = gasFit.freq(2)+freqChange;
                gas_fit_guess(2,2) = gasFit.freq(2);
            end
        end
    end
    
    % Refit
    gasFit= NMR_TimeFit(gas_fid,t,...
        gas_fit_guess(:,1),gas_fit_guess(:,2),...
        gas_fit_guess(:,3),gas_fit_guess(:,4),...
        linebroadening,[]);
    gasFit= gasFit.fitTimeDomainSignal();
    
    freqChange = freqChange-5;
    if(freqChange <= 0)
        gasFit= gasFit.fitTimeDomainSignal();
        disp('Crap, cant find 2 gas-phase peaks!')
        break
%         error('Crap, cant find 2 gas-phase peaks!');
    end
    
    % Look for bad gas freqs
    badGasFreq_low = gasFit.freq < -200;
    badGasFreq_high =  gasFit.freq > 100;
end

% Put frequencies into ppm units
ref_freq_gas = gasFit.freq(1) + tot_shift_hz;
ref_freq = gasFit.freq(1) - 7430 + tot_shift_hz; %7420 freq offset for later subjs %7380 = freq offset for 94B
gasFit_ppm.freq = (gasFit.freq - ref_freq_gas)/xeFreqMHz;
gasFit_ppm.fwhm = gasFit.fwhm/xeFreqMHz;

% Report gas-phase fits
% fprintf ('\n*** 3-Peak fit of gas-phase signal ***\n');
% displayPPMFits(gasFit);

%%
% gasFit = NMR_TimeFit(gas_pfile.data, t,1,-3832,20,0, [],[]);
% gasFit = gasFit.fitTimeDomainSignal();

% fids = fids(:,1:end-22); % Added this line to run dynamics on cal
avgdata = mean(fids(:,150:250),2);
nmrFit = NMR_TimeFit(avgdata, t, area_orig,freq_orig,fwhm_orig,phase_orig,[],[]);
nmrFit = nmrFit.fitTimeDomainSignal();

% Create guesses from fitted data
area = nmrFit.area;
freq = nmrFit.freq;
fwhm = nmrFit.fwhm;
phase = nmrFit.phase; %phase(2) = -160;

ref_freq = nmrFit.freq(end);

% Fit spectra
nComp = length(area);
startingTimePoints = 1:skipSize:(nFrames-nToAvg);
nTimePoints = length(startingTimePoints)-22;
           
area_dyn = zeros(nTimePoints,nComp);
freq_dyn = zeros(nTimePoints,nComp);
fwhm_dyn = zeros(nTimePoints,nComp);
phase_dyn = zeros(nTimePoints,nComp);
t_dyn = zeros(nTimePoints,nComp);

residualSpectrum = zeros(npts,nTimePoints);
fitSignal = zeros(512,1);
tzero = min(t(:));

tic

% phase(2) = mean(;
% nb = [inf inf inf inf];
% phase_ub = [inf, -100, inf, inf];
% phase_lb = [-inf, -105, -inf, -inf];

for k = 1:peaks
    disp([sprintf('%8.3f', nmrFit.area(k)/sum(nmrFit.area(end))), ' ' ...
        sprintf('%8.1f',(nmrFit.freq(k)-ref_freq)/xeFreqMHz),  '  ' ...
        sprintf('%8.1f',nmrFit.fwhm(k)/xeFreqMHz), '  ' ...
        sprintf('%9.1f',nmrFit.phase(k))]);
end

parfor iTimePoint = 1:nTimePoints
    startIdx = startingTimePoints(iTimePoint);
	avgdata = mean(fids(:,startIdx + (1:nToAvg) - 1),2);  
    avgdata_orig = mean(raw_fids(:,startIdx + (1:nToAvg) - 1),2); 
        
    nmrFit = NMR_TimeFit(avgdata, t, area,freq,fwhm,phase,[],[]);
%     nmrFit = nmrFit.setBounds( area_lb, area_ub, freq_lb, freq_ub,...
%         fwhm_lb, fwhm_ub, phase_lb, phase_ub);
%     nmrFit = nmrFit.setBounds(-nb,nb,freq_lb,freq_ub,fwhm_lb,fwhm_ub,phase_lb,phase_ub);
    nmrFit = nmrFit.fitTimeDomainSignal();
%     nmrFit.plotFit()

    area_dyn(iTimePoint,:) = nmrFit.area(:);
    freq_dyn(iTimePoint,:) = (nmrFit.freq(:)-ref_freq)/xeFreqMHz; %ppm units
    fwhm_dyn(iTimePoint,:) = nmrFit.fwhm(:)/xeFreqMHz; %ppm units
    phase_dyn(iTimePoint,:) = nmrFit.phase(:);
    
    startT = t_tr(startIdx);
    stuffThis = repmat(startT + floor(0.5*nToAvg*dwell_time),[1 nComp]);
    t_dyn(iTimePoint,:) = stuffThis;
    
    zeroPaddedTime = tzero + dwell_time*((1:nmrFit.zeroPadSize)-1)';
    fittedSignal = nmrFit.calcTimeDomainSignal(zeroPaddedTime);
    fittedSpectrum = dwell_time*fftshift(fft(fittedSignal));
    residualSpectrum(:,iTimePoint) = (nmrFit.spectralDomainSignal - fittedSpectrum);
    snrs(iTimePoint) = snr(fittedSignal,avgdata_orig-fittedSignal);     
end
toc

% Create Structure for Easy Saving of Data
fnames = {'area','freq','fwhm','phase','t','resSpec','ref_freq','snrs'};

dyn.(fnames{1}) = area_dyn;
dyn.(fnames{2}) = freq_dyn;
dyn.(fnames{3}) = fwhm_dyn;
dyn.(fnames{4}) = phase_dyn;
dyn.(fnames{5}) = t_dyn;
dyn.(fnames{6}) = residualSpectrum;
dyn.(fnames{7}) = ref_freq;
dyn.(fnames{8}) = snrs;

dyn.processed_date = datestr(now);
dyn.raw_file = twix_path;

 
