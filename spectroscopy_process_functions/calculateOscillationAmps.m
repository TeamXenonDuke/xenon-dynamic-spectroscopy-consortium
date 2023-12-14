function [amp, detrend, fitted, area_fit, rbcFit, gof] = calculateOscillationAmps(dyn,BHs,b,fit_type,varargin)

% calculateOscillationAmps returns the oscillation amplitude of all RBC spectral
% parameters (area, freq, fwhm, fwhmG, phase), along with the amplitude of
% the barrier oscillations (b), and the hr. The signal is first run though
% a high pass filter, and then fit to a sine function to determine the
% amplitude.
%
% INPUTS 
% dyn: structure containing dynamic spectroscopy fit information
% BHs: vector with time of BH start and BH end
% b: high pass filter
% 
% Optional: 
%
% getOscillationAmps(dyn,BHs,b,'oscType',method)
% available methods are:
% 'sine' - sine wave identification of amplitude
% 'peaks' - peak finding identification of amplitude
% 'peaks_lu' - peak finding identification of amplitude guided by sine fit

p = inputParser;
addParameter(p,'oscType',fit_type,...
   @(x) any(validatestring(x,{'sine','peaks', 'peaks_lu'})));
parse(p,varargin{:});
   
if isfield(dyn, 'fwhmL')
    dyn.fwhm = dyn.fwhmL;
    dyn.fwhmG(1,:) = [];
end 

% Remove first data point to avoid ringing when filtering 
dyn.t(1,:) = []; dyn.area(1,:) = []; dyn.freq(1,:) = []; dyn.fwhm(1,:) = [];
dyn.phase(1,:) = []; %dyn.snrs(1) = [];

if BHs(2) > dyn.t(end)
    BHs(2) = dyn.t(end);
end 

[BHstart, BHend] = findBHs(dyn.t(:,1), BHs);


%% Oscillation Information

clear opts
ft = fittype( 'sin1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Lower = [0 0 -Inf];
opts.Upper = [Inf Inf Inf];
xData = dyn.t(BHstart:BHend)';

% [rbcFit, gof] = fit(dyn.t(BHstart:BHend,1),dyn.area(BHstart:BHend,1),'exp2'); % biexponential decay
[rbcFit, gof] = fit(dyn.t(BHstart:BHend,1),dyn.area(BHstart:BHend,1),'exp1');
rbcNorm = rbcFit(dyn.t(:,1));

area_detrend = filtfilt(b,1,(dyn.area(:,1)-rbcNorm)./rbcNorm);
area_detrend = area_detrend(BHstart:BHend,1);

[area_fit, area_gof] = fit( xData, area_detrend, ft, opts );
areaVal = area_fit(dyn.t(BHstart:BHend,1));
amp_area = abs(area_fit.a1*2);
hr = area_fit.b1/(2*pi)*60;

freq_detrend = filtfilt(b,1,dyn.freq(:,1)-mean(dyn.freq(BHstart:BHend,1)));
freq_detrend = freq_detrend(BHstart:BHend,1);
    
fwhm_detrend = filtfilt(b,1,dyn.fwhm(:,1)-mean(dyn.fwhm(BHstart:BHend,1)));
fwhm_detrend = fwhm_detrend(BHstart:BHend,1);

dyn.phase(:,1) = rad2deg(unwrap(deg2rad(dyn.phase(:,1))));
phase_detrend = filtfilt(b,1,dyn.phase(:,1)-mean(dyn.phase(BHstart:BHend,1)));
phase_detrend = phase_detrend(BHstart:BHend,1);

% barrier fitting
bFit = fit(dyn.t(BHstart:BHend,1),dyn.area(BHstart:BHend,2),'exp1');
bNorm = bFit(dyn.t(:,1));
b_detrend = filtfilt(b,1,(dyn.area(:,2)-bNorm)./bNorm);
b_detrend = b_detrend(BHstart:BHend);
b_fit = fit( xData, b_detrend, ft, opts );
b_area = abs(b_fit.a1*2);

switch p.Results.oscType   
case 'sine' 
% use sine wave fitting to determine oscillation amplitudes

    % Set frequency to hr from amplitude
    opts.Lower = [0 area_fit.b1 -Inf];
    opts.Upper = [Inf area_fit.b1 Inf];

    [freq_fit, freq_gof] = fit( xData, freq_detrend, ft, opts );
    freqVal = freq_fit(dyn.t(BHstart:BHend,1));
    amp_freq = freq_fit.a1*2;

    fwhm_fit = fit( xData, fwhm_detrend, ft, opts );
    fwhmVal = fwhm_fit(dyn.t(BHstart:BHend,1));
    amp_fwhm = fwhm_fit.a1*2;

    phase_fit = fit( xData, phase_detrend, ft, opts );
    phaseVal = phase_fit(dyn.t(BHstart:BHend,1));
    amp_phase = phase_fit.a1*2;
    
    % Create Structure for Easy Saving of Data
    fnamesF = {'area','freq','fwhm','phase'};

    fitted.(fnamesF{1}) = areaVal;
    fitted.(fnamesF{2}) = freqVal;
    fitted.(fnamesF{3}) = fwhmVal;
    fitted.(fnamesF{4}) = phaseVal;
    
case 'peaks'
    % use peak identification to determine oscillation amplitudes
    
    tr = dyn.t(2) - dyn.t(1);
    
    [~, ay] = findpeaks(abs(smooth(area_detrend)), 'MinPeakDistance', 68/60/6*1/tr);
    area_max = area_detrend(ay); ay = dyn.t(ay+BHstart-1,1);
    % This section of the code will force only one point from peak and
    % trough, considering there are maximum two peaks and troughs available
    % per period.
    area_max2 = area_max;
    area_max_sign = sign(area_max);
    for i = 1:size(area_max2)   
        if i == size(area_max2,1)
            break
        elseif area_max_sign(i) == area_max_sign(i+1)

            if area_max(i)> abs(area_max(i+1))
                area_max2(i+1) = -1;
            else
                area_max2(i) = -1;
            end
        end

        mask = area_max2~=-1.0;
    end
    
    area_max = area_max2.* mask; 
    area_max(area_max==0) = [];
    ay = ay.* mask; ay(ay==0) = [];

    amp_area = median(area_max(area_max>0))-median(area_max(area_max<0));
    
    [~, fy] = findpeaks(abs(smooth(freq_detrend)),'MinPeakDistance', 68/60/6*1/tr);
    freq_max = freq_detrend(fy); fy = dyn.t(fy+BHstart-1,1);
    amp_freq = median(freq_max(freq_max>0))-median(freq_max(freq_max<0));
    
    [~, fwy] = findpeaks(abs(smooth(fwhm_detrend)),'MinPeakDistance', 68/60/6*1/tr);
    fwhm_max = fwhm_detrend(fwy); fwy = dyn.t(fwy+BHstart-1,1);
    amp_fwhm = median(fwhm_max(fwhm_max>0))-median(fwhm_max(fwhm_max<0));
    
    [~, py] = findpeaks(abs(smooth(phase_detrend)),'MinPeakDistance', 68/60/6*1/tr);
    phase_max = phase_detrend(py); py = dyn.t(py+BHstart-1,1);
    amp_phase = median(phase_max(phase_max>0))-median(phase_max(phase_max<0));
    
   
%     % Plot peak identification results
%     figure(11), clf 
%     subplot(4,1,1), set(gcf, 'Position',[875    -200    550    725]), hold on
%     plot(dyn.t(BHstart:BHend,1),area_detrend*100,'k--')
%     plot(dyn.t(BHstart:BHend,1),areaVal*100,'color',[0, 0.4470, 0.7410],'Linewidth',1.5)
%     plot(dyn.t(BHstart:BHend,1),median(area_max(area_max>0))*100*ones(1,length(BHstart:BHend)),'--','Linewidth',.25,'Color',[0, 0.4470, 0.7410])
%     plot(dyn.t(BHstart:BHend,1),median(area_max(area_max<0))*100*ones(1,length(BHstart:BHend)),'--','Linewidth',.25,'Color',[0, 0.4470, 0.7410])
%     plot(ay,area_max*100,'.','color',[0.8500    0.3250    0.0980] ,'MarkerSize',30);
%     plot(dyn.t(BHstart:BHend,1),zeros(1,length(BHstart:BHend)),'Linewidth',.25,'Color',[.5 .5 .5])
%     xlim(BHs), ylabel('Amplitude (%)')
%     subplot(4,1,2), hold on
%     plot(dyn.t(BHstart:BHend,1),freq_detrend,'k--')
%     plot(dyn.t(BHstart:BHend,1),smooth(freq_detrend,20),'Linewidth',2)
%     plot(fy,freq_max,'.','color',[0.8500    0.3250    0.0980] ,'MarkerSize',30);
%     plot(dyn.t(BHstart:BHend,1),zeros(1,length(BHstart:BHend)),'Linewidth',.25,'Color',[.5 .5 .5])
%     xlim(BHs), ylabel('Freq (ppm)')
%     subplot(4,1,3), hold on
%     plot(dyn.t(BHstart:BHend,1),fwhm_detrend,'k--')
%     plot(dyn.t(BHstart:BHend,1),smooth(fwhm_detrend,20),'Linewidth',2)
%     plot(fwy,fwhm_max,'.','color',[0.8500    0.3250    0.0980] ,'MarkerSize',30);
%     plot(dyn.t(BHstart:BHend,1),zeros(1,length(BHstart:BHend)),'Linewidth',.25,'Color',[.5 .5 .5])
%     xlim(BHs), ylabel('FWHM (ppm)')
%     subplot(4,1,4), hold on
%     plot(dyn.t(BHstart:BHend,1),phase_detrend,'k--')
%     plot(dyn.t(BHstart:BHend,1),smooth(phase_detrend,20),'Linewidth',2)
%     plot(py,phase_max,'.','color',[0.8500    0.3250    0.0980] ,'MarkerSize',30); 
%     plot(dyn.t(BHstart:BHend,1),zeros(1,length(BHstart:BHend)),'Linewidth',.25,'Color',[.5 .5 .5])
%     xlim(BHs), xlabel('Time (s)'), ylabel('Phase ({\circ})')
    
    % Create Structure for Easy Saving of Data
    fnamesF = {'area','at','freq','ft','fwhm','fwt','phase','pt','area_max','area_min',...
            'freq_max','freq_min','fwhm_max','fwhm_min','phase_max','phase_min','area_sine'};

    fitted.(fnamesF{1}) = area_max; fitted.(fnamesF{2}) = ay;
    fitted.(fnamesF{3}) = freq_max; fitted.(fnamesF{4}) = fy;
    fitted.(fnamesF{5}) = fwhm_max; fitted.(fnamesF{6}) = fwy;
    fitted.(fnamesF{7}) = phase_max; fitted.(fnamesF{8}) = py;
    fitted.(fnamesF{9}) = median(area_max(area_max>0)); fitted.(fnamesF{10}) = median(area_max(area_max<0));
    fitted.(fnamesF{11}) = median(freq_max(freq_max>0)); fitted.(fnamesF{12}) = median(freq_max(freq_max<0));
    fitted.(fnamesF{13}) = median(fwhm_max(fwhm_max>0)); fitted.(fnamesF{14}) = median(fwhm_max(fwhm_max<0));
    fitted.(fnamesF{15}) = median(phase_max(phase_max>0)); fitted.(fnamesF{16}) = median(phase_max(phase_max<0));
    fitted.(fnamesF{17}) = areaVal;
case 'peaks_lu'
    % custom peak finding using sine fit to guide peak finding
    tr = dyn.t(2) - dyn.t(1);
    [~, ay_sine] = findpeaks(abs(areaVal), 'MinPeakDistance', 68/60/6*1/tr);
    % initialize area_max, freq_max, fwhm_max, and phase_max
    area_max = zeros(length(ay_sine), 1);
    freq_max = zeros(length(ay_sine), 1);
    fwhm_max = zeros(length(ay_sine), 1);
    phase_max = zeros(length(ay_sine), 1);
    % initialize ay, fy, fwy, and py
    ay = zeros(length(ay_sine), 1);
    fy = zeros(length(ay_sine), 1);
    fwy = zeros(length(ay_sine), 1);
    py = zeros(length(ay_sine), 1);
    % for loop to iterate over each peak in the sine fit
    for i = 1:length(ay_sine)
        ind_start = max(ay_sine(i)-floorDiv(median(diff(ay_sine)), 2), 1);
        ind_end = min(ay_sine(i)+floorDiv(median(diff(ay_sine)), 2), length(area_detrend));
        partial_area = zeros(length(area_detrend), 1);
        partial_freq = zeros(length(area_detrend), 1);
        partial_fwhm = zeros(length(area_detrend), 1);
        partial_phase = zeros(length(area_detrend), 1);
        % find only the positive peak if sine is positive else find
        % negative peak
        if areaVal(ay_sine(i)) >=0
            % smooth the area_detrend
            signal_area = smooth(area_detrend);
            signal_freq = freq_detrend;
            signal_fwhm = fwhm_detrend;
            signal_phase = phase_detrend;
        else
            signal_area = -smooth(area_detrend);
            signal_freq = -freq_detrend;
            signal_fwhm = -fwhm_detrend;
            signal_phase = -phase_detrend;
        end
        % loop update
        partial_area(ind_start:ind_end) = signal_area(ind_start:ind_end);
        [~, ay_single] = findpeaks(partial_area, 'NPeaks', 1, 'SortStr', 'Descend');
        % if a peak is found, update. peak will be found most of the time.
        if ~isempty(ay_single)
            area_max(i) = area_detrend(ay_single);
            ay(i) = dyn.t(ay_single+BHstart-1,1);
        end
        partial_freq(ind_start:ind_end) = signal_freq(ind_start:ind_end);
        [~, fy_single] = findpeaks(partial_freq, 'NPeaks', 1, 'SortStr', 'Descend');
        if ~isempty(fy_single)
            freq_max(i) = freq_detrend(fy_single);
            fy(i) = dyn.t(fy_single+BHstart-1,1);
        end
        partial_fwhm(ind_start:ind_end) = signal_fwhm(ind_start:ind_end);
        [~, fwy_single] = findpeaks(partial_fwhm, 'NPeaks', 1, 'SortStr', 'Descend');
        if ~isempty(fwy_single)
            fwhm_max(i) = fwhm_detrend(fwy_single);
            fwy(i) = dyn.t(fwy_single+BHstart-1,1);
        end
        partial_phase(ind_start:ind_end) = signal_phase(ind_start:ind_end);
        [~, py_single] = findpeaks(partial_phase, 'NPeaks', 1, 'SortStr', 'Descend');
        if ~isempty(py_single)
            phase_max(i) = phase_detrend(py_single);
            py(i) = dyn.t(py_single+BHstart-1,1);
        end
    end
    % remove entries that are zeros, which were not updated since no peak
    % found.
    area_max = nonzeros(area_max);
    ay = nonzeros(ay);
    freq_max = nonzeros(freq_max);
    fy = nonzeros(fy);
    fwhm_max = nonzeros(fwhm_max);
    fwy = nonzeros(fwy);
    phase_max = nonzeros(phase_max);
    py = nonzeros(py);
    % calculate the median amp_area, amp_freq, amp_fwhm, and amp_phase
    amp_area = median(area_max(area_max>0))-median(area_max(area_max<0));
    amp_freq = median(freq_max(freq_max>0))-median(freq_max(freq_max<0));
    amp_fwhm = median(fwhm_max(fwhm_max>0))-median(fwhm_max(fwhm_max<0));
    amp_phase = median(phase_max(phase_max>0))-median(phase_max(phase_max<0));
    % Create Structure for Easy Saving of Data
    fnamesF = {'area','at','freq','ft','fwhm','fwt','phase','pt','area_max','area_min',...
            'freq_max','freq_min','fwhm_max','fwhm_min','phase_max','phase_min','area_sine'};

    fitted.(fnamesF{1}) = area_max; fitted.(fnamesF{2}) = ay;
    fitted.(fnamesF{3}) = freq_max; fitted.(fnamesF{4}) = fy;
    fitted.(fnamesF{5}) = fwhm_max; fitted.(fnamesF{6}) = fwy;
    fitted.(fnamesF{7}) = phase_max; fitted.(fnamesF{8}) = py;
    fitted.(fnamesF{9}) = median(area_max(area_max>0)); fitted.(fnamesF{10}) = median(area_max(area_max<0));
    fitted.(fnamesF{11}) = median(freq_max(freq_max>0)); fitted.(fnamesF{12}) = median(freq_max(freq_max<0));
    fitted.(fnamesF{13}) = median(fwhm_max(fwhm_max>0)); fitted.(fnamesF{14}) = median(fwhm_max(fwhm_max<0));
    fitted.(fnamesF{15}) = median(phase_max(phase_max>0)); fitted.(fnamesF{16}) = median(phase_max(phase_max<0));
    fitted.(fnamesF{17}) = areaVal;
end 
    
% Create Structure for Easy Saving of Data
switch p.Results.oscType   
    case 'sine'     
        fnames = {'area','freq','fwhm','phase','hr','b','area_gof', 'freq_gof'}; 
        amp.(fnames{1}) = amp_area;
        amp.(fnames{2}) = amp_freq;
        amp.(fnames{3}) = amp_fwhm;
        amp.(fnames{4}) = amp_phase;
        amp.(fnames{5}) = hr;
        amp.(fnames{6}) = b_area;
        amp.(fnames{7}) = area_gof;
        amp.(fnames{8}) = freq_gof;
    case 'peaks'
        fnames = {'area','freq','fwhm','phase','hr','b','area_gof'};
        amp.(fnames{1}) = amp_area;
        amp.(fnames{2}) = amp_freq;
        amp.(fnames{3}) = amp_fwhm;
        amp.(fnames{4}) = amp_phase;
        amp.(fnames{5}) = hr;
        amp.(fnames{6}) = b_area;
        amp.(fnames{7}) = area_gof;
    case 'peaks_lu'
        fnames = {'area','freq','fwhm','phase','hr','b','area_gof'};
        amp.(fnames{1}) = amp_area;
        amp.(fnames{2}) = amp_freq;
        amp.(fnames{3}) = amp_fwhm;
        amp.(fnames{4}) = amp_phase;
        amp.(fnames{5}) = hr;
        amp.(fnames{6}) = b_area;
        amp.(fnames{7}) = area_gof;
end 

% Create Structure for Easy Saving of Data
fnamesD = {'area','freq','fwhm','phase'};

detrend.(fnamesD{1}) = area_detrend;
detrend.(fnamesD{2}) = freq_detrend;
detrend.(fnamesD{3}) = fwhm_detrend;
detrend.(fnamesD{4}) = phase_detrend;

end 
 