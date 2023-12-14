function [amp, nmrFit_ppm] = dynamicRBCtoBarrier(raw_path,dyn_path,BHs,save_fig_flag,figname,suppressOutput)
% raw_path: raw data file either .dat or .7
% dyn_path: dynamic spectroscopy matlab structure or structure location
% BHs: either a vector containing the time of inhale/exhale or the filepath
%      for a saved matlab vector with the BH information
% save_fig_flag: logical input for saving figures. Figures are saved in the
%                same folder as dyn_path.
% figname: adds text after default figure name (e.g. '_2' -> dynV_2)

%% Example inputs
% raw_path = 'D:\Elly\Documents\Duke\CIVM Research\Siemens Data\1 - Healthy Volunteers\18-10-01 SUBJECT 000-001C\2 - Raw\meas_MID00026_FID69955_Xe_fid_DynamicSpec_high_flip.dat';
%
% 
% dyn_name = 'dynV5_2';
% 
% dyn_path = [path(1:end-7),'4 - Spectra\',dyn_name,'.mat'];
%
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
    [dyn_folder, ~, scanner] = fileparts(raw_path);
else 
    [dyn_folder, ~, scanner] = fileparts(raw_path);
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

if isempty(BHs)
    BHs = [2 10];
elseif ischar(BHs)
    load(BHs);
elseif isvector(BHs)
end 

interest = [0 dyn.t(end)];
[~, nComp] = size(dyn.t);
peaks = nComp;
    
BHstart = find(round(dyn.t(:,1),2) == BHs(1));
BHend = find(round(dyn.t(:,1),2) == BHs(2));

if ~exist('BHend','var')
    BHend = length(dyn.t);
end 
m3rd = BHstart:BHend;

if ~isfield(dyn,'fwhmG') % if the data is not Voigt, set all FWHM_G to 0
    dyn.fwhmG = zeros(size(dyn.t));
    dyn.fwhmL = dyn.fwhm;
    fitType = 'L';
else
    fitType = 'V';
end 

%% Create Dynamic By Resonance Figure 

% Colors
colors = [0.8500    0.3250    0.0980 
          0.4660    0.6740    0.1880
          0         0.4470    0.7410];

plotlim = [2 5 1.5; 1.5 .5 .5; 5 5 .5; 30 6 12]; % area, freq, fwhm, phase
         
vSpace = .03;
wSpace = .05;
leftSpace = .06;
rightSpace = 0.075;
topSpace = .1;
bottomSpace = .12;
plotHeight = ((1-topSpace)-3*vSpace-bottomSpace)/4;
plotWidth = ((1-rightSpace)-2*wSpace-wSpace)/3;

linewidth = 1.5;
font_label = 12;
font_axis = 12;

t_plot = dyn.t(:,1);
startInhale = 1;

if suppressOutput ~= 1
    
figure(1), clf
set(gcf,'Position',[20 -200 1200 675])
axHandles = zeros(1,5*nComp);

for iComp = 1:nComp
       
    % Area
    axHandles(5*(iComp-1)+1) = subplot(4,nComp,iComp); hold on;
    plot(t_plot(1:end-5),dyn.area(startInhale:end-5,iComp)/max(dyn.area(50:end,end)),'-','LineWidth',linewidth,'color',colors(iComp,:));
    set(gca,'Xtick','','FontSize',font_axis)
    set(gca,'Position',[leftSpace+(iComp-1)*plotWidth+(iComp-1)*wSpace bottomSpace+3*vSpace+3*plotHeight plotWidth plotHeight])
    if iComp == 1
        ylabel('Amplitude','FontSize',font_label);
    end
    ylim([0 max(dyn.area(10:end,iComp))/max(dyn.area(50:end,end))])
%     ylim([0 plotlim(1,iComp)])
%     title(titles(iComp),'FontSize',20);
    
    for idx = 1:length(BHs)/2  
        fill([BHs(2*idx-1) BHs(2*idx-1) BHs(2*idx) BHs(2*idx)],...
            [-500 500 500 -500],[.75 .75 .75],'LineStyle','none','FaceAlpha',.5)
    end   
    set(gca,'Layer','top')

    % Frequency
    axHandles(5*(iComp-1)+2) = subplot(4,nComp,nComp+iComp); hold on;
    plot(t_plot,dyn.freq(startInhale:end,iComp),'-','LineWidth',linewidth,'color',colors(iComp,:));
    set(gca,'Xtick','','FontSize',font_axis)
    set(gca,'Position',[leftSpace+(iComp-1)*plotWidth+(iComp-1)*wSpace bottomSpace+2*vSpace+2*plotHeight plotWidth plotHeight])
    if iComp == 1
        ylabel('Shift (ppm)','FontSize',font_label); %ylabel('Frequency (Hz)');
    end 
    ylim([mean(dyn.freq(m3rd,iComp))-plotlim(2,iComp),mean(dyn.freq(m3rd,iComp))+plotlim(2,iComp)])

    for idx = 1:length(BHs)/2  
        fill([BHs(2*idx-1) BHs(2*idx-1) BHs(2*idx) BHs(2*idx)],...
            [-500 500 500 -500],[.75 .75 .75],'LineStyle','none','FaceAlpha',.5)
    end   
    set(gca,'Layer','top')
    
    
    % Linewidth
    axHandles(5*(iComp-1)+3) = subplot(4,nComp,2*nComp+iComp); hold on;
    
    if iComp == 2 && sum(dyn.fwhmG(:)) > 0
        h(1,:) = plot(dyn.t,dyn.fwhmL(1:end,iComp),'-','LineWidth',1,'color',[0.4660    0.6740    0.1880]);
        h(2,:) = plot(dyn.t,dyn.fwhmG(1:end,iComp),':','LineWidth',2,'color',[0.4660    0.6740    0.1880]*.7);
        ymin = min(mean(dyn.fwhmL(m3rd,iComp))-plotlim(3,iComp),mean(dyn.fwhmG(m3rd,iComp))-plotlim(3,iComp));
        ymax = max(mean(dyn.fwhmL(m3rd,iComp))+plotlim(3,iComp),mean(dyn.fwhmG(m3rd,iComp))+plotlim(3,iComp));
        ylim([ymin ymax])
    else
    plot(dyn.t,dyn.fwhmL(1:end,iComp),'-','LineWidth',linewidth,'color',colors(iComp,:));   
    ylim([mean(dyn.fwhmL(m3rd,iComp))-plotlim(3,iComp),mean(dyn.fwhmL(m3rd,iComp))+plotlim(3,iComp)]) 
    end 
    
    set(gca,'Xtick','','FontSize',font_axis)
    set(gca,'Position',[leftSpace+(iComp-1)*plotWidth+(iComp-1)*wSpace bottomSpace+vSpace+plotHeight plotWidth plotHeight])
    
    if iComp == 1
        ylabel('FWHM (ppm)','FontSize',font_label);
    end 
    
    for idx = 1:length(BHs)/2  
        fill([BHs(2*idx-1) BHs(2*idx-1) BHs(2*idx) BHs(2*idx)],...
            [-500 500 500 -500],[.75 .75 .75],'LineStyle','none','FaceAlpha',.5)
    end   
    set(gca,'Layer','top')
    
    if iComp == 2 && sum(dyn.fwhmG(:)) > 0
        legend(h(1:2),{'FWHM','FWHM_G'},'fontSize',8,'Location','SouthEast')
    end 
    
    % Phase
    axHandles(5*(iComp-1)+5) = subplot(5,nComp,4*nComp + iComp);
    hold on;
%     dyn.phase(startInhale:end,iComp) = rad2deg(unwrap(deg2rad(dyn.phase(startInhale:end,iComp))));
%     dyn.phase(:,iComp) = rad2deg(unwrap(deg2rad(dyn.phase(:,iComp))));
    plot(t_plot,dyn.phase(startInhale:end,iComp),'-','LineWidth',linewidth,'color',colors(iComp,:));
    set(gca,'FontSize',font_axis)
    set(gca,'Position',[leftSpace+(iComp-1)*plotWidth+(iComp-1)*wSpace bottomSpace plotWidth plotHeight])
    if iComp == 1
        ylabel('Phase','FontSize',font_label);
    end 
    xlabel('Time (sec)','FontSize',font_label);
    ylim([mean(dyn.phase(m3rd,iComp))-plotlim(4,iComp),mean(dyn.phase(m3rd,iComp))+plotlim(4,iComp)])
    
    for idx = 1:length(BHs)/2  
        fill([BHs(2*idx-1) BHs(2*idx-1) BHs(2*idx) BHs(2*idx)],...
            [-500 500 500 -500],[.75 .75 .75],'LineStyle','none','FaceAlpha',.5)
    end  
    set(gca,'Layer','top')
end

warning('off','MATLAB:linkaxes:RequireDataAxes')
linkaxes(axHandles,'x'); 
xlim([0 dyn.t(end)])

if save_fig_flag == 1
    options.Format = 'tiff';
    hgexport(gcf,[dyn_folder,'\dyn',fitType,figname,'.tif'],options)   
end

%% Plot Detrended RBC Oscillation and Calculate Amplitudes
% rbc_signal_max = max(dyn.area(:,1));
% first_real_signal = find(dyn.area(:,1)> max(dyn.area(:,1)/10),1);
% last_real_signal = find(dyn.area(:,1)> max(dyn.area(:,1)/10),1,'last');
% b = highpassfilter(last_real_signal-first_real_signal);
b = highpassfilter(length(dyn.area(:,1)));
[amp, detrend, fitted] = getOscillationInfo(dyn,BHs,b);
% if exist('dyn.snrs')
    amp.snr = mean(dyn.snrs(BHstart:BHend));
% else 
%     amp.snr = zeros(1,length(dyn.area));
%     dyn.snrs = zeros(1,length(dyn.area));
% end 
figure(2), clf
set(gcf, 'Position',[875    -200    550    725])
axFontFs = 13;

subplot(4,1,1), hold on
plot(dyn.t(BHstart:BHend),detrend.area*100,'k--')
plot(dyn.t(BHstart:BHend),fitted.area*100,'Linewidth',2,'color',colors(1,:)); 
xlim([dyn.t(BHstart) dyn.t(BHend)])
xlabel('Time (s)'), ylabel('Amplitude')
set(gca,'FontSize',axFontFs)
ylim([-10 10])

subplot(4,1,2), hold on
plot(dyn.t(BHstart:BHend),detrend.freq,'k--')
plot(dyn.t(BHstart:BHend),fitted.freq,'Linewidth',2,'color',colors(1,:)); 
xlim([dyn.t(BHstart) dyn.t(BHend)])
xlabel('Time (s)'), ylabel('Chemical Shift (ppm)')
set(gca,'FontSize',axFontFs)
ylim([-.5 .5])

subplot(4,1,3), hold on
plot(dyn.t(BHstart:BHend),detrend.fwhm,'k--')
plot(dyn.t(BHstart:BHend),fitted.fwhm,'Linewidth',2,'color',colors(1,:)); 
xlim([dyn.t(BHstart) dyn.t(BHend)])
xlabel('Time (s)'), ylabel('FWHM (ppm)')
set(gca,'FontSize',axFontFs)
ylim([-.5 .5])

subplot(4,1,4), hold on
plot(dyn.t(BHstart:BHend),detrend.phase,'k--')
plot(dyn.t(BHstart:BHend),fitted.phase,'Linewidth',2,'color',colors(1,:)); 
xlim([dyn.t(BHstart) dyn.t(BHend)])
xlabel('Time (s)'), ylabel('Phase (degrees)')
set(gca,'FontSize',axFontFs)
ylim([-10 10])

if save_fig_flag == 1
    options.Format = 'tiff';
    hgexport(gcf,[dyn_folder,'\dyn',fitType,'_oscs_new',figname,'.tif'],options)
end  
%% Plot SNR of Dynamic Spectroscopy

figure(3), clf
set(gcf,'Position',[275 -100 450 200]), hold on
plot(dyn.t,dyn.snrs(1:end),'-','LineWidth',linewidth,'color',colors(iComp,:));
for idx = 1:length(BHs)/2  
    fill([BHs(2*idx-1) BHs(2*idx-1) BHs(2*idx) BHs(2*idx)],...
        [-500 500 500 -500],[.75 .75 .75],'LineStyle','none','FaceAlpha',.5)
end   
set(gca,'Layer','top')
xlim(interest)
ylim([10 30])
xlabel('Time (s)'), ylabel('SNR')

if save_fig_flag == 1
    options.Format = 'tiff';
    hgexport(gcf,[dyn_folder,'\dyn',fitType,'_error',figname,'.tif'],options)   
end

%% Display oscillation results
disp([10 'file: ',raw_path])
disp([10 'RBC Dynamic Oscillation Amplitudes:'])
disp('   Area      Freq      FWHM      Phase      HR      SNR');
disp([sprintf('%8.1f', amp.area*100), ' ' ...
    sprintf('%8.2f',amp.freq),  '  ' ...
    sprintf('%8.2f',amp.fwhm), '  ' ...
    sprintf('%8.1f',amp.phase), '   '...
    sprintf('%6.0f',amp.hr),'  ',...
    sprintf('%7.1f',amp.snr)]);

%% Calculate Static Spectral Parameters
 
% Read in twix or P file and define associated variables
if strcmp(scanner,'.dat')
    % Twix file from Siemens
    twix = readtwix(raw_path);
    npts = twix.hdr.Config.RawCol*2;                   % Number of samples                  % Receiver bandwidth (kHz)
    dwell_time = twix.hdr.Config.DwellTime*10^-9;                                 % Time between each sample
    fids = twix.data;
    xeGam = twix.hdr.Meas.lFrequency*10e-7; %34.091516
elseif strcmp(scanner,'.7')
    % Pfile from GE
    pfile = GE.Pfile.read(raw_path);
    MRI.DataProcessing.checkForOverranging(pfile); % Check for overranging
    pfile = MRI.DataProcessing.removeBaselineViews(pfile); % Remove baselines
    bw = 1000*pfile.rdb.rdb_hdr_user12;                    % Receiver bandwidth (kHz)
    dwell_time = 1/(2*bw);                                 % Time between each sample
    dwell_time = Math.nearestMultipleOf(dwell_time,0.000002);
    npts = pfile.rdb.rdb_hdr_frame_size;
    fids = pfile.data;
    xeGam = 17.660445;
else 
    error('Unknown Raw File Type')
end 
t = dwell_time*(0:(npts-1))';

% Constants
xe_shift = (0.3/5.5)*(273/298)*(0.548);
tot_shift = -xe_shift;
tot_shift_hz = tot_shift*xeGam;

% Three peak fit
area_orig = [1 1 1];
freq_orig = [0 -12 -220]*xeGam;
% freq_orig = [0 -12 -218]*xeGam+7500; For 90 degree spoiled
fwhm_orig = [7 6 2]*xeGam;
fwhmL_orig = [7 4 2]*xeGam;
fwhmG_orig = [0 5 0]*xeGam;
phase_orig = [0 0 0];
ngas = 3;
nbar = 2;

% Referance values from healthy population
area_ref = [0.3598    0.6131    0.0565];
freq_ref = [216.0468  197.3733   -0.3411];
fwhm_ref = [9.9863    8.6263    2.0852];
phase_ref = [-48.2896 -140.9711  -66.7937];

ymax = 1;

avgdata = mean(fids(:,BHstart:BHstart+50),2);

if sum(dyn.fwhmG(:)) > 0  
    nmrFit = NMR_TimeFit_v(avgdata, t, area_orig,freq_orig,fwhmL_orig,fwhmG_orig,phase_orig,[],[]);
else 
    nmrFit = NMR_TimeFit(avgdata, t, area_orig,freq_orig,fwhm_orig,phase_orig,[],[]);
end 
nmrFit = nmrFit.fitTimeDomainSignal();

nmrFit_ref = NMR_TimeFit(avgdata, t, area_ref,freq_ref*xeGam+nmrFit.freq(ngas),fwhm_ref*xeGam,phase_ref,[],[]);
%%
ref_freq = nmrFit.freq(end);
f = xeGam*150+ref_freq:xeGam*250+ref_freq;
lorentz_ref = calcComponentLorentzianCurves(nmrFit_ref,f);
ppm_shift = (nmrFit.freq(ngas)-ref_freq)/xeGam; % shift reference spectrum to match calculated
ppm = (f+tot_shift_hz-ref_freq)/xeGam;

% Calculate fitted and residual spectrums usign time domain
% signal so that amplitudes are correct even if signal is
% truncated at the end

% Calculate spectrum
fittedFID = nmrFit.calcTimeDomainSignal(nmrFit.t);
fittedSpectrum = dwell_time*fftshift(fft(fittedFID));
residualSpectrum = nmrFit.spectralDomainSignal - fittedSpectrum;

F = linspace(-0.5,0.5,length(nmrFit.t)*2+1)/dwell_time;
F = F(1:(end-1)); % Take off last sample to have nSamples
PPM = (F+tot_shift_hz-ref_freq)/xeGam;

tplot = linspace(nmrFit.t(1),2*nmrFit.t(end),length(F));
fittedTime = nmrFit.calcComponentTimeDomainSignal(tplot);
fittedSpectrumComp = dwell_time*fftshift(fft(fittedTime));

% SNR = snr(fittedSpectrum,residualSpectrum);
avgdata = mean(fids(:,BHstart:BHstart + 4),2);

if sum(dyn.fwhmG(:)) > 0  
    nmrFit_SNR = NMR_TimeFit_v(avgdata, t, area_orig,freq_orig,fwhmL_orig,fwhmG_orig,phase_orig,[],[]);
else 
    nmrFit_SNR = NMR_TimeFit(avgdata, t, area_orig,freq_orig,fwhm_orig,phase_orig,[],[]);
end 
nmrFit_SNR = nmrFit_SNR.fitTimeDomainSignal();

timeFit=nmrFit.calcTimeDomainSignal(nmrFit_SNR.t);  % calculate fitted data in time domain
timeRes=timeFit - nmrFit_SNR.timeDomainSignal;
n25pct = round(length(timeRes)/4);
std25 = std(timeRes(end-n25pct:end)); % take standard deviation of tail of residuals
SNR_dis=nmrFit_SNR.area/std25; % will be a matrix for multiple peak fit

%%
figure(4), clf
ax2 = subplot(3,1,1); hold on
set(gcf,'Position',[1275 -100 600 600]); 
plot(ppm-ppm_shift,lorentz_ref/sum(max(lorentz_ref(:,nbar))),':k','Linewidth',2)

plotorder = [2 3 1];
for i = 1:peaks
    plot(PPM,abs(fittedSpectrumComp(:,plotorder(i)))/abs(max(fittedSpectrumComp(:,3))),'Color',colors(i,:),'Linewidth',3)  
end 
xticks(150:20:250), set(gca,'XDir','reverse'), box on
ylabel('Component Intensity'), xlabel('Chemical Shift (ppm)');
ylim([0 ymax]), xlim([150 250])

ax4 = subplot(3,1,2); hold on
plot((nmrFit.f-ref_freq)/xeGam,real(nmrFit.spectralDomainSignal),'.k','markersize',16);
plot((nmrFit.f-ref_freq)/xeGam,real(fittedSpectrum),'-g','Linewidth',2);
plot((nmrFit.f-ref_freq)/xeGam,real(residualSpectrum),'.r','markersize',8);
xticks(150:20:250), set(ax4,'XDir','reverse'), box on
ylabel('Real'), xlabel('Chemical Shift (ppm)');
legend({'Measured','Fitted','Residual'},'Location','SouthEast');

ax5 = subplot(3,1,3); hold on
plot((nmrFit.f-ref_freq)/xeGam,imag(nmrFit.spectralDomainSignal),'.k','markersize',16);
plot((nmrFit.f-ref_freq)/xeGam,imag(fittedSpectrum),'-g','Linewidth',2);
plot((nmrFit.f-ref_freq)/xeGam,imag(residualSpectrum),'.r','markersize',8);
xlabel('Chemical Shift (ppm)');
ylabel('Imaginary');

% Keep all x axes in sinc
linkaxes([ax2,ax4 ax5],'x');   
xlim([150 250])

if save_fig_flag == 1
    options.Format = 'tiff';
    hgexport(gcf,[dyn_folder,'\static',fitType,figname,'.tif'],options)   
end  


positive_phase = nmrFit.phase - nmrFit.phase(2);
positive_phase(positive_phase < 0) = positive_phase(positive_phase < 0) + 360;

nmrFit_ppm.area = nmrFit.area/sum(nmrFit.area(nbar));
nmrFit_ppm.freq = (nmrFit.freq-ref_freq)/xeGam;
nmrFit_ppm.fwhm = nmrFit.fwhm/xeGam;
nmrFit_ppm.phase = positive_phase;
nmrFit_ppm.SNR_dis = SNR_dis;

if sum(dyn.fwhmG(:)) > 0
    disp([10 'Static Parameters:'])
    disp('    Area     Freq      FWHM      FWHMg     Phase       SNR');
    nmrFit.fwhmG = abs(nmrFit.fwhmG);
    for k = 1:peaks
        disp([sprintf('%8.2f', nmrFit.area(k)/sum(nmrFit.area(nbar))), ' ' ...
            sprintf('%8.1f',(nmrFit.freq(k)-ref_freq)/xeGam),  '  ' ...
            sprintf('%8.1f',nmrFit.fwhm(k)/xeGam), '  ' ...
            sprintf('%8.1f',nmrFit.fwhmG(k)/xeGam), '  ' ...
            sprintf('%9.1f',positive_phase(k)), '  '...
            sprintf('%9.1f',SNR_dis(k))]);
    end
    nmrFit_ppm.fwhmG = nmrFit.fwhmG/xeGam;
else 
    disp([10 'Static Parameters:'])
    disp('    Area     Freq      FWHM       Phase       SNR');
    nmrFit.phase(nmrFit.phase < 0) = nmrFit.phase(nmrFit.phase < 0)+360;
    for k = 1:peaks
        disp([sprintf('%8.2f', nmrFit.area(k)/sum(nmrFit.area(nbar))), ' ' ...
            sprintf('%8.1f',(nmrFit.freq(k)-ref_freq)/xeGam),  '  ' ...
            sprintf('%8.1f',nmrFit.fwhm(k)/xeGam), '  ' ...
            sprintf('%9.1f',positive_phase(k)), '  '...
            sprintf('%9.1f',SNR_dis(k))]);
    end
end 

%% RBC:bar plot
       
figure(5), clf
set(gcf,'Position',[247         347        1180         556])
hold on

% Area
plot(t_plot(1:end-5),dyn.area(startInhale:end-5,1)./dyn.area(startInhale:end-5,2),'-','LineWidth',3,'color',colors(3,:));
plot(t_plot(1:end-5),dyn.fwhmL(startInhale:end-5,3)./max(dyn.fwhmL(BHstart:end,3)),'-','LineWidth',2,'color',colors(2,:));
plot(t_plot(1:end-5),dyn.freq(startInhale:end-5,3)./max(dyn.freq(BHstart:end,3)+0.5)+0.5,'-','LineWidth',2,'color',colors(1,:));
set(gca,'FontSize',font_axis)
ylabel('RBC:barrier','FontSize',font_label);
xlabel('time (s)')
%     ylim([0 plotlim(1,iComp)])
%     title(titles(iComp),'FontSize',20);

for idx = 1:length(BHs)/2  
    fill([BHs(2*idx-1) BHs(2*idx-1) BHs(2*idx) BHs(2*idx)],...
        [-500 500 500 -500],[.75 .75 .75],'LineStyle','none','FaceAlpha',.5)
end   
set(gca,'Layer','top')
ylim([0.2 1])

if save_fig_flag == 1
    options.Format = 'tiff';
    hgexport(gcf,[dyn_folder,'\RBC_to_Bar',fitType,figname,'.tif'],options) 
    amp = 0; nmrFit_ppm = 0;
end  

else
%% RBC:bar plot
       
figure(5), clf
set(gcf,'Position',[247         347        1180         556])
hold on

% Area
plot(t_plot(1:end-5),dyn.area(startInhale:end-5,1)./dyn.area(startInhale:end-5,2),'-','LineWidth',3,'color',colors(3,:));
% plot(t_plot(1:end-5),dyn.fwhmL(startInhale:end-5,3)./max(dyn.fwhmL(BHstart:end,3)),'-','LineWidth',2,'color',colors(2,:));
% plot(t_plot(1:end-5),dyn.freq(startInhale:end-5,3)./max(dyn.freq(BHstart:end,3)+0.5)+0.5,'-','LineWidth',2,'color',colors(1,:));
set(gca,'FontSize',font_axis)
ylabel('RBC:Barrier','FontSize',font_label);
xlabel('Time (s)')
%     ylim([0 plotlim(1,iComp)])
%     title(titles(iComp),'FontSize',20);

for idx = 1:length(BHs)/2  
    fill([BHs(2*idx-1) BHs(2*idx-1) BHs(2*idx) BHs(2*idx)],...
        [-500 500 500 -500],[.75 .75 .75],'LineStyle','none','FaceAlpha',.5)
end   
set(gca,'Layer','top')
ylim([0.2 1])

if save_fig_flag == 2
    options.Format = 'tiff';
    hgexport(gcf,[dyn_folder,'\RBC_to_Bar',fitType,figname,'.tif'],options) 
    amp = 0; nmrFit_ppm = 0;
end  

end


end 
