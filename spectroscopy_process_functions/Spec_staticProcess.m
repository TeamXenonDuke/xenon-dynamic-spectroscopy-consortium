%% Static Spectroscopy Processing Using the Functions
%% Static Spectroscopy Processing 
% - readRawDyn
% - calculate Static Spectroscopy
% - create a table.
% - prepare the report

%%
clear all

disp("Select Dixon File ...")
[file, path] = uigetfile('*.*', 'Select file');  % starts in current dir 
file_with_path = strcat(path, file);  % join path and filename to open

% Twix file from Siemens
twix = mapVBVD(file_with_path);

fids = double(squeeze(twix.image.unsorted()));
dwell_time = twix.hdr.Config.DwellTime*10^-9;  

xeGam = twix.hdr.Dicom.lFrequency*10e-7; 
tr = twix.hdr.Config.TR(1)*1E-6;  
npts = twix.hdr.Config.RawCol*2;
nAvg = ceil(tr);

%[fids, dwell_time, npts, tr, xeGam] = readRawDyn(raw_path);
nFrames = size(fids,2);             % Number of Dis. Frames
t_tr = tr*(1:nFrames);

BHs = [2, 10];
[BHstart, ~] = findBHs(t_tr, BHs);

% Calculate Static Spectral Parameters
t = dwell_time*(0:(npts-1))';

% Three peak fit
area_orig = [1 1 1];
freq_orig = [0 -20.7 -218.4]*xeGam;
% freq_orig = [0 -12 -218]*xeGam+7500; For 90 degree spoiled
fwhm_orig = [7 6 2]*xeGam;
fwhmL_orig = [8.8 5.0 2]*xeGam;
fwhmG_orig = [0 6.1 0]*xeGam;
phase_orig = [0 0 0];
peaks = 3; ngas = 3; nbar = 2;

avgdata = mean(fids(:,BHstart:BHstart+nAvg),2);

nmrFit = NMR_TimeFit_v(avgdata, t, area_orig,freq_orig,fwhmL_orig,fwhmG_orig,phase_orig,[],[]);
nmrFit = nmrFit.fitTimeDomainSignal();
ref_freq = nmrFit.freq(ngas);

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

positive_phase = nmrFit.phase - nmrFit.phase(2);
positive_phase(positive_phase < 0) = positive_phase(positive_phase < 0) + 360;

nmrFit_ppm.area = nmrFit.area/sum(nmrFit.area(nbar));
nmrFit_ppm.freq = (nmrFit.freq-ref_freq)/xeGam;
nmrFit_ppm.fwhm = nmrFit.fwhm/xeGam;
nmrFit_ppm.fwhmG = nmrFit.fwhmG/xeGam;
nmrFit_ppm.phase = positive_phase;
nmrFit_ppm.SNRresid = SNRresid;
nmrFit_ppm.SNRsnf = SNRsnf;

disp([10 'Static Parameters:'])
disp('    Area     Freq      FWHM      FWHMg     Phase     SNRresid   SNRsnf');
for k = 1:peaks
    disp([sprintf('%8.2f', nmrFit.area(k)/sum(nmrFit.area(nbar))), ' ' ...
        sprintf('%8.1f',(nmrFit.freq(k)-ref_freq)/xeGam),  '  ' ...
        sprintf('%8.1f',nmrFit.fwhm(k)/xeGam), '  ' ...
        sprintf('%8.1f',nmrFit.fwhmG(k)/xeGam), '  ' ...
        sprintf('%9.1f',positive_phase(k)), '  '...
        sprintf('%9.1f',SNRresid(k)), '  '...
        sprintf('%9.1f',SNRsnf(k))]);
end


%%
import mlreportgen.ppt.*;

slidesFile = ' Spectroscopy Summary.pptx';
slides = Presentation(slidesFile,[extractBefore(which('dynamicSummaryPPT'),'dynamicSummaryPPT'),'spectSummary.pptx']);

%%
function slides = dynamicSummaryPPTall(slides,subj,folder,imgName,dates,amp,nmrFit_ppm)

import mlreportgen.ppt.*;

for slideNum = 1:length(amp)
    % Add a summary slide
    summarySlide = add(slides,'SpectralSummary');
    replace(summarySlide,'Static Text','Static Spectroscopy');
    replace(summarySlide,'StaticValues Text','Static Spectroscopy (Barrier Voigt)');
    replace(summarySlide,'norm Text','*normalized to barrier peak');
    replace(summarySlide,'box1',' ');
    replace(summarySlide,'box2',' ');
    replace(summarySlide,'box3',' ');
    replace(summarySlide,'scanDate Text','Scan Date:');
    
    % Add Static table
    static_params = {''; ['Intensity', 10, 'Ratio*']; ['Shift', 10, '(ppm)'];...
        ['FWHM', 10, '(ppm)']; ['FWHMg', 10, '(ppm)'];...
        ['Phase', 10, '(degrees)']; ['5 FID SNR', 10, '(amp/noise)']};

    if sum(nmrFit_ppm(slideNum).fwhmG(:)) > 0
        fitType = 'V';
        replace(summarySlide,'Title 1',['Subject ',subj,' (Voigt)']);
            rbc_static = {''; sprintf('%3.2f', nmrFit_ppm(slideNum).area(1));...
            sprintf('%4.1f',nmrFit_ppm(slideNum).freq(1));
            sprintf('%2.1f',nmrFit_ppm(slideNum).fwhm(1)); '-----';
            sprintf('%3.1f',nmrFit_ppm(slideNum).phase(1)); ...
            sprintf('%3.1f',nmrFit_ppm(slideNum).SNR_dis(1))};

        rbc_ref = {'Ref. Values'; ['(0.59 ',177,' 0.12)'];
            ['(218.4 ',177,' 0.4)']; ['(8.7 ',177,' 0.3)'];...
            '-----'; ['(81.9 ',177,' 3.6)']; ''};

        bar_static = {''; '1.00';
            sprintf('%4.1f',nmrFit_ppm(slideNum).freq(2)); 
            sprintf('%2.1f',nmrFit_ppm(slideNum).fwhm(2)); 
            sprintf('%3.2f',nmrFit_ppm(slideNum).fwhmG(2));
            '0.0'; 
            sprintf('%3.1f',nmrFit_ppm(slideNum).SNR_dis(2))};

        bar_ref = {'Ref. Values'; '(1.0)';...
            ['(197.7 ',177,' 0.3)']; ['(5.0 ',177,' 0.3)'];...
            ['(6.1 ',177,' 0.3)']; '(0.0)'; ''};
    else 
        
    end

    staticTable = Table([static_params, rbc_static, rbc_ref, bar_static, bar_ref]);
    colSpecsAmp(1) = ColSpec('1in'); 
    colSpecsAmp(2) = ColSpec('0.75in');
    colSpecsAmp(3) = ColSpec('1in');
    colSpecsAmp(4) = ColSpec('0.7in');
    colSpecsAmp(5) = ColSpec('1in');

    staticTable.ColSpecs = colSpecsAmp;
    staticTable.Height = '2.5in';

    staticTable.ColSpecs(1).HAlign = 'Center'; staticTable.ColSpecs(1).VAlign = 'Middle';
    staticTable.ColSpecs(1).FontSize = '11pt'; 
    staticTable.ColSpecs(2).HAlign = 'Center'; staticTable.ColSpecs(2).VAlign = 'Middle';
    staticTable.ColSpecs(2).FontSize = '12pt';
    staticTable.ColSpecs(3).HAlign = 'Center'; staticTable.ColSpecs(3).VAlign = 'Middle';
    staticTable.ColSpecs(3).FontSize = '11pt'; staticTable.ColSpecs(3).FontColor = '#767171';
    staticTable.ColSpecs(4).HAlign = 'Center'; staticTable.ColSpecs(4).VAlign = 'Middle';
    staticTable.ColSpecs(4).FontSize = '12pt';
    staticTable.ColSpecs(5).HAlign = 'Center'; staticTable.ColSpecs(5).VAlign = 'Middle';
    staticTable.ColSpecs(5).FontSize = '11pt'; staticTable.ColSpecs(5).FontColor = '#767171';

    staticTable.row(1).Style = {VAlign('bottom')}; 
    staticTable.row(1).FontSize = '12pt';
    staticTable.row(7).Style = {VAlign('Top')}; 

    staticTable.StyleName = 'No Style, No Grid';
    staticTable.X = '8.79in'; hrTable.Y = '4.32in';
    replace(summarySlide,'static table',staticTable);

    p1 = mlreportgen.ppt.Picture([folder,'\dyn',fitType,imgName{slideNum},'.tif']);
    replace(summarySlide,'dynPicture',p1);
    p2 = mlreportgen.ppt.Picture([folder,'\dyn',fitType,'_oscs_new',imgName{slideNum},'.tif']);
    replace(summarySlide,'oscillationPicture',p2);
    p3 = mlreportgen.ppt.Picture([folder,'\static',fitType,imgName{slideNum},'.tif']);
    replace(summarySlide,'staticPicture',p3);
    p4 = mlreportgen.ppt.Picture([folder,'\dyn',fitType,'_error',imgName{slideNum},'.tif']);
    replace(summarySlide,'snrPicture',p4);
    
end 

    % Add RBC_to_barrier slide
%     rbc_to_barSlide = add(slides,'rbc_to_bar');
%     replace(rbc_to_barSlide,'Title 1',['Subject ',subj,' (Voigt)']);
%     replace(rbc_to_barSlide,'Osc Text','RBC:Barrier');
%     p5 = mlreportgen.ppt.Picture([folder,'\rbc_to_bar',fitType,imgName{slideNum},'.tif']);
%     replace(rbc_to_barSlide,'oscillationPicture',p5);
end 

%%
%plotStaticSpectroscopy(file_with_path, nmrFit)

%%
function plotStaticSpectroscopy(raw_path, nmrFit)

%save_fig_flag, save_fig_path
%[~, dwell_time, ~, ~, xeGam] = readRawDyn(raw_path);

peaks = 3;
nbar = 2; ngas = 3;
ref_freq = nmrFit.freq(end);

% Constants
xe_shift = (0.3/5.5)*(273/298)*(0.548);
tot_shift = -xe_shift;
tot_shift_hz = tot_shift*xeGam;

%% Calculations

% Calculate fitted and residual spectrums using time domain
% signal so that amplitudes are correct even if signal is
% truncated at the end

% Calculate spectrum
fittedFID = nmrFit.calcTimeDomainSignal(nmrFit.t);
fittedSpectrum = dwell_time*fftshift(fft(fittedFID));
residualSpectrum = nmrFit.spectralDomainSignal - fittedSpectrum;

f = linspace(-0.5,0.5,length(nmrFit.t)*2+1)/dwell_time;
f = f(1:(end-1)); % Take off last sample to have nSamples
ppm = (f+tot_shift_hz-ref_freq)/xeGam;
ppm_shift = (nmrFit.freq(ngas)-ref_freq)/xeGam; % shift reference spectrum to match calculated

tplot = linspace(nmrFit.t(1),2*nmrFit.t(end),length(f));
fittedTime = nmrFit.calcComponentTimeDomainSignal(tplot);
fittedSpectrumComp = dwell_time*fftshift(fft(fittedTime));

%%
if isprop(nmrFit,'fwhmG')
    % Referance values from healthy population for voigt
    area_ref = [0.59    1.0    0.17];
    freq_ref = [218.4  197.7   0];
    fwhm_ref = [8.72    4.97    1.45];
    fwhmG_ref = [0    6.1    0];
    phase_ref = [81.9 0  248.1]+nmrFit.phase(2);

    nmrFit_ref = NMR_TimeFit_v(ones(length(nmrFit.t),1), nmrFit.t,...
        area_ref,freq_ref*xeGam+nmrFit.freq(ngas),fwhm_ref*xeGam,...
        fwhmG_ref*xeGam, phase_ref,[],[]);
    
    fitted_ref = nmrFit_ref.calcComponentTimeDomainSignal(tplot);
    spectrum_ref = fftshift(fft(fitted_ref)); 
else 
    % Referance values from healthy population for Lorentzian 
    area_ref = [0.3598    0.6131    0.0565];
    freq_ref = [216.0468  197.3733   -0.3411];
    fwhm_ref = [9.9863    8.6263    2.0852];
    phase_ref = [-48.2896 -140.9711  -66.7937];

    nmrFit_ref = NMR_TimeFit(ones(length(nmrFit.t),1), nmrFit.t,...
        area_ref,freq_ref*xeGam+nmrFit.freq(ngas),fwhm_ref*xeGam,...
        phase_ref,[],[]);
    
    spectrum_ref = calcComponentLorentzianCurves(nmrFit_ref,f);
end 

%% Plotting

% Colors
colors = [0.8500    0.3250    0.0980 
          0.4660    0.6740    0.1880
          0         0.4470    0.7410];
      
figure(4), clf
ax2 = subplot(3,1,1); hold on
set(gcf,'Position',[1275 -100 600 600]); 
plot(ppm-ppm_shift, abs(spectrum_ref)/sum(max(abs(spectrum_ref(:,3)))),':k','Linewidth',2)

plotorder = [2 3 1];
for i = 1:peaks
    plot(ppm,abs(fittedSpectrumComp(:,plotorder(i)))/abs(max(fittedSpectrumComp(:,3))),'Color',colors(i,:),'Linewidth',3)  
end 
xticks(150:20:250), set(gca,'XDir','reverse'), box on
ylabel('Component Intensity'), xlabel('Chemical Shift (ppm)');
ylim([0 1]), xlim([150 250])

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
    hgexport(gcf,save_fig_path,options)   
end
end

