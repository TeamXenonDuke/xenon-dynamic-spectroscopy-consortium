function plotStaticSpectroscopy(raw_path, nmrFit, save_fig_flag, save_fig_path)

[~, dwell_time, ~, ~, xeGam] = readRawDyn(raw_path);

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

if isprop(nmrFit,'fwhmG')
    % Reference values from healthy population for voigt
    % BD 7/25/22 - these values don't affect report - just fit guesses?
    area_ref = [0.60    1.0    0.17]; % Ari Hgb corrected TR=15ms
    freq_ref = [218.2  197.7   0]; % Bas curated population
    fwhm_ref = [8.72    4.97    1.45];
    fwhmG_ref = [0    6.1    0];
    phase_ref = [81.9 0  248.1]+nmrFit.phase(2);

    nmrFit_ref = NMR_TimeFit_v(ones(length(nmrFit.t),1), nmrFit.t,...
        area_ref,freq_ref*xeGam+nmrFit.freq(ngas),fwhm_ref*xeGam,...
        fwhmG_ref*xeGam, phase_ref,[],[]);
    
    fitted_ref = nmrFit_ref.calcComponentTimeDomainSignal(tplot);
    spectrum_ref = fftshift(fft(fitted_ref)); 
else 
    % Reference values from healthy population for Lorentzian 
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
plot(ppm-ppm_shift,abs(spectrum_ref)/sum(max(abs(spectrum_ref(:,3)))),':k','Linewidth',2)

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