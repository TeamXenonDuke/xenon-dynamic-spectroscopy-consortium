% Create Bland-Altman Plots

spectData = readtable('D:\OneDrive\Documents\Duke\CIVM Research\Siemens Data\Processed Data\spectroscopy_stats_master.xlsx','Sheet','Pairs');
% spectData = readtable('D:\OneDrive\Documents\Duke\CIVM Research\Siemens Data\Processed Data Peaks\all_stats_peaks19-07-26_removed.xlsx','Sheet','Pairs');
% spectData = readtable('D:\OneDrive\Documents\Duke\CIVM Research\Siemens Data\Processed Data GOF\all_stats19-10-21.xlsx');
% Amp_Osc_mean = (spectData{1:2:end,5} + spectData{:,24})/2;
% Amp_Osc_diff = (spectData{:,5} - spectData{:,24});
% 
% 
% Amp_Osc_mean = (spectData{:,5} + spectData{:,24})/2;
% Amp_Osc_diff = (spectData{:,5} - spectData{:,24});

labels = zeros(size(spectData,1),1);
labels(strcmp(spectData.Label, 'Healthy')) = 1;
labels(strcmp(spectData.Label, 'IPF')) = 2;
labels(strcmp(spectData.Label, 'NSIP')) = 3;
labels(strcmp(spectData.Label, 'Alpha 1')) = 4;
labels(strcmp(spectData.Label, 'RT')) = 5;
labels(strcmp(spectData.Label, 'COPD')) = 6;
labels(labels == 0) = 7;

figure(3); clf
subplot(2,2,1), box on
blandAltman(spectData.Amp_Osc,spectData.Amp_Osc2, labels)
title('RBC Amplitude Oscillations')
subplot(2,2,2), box on
blandAltman(spectData.Shift_Osc,spectData.Shift_Osc2, labels)
title('RBC Chemical Shift Oscillations')
legend('Healthy','IPF','NSIP','Alpha1','RT','COPD','Other')
subplot(2,2,3), box on
blandAltman(spectData.FWHM_Osc,spectData.FWHM_Osc2, labels)
title('RBC FWHM Oscillations')
subplot(2,2,4), box on
blandAltman(spectData.Phase_Osc,spectData.Phase_Osc2, labels)
title('RBC Phase Oscillations')


figure(6); clf
subplot(2,2,1), box on
blandAltman(spectData.RBC_Bar_Ratio,spectData.RBC_Bar_Ratio2, labels)
title('RBC:Barrier')
subplot(2,2,2), box on
blandAltman(spectData.RBC_Shift,spectData.RBC_Shift2, labels)
title('RBC Chemical Shift')
legend('Healthy','IPF','NSIP','Alpha1','RT','COPD','Other')
subplot(2,2,3), box on
blandAltman(spectData.RBC_FWHM,spectData.RBC_FWHM2, labels)
title('RBC FWHM')
subplot(2,2,4), box on
blandAltman(spectData.RBC_Phase,spectData.RBC_Phase2, labels)
title('RBC Phase')

figure(7); clf
subplot(2,2,2), box on
blandAltman(spectData.Bar_Shift,spectData.Bar_Shift2, labels)
title('Barrier Chemical Shift')
legend('Healthy','IPF','NSIP','Alpha1','RT','COPD','Other')
subplot(2,2,3), box on
blandAltman(spectData.Bar_FWHM_L,spectData.Bar_FWHM_L2, labels)
title('Barrier FWHM L')
subplot(2,2,4), box on
blandAltman(spectData.Bar_FWHM_G,spectData.Bar_FWHM_G2, labels)
title('Barrier FWHM G')