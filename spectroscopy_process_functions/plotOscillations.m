function plotOscillations(dyn, BHs, detrend, fitted, rbc_axis_lim, save_fig_flag, save_fig_path)
% plotOscillations: plot RBC amplitude, chem shift, linewidth, and phase
%   oscillations
%
%   dyn: dynamic spectroscopy matlab structure or structure location
%   BHs: either a vector containing the time of inhale/exhale or the filepath
%       for a saved matlab vector with the BH information. Default
%       values are [2 10].
%   detrend: structure containing detrended oscillations
%   fitted: structure containing fits to detrended oscillations
%   rbc_axis_lim: sets axis limits on rbc oscillation plots to +-rbc_axis_lim
%       if zero, automatic limits are used
%   save_fig_flag: logical input for saving figures
%   save_fig_path: path to save figures to


[BHstart, BHend] = findBHs(dyn.t(:,1), BHs);

if length(fitted.freq) == length(dyn.t(BHstart:BHend))
    fitType = 'sine';
else 
    fitType = 'peaks';
end 

% Colors
colors = [0.8500    0.3250    0.0980 
          0.4660    0.6740    0.1880
          0         0.4470    0.7410];

switch fitType
case 'sine'
    figure(2), clf
    set(gcf, 'Position',[875    -200    550    725])
    axFontFs = 13;

    subplot(4,1,1), hold on
    plot(dyn.t(BHstart:BHend),detrend.area*100,'k--')
    plot(dyn.t(BHstart:BHend),fitted.area*100,'Linewidth',2,'color',colors(1,:)); 
    plot(dyn.t(BHstart:BHend,1),zeros(1,length(BHstart:BHend)),'Linewidth',.25,'Color',[.5 .5 .5])
    xlim([dyn.t(BHstart) dyn.t(BHend)])
    xlabel('Time (s)'), ylabel('Amplitude')
    set(gca,'FontSize',axFontFs)
    if rbc_axis_lim~=0
        if rbc_axis_lim>0
            ylim([-rbc_axis_lim rbc_axis_lim])
        else
            ylim([rbc_axis_lim -rbc_axis_lim])
        end
    end

    subplot(4,1,2), hold on
    plot(dyn.t(BHstart:BHend),detrend.freq,'k--')
    plot(dyn.t(BHstart:BHend),fitted.freq,'Linewidth',2,'color',colors(1,:)); 
    plot(dyn.t(BHstart:BHend,1),zeros(1,length(BHstart:BHend)),'Linewidth',.25,'Color',[.5 .5 .5])
    xlim([dyn.t(BHstart) dyn.t(BHend)])
    xlabel('Time (s)'), ylabel('Chemical Shift (ppm)')
    set(gca,'FontSize',axFontFs)
    ylim([-.5 .5])

    subplot(4,1,3), hold on
    plot(dyn.t(BHstart:BHend),detrend.fwhm,'k--')
    plot(dyn.t(BHstart:BHend),fitted.fwhm,'Linewidth',2,'color',colors(1,:)); 
    plot(dyn.t(BHstart:BHend,1),zeros(1,length(BHstart:BHend)),'Linewidth',.25,'Color',[.5 .5 .5])
    xlim([dyn.t(BHstart) dyn.t(BHend)])
    xlabel('Time (s)'), ylabel('FWHM (ppm)')
    set(gca,'FontSize',axFontFs)
    ylim([-.5 .5])

    subplot(4,1,4), hold on
    plot(dyn.t(BHstart:BHend),detrend.phase,'k--')
    plot(dyn.t(BHstart:BHend),fitted.phase,'Linewidth',2,'color',colors(1,:)); 
    plot(dyn.t(BHstart:BHend,1),zeros(1,length(BHstart:BHend)),'Linewidth',.25,'Color',[.5 .5 .5])
    xlim([dyn.t(BHstart) dyn.t(BHend)])
    xlabel('Time (s)'), ylabel('Phase (degrees)')
    set(gca,'FontSize',axFontFs)
    ylim([-10 10])
case 'peaks'
    figure(2), clf
    set(gcf, 'Position',[875    -200    550    725])
    axFontFs = 13;

    subplot(4,1,1), hold on
    plot(dyn.t(BHstart:BHend,1),detrend.area*100,'k--')
    plot(fitted.at,fitted.area*100,'.','color',[0.8500    0.3250    0.0980] ,'MarkerSize',30);
    plot(dyn.t(BHstart:BHend,1),zeros(1,length(BHstart:BHend)),'Linewidth',.25,'Color',[.5 .5 .5])
    plot(dyn.t(BHstart:BHend,1),fitted.area_sine*100,'color',[0, 0.4470, 0.7410],'Linewidth',1.5)
    plot(dyn.t(BHstart:BHend,1),fitted.area_max*100*ones(1,length(BHstart:BHend)),'--','Linewidth',.25,'Color',[0, 0.4470, 0.7410])
    plot(dyn.t(BHstart:BHend,1),fitted.area_min*100*ones(1,length(BHstart:BHend)),'--','Linewidth',.25,'Color',[0, 0.4470, 0.7410])
    xlim(BHs), ylabel('Amplitude (%)')
    if rbc_axis_lim~=0
        if rbc_axis_lim>0
            ylim([-rbc_axis_lim rbc_axis_lim])
        else
            ylim([rbc_axis_lim -rbc_axis_lim])
        end
    end
    
    subplot(4,1,2), hold on
    plot(dyn.t(BHstart:BHend,1),detrend.freq,'k--')
%     plot(dyn.t(BHstart:BHend,1),smooth(freq_detrend,20),'Linewidth',2)
    plot(fitted.ft,fitted.freq,'.','color',[0.8500    0.3250    0.0980] ,'MarkerSize',30);
    plot(dyn.t(BHstart:BHend,1),zeros(1,length(BHstart:BHend)),'Linewidth',.25,'Color',[.5 .5 .5])
    plot(dyn.t(BHstart:BHend,1),fitted.freq_max*ones(1,length(BHstart:BHend)),'--','Linewidth',.25,'Color',[0, 0.4470, 0.7410])
    plot(dyn.t(BHstart:BHend,1),fitted.freq_min*ones(1,length(BHstart:BHend)),'--','Linewidth',.25,'Color',[0, 0.4470, 0.7410])
    xlim(BHs), ylabel('Freq (ppm)')
    
    subplot(4,1,3), hold on
    plot(dyn.t(BHstart:BHend,1),detrend.fwhm,'k--')
%     plot(dyn.t(BHstart:BHend,1),smooth(fwhm_detrend,20),'Linewidth',2)
    plot(fitted.fwt,fitted.fwhm,'.','color',[0.8500    0.3250    0.0980] ,'MarkerSize',30);
    plot(dyn.t(BHstart:BHend,1),zeros(1,length(BHstart:BHend)),'Linewidth',.25,'Color',[.5 .5 .5])
    plot(dyn.t(BHstart:BHend,1),fitted.fwhm_max*ones(1,length(BHstart:BHend)),'--','Linewidth',.25,'Color',[0, 0.4470, 0.7410])
    plot(dyn.t(BHstart:BHend,1),fitted.fwhm_min*ones(1,length(BHstart:BHend)),'--','Linewidth',.25,'Color',[0, 0.4470, 0.7410])
    xlim(BHs), ylabel('FWHM (ppm)')
    
    subplot(4,1,4), hold on
    plot(dyn.t(BHstart:BHend,1),detrend.phase,'k--')
%     plot(dyn.t(BHstart:BHend,1),smooth(phase_detrend,20),'Linewidth',2)
    plot(fitted.pt,fitted.phase,'.','color',[0.8500    0.3250    0.0980] ,'MarkerSize',30);   
    plot(dyn.t(BHstart:BHend,1),zeros(1,length(BHstart:BHend)),'Linewidth',.25,'Color',[.5 .5 .5])
    plot(dyn.t(BHstart:BHend,1),fitted.phase_max*ones(1,length(BHstart:BHend)),'--','Linewidth',.25,'Color',[0, 0.4470, 0.7410])
    plot(dyn.t(BHstart:BHend,1),fitted.phase_min*ones(1,length(BHstart:BHend)),'--','Linewidth',.25,'Color',[0, 0.4470, 0.7410])
    xlim(BHs), xlabel('Time (s)'), ylabel('Phase ({\circ})')
    
    set(gca,'FontSize',axFontFs)
    ylim([-10 10])    
end 

if save_fig_flag == 1
    options.Format = 'tiff';
    hgexport(gcf,save_fig_path,options)
end  
