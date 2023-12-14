function plotDynamics(dyn, BHs, save_fig_flag, save_fig_path)
% Create Dynamic By Resonance Figure 

[~, nComp] = size(dyn.t);

[BHstart, BHend] = findBHs(dyn.t(:,1), BHs);
m3rd = BHstart:BHend;

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

figure(1); clf
set(gcf,'Position',[20 -200 1200 675])
axHandles = zeros(1,5*nComp);
t_plot = dyn.t(:,1);
startInhale = 1;

for iComp = 1:nComp
       
    % Area
    axHandles(5*(iComp-1)+1) = subplot(4,nComp,iComp); hold on;
    plot(t_plot(1:end-5),dyn.area(startInhale:end-5,iComp)/max(dyn.area(50:end,end)),'-','LineWidth',linewidth,'color',colors(iComp,:));
    set(gca,'Xtick','','FontSize',font_axis)
    set(gca,'Position',[leftSpace+(iComp-1)*plotWidth+(iComp-1)*wSpace bottomSpace+3*vSpace+3*plotHeight plotWidth plotHeight])
    if iComp == 1
        ylabel('Amplitude','FontSize',font_label);
    end
    ylim([0 max(dyn.area(startInhale:end,iComp))/max(dyn.area(50:end,end))])
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
    hgexport(gcf,save_fig_path,options)   
end
