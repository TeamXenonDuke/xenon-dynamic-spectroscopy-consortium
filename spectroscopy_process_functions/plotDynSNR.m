function plotDynSNR(dyn, BHs, save_fig_flag, save_fig_path) 

if ~isfield(dyn,'snrs')
    dyn.snrs = zeros(1,length(dyn.area));
end 

figure(3), clf
set(gcf,'Position',[275 -100 450 200]), hold on
plot(dyn.t,dyn.snrs(1:end),'-','LineWidth',2,'color',[0         0.4470    0.7410]);
for idx = 1:length(BHs)/2  
    fill([BHs(2*idx-1) BHs(2*idx-1) BHs(2*idx) BHs(2*idx)],...
        [-500 500 500 -500],[.75 .75 .75],'LineStyle','none','FaceAlpha',.5)
end   
set(gca,'Layer','top')
ylim([10 30])
xlim([dyn.t(1), dyn.t(end)])
xlabel('Time (s)'), ylabel('SNR')

if save_fig_flag == 1
    options.Format = 'tiff';
    hgexport(gcf,save_fig_path,options)   
end