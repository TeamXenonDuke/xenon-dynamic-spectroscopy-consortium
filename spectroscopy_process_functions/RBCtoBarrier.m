rawData = pfiles_all{find(contains(pfiles_all,"93\s1"))}; [dyn_loc, BH_loc] = locateDynfromRaw(rawData);

font_label = 12;
font_axis = 12;

load(dyn_loc);
if isempty(BH_loc)
    BHs = [2 10];
elseif ischar(BH_loc)
    load(BH_loc);
elseif isvector(BH_loc)
    BH_loc = BHs;
end 

% BHs = [12 17]; 
BHstart = find(round(dyn.t(:,1),2) == BHs(1)); BHend = find(round(dyn.t(:,1),2) == BHs(2));


figure(5), clf
set(gcf,'Position',[247         347        1180         556])
hold on

% Area
plot(dyn.t(1:end-5,1),dyn.area(1:end-5,1)./dyn.area(1:end-5,2),'-','LineWidth',3,'color',[0         0.4470    0.7410]);
plot(dyn.t(1:end-5,1),dyn.fwhmL(1:end-5,3)./max(dyn.fwhmL(BHstart:end,3)),'-','LineWidth',2,'color',[0.4660    0.6740    0.1880]);
plot(dyn.t(1:end-5,1),dyn.freq(1:end-5,3)./max(dyn.freq(BHstart:end,3)+0.5)+0.5,'-','LineWidth',2,'color',[0.8500    0.3250    0.0980]);
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

%%
b = highpassfilter(length(dyn.area(:,1)));
[amp, detrend, fitted] = getOscillationInfo(dyn,BHs,b);
cardiacL = round(amp.hr/60*50);

figure(6), clf
set(gcf,'Position',[247         347        1180         556])
hold on

% Area
rbcBar = dyn.area(1:end-5,1)./dyn.area(1:end-5,2);
rbcBar_smooth = smooth(rbcBar,cardiacL);
plot(dyn.t(1:end-5,1),rbcBar_smooth,'-','LineWidth',3,'color',[0         0.4470    0.7410]);
% plot(dyn.t(1:end-5,1),dyn.fwhmL(1:end-5,3)./max(dyn.fwhmL(BHstart:end,3)),'-','LineWidth',2,'color',[0.4660    0.6740    0.1880]);
% plot(dyn.t(1:end-5,1),dyn.freq(1:end-5,3)./max(dyn.freq(BHstart:end,3)+0.5)+1,'-','LineWidth',2,'color',[0.8500    0.3250    0.0980]);
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

rbc2bar_BH = mean(rbcBar(BHstart:BHend));
rbc2bar_exhale = max(rbcBar(BHend:end));
