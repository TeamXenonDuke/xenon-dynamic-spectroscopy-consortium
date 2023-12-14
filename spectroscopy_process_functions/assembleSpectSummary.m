clear all_amp
clear all_fit
clear dates
filelist_twix_dyn;
% filelist_pfiles_dyn;
% all_twix_dyn = all_pfiles_dyn;

subj = "008-004/";  

% filelist_twix_dyn;
subjNums = find(contains(all_twix_dyn,subj));
[~, ~, subjName] = locateDynfromRaw(all_twix_dyn{subjNums(1)})
figname = {'','','_2','_2','_3','_3'};

for idx = 1:2:2*length(subjNums)
    k = subjNums((idx+1)/2);
    [dyn_loc, BHs] = locateDynfromRaw(all_twix_dyn{k})
    [dyn_locL, BHs_L] = locateDynfromRaw(all_twix_dyn{k},'L'); 
    [ampL, nmrFit_ppmL] = dynamicSummary(all_twix_dyn{k},dyn_locL,BHs_L,1,figname{idx});           
    all_amp(idx) = ampL; nmrFit_ppmL.fwhmG = [0 0 0]; all_fit(idx) = nmrFit_ppmL;
    [amp, nmrFit_ppm] = dynamicSummary(all_twix_dyn{k},dyn_loc,BHs,1,figname{idx});
    all_amp(idx+1) = amp; all_fit(idx+1) = nmrFit_ppm;
    dates(idx) = getDynDatesforPPT(all_twix_dyn{k},dyn_locL);
    dates(idx+1) = getDynDatesforPPT(all_twix_dyn{k},dyn_loc); 
    clear nmrFit_ppmL
end

dynamicSummaryPPT(subjName,fileparts(dyn_loc),figname,dates,all_amp,all_fit)
