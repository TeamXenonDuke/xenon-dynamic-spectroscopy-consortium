function assembleSpectSummaryVoigt(subjName, rawfilepath, dynfilepath)

% we currently do a 12s breath hold. Ignore the first 2s (100 frames)
% due to downstream magnetization
BHs = [2 7]; 
save_figs = 1;
figName = '';

[amp, nmrFit_ppm] = dynamicSummary(rawfilepath,dynfilepath,BHs,save_figs,figName);
dates = getDynDatesforPPT(rawfilepath,dynfilepath); 

dynfolder = fileparts(dynfilepath);
imgName = cell(1);

dynamicSummaryPPT(subjName,dynfolder,imgName,dates,amp,nmrFit_ppm)
