% Demo dynamics with peak identification
subjName = '008-004';
rawfilepath = all_twix_dyn{88};
[dynfilepath, ~, subjName] = locateDynfromRaw(rawfilepath);

BHs = [2 7]; 
save_figs = 1;
figName = '_peaks';

[amp, nmrFit_ppm] = dynamicSummary(rawfilepath,dynfilepath,BHs,save_figs,figName,'OscType','peaks');
dates = getDynDatesforPPT(rawfilepath,dynfilepath); 

dynfolder = fileparts(dynfilepath);
imgName = {'_peaks'};

dynamicSummaryPPT([subjName,'_peaks'],dynfolder,imgName,dates,amp,nmrFit_ppm)
