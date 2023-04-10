% Create full summary 
clear all;
filelist_twix_dyn;

import mlreportgen.ppt.*;
subjectType = 'Repeatability';
powerpoint_file_name = 'repeatability_summaries_peaks_v2_20220216.pptx'; 
powerpoint_save_loc = '/Users/eab78/Library/CloudStorage/OneDrive-Personal/Documents/Duke/CIVM Research/Repeatability/';
slidesFile = [powerpoint_save_loc,powerpoint_file_name]; % name and location of output summary
slides = Presentation(slidesFile,[extractBefore(which('dynamicSummaryPPT'),'dynamicSummaryPPT'),'spectSummary.pptx']); % location of Spect Summary Template

% Add a title slide
presentationTitleSlide = add(slides,'Title Slide');
replace(presentationTitleSlide,'Title 1',['Spectroscopy Summaries ',subjectType]);
replace(presentationTitleSlide,'Subtitle 2',datestr(now,'mm/dd/yyyy'));

filelist_pairs;
processNum = 1;

for k = pairs_index_rand  
    [dyn_loc, BH_loc, subjName] = locateDynfromRaw(all_twix_dyn{k});
    
    save_path = [loc,'Processed Data Peaks'];
    load([save_path,'/spect_',subjName]);
    dates = getDynDatesforPPT(all_twix_dyn{k},dyn_loc);
        
    if contains(subjName,"s2") || contains(subjName,"_2")
        figname = {'_2'};
    else 
        figname = {''};
    end 
   
    slides = dynamicSummaryPPTall_peaks(slides, subjName,fileparts(dyn_loc),figname,dates,spect.amp,spect.nmrFit);
    disp(['Added:  ',subjName, '   (num = ',num2str(processNum),')'])
    processNum = processNum + 1;
end 


% Generate and open the presentation
close(slides);

if ispc
    winopen(slidesFile);
end
