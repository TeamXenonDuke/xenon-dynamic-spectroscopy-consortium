% Processing GE long format for healthy test sujects 
filelist_pfiles_dyn;

import mlreportgen.ppt.*;
slidesFile = ['D:\OneDrive\Documents\Duke\CIVM Research\Dynamic Spectroscopy\GE\','IPF Spectroscopy Summary 3.pptx'];
slides = Presentation(slidesFile,[extractBefore(which('dynamicSummaryPPT'),'dynamicSummaryPPT'),'spectSummary.pptx']);

% Add a title slide
presentationTitleSlide = add(slides,'Title Slide');
replace(presentationTitleSlide,'Title 1','Spectroscopy Summaries');
replace(presentationTitleSlide,'Subtitle 2',datestr(now,'mm/dd/yyyy'));

subj = {"003\ScanB","006\ScanA","39\ScanB","41\ScanA","42\ScanA",...
    "42\ScanB","49\ScanA","49\ScanD","055","064","065","003_001","93\s1","93\s2",...
    "93\s3","93\s4","93\s5","93A","94\s1","94\s2","94\s3","94\s4",...
    "94\s5","94A","95\s1","95\s2","95\s3","95\s4","95\s5"}

subj = {"046\ScanB","067\ScanA","068\ScanA","069\ScanB","080",...
    "073\ScanA","076\ScanA"}

for k = 1:length(subj)
    subjNums = find(contains(all_pfiles_dyn,subj{k}));
    [~, ~, subjName] = locateDynfromRaw(all_pfiles_dyn{subjNums})
    figname = {'','','',''};

    for idx = 1:length(subjNums)
        k = subjNums(idx);
        [dyn_loc, BH_loc] = locateDynfromRaw(all_pfiles_dyn{k});  
        [amp, nmrFit_ppm] = dynamicRBCtoBarrier(all_pfiles_dyn{k},dyn_loc,BH_loc,1,figname{idx},1);
%         [amp, nmrFit_ppm] = dynamicSummary(all_pfiles_dyn{k},dyn_loc,BH_loc,1,figname{idx});
        all_amp(idx) = amp; all_fit(idx) = nmrFit_ppm;
        dates(idx) = getDynDatesforPPT(all_pfiles_dyn{k},dyn_locL);
        clear nmrFit_ppmL
    end

    slides = dynamicSummaryPPTall(slides, subjName,fileparts(dyn_loc),figname,dates,all_amp,all_fit)
    clear all_amp
    clear all_fit
end 

% Generate and open the presentation
close(slides);

if ispc
    winopen(slidesFile);
end
