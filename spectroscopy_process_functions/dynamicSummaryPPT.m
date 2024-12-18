function dynamicSummaryPPT(subject_id, path, imgName, dates, amp_all, nmrFit_ppm, rf_excitation)
% TODO add docstring
% 7/23/22 - BD added updated peaks-specific healthy reference values

import mlreportgen.ppt.*;
ppt_file = ['Subject ', subject_id, '_', num2str(round(rf_excitation)), 'ppm Spectroscopy Summary.pptx'];
slidesFile = fullfile(path, ppt_file);
slides = Presentation(slidesFile, [extractBefore(which('dynamicSummaryPPT'), 'dynamicSummaryPPT'), 'spectSummary.pptx']);

% Add a title slide
presentationTitleSlide = add(slides, 'Title Slide');
replace(presentationTitleSlide, 'Title 1', ['Subject ', subject_id, '_', num2str(rf_excitation), 'ppm', 10, 'Spectroscopy Summaries']);
replace(presentationTitleSlide, 'Subtitle 2', datestr(now, 'mm/dd/yyyy'));

scan_strV = 1;
for slideNum = 1:length(imgName)

    % This if else block will make sure, sine fitted report in slide 1 and
    % peaks fitted figures and number are on slide 2
    if slideNum == 1
        amp = amp_all.sine;
        oscs_ref = {'Healthy Ref. Values'; ['(9.4 ', 177, ' 2.7%)']; ...
            ['(0.06 ', 177, ' 0.05)']; ['(0.19 ', 177, ' 0.09)']; ...
            ['(1.3 ', 177, ' 0.7', 176, ')']};
        RBCSNR = nmrFit_ppm.SNRresid(1);
        BarrSNR = nmrFit_ppm.SNRresid(2);
    elseif slideNum == 2
        amp = amp_all.peaks;
        oscs_ref = {'Healthy Ref. Values'; ['(10.3 ', 177, ' 1.5%)']; ...
            ['(0.10 ', 177, ' 0.03)']; ['(0.15 ', 177, ' 0.07)']; ...
            ['(1.9 ', 177, ' 0.4', 176, ')']};
        RBCSNR = nmrFit_ppm.SNRsnf(1);
        BarrSNR = nmrFit_ppm.SNRsnf(2);
    else
        sprintf('Out of slide Number range, should be 2 for sine and peaks')
    end

    % Add a summary slide
    summarySlide = add(slides, 'SpectralSummary');

    % Edit text boxes to fit the word "Membrane" without needing two lines
    barrierTextBox = find(summarySlide,'Barrier');
    barrier2TextBox = find(summarySlide,'Barrier2');
    rbc2TextBox = find(summarySlide,'RBC2');
    set(barrierTextBox,'Width','1.4in','X','10.5in');
    set(barrier2TextBox,'FontSize','11pt','VAlign','middle','Y','1.15in');
    set(rbc2TextBox,'FontSize','11pt','VAlign','middle','Y','1.15in');

    % Add all text elements
    replace(summarySlide, 'Dyn Text', 'Dynamics by Resonance');
    replace(summarySlide, 'RBC', 'RBC');
    replace(summarySlide, 'Barrier', 'Membrane');
    replace(summarySlide, 'Gas', 'Gas');
    replace(summarySlide, 'Osc Text', 'Detrended RBC Oscillations');
    replace(summarySlide, 'Amp Text', 'RBC Oscillation Analysis');
    replace(summarySlide, 'Pkpk Text', '*Peak-to-Peak');
    replace(summarySlide, 'Static Text', 'Static Spectroscopy');
    replace(summarySlide, 'StaticValues Text', 'Static Spectroscopy');
    replace(summarySlide, 'dynSNR Text', 'Dynamic SNR');
    replace(summarySlide, 'norm Text', '*normalized to membrane peak');
    replace(summarySlide, 'box1', ' ');
    replace(summarySlide, 'box2', ' ');
    replace(summarySlide, 'box3', ' ');
    replace(summarySlide, 'RBC2', 'RBC');
    replace(summarySlide, 'Barrier2', 'Membrane');
    replace(summarySlide, 'scanDate Text', 'Scan Date:');
    replace(summarySlide, 'processingDate Text', 'Dynamic Fit Date:');
    replace(summarySlide, 'ampDate Text', 'RBC Fit Date:');

    % Add dates
    replace(summarySlide, 'scanDate', dates.scan);
    replace(summarySlide, 'processingDate', dates.dyn_fit);
    replace(summarySlide, 'ampDate', dates.rbc_fit);

    % Add hr table
    hrTable = Table({'Heart Rate:', sprintf('%2.0f', amp.hr); 'Mean SNR:', sprintf('%2.1f', amp.snr)});
    colSpecs(1) = ColSpec('1.41in');
    colSpecs(2) = ColSpec('0.6in');
    hrTable.ColSpecs = colSpecs;

    hrTable.ColSpecs(1).Bold = 1;
    hrTable.ColSpecs(1).HAlign = 'Right';
    hrTable.ColSpecs(2).HAlign = 'Left';
    hrTable.FontSize = '14pt';
    hrTable.StyleName = 'No Style, No Grid';
    hrTable.Border = 'none';
    hrTable.ColSep = 'none';
    hrTable.RowSep = 'none';
    hrTable.X = '6.55in';
    hrTable.Y = '6.41in';

    if amp.snr <= 20
        hrTable.row(2).Style = {BackgroundColor('yellow')};
    end

    replace(summarySlide, 'HR_SNR', hrTable);

    % Add Oscillation table
    oscs_params = {''; 'Amplitude (%):'; ['Chemical Shift', 10, '(ppm):']; ...
        ['Linewidth', 10, '(ppm):']; 'Phase (degrees):'};
    oscs = {''; sprintf('%3.1f', amp.area*100); sprintf('%3.2f', amp.freq); ...
        sprintf('%3.2f', amp.fwhm); sprintf('%3.1f', amp.phase)};
    % temp - moved oscs_ref to be dependent on sine vs peaks


    ampTable = Table([oscs_params, oscs, oscs_ref]);
    colSpecsAmp(1) = ColSpec('1.6in');
    colSpecsAmp(2) = ColSpec('1in');
    colSpecsAmp(3) = ColSpec('1.43in');
    ampTable.ColSpecs = colSpecsAmp;
    ampTable.Height = '2.2in';

    ampTable.ColSpecs(1).HAlign = 'Right';
    ampTable.ColSpecs(1).VAlign = 'Middle';
    ampTable.ColSpecs(1).FontSize = '14pt';
    ampTable.ColSpecs(2).HAlign = 'Center';
    ampTable.ColSpecs(2).VAlign = 'Middle';
    ampTable.ColSpecs(2).FontSize = '16pt';
    ampTable.ColSpecs(3).HAlign = 'Left';
    ampTable.ColSpecs(3).VAlign = 'Middle';
    ampTable.ColSpecs(3).FontSize = '16pt';
    ampTable.ColSpecs(3).FontColor = '#767171';
    ampTable.row(1).Style = {VAlign('bottom')};
    ampTable.row(1).FontSize = '11pt';
    ampTable.row(5).Style = {VAlign('Middle')};

    ampTable.StyleName = 'No Style, No Grid';
    ampTable.Border = 'none';
    ampTable.ColSep = 'none';
    ampTable.RowSep = 'none';
    ampTable.X = '9.05in';
    hrTable.Y = '1.0in';

    replace(summarySlide, 'amp table', ampTable);

    % Add Static table
    static_params = {''; ['Intensity', 10, 'Ratio*']; ['Shift', 10, '(ppm)']; ...
        ['FWHM', 10, '(ppm)']; ['FWHMg', 10, '(ppm)']; ...
        ['Phase', 10, '(degrees)']; ['1 sec. SNR']};

    if sum(nmrFit_ppm.fwhmG(:)) > 0
        fitType = 'V';
        %scan_str = (slideNum+1)/2;
        if slideNum == 1
            fit_type = ' Sine';
        else
            fit_type = ' Peaks';
        end
        replace(summarySlide, 'Title 1', ['Subject ', subject_id, '_', num2str(rf_excitation), 'ppm (Voigt/', fit_type, ')']);
        rbc_static = {''; sprintf('%3.3f', nmrFit_ppm.area(1)); ...
            sprintf('%4.1f', nmrFit_ppm.freq(1)); ...
            sprintf('%2.1f', nmrFit_ppm.fwhm(1)); '-----'; ...
            sprintf('%3.1f', nmrFit_ppm.phase(1)); ...
            sprintf('%3.1f', RBCSNR)};

        % Ari Hgb corrected RBC/M and Bas curated RBC shift
        rbc_ref = {'Ref. Values'; ['(0.60 ', 177, ' 0.08)']; ...
            ['(218.2 ', 177, ' 0.4)']; ['(8.7 ', 177, ' 0.3)']; ...
            '-----'; ['(81.9 ', 177, ' 3.6)']; ''};

        bar_static = {''; '1.00'; ...
            sprintf('%4.1f', nmrFit_ppm.freq(2)); ...
            sprintf('%2.1f', nmrFit_ppm.fwhm(2)); ...
            sprintf('%3.2f', nmrFit_ppm.fwhmG(2)); ...
            '0.0'; ...
            sprintf('%3.1f', BarrSNR)};

        bar_ref = {'Ref. Values'; '(1.0)'; ...
            ['(197.7 ', 177, ' 0.3)']; ['(5.0 ', 177, ' 0.3)']; ...
            ['(6.1 ', 177, ' 0.3)']; '(0.0)'; ''};

        scan_strV = scan_strV + 1;
    else
        fitType = 'L';
        scan_str = (slideNum + 1) / 2;
        replace(summarySlide, 'Title 1', ['Subject ', subject_id, ' s', num2str(scan_str)]);
        rbc_static = {''; sprintf('%3.3f', nmrFit_ppm.area(1)); ...
            sprintf('%4.1f', nmrFit_ppm.freq(1)); ...
            sprintf('%2.1f', nmrFit_ppm.fwhm(1)); '-----'; ...
            sprintf('%3.1f', nmrFit_ppm.phase(1)); ...
            sprintf('%3.1f', RBCSNR)};

        rbc_ref = {'Ref. Values'; ['(0.59 ', 177, ' 0.14)']; ...
            ['(216.0 ', 177, ' 0.7)']; ['(10.0 ', 177, ' 0.4)']; ...
            '-----'; ['(27.3 ', 177, ' 3.6)']; ''};

        bar_static = {''; '1.00'; ...
            sprintf('%4.1f', nmrFit_ppm.freq(2)); ...
            sprintf('%2.1f', nmrFit_ppm.fwhm(2)); '-----'; ...
            '0.0'; ...
            sprintf('%3.1f', BarrSNR)};

        bar_ref = {'Ref. Values'; '(1.0)'; ...
            ['(197.7 ', 177, ' 0.3)']; ['(5.0 ', 177, ' 0.3)']; ...
            '-----'; '(0.0)'; ''};
    end

    staticTable = Table([static_params, rbc_static, rbc_ref, bar_static, bar_ref]);
    colSpecsAmp(1) = ColSpec('1in');
    colSpecsAmp(2) = ColSpec('0.75in');
    colSpecsAmp(3) = ColSpec('1in');
    colSpecsAmp(4) = ColSpec('0.7in');
    colSpecsAmp(5) = ColSpec('1in');

    staticTable.ColSpecs = colSpecsAmp;
    staticTable.Height = '2.5in';

    staticTable.ColSpecs(1).HAlign = 'Center';
    staticTable.ColSpecs(1).VAlign = 'Middle';
    staticTable.ColSpecs(1).FontSize = '11pt';
    staticTable.ColSpecs(2).HAlign = 'Center';
    staticTable.ColSpecs(2).VAlign = 'Middle';
    staticTable.ColSpecs(2).FontSize = '12pt';
    staticTable.ColSpecs(3).HAlign = 'Center';
    staticTable.ColSpecs(3).VAlign = 'Middle';
    staticTable.ColSpecs(3).FontSize = '11pt';
    staticTable.ColSpecs(3).FontColor = '#767171';
    staticTable.ColSpecs(4).HAlign = 'Center';
    staticTable.ColSpecs(4).VAlign = 'Middle';
    staticTable.ColSpecs(4).FontSize = '12pt';
    staticTable.ColSpecs(5).HAlign = 'Center';
    staticTable.ColSpecs(5).VAlign = 'Middle';
    staticTable.ColSpecs(5).FontSize = '11pt';
    staticTable.ColSpecs(5).FontColor = '#767171';

    staticTable.row(1).Style = {VAlign('bottom')};
    staticTable.row(1).FontSize = '11pt';
    staticTable.row(7).Style = {VAlign('Top')};

    staticTable.StyleName = 'No Style, No Grid';
    % BD no longer commenting out .Border to get things running
    staticTable.Border = 'none';
    staticTable.ColSep = 'none';
    staticTable.RowSep = 'none';
    staticTable.X = '0.13in';
    hrTable.Y = '1.08in';
    replace(summarySlide, 'static table', staticTable);

    % Add figures
    p1 = mlreportgen.ppt.Picture(fullfile(path, ['dyn', fitType, imgName{slideNum}, '.tif']));
    replace(summarySlide, 'dynPicture', p1);
    p2 = mlreportgen.ppt.Picture(fullfile(path, ['dyn', fitType, '_oscs_new', imgName{slideNum}, '.tif']));
    replace(summarySlide, 'oscillationPicture', p2);
    p3 = mlreportgen.ppt.Picture(fullfile(path, ['static', fitType, imgName{slideNum}, '.tif']));
    replace(summarySlide, 'staticPicture', p3);
    p4 = mlreportgen.ppt.Picture(fullfile(path, ['dyn', fitType, '_error', imgName{slideNum}, '.tif']));
    replace(summarySlide, 'snrPicture', p4);
end

% Generate and open the presentation
close(slides);

if ispc
    winopen(slidesFile);
end

end