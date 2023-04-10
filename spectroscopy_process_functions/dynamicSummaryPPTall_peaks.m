function slides = dynamicSummaryPPTall_peaks(slides,subj,folder,imgName,dates,amp,nmrFit_ppm)

import mlreportgen.ppt.*;

for slideNum = 1:length(amp)
    % Add a summary slide
    summarySlide = add(slides,'SpectralSummary');
    
    % Add all text elements
    replace(summarySlide,'Dyn Text','Dynamics by Resonance');
    replace(summarySlide,'RBC','RBC');
    replace(summarySlide,'Barrier','Barrier');
    replace(summarySlide,'Gas','Gas');
    replace(summarySlide,'Osc Text','Detrended RBC Oscillations');
    replace(summarySlide,'Amp Text','RBC Oscillation Amplitude*');
    replace(summarySlide,'Pkpk Text','*Peak-to-Peak');
    replace(summarySlide,'Static Text','Static Spectroscopy');
    replace(summarySlide,'StaticValues Text','Static Spectroscopy (Barrier Voigt)');
    replace(summarySlide,'dynSNR Text','Dynamic SNR');
    replace(summarySlide,'norm Text','*normalized to barrier peak');
    replace(summarySlide,'box1',' ');
    replace(summarySlide,'box2',' ');
    replace(summarySlide,'box3',' ');
    replace(summarySlide,'RBC2','RBC');
    replace(summarySlide,'Barrier2','Barrier');
    replace(summarySlide,'scanDate Text','Scan Date:');
    replace(summarySlide,'processingDate Text','Dynamic Fit Date:');
    replace(summarySlide,'ampDate Text','RBC Fit Date:');
    
    % Add dates
    replace(summarySlide,'scanDate',dates(slideNum).scan);
    replace(summarySlide,'processingDate',dates(slideNum).dyn_fit);
    replace(summarySlide,'ampDate',dates(slideNum).rbc_fit);

    % Add hr table
    hrTable = Table({'Heart Rate:',sprintf('%2.0f',amp(slideNum).hr); 'Mean SNR:', sprintf('%2.1f',amp(slideNum).snr)});
    colSpecs(1) = ColSpec('1.41in'); 
    colSpecs(2) = ColSpec('0.6in');
    hrTable.ColSpecs = colSpecs;

    hrTable.ColSpecs(1).Bold = 1;
    hrTable.ColSpecs(1).HAlign = 'Right';
    hrTable.ColSpecs(2).HAlign = 'Left';
    hrTable.FontSize = '14pt';
    hrTable.StyleName = 'No Style, No Grid';
    hrTable.Border = 'none'; hrTable.ColSep = 'none'; hrTable.RowSep = 'none';
    hrTable.X = '6.55in'; hrTable.Y = '6.41in';
    
    if amp(slideNum).snr <= 20
        hrTable.row(2).Style = {BackgroundColor('yellow')};
    end 

    replace(summarySlide,'HR_SNR',hrTable);

    % Add Oscillation table
    oscs_params = {''; 'Amplitude (%):'; ['Chemical Shift',10,'(ppm):']; ...
        ['Linewidth', 10, '(ppm):']; 'Phase (degrees):'};
    oscs = {''; sprintf('%3.1f', amp(slideNum).area); sprintf('%3.2f',amp(slideNum).freq); ...
        sprintf('%3.2f',amp(slideNum).fwhm); sprintf('%3.1f',amp(slideNum).phase)};
    oscs_ref = {'Healthy Ref. Values'; ['(9.4 ',177,' 2.7%)'];...
        ['(0.06 ',177,' 0.05)']; ['(0.19 ',177,' 0.09)'];...
        ['(1.3 ',177,' 0.7',176,')']};

    ampTable = Table([oscs_params, oscs, oscs_ref]);
    colSpecsAmp(1) = ColSpec('1.6in'); 
    colSpecsAmp(2) = ColSpec('1in');
    colSpecsAmp(3) = ColSpec('1.43in');
    ampTable.ColSpecs = colSpecsAmp;
    ampTable.Height = '2.2in';

    ampTable.ColSpecs(1).HAlign = 'Right'; ampTable.ColSpecs(1).VAlign = 'Middle';
    ampTable.ColSpecs(1).FontSize = '14pt'; 
    ampTable.ColSpecs(2).HAlign = 'Center'; ampTable.ColSpecs(2).VAlign = 'Middle';
    ampTable.ColSpecs(2).FontSize = '16pt';
    ampTable.ColSpecs(3).HAlign = 'Left'; ampTable.ColSpecs(3).VAlign = 'Middle';
    ampTable.ColSpecs(3).FontSize = '16pt'; ampTable.ColSpecs(3).FontColor = '#767171'; 
    ampTable.row(1).Style = {VAlign('bottom')}; 
    ampTable.row(1).FontSize = '11pt'; 
    ampTable.row(5).Style = {VAlign('Middle')}; 

    ampTable.StyleName = 'No Style, No Grid';
    ampTable.Border = 'none'; ampTable.ColSep = 'none'; ampTable.RowSep = 'none';
    ampTable.X = '9.05in'; hrTable.Y = '1.0in';

    replace(summarySlide,'amp table',ampTable);

    % Add Static table
    static_params = {''; ['Intensity', 10, 'Ratio*']; ['Shift', 10, '(ppm)'];...
        ['FWHM', 10, '(ppm)']; ['FWHMg', 10, '(ppm)'];...
        ['Phase', 10, '(degrees)']; ['5 FID SNR', 10, '(amp/noise)']};

    if sum(nmrFit_ppm(slideNum).fwhmG(:)) > 0
        fitType = 'V';
        replace(summarySlide,'Title 1',['Subject ',subj,' (Voigt)']);
            rbc_static = {''; sprintf('%3.2f', nmrFit_ppm(slideNum).area(1));...
            sprintf('%4.1f',nmrFit_ppm(slideNum).freq(1));
            sprintf('%2.1f',nmrFit_ppm(slideNum).fwhm(1)); '-----';
            sprintf('%3.1f',nmrFit_ppm(slideNum).phase(1)); ...
            sprintf('%3.1f',nmrFit_ppm(slideNum).SNR_dis(1))};

        rbc_ref = {'Ref. Values'; ['(0.59 ',177,' 0.12)'];
            ['(218.4 ',177,' 0.4)']; ['(8.7 ',177,' 0.3)'];...
            '-----'; ['(81.9 ',177,' 3.6)']; ''};

        bar_static = {''; '1.00';
            sprintf('%4.1f',nmrFit_ppm(slideNum).freq(2)); 
            sprintf('%2.1f',nmrFit_ppm(slideNum).fwhm(2)); 
            sprintf('%3.2f',nmrFit_ppm(slideNum).fwhmG(2));
            '0.0'; 
            sprintf('%3.1f',nmrFit_ppm(slideNum).SNR_dis(2))};

        bar_ref = {'Ref. Values'; '(1.0)';...
            ['(197.7 ',177,' 0.3)']; ['(5.0 ',177,' 0.3)'];...
            ['(6.1 ',177,' 0.3)']; '(0.0)'; ''};
    else 
        fitType = 'L';
        scan_str = (slideNum+1)/2;
        replace(summarySlide,'Title 1',['Subject ',subj]);
        rbc_static = {''; sprintf('%3.2f', nmrFit_ppm(slideNum).area(1));
        sprintf('%4.1f',nmrFit_ppm(slideNum).freq(1));
        sprintf('%2.1f',nmrFit_ppm(slideNum).fwhm(1)); '-----';
        sprintf('%3.1f',nmrFit_ppm(slideNum).phase(1));
        sprintf('%3.1f',nmrFit_ppm(slideNum).SNR_dis(1))};

        rbc_ref = {'Ref. Values'; ['(0.59 ',177,' 0.14)'];
            ['(216.0 ',177,' 0.7)']; ['(10.0 ',177,' 0.4)'];
            '-----'; ['(27.3 ',177,' 3.6)']; ''};

        bar_static = {''; '1.00'; ...
            sprintf('%4.1f',nmrFit_ppm(slideNum).freq(2)); 
            sprintf('%2.1f',nmrFit_ppm(slideNum).fwhm(2)); '-----';
            '0.0'; 
            sprintf('%3.1f',nmrFit_ppm(slideNum).SNR_dis(2))};

        bar_ref = {'Ref. Values'; '(1.0)';...
            ['(197.7 ',177,' 0.3)']; ['(5.0 ',177,' 0.3)'];...
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

    staticTable.ColSpecs(1).HAlign = 'Center'; staticTable.ColSpecs(1).VAlign = 'Middle';
    staticTable.ColSpecs(1).FontSize = '11pt'; 
    staticTable.ColSpecs(2).HAlign = 'Center'; staticTable.ColSpecs(2).VAlign = 'Middle';
    staticTable.ColSpecs(2).FontSize = '12pt';
    staticTable.ColSpecs(3).HAlign = 'Center'; staticTable.ColSpecs(3).VAlign = 'Middle';
    staticTable.ColSpecs(3).FontSize = '11pt'; staticTable.ColSpecs(3).FontColor = '#767171';
    staticTable.ColSpecs(4).HAlign = 'Center'; staticTable.ColSpecs(4).VAlign = 'Middle';
    staticTable.ColSpecs(4).FontSize = '12pt';
    staticTable.ColSpecs(5).HAlign = 'Center'; staticTable.ColSpecs(5).VAlign = 'Middle';
    staticTable.ColSpecs(5).FontSize = '11pt'; staticTable.ColSpecs(5).FontColor = '#767171';

    staticTable.row(1).Style = {VAlign('bottom')}; 
    staticTable.row(1).FontSize = '12pt';
    staticTable.row(7).Style = {VAlign('Top')}; 

    staticTable.StyleName = 'No Style, No Grid';
    staticTable.Border = 'none'; staticTable.ColSep = 'none'; staticTable.RowSep = 'none';
    staticTable.X = '0.13in'; hrTable.Y = '1.08in';
    replace(summarySlide,'static table',staticTable);

    % Add figures
    p1 = mlreportgen.ppt.Picture([folder,'/dyn',fitType,imgName{slideNum},'.tif']);
    replace(summarySlide,'dynPicture',p1);
    p2 = mlreportgen.ppt.Picture([folder,'/dyn',fitType,'_oscs_peaks',imgName{slideNum},'.tif']);
    replace(summarySlide,'oscillationPicture',p2);
    p3 = mlreportgen.ppt.Picture([folder,'/static',fitType,imgName{slideNum},'.tif']);
    replace(summarySlide,'staticPicture',p3);
    p4 = mlreportgen.ppt.Picture([folder,'/dyn',fitType,'_error',imgName{slideNum},'.tif']);
    replace(summarySlide,'snrPicture',p4);
end 


    % Add RBC_to_barrier slide
%     rbc_to_barSlide = add(slides,'rbc_to_bar');
%     replace(rbc_to_barSlide,'Title 1',['Subject ',subj,' (Voigt)']);
%     replace(rbc_to_barSlide,'Osc Text','RBC:Barrier');
%     p5 = mlreportgen.ppt.Picture([folder,'\rbc_to_bar',fitType,imgName{slideNum},'.tif']);
%     replace(rbc_to_barSlide,'oscillationPicture',p5);
end 