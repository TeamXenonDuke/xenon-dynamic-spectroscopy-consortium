filelist_twix_dyn;
% pfiles_all_dynamic;

Nsubj = length(all_twix_dyn);

tic
for idx = 1:Nsubj
    [dyn_loc, BHs, subjName] = locateDynfromRaw(all_twix_dyn{idx});
    [yyyy, mm, dd] = getScanDate(all_twix_dyn{idx});
    scanDate = [yyyy,'-',mm,'-',dd];
    %     save_path = fileparts(dyn_loc{k});
    save_path = '/Users/eab78/OneDrive/Documents/Duke/CIVM Research/Siemens Data/Processed Data';
    if isempty(dyn_loc), continue, end % skip subject if dyn data cannot be found 
   
    if exist([save_path,'/spect_',subjName,'.mat'],'file')
        load([save_path,'/spect_',subjName]);
        save_stats = 0;
        disp(['Already Exists: ',subjName])
    else 
        load(dyn_loc)
        
    if contains(subjName,'s')
        figname = ['_',extractAfter(subjName,'s')];
    else 
        figname = '';
    end 

        [amp, nmrFit] = dynamicSummary(all_twix_dyn{idx},dyn_loc,BHs,1,figname,'oscType','sine');
        save_stats = 1;
    
        amp.area = amp.area*100;

        static.rbc_area = nmrFit.area(1);
        static.rbc_freq = nmrFit.freq(1);
        static.rbc_fwhm = nmrFit.fwhm(1);
        static.rbc_fwhmG = 0;
        static.rbc_phase = nmrFit.phase(1); 
%         static.rbc_ppm.SNR_dis = SNR_dis

        static.bar_area = 1;
        static.bar_freq = nmrFit.freq(2);
        static.bar_fwhm = nmrFit.fwhm(2);
        static.bar_fwhmG = nmrFit.fwhmG(2);
        static.bar_phase = 0;

        static.gas_area = nmrFit.area(3);
        static.gas_freq = nmrFit.freq(3);
        static.gas_fwhm = nmrFit.fwhm(3);
        static.gas_fwhmG = 0;
        static.gas_phase = nmrFit.phase(3);   

        spect.dyn = dyn;
        spect.amp = amp;
        spect.nmrFit = nmrFit;
        spect.static = static;
        spect.scanDate = scanDate;
        spect.processed_date = datestr(now);
        spect.raw_path = all_twix_dyn{idx};
        spect.dyn_path = dyn_loc;
        spect.subjName = subjName;

        saveName = subjName;
        replace(saveName,'-','_');
        save([save_path,'/spect_',saveName],'spect');
        disp(['Saved:  ',saveName,'/spect_',subjName])
    end 
    
end 
toc 

    