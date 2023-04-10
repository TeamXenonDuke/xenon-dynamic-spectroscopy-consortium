clear;

filelist_twix_dyn;

subj = 1:length(all_twix_dyn);
% getPairsIdx;
% subj = pairs;
filename = ['all_stats',datestr(now,'_yy-mm-dd'),'.xlsx'];
Nsubj = length(subj);
% Subject = cell(Nsubj,1); 

for idx = 1:Nsubj
    
    k = subj(idx);
    [dyn_loc, BH_loc, subjName] = locateDynfromRaw(all_twix_dyn{k});
    save_path = '/Users/eab78/OneDrive/Documents/Duke/CIVM Research/Siemens Data/Processed Data/';
    load([save_path,'/spect_',subjName]);
    
    if length(subjName(strfind(subjName,'-'):end)) > 3
        subjID{idx,1} = subjName(1:strfind(subjName,'-') + 3);
    else
        subjID{idx,1} = subjName(1:strfind(subjName,'-') + 2);
    end 
    
    visitScan = extractAfter(subjName,subjID{idx,1});
    
    if contains(subjID{idx,1},'006')
        if ~contains(visitScan,'s'), visit{idx,1} = ''; scan{idx,1} = 'pre'; end
        if contains(visitScan,'s2') || contains(visitScan,'_2'), visit{idx,1} = ''; scan{idx,1} = 'post'; end
        if contains(visitScan,'s3') || contains(visitScan,'_3'), visit{idx,1} = ''; scan{idx,1} = 'final'; end
    elseif contains(visitScan,'s') || exist([save_path,'\spect_',subjName,'s2.mat'],'file') || exist([save_path,'\spect_',subjName,'s3.mat'],'file')
        if ~contains(subjName,'s'), visitScan = [visitScan, 's1']; end
        visit{idx,1} = visitScan(1:end-2);
        visitScan(1:end-2) = [];
        scan{idx,1} = visitScan;
        
    else
        visit{idx,1} = visitScan;
        scan{idx,1} = '';
    end
    
    twix = readtwix(all_twix_dyn{k});
    age{idx,1} = twix.hdr.Meas.flPatientAge;
    tr(idx,1) = twix.hdr.Config.TR(1)*1E-6;  
    
    if isfield(twix.hdr.Phoenix, 'sWiPMemBlock')
        % Duke twix file
        if isfield(twix.hdr.Phoenix.sWiPMemBlock,'adFree')
            voltage{idx,1} = twix.hdr.Phoenix.sWiPMemBlock.adFree{4};
        elseif isfield(twix.hdr.Phoenix.sWiPMemBlock,'alFree')
            % Duke twix file using UVA sequence
            voltage{idx,1} = twix.hdr.Phoenix.sWiPMemBlock.alFree{1};
        else
            voltage{idx,1} = [];
        end 
    elseif isfield(twix.hdr.Phoenix, 'sWipMemBlock')
        % New Duke or UVA twix file
        voltage{idx,1} = twix.hdr.Phoenix.sWipMemBlock.alFree{1};
    else
        voltage{idx,1} = [];
    end 
       
    Subject{idx,1} = spect.subjName;
    Amp_Osc(idx,1) = spect.amp.area;
    Shift_Osc(idx,1) = spect.amp.freq;
    FWHM_Osc(idx,1) = spect.amp.fwhm;
    Phase_Osc(idx,1) = spect.amp.phase;
    HR(idx,1) = spect.amp.hr;
    SNR(idx,1) = spect.amp.snr;
    rbc_rmse(idx,1) = spect.amp.area_gof.rmse;
    freq_rmse(idx,1) = spect.amp.freq_gof.rmse;
    rbc_r2(idx,1) = spect.amp.area_gof.rsquare;
    freq_r2(idx,1) = spect.amp.freq_gof.rsquare;

    RBC_Bar_Ratio(idx,1) = spect.static.rbc_area;
    RBC_Shift(idx,1) = spect.static.rbc_freq;
    RBC_FWHM(idx,1) = spect.static.rbc_fwhm;
    RBC_Phase(idx,1) = spect.static.rbc_phase;
    RBC_SNR(idx,1) = spect.nmrFit.SNR_dis(1);

    Bar_Bar_Ratio(idx,1) = spect.static.bar_area;
    Bar_Shift(idx,1) = spect.static.bar_freq;
    Bar_FWHM_L(idx,1) = spect.static.bar_fwhm;
    Bar_FWHM_G(idx,1) = spect.static.bar_fwhmG;
    Bar_Phase(idx,1) = spect.static.bar_phase;
    Bar_SNR(idx,1) = spect.nmrFit.SNR_dis(2);

    Gas_Bar_Ratio(idx,1) = spect.static.gas_area;
    Gas_Shift(idx,1) = spect.static.gas_freq;
    Gas_FWHM(idx,1) = spect.static.gas_fwhm;
    Gas_Phase(idx,1) = spect.static.gas_phase;   
    Gas_SNR(idx,1) = spect.nmrFit.SNR_dis(3);
    
    if isfield(spect,'scanDate')
        scan_date{idx,1} = spect.scanDate;
    else
        [yyyy, mm, dd] = getScanDate(all_twix_dyn{idx});
        scan_date{idx,1} = [yyyy,'-',mm,'-',dd];
    end
    disp(['Finished:  ',subjName])
end 
   
% output for peaks without statics
% T = table(subjID,visit,scan,repmat("",[length(subjID),1]), scan_date,...
%     repmat("",[length(subjID),1]),repmat("",[length(subjID),1]),...
%     age,tr,repmat("",[length(subjID),1]),voltage,...
%     Amp_Osc,Shift_Osc,FWHM_Osc,Phase_Osc,HR,SNR);

T = table(subjID,visit,scan,repmat("",[length(subjID),1]), scan_date,...
    repmat("",[length(subjID),1]),repmat("",[length(subjID),1]),...
    age,tr,repmat("",[length(subjID),1]),voltage,...
    Amp_Osc,Shift_Osc,FWHM_Osc,Phase_Osc,HR,SNR,...
    rbc_rmse,rbc_r2,freq_rmse,freq_r2,...
    RBC_Bar_Ratio,RBC_Shift,RBC_FWHM,RBC_Phase,RBC_SNR,...
    Bar_Bar_Ratio,Bar_Shift,Bar_FWHM_L,Bar_FWHM_G,Bar_Phase,Bar_SNR,...
    Gas_Bar_Ratio,Gas_Shift,Gas_FWHM,Gas_Phase,Gas_SNR);

writetable(T,[save_path,filename]);
  
    