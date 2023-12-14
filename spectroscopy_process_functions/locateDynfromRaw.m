function [dyn_loc, BHs, subjName] = locateDynfromRaw(raw_path,varargin)

[raw_folder, ~, scanner] = fileparts(raw_path);
filename = 'dynV5';
default_BHs = [2 7];

% for Siemens data
if strcmp(scanner,'.dat')
    
    % Find end of dyn variable if more than one flavor of dynamic scan
    [~, dynName] = fileparts(raw_path);
    dynName = extractAfter(dynName,'Dynamic_');
    if contains(dynName,"Spec_170118")
        dynName(end-10:end) = [];
    elseif contains(dynName,"Spec_")
        dynName = dynName(6:end);
    elseif contains(dynName,"Cali_combo")
        dynName = [];
    elseif strcmp(dynName,"_Cali_15TR_7430HZ")
        dynName = [];
    elseif contains(dynName,"Cali.")
        dynName = [];
    elseif contains(dynName,"Cali_2")
        dynName = dynName(end-1:end);
    elseif contains(dynName,"Cali_3")
        dynName = dynName(end-1:end);
    end
    
    if exist([raw_folder,'/dynV','.mat'],'file') && exist([raw_folder,'/dynV_Cali_combo','.mat'],'file')
            if contains(raw_path,'Dynamic_Cali_combo')
                dyn_loc = [raw_folder,'/dynV_Cali_combo','.mat'];
                BH_loc = [raw_folder,'/BHs','.mat'];
                subjName = extractAfter(extractAfter(raw_folder,'Siemens Data/'),'/');
            else 
                dyn_loc = [raw_folder,'/dynV','.mat'];
                BH_loc = [raw_folder,'/BHs','.mat'];
                subjName = [extractAfter(extractAfter(raw_folder,'Siemens Data/'),'/'),'s2'];               
            end 
    elseif exist([raw_folder,'/dynV','.mat'],'file')
        dyn_loc = [raw_folder,'/dynV','.mat'];
        BH_loc = [raw_folder,'/BHs','.mat'];
        subjName = extractAfter(extractAfter(raw_folder,'Siemens Data/'),'/');  
        % Check to see if raw path corresponds with first or second
        % dynamic scan
        [fpath, fname, fext] = fileparts(raw_path);
        f_all = dir([fullfile(fpath,'*.dat')]);
        f_a = contains({f_all.name},'_Dynamic') | contains({f_all.name},'_cali') | contains({f_all.name},'_Cali');
        f_b = -1*contains({f_all.name},'fid_cali_'); % remove calibration only scans
        f_c = f_a + f_b;
        f_all(f_c == 0) = [];
        
        if length(f_all) > 1
            fnum = str2num(extractBefore(extractAfter(raw_path,'meas_MID'),'_'));
            for k = 1:length(f_all)
                   f(k) = str2num(extractBefore(extractAfter(f_all(k).name,'meas_MID'),'_'));
            end 
            if find(f == fnum) ~= 1
                add_dyn = ['_',num2str(find(f == fnum,1))];
                add_name = ['s',num2str(find(f == fnum,1))];
            else 
                add_dyn = '';
                add_name = '';
            end
        else 
            add_dyn = '';
            add_name = '';
        end
        dyn_loc = [raw_folder,'/dynV',add_dyn,'.mat'];
        BH_loc = [raw_folder,'/BHs','.mat'];
        subjName = [extractAfter(extractAfter(raw_folder,'Siemens Data/'),'/'),add_name]; 
    elseif exist([raw_folder,'/dynV_Cali_combo','.mat'],'file')
            % Check to see if raw path corresponds with first or second
            % cali_combo scan
            [fpath, fname, fext] = fileparts(raw_path);
            f_all = dir(fullfile(fpath, '*_Dynamic_Cali*'));
            if length(f_all) > 1
                fnum = str2num(extractBefore(extractAfter(raw_path,'meas_MID'),'_'));
                for k = 1:length(f_all)
                       f(k) = str2num(extractBefore(extractAfter(f_all(k).name,'meas_MID'),'_'));
                end 
                if find(f == fnum) ~= 1
                    add_dyn = ['_',num2str(find(f == fnum))];
                    add_name = ['s',num2str(find(f == fnum))];
                else 
                    add_dyn = '';
                    add_name = '';
                end
            else 
                add_dyn = '';
                add_name = '';
            end
        dyn_loc = [raw_folder,'/dynV_Cali_combo',add_dyn,'.mat'];
        BH_loc = [raw_folder,'/BHs','.mat'];
        subjName = [extractAfter(extractAfter(raw_folder,'Siemens Data/'),'/'),add_name];  
        if ~exist([raw_folder,'/dynV_Cali_combo',add_dyn,'.mat'],'file') && exist([raw_folder,'/dynV_',dynName,'.mat'],'file')
            dyn_loc = [raw_folder,'/dynV_',dynName,'.mat']; 
        end
    elseif exist([raw_folder,'/dynV_',dynName,'.mat'],'file')
        dyn_loc = [raw_folder,'/dynV_',dynName,'.mat'];
        BH_loc = [raw_folder,'/BHs','.mat'];
        subjName = [extractAfter(extractAfter(raw_folder,'Siemens Data/'),'/')];
        % Check to see if raw path corresponds with first or second
        % scan
        [fpath, fname, fext] = fileparts(raw_path);
        f_all = dir(fullfile(fpath, '*_Dynamic_*'));
        if length(f_all) > 1
            fnum = str2num(extractBefore(extractAfter(raw_path,'meas_MID'),'_'));
            for k = 1:length(f_all)
                   f(k) = str2num(extractBefore(extractAfter(f_all(k).name,'meas_MID'),'_'));
            end 
            if find(f == fnum) ~= 1
                add_dyn = ['_',num2str(find(f == fnum))];
                add_name = ['s',num2str(find(f == fnum))];
            else 
                add_dyn = '';
                add_name = '';
            end
        else 
            add_dyn = '';
            add_name = '';
        end
        dyn_loc = [raw_folder,'/dynV_',dynName,add_dyn,'.mat'];
        BH_loc = [raw_folder,'/BHs',dynName,add_dyn,'.mat'];
        subjName = [extractAfter(extractAfter(raw_folder,'Siemens Data/'),'/'),add_name];
    elseif exist([raw_folder,'/dynV_Cali','.mat'],'file')
            % Check to see if raw path corresponds with first or second
            % cali_combo scan
            [fpath, fname, fext] = fileparts(raw_path);
            f_all = dir(fullfile(fpath, '*_Dynamic_Cali*'));
            if length(f_all) > 1
                fnum = str2num(extractBefore(extractAfter(raw_path,'meas_MID'),'_'));
                for k = 1:length(f_all)
                       f(k) = str2num(extractBefore(extractAfter(f_all(k).name,'meas_MID'),'_'));
                end 
                if find(f == fnum) ~= 1
                    add_dyn = ['_',num2str(find(f == fnum))];
                    add_name = ['s',num2str(find(f == fnum))];
                else 
                    add_dyn = '';
                    add_name = '';
                end
            else 
                add_dyn = '';
                add_name = '';
            end
        dyn_loc = [raw_folder,'/dynV_Cali',add_dyn,'.mat'];
        BH_loc = [raw_folder,'/BHs','.mat'];
        subjName = [extractAfter(extractAfter(raw_folder,'Siemens Data/'),'/'),add_name];
    else
        dyn_loc = ''; BH_loc = ''; subjName = '';
        disp(['Could not locate processed dynamic data for: ',raw_path]) 
    end 
    
% for GE data
elseif strcmp(scanner,'.7')
    subjName = extractAfter(raw_folder,'Subject_');
    if contains(subjName, "002_")
        subjName = erase(subjName,'002_');
        subjName = erase(subjName,'/Scan'); subjName = erase(subjName,'/');
        subjName = erase(subjName,'Dynamic_Spec');
        if contains(subjName, "Series")
            subjName = extractBefore(subjName,'Series');
        end 
        while subjName(1) == '0'
            subjName(1) = '';
        end 
        
        dyn_loc = strcat('D:/OneDrive/Documents/Duke/CIVM Research/MATLAB/ReadingGE/dynData/',subjName,'/dynV5.mat');
        BH_loc = strcat('D:/OneDrive/Documents/Duke/CIVM Research/MATLAB/ReadingGE/dynData/',subjName,'/BHs.mat');        
    elseif contains(raw_folder, '003_')
        subjName = strrep(subjName,'_','-');
        
        dyn_loc = strcat('D:/OneDrive/Documents/Duke/CIVM Research/MATLAB/ReadingGE/dynData/',subjName,'/dynV5.mat');
        BH_loc = strcat('D:/OneDrive/Documents/Duke/CIVM Research/MATLAB/ReadingGE/dynData/',subjName,'/BHs.mat');
    else
        dyn_loc = ''; BH_loc = ''; subjName = '';
        disp(['Could determine subject name processed dynamic data for: ',raw_path]);   
    end 
else
    dyn_loc = ''; BH_loc = ''; subjName = '';
    disp(['Could not locate processed dynamic data for: ',raw_path])      
end 
    
if exist(BH_loc,'file')
    load(BH_loc);
else
    BHs = default_BHs;
end 

switch nargin
    case 2
        dyn_loc = strrep(dyn_loc,'dynV5','dyn');
        dyn_loc = strrep(dyn_loc,'dynV','dyn');
end

% If raw file or dyn has extra "_text" take off when reporting subjName
% if contains(subjName,"__")
%     subjName = extractBefore(subjName,'__');
% end
    

