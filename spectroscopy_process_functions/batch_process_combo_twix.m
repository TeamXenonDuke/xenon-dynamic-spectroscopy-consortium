clear all_twix_dyn;
f1 = filepath('D:\OneDrive\Documents\Duke\CIVM Research\Siemens Data\'); 
[fpath, fname, fext] = fileparts(f1);
fall = dir(fullfile(fpath, '*_Dynamic_Cali_combo*'));
for nscan = 1:size(fall,1)
    all_twix_dyn{nscan,1} = [fall(nscan).folder,'\',fall(nscan).name];
end 

num = 2;
for k = 1:length(all_twix_dyn)
%     raw_filepath = filepath();
    raw_filepath = all_twix_dyn{k}
    [twix_path, raw_file] = fileparts(raw_filepath);
    
    dynName = extractAfter(raw_file,'Dynamic');
    if strcmp(dynName,"Spec_170118")
        dynName(end-10:end) = [];
    end
    if strcmp(dynName,"_170118")
        dynName(end-6:end) = [];
    end
    
    dynName2 = dynName;
    while exist([twix_path,'\dynV',dynName2,'.mat'],'file') ~= 0
        dynName2 = [dynName,'_',num2str(num)];
        num = num + 1;
    end 
    
    [twix_path,'\dynV',dynName2] 
    dyn = fitDynamicSpec(raw_filepath);  
  
    save([twix_path,'\dynV',dynName2],'dyn')
  
    dyn = fitFBtwix(raw_filepath, 3, 5);   
    save([twix_path,'\dyn',dynName2],'dyn')
end    

% dynamicSummary(raw_path,dyn_path,BHs,save_fig_flag,figname)