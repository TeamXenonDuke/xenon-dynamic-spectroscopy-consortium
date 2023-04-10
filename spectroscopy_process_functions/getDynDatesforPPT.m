function dates = getDynDatesforPPT(varargin)

if nargin == 2
    inputType = 'raw_dyn_path';
    raw_path = varargin{1};
    dyn_path = varargin{2};
elseif nargin == 1 && contains(contains(varargin{1},'.dat')) || contains(contains(varargin{1},'.7'))
    inputType = 'raw_dyn_path';
    raw_path = contains(varargin{1});
    dyn_path = locateDynfromRaw(raw_path);
elseif nargin == 1 && contains(contains(varargin{1},'.mat'))
    inputType = 'struct';
    raw_path = spect.dyn.raw_file;
else 
    disp('data type not supported')
end 
    
[yyyy, mm, dd] = getScanDate(raw_path);
dates.scan = [yyyy,'-',mm,'-',dd];

switch inputType
    case 'raw_dyn_path'
        load(dyn_path);

        if isfield(dyn,'processed_date')
            dates.dyn_fit = datestr(dyn.processed_date,'yyyy-mm-dd');
        else
            file = dir(dyn_path);
            dates.dyn_fit = [datestr(file.date,'yyyy-mm-dd'),'*'];
        end 
        dates.rbc_fit = datestr(now,'yyyy-mm-dd');
        
    case 'struct'
        if isfield(spect.dyn,'processed_date')
            dates.dyn_fit = datestr(spect.dyn.processed_date,'yyyy-mm-dd');
        else       
            file = dir(locateDynfromRaw(raw_path));
            dates.dyn_fit = [datestr(file.date,'yyyy-mm-dd'),'*'];
        end 
        
        dates.rbc_fit = datestr(spect.processed_date,'yyyy-mm-dd');
end
