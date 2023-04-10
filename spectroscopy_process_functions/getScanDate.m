function [yyyy, mm, dd] = getScanDate(raw_path)

[~, ~, scanner] = fileparts(raw_path);

% Read in twix or P file and get date of scan aquisition
if strcmp(scanner,'.dat')
    twix = mapVBVD(raw_path);
    scanDate = twix.hdr.Phoenix.tReferenceImage0;   
    scanDate = strsplit(scanDate,'.');
    scanDate = scanDate{end};
    yyyy = scanDate(1:4);
    mm = scanDate(5:6);
    dd = scanDate(7:8);
elseif strcmp(scanner,'.7')
    pfile = GE.Pfile.read(raw_path);
    if(pfile.exam.ex_datetime == 0)
        %Create exam and series dates in YYYY_MM_DD format
        yyyy = ['20',pfile.rdb.rdb_hdr_scan_date(8:9)'];
        mm = pfile.rdb.rdb_hdr_scan_date(1:2)';
        dd = pfile.rdb.rdb_hdr_scan_date(4:5)';
    else
        date_number = pfile.exam.ex_datetime/86400 + datenum(1970,1,1);
        scanDate = datestr(date_number,'yyyymmdd');
        yyyy = scanDate(1:4);
        mm = scanDate(5:6);
        dd = scanDate(7:8);
    end
elseif strcmp(scanner,'.h5')
    % Dataset header variables are in the xml field - gives string 
    dataset_header = h5read(raw_path, '/dataset/xml');
    xml_struct = xml2struct(dataset_header); % the xml2struct converts the string
    % if there is no scanDate, set it to today
    try
        scanDate = xml_struct.ismrmrdHeader.studyInformation.studyDate.Text;
    catch
        scanDate = num2str(yyyymmdd(datetime('today')));
    end

    yyyy = scanDate(1:4);
    mm = scanDate(5:6);
    dd = scanDate(7:8);
end 