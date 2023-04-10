function [fids, dwell_time, npts, tr, xeFreqMHz, rf_excitation] = readRawDyn(raw_path)

[~, ~, scanner] = fileparts(raw_path);

% Read in twix or P file and define associated variables
if strcmp(scanner, '.dat')
    % Twix file from Siemens
    twix = readtwix(raw_path);
    npts = size(twix.data, 1); % Number of samples per FID
    dwell_time = twix.hdr.MeasYaps.sRXSPEC.alDwellTime{1, 1} * 10^-9; % Receiver bandwidth (kHz); works for both calibration types
    tr = twix.hdr.Config.TR(1) * 1E-6; % Time between each sample
    fids = twix.data;
    if isfield(twix.hdr.Config, 'Frequency')
        % UVA Siemens File
        xeFreqMHz = twix.hdr.Config.Frequency * 10e-7; %34.093484
    elseif isfield(twix.hdr.Meas, 'lFrequency')
        % Duke Siemens File
        xeFreqMHz = twix.hdr.Meas.lFrequency * 10e-7; %34.091516
    end

    % Magnetic Field Strength
    mag_fstregth = twix.hdr.Dicom.flMagneticFieldStrength;

    % This if, else block will read the exitation from different fields,
    % depending on our scan version

    if isfield(twix.hdr.Phoenix, 'sWiPMemBlock')
        if isfield(twix.hdr.Phoenix.sWiPMemBlock, 'adFree')
            excitation = twix.hdr.Phoenix.sWiPMemBlock.adFree{1, 9};
        elseif isfield(twix.hdr.Phoenix.sWiPMemBlock, 'alFree') % BD added 7/5/22 per Ari to make cali 2104 files work
            excitation = twix.hdr.Phoenix.sWiPMemBlock.alFree{1, 5};
        end
    elseif isfield(twix.hdr.Phoenix, 'sWipMemBlock')
        if isfield(twix.hdr.Phoenix.sWipMemBlock, 'alFree')
            excitation = twix.hdr.Phoenix.sWipMemBlock.alFree{1, 5};
        end
    end

    gyro_ratio = 11.777; % gyromagnetic ratio of 129Xe in MHz/Tesla

    % RF excitation will be in ppm, either 218ppm or 208 ppm
    rf_excitation = round(excitation/(gyro_ratio * mag_fstregth));

elseif strcmp(scanner, '.7')
    % Pfile from GE
    pfile = GE.Pfile.read(raw_path);
    MRI.DataProcessing.checkForOverranging(pfile); % Check for overranging
    pfile = MRI.DataProcessing.removeBaselineViews(pfile); % Remove baselines
    bw = 1000 * pfile.rdb.rdb_hdr_user12; % Receiver bandwidth (kHz)
    dwell_time = 1 / (2 * bw); % Time between each sample
    dwell_time = Math.nearestMultipleOf(dwell_time, 0.000002);
    npts = pfile.rdb.rdb_hdr_frame_size;
    tr = pfile.image.tr * 1E-6;
    fids = pfile.data;
    xeFreqMHz = 17.660445;

elseif strcmp(scanner, '.h5')
    % mrd file
    % Reading the k-space data
    dataset_data = h5read(raw_path, '/dataset/data');
    all_kspace_data = dataset_data.data;

    %Reshaping the k-space
    npts = size(all_kspace_data{1}, 1) / 2; % 512 for duke
    nfids = size(all_kspace_data, 1); % 520 for duke

    fids = [];
    for ii = 1:nfids
        fid = all_kspace_data{ii};
        i = 1;
        for j = 1:npts
            fids(j, ii) = double(complex(fid(i), fid(i+1)));
            i = i + 2;
        end
    end
    % Reading dwell_time from a k-space fid
    cali_data_head = dataset_data.head;
    dwell_time_all = cali_data_head.sample_time_us;
    dwell_time = double(dwell_time_all(1)) * 10^-6;

    % Dataset header variables are in the xml field - gives string
    dataset_header = h5read(raw_path, '/dataset/xml');

    xml_struct = xml2struct(dataset_header); % the xml2struct converts the string
    tr_mrd = str2double(xml_struct.ismrmrdHeader.sequenceParameters.TR.Text);

    % Threshold to not scale by 1E-3 some tr_mrd are in seconds
    if tr_mrd < 1E-1
        tr = tr_mrd; % this is 0.015 for duke
    else
        tr = tr_mrd * 1E-3;
    end

    site_name = xml_struct.ismrmrdHeader.acquisitionSystemInformation.systemVendor.Text;

    % IOWA GE data
    switch site_name

        case 'GE'

            % conjugating the data solved the issues but we don't exactly know
            % the reason

            fids = conj(fids);
            % remove the first two points due to initial recovery time
            % and zero pad the end
            fids = fids(3:end, :);
            fids(size(fids, 1):size(fids, 1)+2, :) = 0;

            fids = fids / max(abs(fids(:)), [], 'all'); % max Normalizing
            xeFreqMHz_mrd = str2double(xml_struct.ismrmrdHeader.encoding.trajectoryDescription.userParameterDouble{1, 2}.value.Text);
            xeFreqMHz = xeFreqMHz_mrd * 1E-6; % 34.0923 MHz for duke
            mag_fstrength = str2double(xml_struct.ismrmrdHeader.acquisitionSystemInformation.systemFieldStrengthu_T.Text);
            % see if the excitation is in the mrd file
            try
                excitation = str2double(xml_struct.ismrmrdHeader.encoding.trajectoryDescription.userParameterDouble{1, 3}.value.Text);
                % RF excitation will be in ppm, either 218ppm or 208 ppm
                rf_excitation = round(abs(excitation-xeFreqMHz*1E6)/(xeFreqMHz), 1);
            catch
                rf_excitation = 218;
            end
            % Philips - Cincinnati Children's Hospital Medical Center Data
        case 'Philips'
            xeFreqMHz = str2double(xml_struct.ismrmrdHeader.userParameters.userParameterLong.value.Text) * 1E-6;
            mag_fstrength = xml_struct.ismrmrdHeader.acquisitionSystemInformation.systemFieldStrengthu_T.Text;
            rf_excitation = str2double(xml_struct.ismrmrdHeader.userParameters.userParameterDouble.value.Text);

        case 'SIEMENS'
            xeFreqMHz_mrd = str2double(xml_struct.ismrmrdHeader.encoding.trajectoryDescription.userParameterDouble{1, 2}.value.Text);
            xeFreqMHz = xeFreqMHz_mrd * 1E-6;
            dwell_time = str2double(xml_struct.ismrmrdHeader.encoding.trajectoryDescription.userParameterDouble{1}.value.Text);
            % see if the excitation is in the mrd file
            try
                excitation = str2double(xml_struct.ismrmrdHeader.encoding.trajectoryDescription.userParameterDouble{1, 3}.value.Text);
                % RF excitation will be in ppm, either 218ppm or 208 ppm
                rf_excitation = round(abs(excitation-xeFreqMHz*1E6)/(xeFreqMHz), 1);
            catch
                rf_excitation = 218;
            end
        otherwise
            error('Not a valid site name')
    end

else
    error('Unknown Raw File Type')
end
