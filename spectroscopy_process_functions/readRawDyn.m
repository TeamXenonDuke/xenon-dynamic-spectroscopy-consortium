function [fids, dwell_time, npts, tr, xeFreqMHz, rf_excitation] = readRawDyn(raw_path)
    
    global droppt_N;
    
    if isempty(droppt_N)
        droppt_N = 0; % default if not set in main
    end

    [~, twix_filename, scanner] = fileparts(raw_path);
    
    gyro_ratio = 11.777; % gyromagnetic ratio of 129Xe in MHz/Tesla
    
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
                if contains(twix_filename, '2002')
                    excitation = twix.hdr.Phoenix.sWipMemBlock.alFree{1, 4};
                else
                    excitation = twix.hdr.Phoenix.sWipMemBlock.alFree{1, 5};
                end
            end
        end
    
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
    
        % read in mrd file dataset and ismrmrdHeader
        dataset = ismrmrd.Dataset(raw_path, 'dataset');
        ismrmrd_header = ismrmrd.xml.deserialize(dataset.readxml);
    
        % convert user parameter fields to maps for easy query
        general_user_params_long = containers.Map();
        data_struct = ismrmrd_header.userParameters.userParameterLong;
        for i = 1:numel(data_struct)
            general_user_params_long(data_struct(i).name) = data_struct(i).value;
        end
    
        % read in key variables
        vendor = ismrmrd_header.acquisitionSystemInformation.systemVendor;
        dwell_time = double(dataset.readAcquisition().head.sample_time_us(1)) * 1e-6; % in s
        xeFreqMHz = general_user_params_long("xe_center_frequency") * 1e-6; % gas excitation frequency in MHz
        tr = ismrmrd_header.sequenceParameters.TR(2) * 1e-3; % dissolved TR in s
        field_strength = ismrmrd_header.acquisitionSystemInformation.systemFieldStrength_T; % in T
        freq_dis_excitation_hz = general_user_params_long("xe_dissolved_offset_frequency"); % in Hz
    
        % calculate rf excitation in ppm
        rf_excitation = round(freq_dis_excitation_hz/(gyro_ratio * field_strength));
    
        % read k-space data
        npts = size(dataset.readAcquisition(1).data{1},1);
        nfids = dataset.getNumberOfAcquisitions;
        fids_cell = dataset.readAcquisition().data;
        fids = zeros(npts, nfids);
        for i=1:nfids
            fids(:,i) = transpose(double(fids_cell{i}(:,1)));
        end

    else
        error('Unknown Raw File Type')
    end

    if droppt_N > 0
        fids = fids(droppt_N + 1:end, :);
        fids(size(fids, 1):size(fids, 1) + droppt_N, :) = 0;
    end
    % if data from GE scanner, take complex conjugate
    if strcmp(scanner, '.h5') && strcmpi(vendor,'ge')
        fids = conj(fids);
    end
end
