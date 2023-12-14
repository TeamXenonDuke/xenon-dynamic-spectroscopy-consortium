function siftedFIDs = SIFT(fids,bw,tr)
% Assumes the last FID is a gas excitation

display_flag = 0;
% bw = 1000*pfile.rdb.rdb_hdr_user12;   % Receiver bandwidth (kHz)
% tr = pfile.image.tr*1E-6;

% Create arrays of sample times and frequencies
npts = size(fids,1);                    % Number of samples
dwell_time = 1/(2*bw);                  % Time between each sample
t = dwell_time*(0:(npts-1))';           % Direct time domain vector
nFrames = size(fids,2)-1;               % Number of fids
t_tr = tr*((1:nFrames)-1);              % Indirect time domain vector
f = linspace(-0.5,0.5,nFrames)/tr;      % Indirect frequency vector


% Take indirect Fourier transform of FIDs
dis_pfile = fids(:,1:end-1); % change pfile to fid
sift_fft = fftshift(fft(dis_pfile,[],2),2);
thresh_sift_fft = sift_fft;

% Find baseline noise for each indirect frequency. Set threshold for each
% indirect frequncy to be 2*std of noise + the average value of the noise.
noise = sift_fft(:,end-round(size(sift_fft,2)/6):end); % 
thresh = 2*std(abs(noise'))+mean(abs(noise'));

% Set all points below the threshold in indirect frequency spectrum to zero
for k = 1:npts
    change = abs(sift_fft(k,:)) < thresh(k);
    thresh_sift_fft(k,change) = 0;
end 

% Set any points in indirect frequency spectrum that are isolated spikes
% (don't have any components on either side) to zero.
ts = find(abs(thresh_sift_fft) > 0);
ts(ts > length(thresh_sift_fft)) = [];
change2 = thresh_sift_fft(ts+1) == 0;
thresh_sift_fft(ts(change2)) = 0;

% Take inverse Fourier transform with respect to indirect time domain
siftedFIDs = ifft(ifftshift(thresh_sift_fft,2),'',2);

% Compare SIFTed FIDs to orginal and calculate residuals. Replace any
% sifted FIDs with residuals greater than 2*std(residual)
resid = real(fids(:,1:end-1)-siftedFIDs);
ch = (sum(abs(resid)) >= (2*std(resid(:)))*size(resid,2));
siftedFIDs(ch) = dis_pfile(ch);
siftedFIDs(end-10:end,:) = siftedFIDs(end-20:end-10,:);

siftedFIDs(:,end+1) = fids(:,end);

if display_flag == 1
    figure(1); surf(abs(fids)), shading interp
    figure(2); surf(abs(sift_fft)), shading interp
    figure(3); surf(abs(thresh_sift_fft)), shading interp
    figure(4); surf(abs(siftedFIDs(:,1:end-1))), shading interp
    figure(5); surf((real(resid))), shading interp
end 

end
