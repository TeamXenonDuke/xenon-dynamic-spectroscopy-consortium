function b = highpassfilter(lengthdata)
%TEMP Returns a discrete-time filter object. With order N.

% MATLAB Code
% Generated by MATLAB(R) 9.2 and the Signal Processing Toolbox 7.4.
% Generated on: 19-Jul-2017 11:45:18

% FIR Window Highpass filter designed using the FIR1 function.

% All frequency values are in Hz.
Fs = 50;  % Sampling Frequency - look at this

N     = 2*floor(lengthdata/2/3)-2;      % Order
Fc    = 0.5;      % Cutoff Frequency
flag  = 'scale';  % Sampling Flag
Alpha = 2.5;      % Window Parameter

% Create the window vector for the design algorithm.
win = gausswin(N+1, Alpha);

% Calculate the coefficients using the FIR1 function.
b  = fir1(N, Fc/(Fs/2), 'high', win, flag);
Hd = dfilt.dffir(b);

