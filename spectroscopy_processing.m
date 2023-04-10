clc;
clear;
close all;

% Ensure spectroscopy processing functions are in the matlab path
% Find the directory with the script
[filepath, name, ext] = fileparts(which('spectroscopy_processing.m'));
project_dir = fullfile(filepath, 'spectroscopy_process_functions');
addpath(project_dir); % add the functions to the matlab path. Ok if already there.

% Prompt user to select file for processing from current starting directory
disp("Select Calibration Twix File (.dat) or ismrmd file (.h5) ...");

[file, path] = uigetfile('*.*', 'Select file');
file_with_path = fullfile(path, file);
output_dir = fullfile(path, 'Spectroscopy');
% Create the Spectroscopy folder
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end
% Fit dynamic spectra -- this should create dynV.mat file
dyn = fitDynamicSpec(file_with_path, 'dynV');
% Get the mat file path
dyn_filepath = fullfile(path, 'Spectroscopy', 'dynV.mat');
% Get the subject ID from folder containing twix file
if contains(path, '\') % windows path
    path_cell = strsplit(path, '\');
    subject_id = path_cell{end-1};
elseif contains(path, '/') % linux/macos path
    path_cell = strsplit(path, '/');
    subject_id = path_cell{end-1};
else
    subject_id = path;
end

%% Need RF Excitation to update on the report
[fids, dwell_time, npts, tr, xeFreqMHz, rf_excitation] = readRawDyn(file_with_path);

%% Preparing the ppt summary

% Turn off figure visibility to suppress pop ups
set(groot, 'defaultFigureVisible', 'off');

% we currently do a 12s breath hold. Ignore the first 2s (100 frames)
% due to downstream magnetization
BHs = [2 7]; 

%set figure saving flag
save_figs = 1;

%set axis limits on rbc oscillation plots. if zero, auto limits used
rbc_axis_lim = 0; %limits are set as +-rbc_axis_lim (percent)

% Dynamic Summary after fitting Sine and Peaks
[amp_all, nmrFit, nmrFit_ppm, dyn, detrend_sine, fitted_sine, detrend_peaks, fitted_peaks] = ...
    dynamicSummary(file_with_path, dyn_filepath, BHs, save_figs, rbc_axis_lim);
dates = getDynDatesforPPT(file_with_path, dyn_filepath); % gives scan date, dynamic fit date, etc.

% Sine and Peaks fitting respectively, the figures will be fetechted
% accordingly
imgName{1} = 'sine';
imgName{2} = 'peaks';

%% Creating the pptx, should have two slides with sine and peaks fitted data
dynamicSummaryPPT(subject_id, output_dir, imgName, dates, amp_all, nmrFit_ppm, rf_excitation);
dynamicSummaryCSV(subject_id, output_dir, dates, amp_all, nmrFit_ppm, rf_excitation);
% Changing to the home directory
% cd(current_dir)

% close all

% Turn default figure visibility back on
set(groot, 'defaultFigureVisible', 'on');
close all;