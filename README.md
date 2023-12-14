# Spectroscopy_Processing_Production

This MATLAB repository will process our calibration scan and produce a pptx file. User should dowonload/clone this repository and open 'spectroscopy_processing.m' file on MATLAB (MATLAB Current Folder must contain the "spectroscopy_process_functions" folder). In order to read MRD files you will need to download/clone the ismrmrd repository and add it to your MATLAB path: https://github.com/ismrmrd/ismrmrd

To process a subject:

1. Run the 'spectroscopy_processing.m' script.
2. User will be asked to select the calibration twix (.dat) file. Select the file. 

The code will run and create figures that will also pop up in the screen. Afterwards, figures will be closed automatically. Finally, the report will be available in the same directory where the twix file is. 

### NOTE to the Analyst:
Analyst will be responsible to run this script and prepare the spectroscopy report. Before uploading the report on the corresponding slack thread, analyst MUST open the pptx and check all the figures are available and nothing is missing. Sometimes the Static Spectroscopy figure might be missing. In that case, analyst should rerun the code and see the figures. 
