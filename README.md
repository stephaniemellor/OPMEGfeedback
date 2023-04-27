# OPMEGfeedback
Analysis scripts to reproduce the results from *Real-time, model-based magnetic field correction for moving, wearable MEG*. 
The corresponding data can be found here: 

## Dependencies
- OPM: https://github.com/stephaniemellor/OPM
- SPM12: https://github.com/spm/spm12
- FIL-OPMEG/optitrack: https://github.com/FIL-OPMEG/optitrack
- ICP: https://uk.mathworks.com/matlabcentral/fileexchange/27804-iterative-closest-point
- linspecer: https://uk.mathworks.com/matlabcentral/fileexchange/42673-beautiful-and-distinguishable-line-colors-colormap (not necessary, only changes colours)
- export_fig: https://github.com/altmany/export_fig (not necessary, only used for saving figures)
- BrewerMap: https://github.com/DrosteEffect/BrewerMap

## Main scripts
1. Background noise test, background_noise_feedback_tests_analysis.m
2. Auditory evoked response, moving_MMN_analysis.m
3. External coil experiments, coilSweepAnalysis.m
