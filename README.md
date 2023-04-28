# OPMEGfeedback
Analysis scripts to reproduce the figures from *Real-time, model-based magnetic field correction for moving, wearable MEG*. 

## Dependencies
- OPM: https://github.com/tierneytim/OPM (although need version of spm_opm_epoch_trigger from https://github.com/stephaniemellor/OPM)
- SPM12: https://github.com/spm/spm12
- FIL-OPMEG/optitrack: https://github.com/FIL-OPMEG/optitrack
- ICP: https://uk.mathworks.com/matlabcentral/fileexchange/27804-iterative-closest-point
- linspecer: https://uk.mathworks.com/matlabcentral/fileexchange/42673-beautiful-and-distinguishable-line-colors-colormap (not necessary, changes colours)
- export_fig: https://github.com/altmany/export_fig (not necessary, used for saving figures)
- BrewerMap: https://github.com/DrosteEffect/BrewerMap (not necessary, used for colourmap)

## Main scripts
1. Background noise test, background_noise_feedback_tests_analysis.m
2. Auditory evoked response, moving_MMN_analysis.m
3. External coil experiments, coilSweepAnalysis.m
