#!/bin/env bash

dir=/project/4180000.36/awake/mouse_scripts/

file="/project/4180000.36/awake/mouse_scripts/sub-0100100_ses-1_run-1_task-rest_bold.sh"

#remove module load 
sed -i 's/module load rabies//g' $file
#update rabies version
sed -i 's/0.5.1/0.6.0/g' $file
#remove --label in preprocess
sed -i 's/--labels \/home\/traaffneu\/joagra\/code\/awake\/assets\/template\/mouse\/labels.nii.gz //g' $file

#replace --conf_list with --nuisance_regressors in confound_correction
sed -i 's/--conf_list/--nuisance_regressors/g' $file


#add fc and dr analysis
sed -i 's/vpm_r.nii.gz/vpm_r.nii --ROI_labels_file \/home\/traaffneu\/joagra\/code\/awake\/assets\/template\/mouse\/labels.nii.gz --FC_matrix --prior_maps \/home\/traaffneu\/joagra\/code\/awake\/assets\/template\/mouse\/ica.nii.gz --DR_ICA --prior_bold_idx 1 2 --prior_confound_idx 3 4/ g' $file

#make sure we get also clean time series
sed -i 's/frame_censoring_mask//g' $file
