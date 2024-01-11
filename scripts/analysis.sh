#!/bin/bash

#cd /home/traaffneu/margal/awake_code/awake/scripts
# ./analysis.sh

# cd "/opt/singularity/3.10.3/bin/singularity"
module load singularity
module load ANTs
module unload ANTs
# module unload freesurfer
# module unload fsl 

cd /project/4180000.36/AwakeRodent/scratch/

bids_folder='01_mice'
subj_num='100100'
ses_num='1'
seed_mask_list=("/project/4180000.36/AwakeRodent/scratch/template/seed_S1-right_mouse.nii.gz" "/project/4180000.36/AwakeRodent/scratch/template/seed_S1-left_mouse.nii.gz")

orig_bids_dir=/project/4180000.36/AwakeRodent/bids/$bids_folder/
BIDS_input=/project/4180000.36/AwakeRodent/scratch/bids/$bids_folder/sub-0100100
preprocess_outputs=/project/4180000.36/AwakeRodent/scratch/RABIES_preprocess/sub-0100100_ses-1/
confound_correction_outputs=$preprocess_outputs/confound_correction_outputs/
analysis_outputs=$preprocess_outputs/analysis_outputs

singularity run -B $BIDS_input:/BIDS_input:ro \
                     -B $preprocess_outputs:/preprocess_outputs/ \
                     -B $confound_correction_outputs:/confound_correction_outputs/ \
                     -B $analysis_outputs:/analysis_outputs/ \
                     /opt/rabies/0.5.0/rabies.sif -p MultiProc analysis /confound_correction_outputs/ /analysis_outputs/ \
                     --seed_list "${seed_mask_list[@]}"


