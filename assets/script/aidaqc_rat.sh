#!/bin/bash
#SBATCH --job-name=mouseqc
#SBATCH --nodes=1
#SBATCH --time=48:00:00
#SBATCH --mail-type=FAIL
#SBATCH --partition=batch
#SBATCH --mem=128GB

root_dir="/project/4180000.36/awake"
bids_dir="/bids/rat"
aidaqc_dir="aidaqc_rat"

apptainer exec /groupshare/traaffneu/preclinimg/apptainer/aidaqc.sif ParsingData.py --initial_path ${root_dir}/${bids_dir} --output_path ${root_dir}/${aidaqc_dir} --format nifti
