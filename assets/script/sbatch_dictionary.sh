#!/bin/bash
#SBATCH --job-name=dictionary
#SBATCH --nodes=1
#SBATCH --time=72:00:00
#SBATCH --mail-type=FAIL
#SBATCH --partition=batch
#SBATCH --mem=256GB
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=32

apptainer exec /groupshare/traaffneu/preclinimg/apptainer/nvim.sif python3 /home/traaffneu/joagra/code/awake/assets/script/dictionary.py

#### clean up####
 
rm -rf $TMPDIR/*
