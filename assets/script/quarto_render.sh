#!/bin/bash
#SBATCH --job-name=render
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-type=FAIL
#SBATCH --partition=batch
#SBATCH --mem=64GB


apptainer exec /groupshare/traaffneu/preclinimg/apptainer/nvim.sif quarto render /home/traaffneu/joagra/code/awake/analysis/dataset_summary.qmd
apptainer exec /groupshare/traaffneu/preclinimg/apptainer/nvim.sif quarto render /home/traaffneu/joagra/code/awake/analysis/dataset_analysis.qmd
