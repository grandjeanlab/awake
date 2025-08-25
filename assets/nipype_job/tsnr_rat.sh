#!/bin/bash
#SBATCH --job-name=nipype
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=32
#SBATCH --time=01:00:00
#SBATCH --mail-type=FAIL
#SBATCH --partition=batch
#SBATCH --mem=64GB
 
source /home/traaffneu/joagra/python/nipype/.venv/bin/activate
 
python /home/traaffneu/joagra/code/awake/assets/nipype_job/tsnr.py rat
