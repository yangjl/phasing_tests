#!/bin/bash -l
#SBATCH -D /home/jolyang/Documents/Github/phasing_tests
#SBATCH -J phase_mom
#SBATCH -o /home/jolyang/Documents/Github/phasing_tests/slurm/phase_out-%j.txt
#SBATCH -e /home/jolyang/Documents/Github/phasing_tests/slurm/phase_error-%j.txt


#$1: 1-160; row number of the all_files.txt
row=$1
R --no-save "--args $row $SLURM_JOB_ID" < profiling/2.cj_data/2.C.2_cj_phasing.R

#usage: sbatch -p serial slurm-script/run_phasing.sh 1
