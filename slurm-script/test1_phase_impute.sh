#!/bin/bash -l
#SBATCH -D /home/jolyang/Documents/Github/phasing_tests
#SBATCH -J phase_test
#SBATCH -o /home/jolyang/Documents/Github/phasing_tests/slurm/phase_out-%j.txt
#SBATCH -e /home/jolyang/Documents/Github/phasing_tests/slurm/phase_error-%j.txt
#SBATCH --array=1-100

#only argument is mean (of poisson) of number of crossovers
#crossovers=$1
#R --no-save "--args $SLURM_JOB_ID" < profiling/1.test_codes/1.A.2_test_phase_mom.R
R --no-save "--args $SLURM_JOB_ID" < profiling/1.test_codes/1.C.1_imputeKids_test.R

#sbatch --p serial slurm/testrun.sh