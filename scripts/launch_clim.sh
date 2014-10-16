#!/bin/tcsh
#
# LSF batch script to run an MPI application
#
#BSUB -P UCLB0005
#BSUB -W 04:00               # wall-clock time (hrs:mins)
#BSUB -n 4                  # number of tasks in job         
#BSUB -J popclim               # job name
#BSUB -o popclim.%J.out        # output file name in which %J is replaced by the job ID
#BSUB -e popclim.%J.err        # error file name in which %J is replaced by the job ID
#BSUB -q caldera             # queue


module load nco
python climatologies.py
