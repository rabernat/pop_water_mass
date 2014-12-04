#!/bin/tcsh
#
# LSF batch script to run an MPI application
#
#BSUB -P UCLB0005 
#BSUB -W 5:00               # wall-clock time (hrs:mins)
#BSUB -n 25                  # number of tasks in job         
#BSUB -R "span[ptile=25]"    
#BSUB -J popclim               # job name
#BSUB -o trans.%J.out        # output file name in which %J is replaced by the job ID
#BSUB -e trans.%J.err        # error file name in which %J is replaced by the job ID
#BSUB -q geyser             # queue
#BSUB -N                    # mail when done

module load python all-python-libs
ipcluster start -n 24 &
python -u output_dens_flux.py
kill %1

