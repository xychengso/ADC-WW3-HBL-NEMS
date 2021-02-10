#!/bin/bash --login
### This script is used to launch the standalone WW3 simulation  ###
### From the regtests.  X.Chen                                   ###

#SBATCH --constraint=hatteras
#SBATCH --job-name=jnsem_test_atm2ocn
#SBATCH -p batch
#SBATCH --time=12:00:00
#SBATCH -n 24
#SBATCH --mail-user=xychen@ht3.renci.org
#SBATCH --mail-type=ALL
#SBATCH --output=ww3_Inlet_test.out.log
#SBATCH --error=ww3_Inlet_test.err.log

REGT_DIR=/home/xychen/ADC-WW3-NWM-NEMS/WW3/regtests
SRC_DIR=/home/xychen/ADC-WW3-NWM-NEMS/WW3/model
INPT_DIR=ww3_tp2.17_shell
cd ${REGT_DIR}
 ./bin/run_test -c Intel -S -s MPI -s NO_PDLIB -w work_a -g a -f -p    \
       mpirun -n 24 ${SRC_DIR} ${INPT_DIR}

