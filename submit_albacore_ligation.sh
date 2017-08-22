#!/bin/bash
# specify number of nodes and cores per node
#PBS -l nodes=1:ppn=8
# specify the time you expect the job to run hh:mm:ss
#PBS -l walltime=12:00:00
#specify the amount of memory needed
#PBS -l mem=16G
# output and error files
#PBS -o ligation_bc.o$PBS_JOBID
#PBS -e ligation_bc.e$PBS_JOBID

# load paths
source ~/.bashrc
source activate py35

# move to current working directory
cd $PBS_O_WORKDIR


# call albacore
echo "read_fast5_basecaller.py -f $1 -k $2 -i $3 -r -s $4  -t 8 -o fast5,fastq"
read_fast5_basecaller.py -f $1 -k $2 -i $3 -r -s $4 -t 8 -o fast5,fastq

