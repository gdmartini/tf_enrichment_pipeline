#!/bin/bash
# Job name:
#SBATCH --job-name=savio_test
#
# Partition:
#SBATCH --partition=savio2
#
# QoS:
#SBATCH --qos=savio_normal
#
# Account:
#SBATCH --account=fc_mhmm
#
# Request one node:
#SBATCH --nodes=1
#
# Processors per task:
#SBATCH --cpus-per-task=1
#
# Number of Processors per Node:
#SBATCH --ntasks-per-node=24
#
# Wall clock limit:
#SBATCH --time=72:00:00
#
#SBATCH -a 1
## Command(s) to run:
module load matlab
module load hdf5
# Make a temporary scratch directory for storing job
# and task information, to coordinate parallelization.
# This directory can then be referenced by assigning it to
# a 'parcluster.JobStorageLocation' property in your script.
mkdir -p /global/scratch/$USER/$SLURM_JOB_ID
matlab -nodisplay -nodesktop < main05_conduct_hmm_inference_savio.m