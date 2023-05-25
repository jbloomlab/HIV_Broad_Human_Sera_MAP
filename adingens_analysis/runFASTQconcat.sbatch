#!/bin/sh
#SBATCH
#PBS -l walltime=24:00:00
jupyter nbconvert --to notebook --execute --inplace --ExecutePreprocessor.timeout=-1 Concat_FASTQ_manually_200327HiSeq.ipynb