#!/bin/bash
#SBATCH
#PBS -l walltime=48:00:00
jupyter nbconvert --to notebook --execute --inplace --ExecutePreprocessor.timeout=-1 analysis_notebook.ipynb