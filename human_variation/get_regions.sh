#!/bin/bash
#SBATCH --job-name=all_bz
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=20gb
#SBATCH --partition=20
##SBATCH --output all_bz-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hlharris@wi.mit.edu

python get_regions.py
