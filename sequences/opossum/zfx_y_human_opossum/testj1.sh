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
bash all_bz.log
tba '(humanmasked_X humanmasked_Y (opossum))' *.*.maf tba.maf >&tba.log
maf_project tba.maf humanmasked_X '(humanmasked_X humanmasked_Y (opossum))' > human_proj.maf
msa_view -o FASTA human_proj.maf > ZFX_msa.fa