#!/bin/bash
#SBATCH --partition=compute
#SBATCH -t 12:00:00
#SBATCH --mem=100G
#SBATCH -c 10
#SBATCH --job-name=02_iprscan
#SBATCH --output=02_run_iprscan-%j.out
#SBATCH --error=02_run_iprscan-%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=roger.huerlimann@oist.jp
#SBATCH --array=0-26%12

### Environment ###
ml interproscan/5.48.83.0

### Directories and Files ###
DIRin=/malabar/06_gene_models_functional_annotation
files=(*fa)
name=${files[${SLURM_ARRAY_TASK_ID}]};name=${name%.fa}

#File name example: malabar0.fa

### Code ###
cd ${DIRin}

interproscan.sh -i ${files[${SLURM_ARRAY_TASK_ID}]} --disable-precalc --goterms --pathways -f xml --tempdir ${name}_temp > ${name}.log
