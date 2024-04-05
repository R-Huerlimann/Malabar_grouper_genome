#!/bin/bash
#SBATCH --partition=compute
#SBATCH -t 12:00:00
#SBATCH --mem=50G
#SBATCH -c 10
#SBATCH --job-name=03_blast
#SBATCH --output=03_run_blast-%j.out
#SBATCH --error=03_run_blast-%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=roger.huerlimann@oist.jp
#SBATCH --array=0-27%28

### Environment ###
ml ncbi-blast/2.10.0+

### Directories and Files ###
DIRin=/malabar/06_gene_models_functional_annotation
files=(*fa)
name=${files[${SLURM_ARRAY_TASK_ID}]};name=${name%.fa}
database=/blast_NR

### Code ###
cd ${outDir}

blastp 	-db ${database}/nr \
	-query ${files[${SLURM_ARRAY_TASK_ID}]} \
	-out 03_blast4omics_${name}_fmt14.xml \
	-show_gis \
	-num_threads 10 \
	-evalue 1e-5 \
	-word_size 3 \
	-num_alignments 20 \
	-outfmt 14 \
	-max_hsps 20
