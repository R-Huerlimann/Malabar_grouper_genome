#!/bin/bash
#SBATCH -p compute
#SBATCH -t 0-12
#SBATCH -c 24
#SBATCH --mem=96G
#SBATCH --job-name=04_kraken2
#SBATCH --error=04_kraken2filt-%j.err
#SBATCH --mail-user=roger.huerlimann@oist.jp
#SBATCH --mail-type=END,FAIL
#SBATCH --array=0-23%24

### Environment ###
ml kraken2

### Directories ###
files=(/malabar/07_larval_analysis/02_Trimmed_reads*val_1*gz)
DIRout=/malabar/07_larval_analysis/03_kraken_results
db=/databases/kraken2_bacterial_archaeal/4.08.20/

name=$(echo ${files[${SLURM_ARRAY_TASK_ID}]} | sed 's|.*02_Trimmed_reads/||' | sed 's|_val.*||')


### Code ###

kraken2 --db ${db} \
		--paired ${files[${SLURM_ARRAY_TASK_ID}]} ${files[${SLURM_ARRAY_TASK_ID}]/val_1/val_2} \
		--gzip-compressed \
		--confidence 0.3 \
		--threads 24 \
		--use-names \
		--classified-out ${DIRout}/${name}_classified#.fastq \
		--unclassified-out ${DIRout}/${name}_unclassified#.fastq \
		--report ${DIRout}/${name}_kraken_report
