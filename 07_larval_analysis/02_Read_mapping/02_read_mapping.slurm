#!/bin/bash
#SBATCH -p compute
#SBATCH -t 0-12
#SBATCH -c 20
#SBATCH --mem=80G
#SBATCH --job-name=02_STAR_mapping
#SBATCH --output=02_read_mapping-%j.out
#SBATCH --error=02_read_mapping-%j.err
#SBATCH --mail-user=roger.huerlimann@oist.jp
#SBATCH --mail-type=END,FAIL
#SBATCH --array=0-23%24

### Environment ###
ml star/2.7.9a

### Directories and Files ###
DIRgenome=/malabar/07_larval_analysis/02_Read_mapping/01_genome_index
DIRin=/malabar/07_larval_analysis/01_Read_QC/03_kraken_results
DIRout=/malabar/07_larval_analysis/02_Read_mapping/
files=(${DIRin}/*unclassified_1*)
name=${files[${SLURM_ARRAY_TASK_ID}]};name=${name##*/};name=${name%_unclass*}

#File name example
#REMD01-1_S1_unclassified_1.fastq

### Code ###
mkdir -p ${DIRout};cd ${DIRout};mkdir -p ${name};cd ${name}

STAR	--runThreadN 20 \
		--genomeDir ${inDirGen} \
		--quantMode GeneCounts \
		--outSAMtype BAM SortedByCoordinate \
		--outFileNamePrefix ${name} \
		--readFilesIn ${files[${SLURM_ARRAY_TASK_ID}]} ${files[${SLURM_ARRAY_TASK_ID}]/1.fas/2.fas}