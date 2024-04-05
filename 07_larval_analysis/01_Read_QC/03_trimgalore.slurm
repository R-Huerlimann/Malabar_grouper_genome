#!/bin/bash
#SBATCH -p compute
#SBATCH -t 0-12
#SBATCH -c 4
#SBATCH --mem=40G
#SBATCH --job-name=03_trimgalore
#SBATCH --output=03_trimgalore-%j.out
#SBATCH --error=03_trimgalore-%j.err
#SBATCH --mail-user=roger.huerlimann@oist.jp
#SBATCH --mail-type=END,FAIL
#SBATCH --array=0-23%24

### Environment ###
module load cutadapt
module load trimgalore

### Directories and Files ###
files=(/malabar/00_Raw_data/RNAseq_Larvae/ID240_AB/*R1_*gz)
DIRout=/malabar/07_larval_analysis/01_Read_QC/02_Trimmed_reads
name=${files[${SLURM_ARRAY_TASK_ID}]};name=${name##*/};name=${name%_R1*}

### Code ###
mkdir -p ${DIRout}

trim_galore \
	--paired \
	--retain_unpaired \
	--basename ${name} \
	--fastqc \
	--length 30 \
	--o ${outDir} \
	--cores 4 \
	${files[${SLURM_ARRAY_TASK_ID}]} ${files[${SLURM_ARRAY_TASK_ID}]/R1_/R2_}