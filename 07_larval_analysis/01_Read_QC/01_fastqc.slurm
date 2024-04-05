#!/bin/bash
#SBATCH -p compute
#SBATCH -t 0-12
#SBATCH -c 20
#SBATCH --mem=20G
#SBATCH --mail-user=roger.huerlimann@oist.jp

### Environment ###
module load fastqc/0.11.9

### Directories and Files ###
DIRinA=/malabar/00_Raw_data/RNAseq_Larvae/ID240_A
DIRinB=/malabar/00_Raw_data/RNAseq_Larvae/ID240_B

DIRoutA=/malabar/larvae/Raw_fastqc/RunA
DIRoutB=/malabar/larvae/Raw_fastqc/RunB

mkdir -p ${DIRoutA}
mkdir -p ${DIRoutB}

### Code ###
fastqc -q -o ${DIRoutA} ${DIRinA}/*gz
fastqc -q -o ${DIRoutB} ${DIRinB}/*gz