#!/bin/bash
#SBATCH -p compute
#SBATCH -t 0-2
#SBATCH -c 50
#SBATCH --mem=100G
#SBATCH --job-name=02_hisat_indexing
#SBATCH --output=02_RNAseq_hisat_indexing.out
#SBATCH --error=02_RNAseq_hisat_indexing.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=roger.huerlimann@oist.jp

### Environment ###
ml hisat2/2.2.1

### Directories and Files ###
DIRin=/malabar/02_HiC_scaffolded_assembly
FILEgenome=Emal_V1_24Chr_P0_rn.sm.fasta
DIRout=/malabar/05_braker_genome_annotation

### Code ###
mkdir -p ${DIRout};cd ${DIRout}

ln -s ${DIRin}/${FILEgenome} ${FILEgenome}

hisat2-build -p 50 \
	     ${FILEgenome} \
	     ${FILEgenome/sm.fasta/sm}
