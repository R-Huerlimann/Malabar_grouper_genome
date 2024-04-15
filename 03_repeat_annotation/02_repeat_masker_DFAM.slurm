#!/bin/bash -l
#SBATCH -p compute
#SBATCH -t 4-0
#SBATCH -c 40
#SBATCH --mem=150G
#SBATCH --job-name=02_repmask_DFAM
#SBATCH --output=02_repeat_masker_DFAM.out
#SBATCH --error=02_repeat_masker_DFAM.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=roger.huerlimann@oist.jp

### Environment ###
ml repeatMasker/4.1.1

### Directories and Files ###
DIRgenome=/malabar/02_HiC_scaffolded_assembly
FILEgenome=Emal_final_scaffolded_phase1.fasta

DIRout=/malabar/03_repeat_annotation/02_repeat_masker_DFAM

mkdir -p ${DIRout};cd ${DIRout}

### Code ###

#RepeatMasker with DFAM (~ 2 days).
RepeatMasker -pa 40 -gff -dir . -species Vertebrata ${DIRgenome}/${FILEgenome}
