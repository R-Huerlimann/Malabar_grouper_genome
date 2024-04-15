#!/bin/bash -l
#SBATCH -p compute
#SBATCH -t 1-0
#SBATCH -c 40
#SBATCH --mem=100G
#SBATCH --job-name=03_repmask_RepMod
#SBATCH --output=03_repeat_masker_RepMod.out
#SBATCH --error=03_repeat_masker_RepMod.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=roger.huerlimann@oist.jp

### Environment ###
ml repeatMasker/4.1.1

### Directories and Files ###
DIRgenome=/malabar/02_HiC_scaffolded_assembly
FILEgenome=Emal_final_scaffolded_phase1.fasta

DIRrepeats=/malabar/03_repeat_annotation/01_repeat_modeler
FILErepeats=Emal_V1_P01_repeatfamilies.fa

DIRout=/malabar/03_repeat_annotation/02_repeat_masker_RepMod

mkdir -p ${DIRout};cd ${DIRout}

### Code ###

#RepeatMasker with RepeatModeller output.
RepeatMasker -pa 40 -gff -dir . -lib ${DIRrepeats}/${FILErepeats} ${DIRgenome}/${FILEgenome}
