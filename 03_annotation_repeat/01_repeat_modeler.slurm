#!/bin/bash
#SBATCH --partition=compute
#SBATCH -t 72:00:00
#SBATCH --mem=120G
#SBATCH -c 24
#SBATCH --job-name=01_repmod
#SBATCH --output=01_repeat_modeler.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=roger.huerlimann@oist.jp

### Load modules ###
module load repeatmodeler/1.6

### Paths and variables###
species=malabar_grouper

DIRgenome=/malabar/02_HiC_scaffolded_assembly
FILEgenome=Emal_final_scaffolded_phase1.fasta

DIRout=/malabar/03_repeat_annotation/01_repeat_modeler

mkdir -p ${DIRout};cd ${DIRout}

### Code ###
BuildDatabase -name ${species} -engine ncbi ${DIRgenome}/${FILEgenome}

RepeatModeler -engine ncbi -pa 24 -LTRStruct -database ${species}
