#!/bin/bash
#SBATCH --partition=compute
#SBATCH -t 95:00:00
#SBATCH --mem=400G
#SBATCH --job-name=06_diamond
#SBATCH --output=06_run_diamond-%j.out
#SBATCH --error=06_run_diamond-%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=roger.huerlimann@oist.jp

### Environment ###
ml diamond/2.0.9
ml ncbi-blast/2.10.0+

### Directories and Files ###
DIRin=/malabar/05_braker_annotation/emal_braker_RNA_PROT_all
FILEaa=augustus.hints.aa

### Code ###
# use swissprot to quickly assign function
update_blastdb.pl --decompress --blastdb_version 5 swissprot

diamond blastp -d swissprot -q ${FILEaa} -o diamond_matches.tsv
