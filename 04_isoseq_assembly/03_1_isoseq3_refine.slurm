#!/bin/bash
#SBATCH -p compute
#SBATCH -t 0-2
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH --job-name=03_refine
#SBATCH --output=03_1_isoseq3_refine.out
#SBATCH --error=03_1_isoseq3_refine.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=roger.huerlimann@oist.jp

#Use isoseq3 to REFINE - takes full length reads (combined dataset of all FL reads) and removes concatemers and polyA tails (if present) 
#to produce full-length non-concatemer reads (FLNC)
#runs fast (< 30 mins)
#demult.combined.consensusreadset.xml + primers.fasta --> flnc.bam

### Environment ###
ml PB_isoseq3

### Directories ###
DIRin1=/malabar/03_isoseq_assembly/Run01
DIRin2=/malabar/03_isoseq_assembly/Run02

FILExml1=Run01_demult_combined_fl.xml
FILExml2=Run02_demult_combined_fl.xml

indices=/malabar/00_Raw_data/Iso-Seq_96_primer_set.fasta

### Code ###
cd ${DIRin1}
isoseq3 refine --require-polya \
	${FILExml1} \
	${indices} \
	${FILExml1/fl.xlm/flnc.bam} \
	--log-level INFO --log-file ${FILExml1/fl.xlm/flnc.log}

cd ${DIRin2}
isoseq3 refine --require-polya \
        ${FILExml2} \
        ${indices} \
        ${FILExml2/fl.xlm/flnc.bam} \
        --log-level INFO --log-file ${FILExml2/fl.xlm/flnc.bam}

#--require-ploya		  use if transcripts have a polyA tail- filters for FL reads that have a poly(A) tail with at least 20 base pairs and removes identified tail.
#--log-level INFO 		log file (name specified with --log-file) will show you how many full length (fl) reads were input and how many full length non-concatemer (flnc) reads were output. 

#Can convert to fasta format:
#bamtools convert -format fastq -in flnc.bam > flnc.fastq

