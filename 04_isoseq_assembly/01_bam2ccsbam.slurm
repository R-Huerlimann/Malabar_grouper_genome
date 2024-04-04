#!/bin/bash
#SBATCH -p compute
#SBATCH -t 4-0
#SBATCH -c 32
#SBATCH --mem=64G
#SBATCH --job-name=01_subread2ccs
#SBATCH --output=01_bam2cssbam.out
#SBATCH --error=01_bam2cssbam.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=roger.huerlimann@oist.jp

#ccs takes subreads (outputs from sequencer) and compiles them into circular consensus sequence reads
#subreads.bam --> ccs.bam

### Environment ###
ml load PB_isoseq3 

### Directories ###
DIRin1=/malabar/01_Raw_data/PB_IsoSeq_Run01
FILEbam1=Emal_ISOseq_Run01.ccs.bam
FILEreport1=Emal_ISOseq_Run01_ccs_report.txt

DIRin2=/malabar/01_Raw_data/PB_IsoSeq_Run02
FILEbam2=Emal_ISOseq_Run02.ccs.bam
FILEreport2=Emal_ISOseq_Run02_ccs_report.txt

DIRout1=/malabar/04_isoseq_assembly/Run01
DIRout2=/malabar/04_isoseq_assembly/Run02

mkdir -p ${DIRout1};mkdir -p ${DIRout2}

### Code ###
ccs ${DIRin1}/*.subreads.bam \
    ${DIRout1}/${FILEbam1} \
	--report-file ${DIRout1}/${FILEreport1} \
    --min-passes 3 \
    --min-rq 0.9 #default = 0.99, min predicted accuracy of a read = the expected percentage of matches in an alignment of the consensus sequence to the true read

ccs ${DIRin2}/*.subreads.bam \
    ${DIRout1}/${FILEbam2}  \
	--report-file ${DIRout2}/${FILEreport2} \
    --min-passes 3 \
    --min-rq 0.9 #default = 0.99, min predicted accuracy of a read = the expected percentage of matches in an alignment of the consensus sequence to the true read

