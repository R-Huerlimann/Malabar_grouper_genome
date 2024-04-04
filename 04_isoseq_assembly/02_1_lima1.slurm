#!/bin/bash
#SBATCH -p compute
#SBATCH -t 2-0
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH --job-name=02_lima
#SBATCH --output=02_1_lima_01.out
#SBATCH --error=02_1_lima_01.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=roger.huerlimann@oist.jp

#DEMULTIPLEX
#lima takes the ccs reads and identifies the full-length (FL) reads (have both 5' and 3' primers present) and removes the primers to produce FL reads.
#If have barcodes also removes these
#Orients sequences in 5' --> 3' direction
#runs really fast (~5 mins)
#ccs.bam (plus fasta file of primers, and barcodes if applicable) --> demult.bam

### Environment ###
ml PB_isoseq3

### Directories ###
DIRin=/malabar/04_isoseq_assembly/Run01
DIRout=/malabar/04_isoseq_assembly/Run01
indices=/malabar/00_Raw_data/Iso-Seq_96_primer_set.fasta

### Code ###
lima 	--isoseq \
		--dump-clips \
		--dump-removed \
		--peek-guess \
		-j 16 \
		${DIRin}/Emal_ISOseq_Run01.ccs.bam \
		${indices} \
		${DIRout}Emal_ISOseq_Run01.demult.bam 

#--dump-clips	  clipped barcode regions are stored in the file prefix.lima.clips
#--dump-removed   records that did not pass provided thresholds or are without barcodes, are stored in the file prefix.removed.bam
#--no-pbi		  omit the pbi index file, to minimize time to result (only needed if importing to SMRT Link) - not working
#--peek-guess	  removes spurious/false positive matches - only required if supply multiple primer pairs
#-j 24 			  specifies no. of threads

#provide fasta file of primers - index_sequences.fa provided by sequencing service
#e.g. >primer_5p
#     AAGCAGTGGTATCAACGCAGAGTACATGGGG
#     >sample_brain_3p
#     AAGCAGTGGTATCAACGCAGAGTACCACATATCAGAGTGCG
#     >sample_liver_3p
#     AAGCAGTGGTATCAACGCAGAGTACACACACAGACTGTGAG

#If created multiple ccs.bam files can run these in parallel
