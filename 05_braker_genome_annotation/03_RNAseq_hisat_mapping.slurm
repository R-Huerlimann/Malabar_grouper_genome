#!/bin/bash
#SBATCH -p compute
#SBATCH -t 2-0
#SBATCH -c 50
#SBATCH --mem=200G
#SBATCH --job-name=03_hisat_mapping
#SBATCH --output=03_RNAseq_hisat_mapping.out
#SBATCH --error=03_RNAseq_hisat_mapping.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=roger.huerlimann@oist.jp

### Environment ###
ml hisat2/2.2.1

### Directories and Files ###
DIRin=/malabar/07_larval_analysis/01_Read_QC #data taken from trimmed and filtered fastq files for larval analysis
FILEfastq=REMlarvae_trim_filt_R1.fastq
FILEsam=Emal_larvae_mapped.sam
DIRout=/malabar/05_braker_genome_annotation

### Code ###
cd ${DIRout}

hisat2  -p 50 \
	      --dta \
	      -x Emal_V1_24Chr_P0_rn \
	      -1 ${DIRin}/${FILEfastq} \
	      -2 ${DIRin}/${FILEfastq/R1/R2} \
	      -S ${FILEsam}

#sorting the sam file
echo "starting sorting"
sort -k 3,3 -k 4,4n ${FILEsam} > ${FILEsam/.sam/sorted.sam}
echo "done sorting"

echo "compressing sam to bam file"
samtools view -S -b ${FILEsam/.sam/sorted.sam} > ${FILEsam/.sam/sorted.bam}
echo "done compressing"
