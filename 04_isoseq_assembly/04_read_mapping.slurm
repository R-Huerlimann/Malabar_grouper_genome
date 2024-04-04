#!/bin/bash
#SBATCH -p compute
#SBATCH -t 1-0
#SBATCH -c 30
#SBATCH --mem=500G
#SBATCH --job-name=04_mapping
#SBATCH --output=04_read_mapping.out
#SBATCH --error=04_read_mapping.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=roger.huerlimann@oist.jp

### Notes ###
# Mapping isoseq reads to genome to be used in braker2 later

### Environment ###
ml PB_isoseq3
ml minimap2/2.17
ml samtools/1.9

### Directories and Files ###
DIRgenome=/malabar/02_Scaffolding
FILEgenome=Emal_V1_24Chr_P0_rn.fasta

DIRin=/malabar/04_isoseq_assembly
FILEfastq=Emal_all_flnc.fastq
FILEsam=Emal_isoseq_fastq.sam

DIRout=${DIRin}/04_read_mapping

### Code ####
mkdir -p ${DIRout};cd ${DIRout}

#extracting fastq files from bam files and concatenaing into one file
echo "converting first bam file to fastq"
bamtools convert -format fastq -in ${DIRin}/Run01/*flnc.bam > ${FILEfastq}
echo "converting second bam file to fastq"
bamtools convert -format fastq -in ${DIRin}/Run02/*flnc.bam >> ${FILEfastq}
echo "done converting"

#renaming headers in fastq file to prevent future issues
echo "renaming headers"
cat ${FILEfastq} | awk '{print (NR%4 == 1) ? "@" ++i : $0}' > ${FILEfastq/.fastq/rn.fastq}
echo "done renaming"

#mapping to genome
echo "starting minimap"
minimap2 -ax splice \
	 -t 30 \
	 -uf \
	 --secondary=no \
	 -C5 \
	 ${DIRgenome}/${FILEgenome} ${FILEfastq/.fastq/rn.fastq} > ${FILEsam}
echo "done mapping"

#sorting the sam file
echo "starting sorting"
sort -k 3,3 -k 4,4n ${FILEsam} > ${FILEsam/.sam/_sorted.sam}
echo "done sorting"

echo "compressing sam to bam file"
samtools view -S -b ${FILEsam/.sam/_sorted.sam} > ${FILEsam/.sam/_sorted.bam}
echo "done compressing"
