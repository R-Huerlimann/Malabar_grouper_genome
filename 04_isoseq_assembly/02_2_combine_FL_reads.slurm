#!/bin/bash
#SBATCH -p compute
#SBATCH -t 0-20
#SBATCH -c 1
#SBATCH --mem=10G
#SBATCH --job-name=02_combine_FL
#SBATCH --output=02_2_combine_FL_reads.out
#SBATCH --error=02_2_combine_FL_reads.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=roger.huerlimann@oist.jp

#When doing demultiplexing with lima, if there are multiple 5p/3p primers then one 5p--3p.bam, 5p--3p.bam.pbi and 5p--3p.xml  output will be produced for each primer pair
#.xml file has the path to the .bam and the .bam.pbi files so keep in the same location they were created
#Combine all of the .bam files into one dataset so all the demultuplexed full-length reads are combine to increase transcript recovery
#runs fast (~ 5mins)

### Environment ###
ml PB_isoseq3

### Directories ###
DIRin1=/malabar/04_isoseq_assembly/Run01
DIRin2=/malabar/04_isoseq_assembly/Run02

### Code ###
cd ${inDir1}
dataset create --type ConsensusReadSet Run01_demult_combined_fl.xml *3p.bam

cd ${inDir2}
dataset create --type ConsensusReadSet Run02_demult_combined_fl.xml *3p.bam
