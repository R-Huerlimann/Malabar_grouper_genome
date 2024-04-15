#!/bin/bash
#SBATCH -p compute
#SBATCH -t 0-2
#SBATCH -c 128
#SBATCH --mem=300G
#SBATCH --job-name=IPA_phased
#SBATCH --output=01_assembly_IPA_phased.out
#SBATCH --error=01_assembly_IPA_phased.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=roger.huerlimann@oist.jp

### Environment ###
module load PB_IPA/20200928

### Directories ###
DIRin1=/malabar/00_raw_data/PB_Run01
DIRin2=/malabar/00_raw_data/PB_Run02
DIRin3=/malabar/00_raw_data/PB_Run03
DIRout=/malabar/01_assembly/

echo "creating directory"
mkdir -p ${DIRout}

### Code ###

echo "starting assembly"
ipa local       -i ${DIRin1}/m*.ccs.bam \
				-i ${DIRin2}/m*.ccs.bam \
				-i ${DIRin3}/m*.ccs.bam \
				--nthreads 32 \
                --tmp-dir ${DIRout}/tempP \
                --verbose \
                --njobs 4 \
                --run-dir ${DIRout}
