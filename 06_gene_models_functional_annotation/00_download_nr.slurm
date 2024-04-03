#!/bin/bash
#SBATCH --partition=compute
#SBATCH -t 144:00:00
#SBATCH --mem=10G
#SBATCH -c 5
#SBATCH --job-name=00_download_nr
#SBATCH --output=00_download_nr.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=roger.huerlimann@oist.jp

### Environment ###
ml ncbi-blast/2.10.0+

### Data and folders ###
DIRout=/blast_NR

### Code ###
mkdir -p ${outDir};cd ${outDir}
update_blastdb.pl --verbose nr
