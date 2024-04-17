# The transcriptional landscape underlying larval development and metamorphosis in the Malabar grouper (Epinephelus malabaricus)
[![DOI](https://zenodo.org/badge/777631592.svg)](https://zenodo.org/doi/10.5281/zenodo.10972118)
## Abstract
Teleost fishes typically have a bi-phasic life cycle with a transition between larval and juvenile phases called metamorphosis, which is known to be regulated by thyroid hormones (TH). However, other hormonal systems might be involved as it is the case in amphibians in which corticosteroids are interacting with TH pathways to trigger and regulate metamorphosis. Unfortunately, such interplay is poorly understood in teleost fishes. In order to investigate the potential involvement of these two hormonal pathways, we used the Malabar grouper (Epinephelus malabaricus) as a model system. We assembled a chromosome-scale genome and conducted a transcriptomic analysis of nine larval developmental stages. We studied the expression patterns of genes involved in TH and corticoid pathways, as well as four biological processes known to be regulated by TH in other teleost species: ossification, pigmentation, visual perception, and metabolism. Surprisingly, we observed an activation of many of the same pathways involved in metamorphosis at an earlier stage, suggesting an additional implication of these pathways in early larval development. Overall, our data reveal that on a common background (TH controlling metamorphosis) evolution is assembling species-specific peculiarities that allow to precisely align the molecular completion of metamorphosis with the ecological constraints.

Published in eLife: https://doi.org/10.7554/eLife.94573.1


This githup repository contains all scripts associated with the sequencing, annotation and assembly of the Malabar grouper genome, as well as the analysis of the developmental transcriptome from day 1 to juvenile fish (around day 60).

## Genome assembly
The genome assembly was carried out using unprocessed PacBio HiFi reads with the diploid aware Improved Phased Assembler (IPA, V1.3.1, https://github.com/PacificBiosciences/pbipa) using default parameters, which resulted in a primary and alternative phase genome. The two phased genomes were assessed using purge_haplotigs using default parameters to generate a genome-wide read-depth histogram; however, no purging was necessary. Completeness of the final assembly was assessed using BUSCO63 (V4.1.2) with the actinopterygii_odb10 database. Scaffolding and phasing were outsourced to Phase Genomics (See [Phase Genomics HiC scaffolding protocol](02_HiC_scaffolding/PhaseGenomics_protocol.txt))


## Genome annotation
https://github.com/R-Huerlimann/Malabar_grouper_genome/tree/main/05_braker_genome_annotation
Genome annotation was carried out as described Ryu et al (2022). Briefly, repeat content analysis was done in RepeatModeler (V2.0.1), RepeatMasker (V4.1.1), the vertebrata library of Dfam (V3.3), and GenomeTools (V1.6.1). The genome annotation was done using BRAKER2 and associated programs. For this, the ISO-seq data from the adult tissue and RNA-seq data from the larval samples (see below for quality control process) were used together with publicly available protein data. Post-processing was carried out as described by Ryu et al. (2022), also see below, using the Swiss-Prot protein database82 (UniProt) with Diamond74 (V2.0.9) and Pfam domains83 identified by InterProScan84 (V5.48.83.0). Gene model statistics were calculated using the get_general_stats.pl script from the eval package85 (V2.2.8). Finally, functional annotation was carried out with the filtered gene models produced by BRAKER. The amino acid sequences were blasted against the non-redundant protein database (downloaded 15. November 2021) using blastp86 (V2.10.0+; parameters: -show_gis -num_threads 10 -evalue 1e-5 -word_size 3 -num_alignments 20 -outfmt 14 -max_hsps 20). Additionally, protein domains were assigned using InterProScan84 (V5.48.83.0; parameters: --disable-precalc --goterms --pathways -f xml). The blast and interproscan results were then loaded into OmicsBox for post-processing. 

For the repeat annotation, the final repeat outpus from DFAM and the *de novo* repeat models (See [Repeat annotation](03_repeat_annotation/))
 were merged and curated with the following commands:
 
Getting gff3 files
```
rmOutToGFF3.pl Emal_V1_24Chr_P0_rn_denovo.fasta.out  > merged_repeats_1.gff3
        
rmOutToGFF3.pl Emal_V1_24Chr_P0_rn_dfam.fasta.out > merged_repeats_2.gff3 
```
Validate, merge, sort, reformat gff files using genometools/1.6.1
```
gt gff3validator merged_repeats_1.gff3

gt gff3 -sort merged_repeats_1.gff3  > merged_repeats_1.sorted.gff

gt merge merged_repeats_1.sorted.gff merged_repeats_2.sorted.gff > merged_repeats.gff3

gt uniq -o merged_repeats_unique.gff3 merged_repeats.gff3
```

Final filtering
```
grep -v -e "Satellite" -e ")n" -e "-rich" merged_repeats_unique.gff3 > merged_repeats_unique_complex.gff3   # isolate complex repeats

cat merged_repeats_unique_complex.gff3 | \
perl -ane '$id; if(!/^\#/){@F = split(/\t/, $_); chomp $F[-1];$id++; $F[-1] .= "\;ID=$id"; $_ = join("\t", @F)."\n"} print $_' \
> merged_repeats_unique_complex_reformat.gff3
```



### Braker output post processing
The postprocessing of the braker output involved annotating the amino acid sequences produced by braker with Inter Pro Scan and SwissProt. This data, together with the information below was then used to filter the data in R and determine which processing step produced the best result based on BUSCO results. 

After copying the two scripts predictionAnalysis.py and selectSupportedSubsets.py (Tomas Bruna, Copyright 2020, Georgia Institute of Technology, USA) into the folder containing the Braker2 results, the following code was run.

```python selectSupportedSubsets.py --fullSupport FULLSUPPORT --anySupport ANYSUPPORT --noSupport NOSUPPORT augustus.hints.gtf hintsfile.gff```

Then, the FULLSUPPORT, ANYSUPPORT, NOSUPPORT, iprscan_result.tx, and diamond_matches.tsv files were then loaded into RStudio and the following R markdown script was run: 07_braker_output_postprocessing.Rmd


## Functional annotation
The braker post processing determined that the gene models filtered for longest isoform and any support from braker, with gene models that had no support but returned a IPS or blast result added back in, showed the best result in terms of number of genes and BUSCO results.

The functional annotation was carried out on the amino acid sequences of augustus.hints_anysupport_nosupportwithdomaindiamond_longest_isoform.aa. This included annotating the filtered amino acid sequences in interproscan and blast, using a specific output formt, which was then loaded into Omicsbox for merging, filtering, and GO-term assignment. Furthermore, KAAS was used to assing KEGG pathways to each gene.

## RNAseq analysis of developmental larval stages in R
The differential gene expression analysis for the larval developmental stages was carried out on the sequencing data from the whole larval and juvenile fish. Before processing, the data from the two lanes  were merged per sample . Low quality bases and adaptor sequences were filtered using Trim Galore (V0.6.5) and cutadapt (V2.10) using default parameters with the exception of “--length 30”. Kraken2 (V2.0.9-beta) was used to remove bacterial reads using the bacterial and archeal database (V4.08.20) and “--confidence 0.3”. Cleaned reads were mapped using STAR (V2.7.9a) with “--quantMode GeneCounts” and “--outSAMtype BAM SortedByCoordinate”, using the filtered gff file produced by the braker2 annotation outlined above for the genome indexing (--genomeSAindexNbases 13, --sjdbOverhang 149). The unstranded mapped reads were then loaded into Rstudio (V2022.02.4) using R (V3.6.3). DESeq2 (V1.36.0) was used for general data analysis, with coseq (V1.20.0) being used for cluster analysis. The cluster analysis was carried out on differentially expressed genes only, as determined through likelihood ratio test (LRT) analysis (full model: design = ~ dph, reduced model: reduced = ~ 1, adjusted p-value threshold: 0.001) in DESeq2. Normalisation was done in DESeq2, while the following parameters were used for coseq: model = "Normal", transformation = "arcsin", seed = 1234, iter=10000. Specific genes belonging to clusters where D03 and/or Day32 showed upregulation were then re-clustered with the same parameter for visualization. Pairwise analysis of differentially expressed genes between two time points was done using the Wald test in DESeq2 (design = ~ dph, adjusted p-value threshold: 0.01, log2FoldChange = ±0.58). Figures were plotted using ggplot298 (V3.4.1), and the analysis made general use of the tidyverse package (V1.3.2). 



## Data availability

### NCBI BioProject
Umbrella BioProject: PRJNA798702

Principle phase BioProject: PRJNA798188

Alternate phase BioProject PRJNA798189

Raw data BioProject: PRJNA794870


### NCBI BioSamples
Genome sequencing: SAMN24662200

ISO-seq (Tissues): SAMN24664212

RNA-Seq (Larval stages): SAMN24664213 - SAMN24664234 / SAMN32359227 - SAMN32359229


### NCBI Genome
Principle phase: JANUFT000000000

Alternate phase: JANUFU000000000


### Final gene models from BRAKER2 and and final functional annotation of gene models
DOI: 10.6084/m9.figshare.25486387

https://figshare.com/articles/dataset/Emal_final_gene_models_gff3/25486387


## References
Ryu, T. et al. A chromosome-scale genome assembly of the false clownfish, Amphiprion ocellaris. G3 12, jkac074 (2022).


