# JS20.01

TOLLIP prevents lipid accumulation and dampens innate immune responses during prolonged intracellular M. tuberculosis macrophage infection

Samba Venkatasubramanian 1, Courtney Plumlee 2, Kimberly A. Dill-McFarland 1, Gemma Pearson 3, Robyn Pryor 1, Scott S. Soleimanpour 3, Matthew Altman 1, Kevin B. Urdahl 2, Javeed A. Shah 1,4

1 Department of Medicine, University of Washington, Seattle, WA.
2 Seattle Children’s Research Institute, Seattle, WA.
3 Department of Medicine, University of Michigan, Ann Arbor, MI.
4 VA Puget Sound Health Care System, Seattle, WA. 

# Abstract

Effective macrophage responses to *Mycobacterium tuberculosis* (Mtb) are necessary for effective control of tuberculosis (TB). TOLLIP is a ubiquitin binding protein that controls multiple macrophage functions via endoplasmic reticulum transport and autophagy. In this study, we characterized the role of TOLLIP on macrophage function during prolonged Mtb infection using mouse knockout models of Mtb infection. *Tollip-/-* mice were susceptible to Mtb and demonstrated increased numbers of lipid-laden foam cells in lung infiltrates. Despite increased antimicrobial responses, *Tollip-/-* macrophages were preferentially Mtb- infected by 28 days after infection. Global gene expression analysis of sorted, Mtb-infected macrophages identified cellular stress as the major causal network contributor to this phenotype. We induced cellular stress by administering exogenous neutral lipids to *Tollip-/-* macrophages, which was associated with increased lipid accumulation and intracellular Mtb replication. These studies demonstrate an important dual role for TOLLIP in controlling innate immune activation and resolving lipid accumulation in Mtb-infected macrophages. 

Keywords: TOLLIP, tuberculosis, macrophages, foam cells, innate immunity, unfolded protein response, cellular homeostasis, lipid metabolism, autophagy

# Experimental design

Mixed bone marrow chimeric mice were infected with Mtb and 28 days after infection, Mtb-infected and Mtb-uninfected wild type and *Tollip-/-* alveolar macrophages (AM) were sorted and RNA-seq was performed. 

# Data cleaning (1.Shah_RNAseq_data_cleaning)
## Overview

Command line tools

1. Combine .fastq.gz files per read per sample
2. Remove sequencing adapters
3. Quality filter sequences
4. Align to reference genome
5. Quality filter alignments
6. Count reads in genes

R

7. Filter protein coding genes
8. Filter low coverage samples
9. Filter PCA outliers
10. Filter rare genes

## Files

* `data_clean/`
	- Raw sequence counts in genes (`Shah.counts.clean.csv`)
	- Cleaned data cleaning metrics from command line tools (`Shah.data.cleaning.metrics.csv`)
	- Sample metadata (`Shah.metadata.csv`)
	- EList object of quality filtered, voom normalized counts, sample metadata, and gene key (`Shah.clean.RData`)
* `figs/
	- Data cleaning metrics from command line tools (`cleaning/`)
* `results/`
	- Raw data cleaning metrics from command line tools (`results_cleaning/`)
	- FastQC outputs for raw and adapter trimmed sequences (`results_fastqc/`)
* `scripts/`
	- Bash script of command line cleaning steps (`RNAseq_mouse_pipeline.sh`)

# Differential expression (2.Shah_RNAseq_WGCNA_contrasts)
## Overview

1. Define differentially expressed (DE) genes 
2. Cluster DE genes into modules
3. Define DE modules

## Files

* `figs/`
	- Gene expression box plots for all genes (`gene_level_contrast/`)
	- Module expression box plots for all modules (`module_Shah_contrast_deepSplit3_minMod50/`)
	- WGCNA soft threshold cutoff (`module_Shah_contrast_deepSplit3_minMod50/SFT_thresholding_power21.png`)
	- Principle component analysis of gene and module expression (`PCA*.png`)
* `results/`
	- Quality filtered, voom normalized counts table (same as in `data_clean/Shah.clean.RData`) and linear model results for all genes (`gene_level/`)
	- Mean voom normalized counts in modules, linear model results, and lists of genes in modules (`module_Shah_contrast_deepSplit3_minMod50/`)
* `scripts/`
	- Function for generation of expression box plots (`RNAseq_boxplot_fxn.R`) 

# Gene set enrichment analysis (3.GSEA)
## Overview

1. GSEA of DE genes
2. GSEA of genes in each module
3. GSEA of DE fold change groups
	- Up in uninfected only
	- Down in uninfected only
	- Up in infected only
	- Down in infected only
	- Up in both
	- Down in both

## Files

* `data_clean/`
	- DE fold change groups for genes in module 0 (`module_0_sorted.RData`)
* `results/`
	- Hallmark term enrichment in genes, modules, and fold change groups (`GSEA/`)

# Additional files

* `data_clean/`
	- Ensembl to HGNC symbol key (downloaded from https://www.genenames.org/download/custom/)
* `figs/`
	- Figure 6 A-D files (`publication/`)
* `scripts/`
	- Scripts for Figure A-D generation (`Fig6*.R`)

***