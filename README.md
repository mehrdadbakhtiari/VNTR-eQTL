# VNTR-eQTL
Scripts for analyzing VNTR-eQTLs in GTEx cohort for the [paper](http://biorxiv.org) (link to be added soon)

# Inputs of the pipeline
1. VNTR genotypes of the entire GTEx cohort as outputted by [adVNTR](https://github.com/mehrdadbakhtiari/adVNTR) in the default format (one text file for each individual) using WGS data aligned to GRCh38 (available from dbGaP phs000424.v7.p2).
2. Phenotype file for GTEx cohort (ethnicity and sex are used). accession: pht002742.v7.p2 and exact version used in paper: phs000424.v7.pht002742.v7.p2.c1.GTEx_Subject_Phenotypes.GRU.txt
3. SNP genotypes file for GTEx cohort (phg001219.v1.GTEx_v8_WGS.genotype-calls-vcf.c1). This file is used to compare effect of VNTR variants and SNPs, and identifying population structure.
4. RNA-expression data for GTEx cohort in different tissues as a table (phe000020.v1.GTEx_RNAseq.expression-data-matrixfmt.c1).

# Output of the pipeline
Linear regression results showing association of the length of each VNTR with the expression level of its nearest gene in each of the 46 tissues.
The pipeline will create a directory called `regression_results` for the output. It contains 46 subdirectories (one for each tissue) and each of them have a tab separated file for each VNTR showing linear test for the VNTR in the corresponding tissue.
```
regression_results
└───Whole-Blood
│   │   100374.txt
│   │   100436.txt
│   │   ...
│   
└───Lung
│   │   100374.txt
│   │   100436.txt
│   │   ...
│   
└───etc 
```
Each text file contains 6 tab separated values: <br>
```
Genename  Chromosome  VNTR Coordinate Start   Effect Size (B)   P-value   Standard Error (Bse)
```

For example, `regression_results/Whole-Blood/423956.txt` (result of the association test for VNTR 423956 in Whole Blood) has the following content:
```
POMC    chr2    25161573        0.21607193665873414     9.109131954009228e-06   0.048056165330533064
```

# Requirements
1. eigensoft. this can be installed with `apt install eigensoft` on linux or `conda install -c bioconda/label/cf201901 eigensoft` using conda package manager on suported systems.
2. 

# How to run
### Preprocessing
```
python extract_expression_by_tissue.py
```
### Finding population structure
```
./principal_component_identification.sh
```
### Computing PEER factors
```
python peer_factor_identification.py
```
### Running association test
### Identifying significance threshold (5% FDR)
