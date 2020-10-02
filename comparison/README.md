# Performance comparison
Scripts for evaluating the performance of adVNTR-NN, VNTRseek, and GangSTR with simulated reads.

## Simulating reads
We used [ART](niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm) to simulate reads as stated in [Method](https://doi.org/10.1101/2020.05.25.114082) section.

## Inputs
* Output (genotype calls) of each tool for the simulated tests.

> **GangSTR** outputs the genotype results along with a VCF file. We used the genotype output to compare the results.
> **VNTRseek** outputs genotypes as a VCF fille. We used the VCF file to compare the results.
> **adVNTR-NN** outputs genotypes by default. We used the genotype output.

## Output
* calc_accuracy_advntr_vntrseek.py outputs the accuracies of adVNTR-NN and VNTRseek for each scenario given the genotype outputs.
* compare_wiht_gangstr_6_20_vntrs.py outputs the accuracy of adVNTR-NN and GangSTR for each scenario given the genotype outputs.

## Requirements
1. [adVNTR (v1.4.0)](https://github.com/mehrdadbakhtiari/adVNTR)
2. [GangSTR (v2.4.5)](https://github.com/gymreklab/GangSTR)
3. [VNTRseek (v1.10.0)](https://github.com/yzhernand/VNTRseek/)
4. Reference VNTR Databases for adVNTR can be downloaded from [here](https://cseweb.ucsd.edu/~mbakhtia/adVNTR/).
Note that we used HG19 reference for VNTRseek but HG38 reference for GangSTR.
5. To use this script, file name and locations should be changed appropriately.
