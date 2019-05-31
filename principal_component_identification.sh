#!/bin/bash

GTEXDIR=../files/
VCFFILE=${GTEXDIR}/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCaller.vcf.gz
SAMPLEFILE=${GTEXDIR}/phs000424.v7.pht002742.v7.p2.c1.GTEx_Subject_Phenotypes.GRU.txt
PRUNED_PED_FILE=${GTEXDIR}/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/652Ind_filtered_biallelic_snps_0.05maf.ped2
PRUNED_MAP_FILE=${GTEXDIR}/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/652Ind_filtered_biallelic_snps_0.05maf.map

ALLPREFIX=PCA_results/pca

# Convert to eigstrat format
parfile=convertf_parfile_v2.txt
echo "genotypename: " ${PRUNED_PED_FILE} > ${parfile}
echo "snpname: " ${PRUNED_MAP_FILE} >> ${parfile}
echo "indivname: " ${PRUNED_PED_FILE} >> ${parfile}
echo "outputformat: EIGENSTRAT" >> ${parfile}
echo "genotypeoutname: " ${ALLPREFIX}.eigenstratgeno >> ${parfile}
echo "snpoutname: " ${ALLPREFIX}.snp >> ${parfile}
echo "indivoutname: " ${ALLPREFIX}.ind >> ${parfile}
echo "familynames: NO" >> ${parfile}
convertf -p ${parfile}

python generate_ind_file.py ${ALLPREFIX}.ind ${SAMPLEFILE} > ${ALLPREFIX}.ind.poplabels

# Run smartpca.pl - GTEx only
smartpca.perl \
    -i ${ALLPREFIX}.eigenstratgeno \
    -a ${ALLPREFIX}.snp \
    -b ${ALLPREFIX}.ind.poplabels \
    -k 10 \
    -o ${ALLPREFIX}.gtex.pca \
    -e ${ALLPREFIX}.gtex.evals \
    -p ${ALLPREFIX}.gtex.plot \
    -l ${ALLPREFIX}.gtex.log \
    -w gtex_pops.txt \
    -y gtex_pops.txt

