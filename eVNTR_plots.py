import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn.apionly as sns
import numpy as np
import pandas as pd
from math import log10

from advntr.models import load_unique_vntrs_data

from run_anova import load_individual_genotypes


vntr_models_dir = '/home/mehrdad/workspace/adVNTR/vntr_data/hg38_selected_VNTRs_Illumina.db'
ref_vntrs = load_unique_vntrs_data(vntr_models_dir)
reference_vntrs = {}
for ref_vntr in ref_vntrs:
        reference_vntrs[ref_vntr.id] = ref_vntr


def plot_expression_genotype_correlation(vntr_id, tissue_name):
    gene_name = reference_vntrs[vntr_id].gene_name
    data_file = 'genotype_expression_correlation/%s/%s/correlation.txt' % (tissue_name, vntr_id)
    with open(data_file) as infile:
        lines = infile.readlines()
    data = []
    labels = []
    for line in lines:
        line = line.strip().split(' ')
        labels.append(str(line[0]))
        data.append([float(e) for e in line[1].split(',')])

    sns.set(style="whitegrid")

    fig = plt.figure()
    ax = fig.add_subplot(111)
    print(len(data))
    # data = np.array(data).transpose()
    sns.boxplot(data=data, ax=ax, showfliers=False, width=0.3, linewidth=1)
    sns.swarmplot(data=data, ax=ax, color=".2")
    ax.set_xticklabels([float(item) for item in labels], size=12)
    ax.set_xlabel('Average Number of Repeats')
    ax.set_ylabel('Gene Expression')
    ax.set_title('Expression of %s in %s tissue' % (gene_name, tissue_name))
    ax.set_ylim(0)
    fig.savefig('expression_%s_%s.pdf' % (gene_name, tissue_name))
    plt.cla()
    plt.clf()
    plt.close()


def plot_variant_caviar_scores(vntr_id, tissue_name):
    caviar_post_file = 'caviar_inputs/%s/%s/caviar_post' % (tissue_name.replace(' ', '-'), vntr_id)
    post = pd.read_csv(caviar_post_file, sep="\t", header=0)
    post = post.sort_values(post.columns[2], ascending=False)
    post = post.reset_index(drop=True)

    vntr_x = reference_vntrs[vntr_id].start_point
    vntr_x = 104629254 # I should read from hg19 :((
    vntr_y = 0
    x = []
    y = []
    for row in range(0,len(post.index)):
        if '%s' % vntr_id in post.values[row][0]:
            vntr_y = float(post.values[row][2])
        else:
            x.append(int(post.values[row][0].split('_')[1]))
            y.append(post.values[row][2])
    gene_name = reference_vntrs[vntr_id].gene_name

    if vntr_id != 111235:
        vntr_x = sum(x) / len(x)

    fig = plt.figure(figsize=(6, 2))
    ax = fig.add_subplot(111)
    ax.scatter(x, y, marker='.', c='gray')
    ax.scatter(vntr_x, vntr_y, marker='*', c='r')
    ax.set_ylabel('Causality probability')
    fig.savefig('caviar_%s_%s.pdf' % (gene_name, tissue_name))
    plt.cla()
    plt.clf()
    plt.close()


def plot_variant_pvalues(vntr_id, tissue_name):
    pvalue_file = 'pvalues/%s/%s/pvalues.txt' % (tissue_name, vntr_id)
    with open(pvalue_file) as infile:
        lines = infile.readlines()
    lines = [line.strip() for line in lines if line.strip() != '']
    y = []
    x = []
    vntr_x = 104629254 # I should read from hg19 :((
    vntr_y = 0
    for line in lines:
        snp, p = line.split()
        if '%s' % vntr_id in snp:
            vntr_y = -log10(float(p))
        else:
            y.append(-log10(float(p)))
            x.append(float(snp.split('_')[1]))
    if vntr_id != 111235:
        vntr_x = sum(x) / len(x)

    gene_name = reference_vntrs[vntr_id].gene_name

    fig = plt.figure(figsize=(6, 2))
    ax = fig.add_subplot(111)
    ax.scatter(x, y, marker='.', c='gray')
    ax.scatter(vntr_x, vntr_y, marker='*', c='r')
    # ax.set_xlabel('Genomic location')
    ax.set_ylabel('-log10 P-value')
    fig.savefig('pvalues_%s_%s.pdf' % (gene_name, tissue_name))
    plt.cla()
    plt.clf()
    plt.close()


def plot_allele_count_distribution():
    genotypes = load_individual_genotypes(reference_vntrs)
    genotyped_vntr_ids = {_id: set() for _id in range(1000000)}

    for individual_id, result_map in genotypes.items():
        for genotyped_vntr, avg_length in result_map.items():
            if avg_length == 'None' or avg_length is None:
                continue
            genotyped_vntr_ids[genotyped_vntr].add(avg_length)
            #    print(genotyped_vntr_ids)
    vntr_genotypes = {_vntr_id: len(lengths) for _vntr_id, lengths in genotyped_vntr_ids.items() if len(lengths) != 0}

    data = {}
    for vntr_id, number_of_genotypes in vntr_genotypes.items():
        number_of_genotypes = min(number_of_genotypes, 15)
        annotation = reference_vntrs[vntr_id].annotation
        if annotation not in ['UTR', 'Coding', 'Promoter']:
            continue
        if annotation not in data.keys():
            data[annotation] = {}
        if number_of_genotypes not in data[annotation].keys():
            data[annotation][number_of_genotypes] = 0
        data[annotation][number_of_genotypes] += 1

    offset = -0.2
    for key in data.keys():
        plt.bar(np.array(data[key].keys()) + offset, data[key].values(), label=key, width=0.2)
        offset += 0.2
    plt.legend()
    plt.xlabel('Number of alleles')
    plt.ylabel('Number of VNTRs')
    plt.savefig('VNTR_genotype_allele_count.pdf')
    plt.cla()
    plt.clf()
    plt.close()


def plot_vntr_polymorphic_rate_based_on_annotation():
    genotypes = load_individual_genotypes(reference_vntrs)
    genotyped_vntr_ids = {_id: set() for _id in range(1000000)}

    for individual_id, result_map in genotypes.items():
        for genotyped_vntr, avg_length in result_map.items():
            if avg_length == 'None' or avg_length is None:
                continue
            genotyped_vntr_ids[genotyped_vntr].add(avg_length)
            #    print(genotyped_vntr_ids)
    vntr_genotypes = {_vntr_id: len(lengths) for _vntr_id, lengths in genotyped_vntr_ids.items() if len(lengths) != 0}

    poly = {}
    total = {}
    for vntr_id, number_of_genotypes in vntr_genotypes.items():
        annotation = reference_vntrs[vntr_id].annotation
        if annotation not in ['UTR', 'Coding', 'Promoter']:
            continue
        if annotation not in total.keys():
            total[annotation] = 0
            poly[annotation] = 0
        if number_of_genotypes > 1:
            poly[annotation] += 1
        total[annotation] += 1

    rate = {key: float(poly[key]) / total[key] for key in poly.keys()}
    plt.bar([0, 1, 2], rate.values(), width=0.5)
    plt.xticks((0, 1, 2), [item for item in rate.keys()], size=12)
    plt.xlabel('Annotation')
    plt.ylabel('Rate of polymorphic VNTRs')
    plt.savefig('polymorphic_vntrs.pdf')
    plt.cla()
    plt.clf()
    plt.close()


def plot_evntrs_and_number_of_tissues():
    plt.xlabel('Shared across number of tissues')
    plt.ylabel('Number of VNTR-eQTLs')
    plt.savefig('eVNTR_uniqueness.pdf')
    plt.cla()
    plt.clf()
    plt.close()


def plot_evntrs_per_tissue():
    # also show unique and non-unique eVNTRs separately (with filling or color of the box)
    plt.xlabel('Tissue name')
    plt.ylabel('Number of VNTR-eQTLs')
    # plt.savefig('eVNTR_uniqueness.pdf')
    plt.cla()
    plt.clf()
    plt.close()


def plot_difference_from_reference():
    plt.xlabel('Number of VNTRs')
    plt.ylabel('Base-pairs difference from reference')
    # plt.savefig('eVNTR_uniqueness.pdf')
    plt.cla()
    plt.clf()
    plt.close()

def plot_tissue_correlations():
    #matshow or annotated heatmaps in matplotlib
    # or hexbin
    pass


def plot_vntr_location_effect():
    # y: pvalue
    # x: distance to transcription start site, gene start, first exon, ...
    pass


def plot_significant_vntrs_and_tissues():
    # show which tissues are affected for a list of 20-30 VNTRs (end of page11 in Melissa's paper)
    pass


#GWAS result? correlation with traits? like hight

if __name__ == '__main__':
     plot_expression_genotype_correlation(331737, 'Lung')
     plot_expression_genotype_correlation(331737, 'Esophagus')
#    plot_variant_pvalues(101865, 'Heart')
#    plot_variant_pvalues(331737, 'Heart')
#    plot_variant_caviar_scores(111235, 'Esophagus')
#    plot_variant_caviar_scores(331737, 'Esophagus')

    # plot_allele_count_distribution()
    # plot_vntr_polymorphic_rate_based_on_annotation()
