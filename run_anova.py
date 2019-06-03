import glob
import os
import operator
import pandas as pd
import numpy as np

import statsmodels.api as sm
from statsmodels.formula.api import ols

from advntr.models import load_unique_vntrs_data

from load_confounders import load_peer_factors, load_population_pcs
from gene_locations import GeneLocations


significant_vntrs = 0
highest_fs = 0
lowest_p = 1e10
VNTR_genotypes_dir = '../gtex_genotypes/'
#VNTR_genotypes_dir = '/pedigree2/projects/adVNTR/gtex_genotypes/'
wgs_id_gtex_id_file = '../GTEX_sample_id_conversion.txt'
vntr_models_dir = '/pedigree2/projects/adVNTR/vntr_data/hg38_selected_VNTRs_Illumina.db'

snp_directory = '../files/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/'
snp_file = snp_directory + 'subset'
snp_file = snp_directory + 'subset_100kb'
original_snp_file = snp_directory + '652Ind_filtered_biallelic_snps_0.05maf'

gene_locations_obj = GeneLocations()

rpkm_directory = '../Expression_by_Tissue/'

ref_vntrs = load_unique_vntrs_data(vntr_models_dir)
reference_vntrs = {}
for ref_vntr in ref_vntrs:
    reference_vntrs[ref_vntr.id] = ref_vntr


def get_average(lst):
    return sum(lst) / len(lst)


def get_wgs_id_to_individual_id_map():
    with open(wgs_id_gtex_id_file) as infile:
        lines = infile.readlines()
    result = {}
    for line in lines:
        line = line.strip().split()
        result[line[1]] = line[0]
    return result


def load_individual_genotypes():
    res = {}
    # res['GTEX-QWERT'][527655] = 2.5
    wgs_id_to_gtex_id = get_wgs_id_to_individual_id_map()
    genotype_files = glob.glob(VNTR_genotypes_dir + '*.out')
#    genotype_files = glob.glob(VNTR_genotypes_dir + 'out_*')

    for genotype_file in genotype_files:
        wgs_id = os.path.basename(genotype_file).split('.')[0]
#        wgs_id = os.path.basename(genotype_file).split('_')[1]
        individual_id = wgs_id_to_gtex_id[wgs_id]
        res[individual_id] = {}
        with open(genotype_file) as infile:
            lines = infile.readlines()
        _vntr_id = None
        for i, line in enumerate(lines):
            line = line.strip()
            if i % 2 == 0:
                _vntr_id = int(line)
            else:
                if len(reference_vntrs[_vntr_id].pattern) >= 6:
                    res[individual_id][_vntr_id] = get_average([float(e) for e in line.split('/')]) if line != 'None' else None
    
    return res


def load_snp_file(vntr_id):
    gene_name = reference_vntrs[vntr_id].gene_name
    snp_file = snp_directory + 'subsets/subset_%s' % vntr_id
    if not os.path.exists(snp_file + '.ped') or not os.path.exists(snp_file + '.map'):
        chromosome = reference_vntrs[vntr_id].chromosome[3:]
        start, end = gene_locations_obj.get_gene_coordinates(gene_name)
        if start is None:
            return {}, []
        start = (start - 50000) / 1000
        end = (end + 50000) / 1000
        command = 'plink --file %s --out %s --chr %s --from-kb %s --to-kb %s --noweb --recode' % (original_snp_file, snp_file, chromosome, start, end)
        os.system(command)

#    return {}, []
    map_file = snp_file + '.map'
    ped_file = snp_file + '.ped'

    result = {}
    #result[GTEX-QWERT]['chr_position_...'] = 0 or 1
    ancestral = {}
    all_snps = []
    snps = []

    with open(map_file) as infile:
        lines = infile.readlines()
    for line in lines:
        chromosome, snp_id, _, position = line = line.strip().split('\t')
        ancestral[snp_id] = snp_id.split('_')[-3]
        all_snps.append(snp_id)

    with open(ped_file) as infile:
        lines = infile.readlines()
    for line_number, line in enumerate(lines):
        line = line.strip().split()
        individual_id = line[0]
        result[individual_id] = {}
        for i in range(6, len(all_snps) + 6):
            snp_index = i - 6
            snp_id = all_snps[snp_index]
            snp_chr, snp_pos = snp_id.split('_')[0:2]
            if True or snp_chr == '10' and 104629210 - 100000 < int(snp_pos) < 104661655 + 100000:
                if line_number == 0:
                    snps.append(snp_id)
                alleles = line[6+2*snp_index], line[6+2*snp_index + 1]
                alleles = [e == ancestral[snp_id] for e in alleles]
                result[individual_id][snp_id] = sum(alleles)

    return result, snps


def run_anova():
    genotypes = load_individual_genotypes()
    genotyped_vntr_ids = {_id: set() for _id in range(1000000)}
    for individual_id, result_map in genotypes.items():
        for genotyped_vntr, avg_length in result_map.items():
            if avg_length == 'None' or avg_length is None:
                continue
            genotyped_vntr_ids[genotyped_vntr].add(avg_length)
#    print(genotyped_vntr_ids)
    vntr_genotypes = {_vntr_id: len(lengths) for _vntr_id, lengths in genotyped_vntr_ids.items() if len(lengths) != 0}
    vntr_genotypes = sorted(vntr_genotypes.items(), key=operator.itemgetter(1), reverse=True)
    print(vntr_genotypes[0:10], vntr_genotypes[-1])

    tissue_name = 'Blood Vessel'
#    tissue_name = 'Brain'
    tissue_rpkm_file = rpkm_directory + tissue_name + '.rpkm'
    df = pd.read_csv(tissue_rpkm_file, delimiter='\t', header=1)
    run_anova_for_vntr(df, genotypes, 111235)
#    exit(0)
    for vntr_id, number_of_genotypes in vntr_genotypes:
        if reference_vntrs[vntr_id].chromosome[3:] == 'Y':
            continue
        if number_of_genotypes != 3:
            continue
        run_anova_for_vntr(df, genotypes, vntr_id, tissue_name)
#        break


def run_anova_for_vntr(df, genotypes, vntr_id = 527655, tissue_name = 'Blood Vessel'):
    gene_name = reference_vntrs[vntr_id].gene_name.replace('-', '__')
    gene_df = df.loc[df['Description'] == gene_name]
    if gene_df.shape[0] == 0:
        # don't have expression for this gene
        return

    gene_df = gene_df.reset_index()
    gene_df = gene_df.drop(columns=['index', 'Name'])

    columns = gene_df.columns
    to_drop_cols = []
    for col in columns:
        if not col.startswith('GTEX'):
            continue
        if col not in genotypes.keys() or vntr_id not in genotypes[col] or genotypes[col][vntr_id] is None:
            to_drop_cols.append(col)
    gene_df = gene_df.drop(columns=to_drop_cols)

    genotypes_row = []
    for i in range(1, len(gene_df.columns)):
        individual_id = gene_df.columns[i]
        if individual_id in genotypes.keys() and vntr_id in genotypes[individual_id].keys():
            genotypes_row.append(genotypes[individual_id][vntr_id])
        else:
            genotypes_row.append(None)
    found_genotypes = sorted(list(set([e for e in genotypes_row if e is not None])))
    print(found_genotypes)
    vntr_genotype_title = '%s_%s_Genotype' % (gene_name, vntr_id)
    gene_df.loc[1] = [vntr_genotype_title] + genotypes_row
#    print(gene_df)
#    print(gene_df.columns)
#    print('----')

#    print('is it 2? ', gene_df.shape[0])
#    for i in range(0, len(snps)):
#        snp_genotype_row = []
#        for j in range(1, len(gene_df.columns)):
#            individual_id = gene_df.columns[j]
#            snp_genotype_row.append(snp_map[individual_id][snps[i]])
#        snp_titles.append('SNP_%s' % snps[i].split('_')[1])
#        gene_df.loc[i+2] = [snp_titles[-1]] + snp_genotype_row

    # add a row for each peer factor
    start_row = 2
    peer_factor_map, peer_factors = load_peer_factors(tissue_name)
    for i in range(0, len(peer_factors)):
        peer_row = []
        for j in range(1, len(gene_df.columns)):
            individual_id = gene_df.columns[j]
            peer_row.append(peer_factor_map[individual_id][peer_factors[i]])
        gene_df.loc[i+start_row] = [peer_factors[i]] + peer_row

    # add a row for each population pc
    start_row = 2 + len(peer_factors)
    pop_pc_map, pop_pcs = load_population_pcs()
    for i in range(0, len(pop_pcs)):
        pop_row = []
        for j in range(1, len(gene_df.columns)):
            individual_id = gene_df.columns[j]
            if individual_id not in pop_pc_map.keys():
                pop_row.append(0)
            else:
                pop_row.append(pop_pc_map[individual_id][pop_pcs[i]])
        gene_df.loc[i+start_row] = [pop_pcs[i]] + pop_row

    temp = gene_df.set_index('Description').transpose()
    # temp: (except SNPs that will be added later)
    # Description   gene_name   vntr_genotype_title SNP1    SNP2    ... SNPn    Peer1   Peer2   PopPC1  PopPC2
    # GTEX-QWETY    0.6         2.5                 2       1       ... 1       0.5     0.6     0.4     0.8

    groups = [list(temp.loc[temp['%s' % vntr_genotype_title] == repeat_count]['%s' % gene_name]) for repeat_count in found_genotypes]
#    print("groups: ", len(groups))
#    for e in groups:
#        print(len(e), e)
#    if len(groups) != 3:
#        print('len(group) != 3')
#        exit(0)

    anova_target = gene_name
    anova_target = 'residuals'

    # find residuals
    confoudners_str = ' + '.join(pop_pcs + peer_factors)
    confounders_lm = ols('%s ~ %s' % (gene_name, confoudners_str), data=temp).fit()
    temp['residuals'] = confounders_lm.resid
#    print(temp)

    temp['const'] = 1
    print(temp.shape)
#    print(vntr_genotype_title)
    vntr_mod = ols('%s ~ %s' % (anova_target, vntr_genotype_title), data=temp).fit()
#    vntr_mod = sm.OLS(temp[gene_name].astype(float), temp[[vntr_genotype_title, 'const']].astype(float)).fit()
#    print(vntr_mod.summary())
    print('summary printed for ', vntr_mod.fvalue, vntr_mod.f_pvalue)
#    aov_table = sm.stats.anova_lm(vntr_mod)
#    print(aov_table)
    significant_vntr = True

    print('DONE with OLS, starting ols for SNPs')
    if vntr_mod.f_pvalue > 0.05:
        print('not a significant VNTR for sure')
        return

    # add rows for each SNP genotype
    snp_map, snps = load_snp_file(vntr_id)
    print('snps are loaded')
    snp_titles = []
    for i in range(len(snps)):
        snp_column = [snp_map[individual_id][snps[i]] for individual_id in temp.index]
        snp_titles.append('SNP_%s' % snps[i].split('_')[1])
        temp[snp_titles[-1]] = snp_column

    best_snp_p = 1e9
    best_snp_f = 0
    best_snp = None

    for snp_id in snp_titles:
        snp_mod = ols('%s ~ %s' % (anova_target, snp_id), data=temp).fit()
        snp_mod = sm.OLS(temp[anova_target], temp[[snp_id, 'const']]).fit()
        snp_vntr_lm = ols('%s ~ %s + %s' % (anova_target, snp_id, vntr_genotype_title), data=temp).fit()
        snp_vntr_lm = sm.OLS(temp[anova_target], temp[[snp_id, vntr_genotype_title, 'const']]).fit()

        # f-value and p-value of adding VNTR to SNP linear model
        vntr_contribution_result = snp_vntr_lm.compare_f_test(snp_mod)[0:2]
#        print(snp_mod.fvalue, snp_mod.f_pvalue, snp_vntr_lm.fvalue, snp_vntr_lm.f_pvalue, vntr_contribution_result, snp_vntr_lm.compare_f_test(vntr_mod)[0:2])
        if snp_mod.fvalue > best_snp_f:
            best_snp_f = snp_mod.fvalue
            best_snp_p = snp_mod.f_pvalue
            best_snp = snp_id
        if vntr_contribution_result[1] > 0.05:
            print('not significant', vntr_contribution_result, snp_mod.fvalue, snp_mod.f_pvalue, snp_vntr_lm.fvalue, snp_vntr_lm.f_pvalue)
            significant_vntr = False
        anova_results = sm.stats.anova_lm(snp_mod, snp_vntr_lm)
#        print(anova_results["Pr(>F)"].values[1]) # this is same as vntr_contribution_result[1]
#        print('----')

    print('significant_vntr: ', significant_vntr)
    print('best_snp info: ', best_snp_f, best_snp_p)
    fs, pv = vntr_mod.fvalue, vntr_mod.f_pvalue
    global highest_fs
    global lowest_p
    print(highest_fs, fs)
    print(lowest_p, pv)
    if fs > highest_fs:
        highest_fs = fs
    if pv < lowest_p:
        lowest_p = pv
    smallest_group = 1e10
    for g in groups:
        smallest_group = min(smallest_group, len(g))
    if pv < 0.05 and smallest_group >= 5:# and len(groups[0]) >= 5 and len(groups[1]) >= 5 and len(groups[2]) >= 5:
        global significant_vntrs
        if significant_vntr:
            significant_vntrs += 1
        for g in groups:
            print(g)
            if len(g) > 0:
                print(get_average(g))
        print(vntr_id, gene_name)
#        exit(0)


if __name__ == '__main__':
    run_anova()
    print(highest_fs)
    print(lowest_p)
    print(significant_vntrs)

