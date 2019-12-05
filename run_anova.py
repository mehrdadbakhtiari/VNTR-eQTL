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


low_pvalue_vntrs = {}
significant_vntrs = {}
top_p_value_vntrs = {}
beat_top_10_snps_vntrs = {}
beat_top_20_snps_vntrs = {}
beat_top_100_snps_vntrs = {}
best_pvalues = {}
caviar_top_5 = {}
caviar_top_1 = {}

min_individuals_in_group = 4
significance_threshold = 0.0005

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

try:
    gene_locations_obj = GeneLocations()
except:
    pass

rpkm_directory = '../Expression_by_Tissue/'

caviar_result_dir = 'caviar_inputs/'

if __name__ == '__main__':
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


def load_individual_genotypes(reference_vntrs, average=True):
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
                    avg_length = get_average([float(e) for e in line.split('/')]) if line != 'None' else None
                    diploid_alleles = [float(e) for e in line.split('/')] if line != 'None' else (None, None)
                    res[individual_id][_vntr_id] = avg_length if average else diploid_alleles
    
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
            if line_number == 0:
                snps.append(snp_id)
            alleles = line[6+2*snp_index], line[6+2*snp_index + 1]
            alleles = [e == ancestral[snp_id] for e in alleles]
            result[individual_id][snp_id] = sum(alleles)

    return result, snps


def add_tissue(result_dict, vntr_id, tissue_name):
    if vntr_id not in result_dict.keys():
        result_dict[vntr_id] = set([])
    result_dict[vntr_id].add(tissue_name.replace(' ', '-'))


def get_caviar_zscore(linear_model, variant_title):
    zscore = linear_model.params[variant_title] / linear_model.bse[variant_title]
    return zscore


def get_caviar_ld_matrix(df, variant_titles):
    result = []
    for i, v1 in enumerate(variant_titles):
        result.append([])
        for v2 in variant_titles:
            X = df[v1].replace('None', np.nan).astype(float)
            Y = df[v2].replace('None', np.nan).astype(float)
            if X.corr(Y) is np.nan:
                result[i].append(0.0)
            else:
                result[i].append(X.corr(Y))
    return pd.DataFrame(result, columns=variant_titles, index=variant_titles)


def lookfor (x, p):
    """look for occurence of x in frame column p and output now numb, ID and score"""
    for row in range(0,len(p.index)):
        if x in p.values[row][0]:
            top = p.values[row][0]
            score = p.values[row][2]
            return(row, top, score)


def run_caviar(caviar_variants, caviar_zscores, gene_df, tissue_name, vntr_id):
    rank = 1e10

    tissue_name = tissue_name.replace(' ', '-')
    temp_dir = caviar_result_dir + '%s/%s/' % (tissue_name, vntr_id)
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)
    ld_file = temp_dir + '/%s_LD' % tissue_name
    z_file = temp_dir + '/%s_Z' % tissue_name
    caviar_output = temp_dir + 'caviar'
    caviar_post_file = caviar_output + '_post'
    if not os.path.exists(caviar_post_file):
        caviar_ld_matrix = get_caviar_ld_matrix(gene_df, caviar_variants)
        caviar_ld_matrix.to_csv(ld_file, sep='\t', header=None, index=None)
        with open(z_file, 'w') as outfile:
            for i in range(len(caviar_variants)):
                outfile.write('%s\t%s\n' % (caviar_variants[i], caviar_zscores[i]))
        caviar_cmd = "CAVIAR -l %s -z %s -o %s -c 1 -f 1 > %s" % (ld_file, z_file, caviar_output, temp_dir+"log")
        os.system(caviar_cmd)

    if not os.path.exists(caviar_post_file):
        print('caviar didnt produce posterior probability file')
    else:
        post = pd.read_csv(caviar_post_file, sep="\t", header=0)
        post = post.sort_values(post.columns[2], ascending=False)
        post = post.reset_index(drop=True)
        I, topvntr, top_vntr_score = lookfor('%s' % vntr_id, post)
        rank = I + 1

    return rank


def get_vntr_alleles(genotypes):
    genotyped_vntr_ids = {_id: set() for _id in range(1000000)}
    for individual_id, result_map in genotypes.items():
        for genotyped_vntr, avg_length in result_map.items():
            if avg_length == 'None' or avg_length is None:
                continue
            genotyped_vntr_ids[genotyped_vntr].add(avg_length)

    vntr_genotypes = {_vntr_id: len(lengths) for _vntr_id, lengths in genotyped_vntr_ids.items() if len(lengths) != 0}
    vntr_genotypes = sorted(vntr_genotypes.items(), key=operator.itemgetter(1), reverse=True)
    return vntr_genotypes


def run_anova():
    genotypes = load_individual_genotypes(reference_vntrs)
    global best_pvalues
    best_pvalues = {_id: 1 for _id in range(1000000)}
#    global beat_top_100_snps_vntrs
#    global top_p_value_vntrs
#    beat_top_100_snps_vntrs = {_id: set() for _id in range(1000000)}
#    top_p_value_vntrs = {_id: set() for _id in range(1000000)}

    vntr_genotypes = get_vntr_alleles(genotypes)
    print(vntr_genotypes[0:10], vntr_genotypes[-1])

    peer_files = glob.glob('PEER_results/*')
    tissue_names = [pf.split('/')[-1][13:-3] for pf in peer_files]
    print(tissue_names)

    with open('important_vntr_ids.txt') as infile:
        lines = infile.readlines()
    important_vntr_ids = set([int(line.strip()) for line in lines if line.strip() != ''])
    print(important_vntr_ids)

    tissue_name = 'Blood Vessel'
    for tissue_name in tissue_names:
        tissue_rpkm_file = rpkm_directory + tissue_name + '.rpkm'
        df = pd.read_csv(tissue_rpkm_file, delimiter='\t', header=1)
        for vntr_id, number_of_genotypes in vntr_genotypes:
            if reference_vntrs[vntr_id].chromosome[3:] == 'Y':
                continue
            if number_of_genotypes <= 1:
                continue
#            if vntr_id not in important_vntr_ids:
#                continue
            run_anova_for_vntr(df, genotypes, vntr_id, tissue_name)
#            exit(0)
#        break

def load_bjarni_genotypes():
#    # res['BJARNI_0'][527655] = 2.5
    genotype_file = 'Bjarni/toMehrdad/VNTR_genotypes'
#    with open(genotype_file) as infile:
#        lines = genotype.readlines()
#        lines = [line.strip().split() for line in lines]
#    for i in range(len(lines)):
#        individual_id = 'Bjarni_%s' % i
#        vntr_id = lines[i][0]
#        genotypes[individual_id] = {}
#        for j in range(1, len(lines[0])):
#            genotypes[individual_id][vntr_id] = lines[i][j]
#    return genotypes
    df = pd.read_csv(genotype_file, delimiter=' ', header=None)
#    df = df.drop(columns=[1, 2, 3])
    df = df.set_index([0])
    df.columns = ['Bjarni_%s' % i for i in range(len(df.columns))]
    df = df.transpose()
    for col in df.columns:
        for row in range(df.shape[0]):
            if df[col][row] not in [None, 'None']:
                df[col][row] = sum([float(e) for e in df[col][row].split('/')])/ 2.0
    return df

def run_anova_for_bjarni():
    genotypes_df = load_bjarni_genotypes()
#    vntr_genotypes = get_vntr_alleles(genotypes)
    expression_file = 'Bjarni/toMehrdad/eMat_cropped'
    df = pd.read_csv(expression_file, delimiter=' ', header=None)
    df = df.drop(columns=[1, 2, 3])
    df = df.set_index([0])
    df.columns = ['Bjarni_%s' % i for i in range(len(df.columns))]
    df = df.transpose()
#    for vntr_id, number_of_alleles in vntr_genotypes:
    pvalues = {}
    for vntr_id in df.columns:
#        if vntr_id != 315000:
#            continue
        if reference_vntrs[vntr_id].chromosome[3:] == 'Y':
            continue
        number_of_alleles = len(set([e for e in genotypes_df[vntr_id] if e not in ['None', None]]))
        if number_of_alleles <= 1:
            continue
        gene_df = df.drop(columns=list([e for e in df.columns if e != vntr_id]))
        anova_target = 'expression_%s' % vntr_id
        gene_df.columns = [anova_target]
#        print(genotypes_df[vntr_id])
        if 'None' in genotypes_df[vntr_id]:
            print('bad')

        gene_name = reference_vntrs[vntr_id].gene_name
        vntr_genotype_title = '%s_%s_Genotype' % (gene_name, vntr_id)
        gene_df[vntr_genotype_title] = genotypes_df[vntr_id]
#        print(gene_df)

        to_remove = []
        for i in range(gene_df.shape[0]):
            if gene_df[vntr_genotype_title][i] in [None, 'None']:
                to_remove.append(i)
 #       print(to_remove)
        gene_df.drop(gene_df.index[to_remove], inplace=True)
        
        vntr_mod = ols('%s ~ %s' % (anova_target, vntr_genotype_title), data=gene_df).fit()
        pvalues[vntr_id] = vntr_mod.f_pvalue
        if vntr_mod.f_pvalue < significance_threshold:
            print(vntr_id)
        continue
        run_anova_for_vntr(df, genotypes, vntr_id)
    print('------')
    print(pvalues[690585]) # not mainly in blood
    print(pvalues[315000])
    print(pvalues[393574])
    print(pvalues[316163])
    print(pvalues[423956])
    print(pvalues[386120])
    print(pvalues[319811])
    print(pvalues[273409])


def run_anova_for_vntr(df, genotypes, vntr_id=527655, tissue_name='Blood Vessel'):
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
        if individual_id in genotypes.keys() and vntr_id in genotypes[individual_id].keys() and genotypes[individual_id][vntr_id] is not None:
            genotypes_row.append(genotypes[individual_id][vntr_id])
        else:
            genotypes_row.append(None)
    vntr_genotype_title = '%s_%s_Genotype' % (gene_name, vntr_id)
    gene_df.loc[1] = [vntr_genotype_title] + genotypes_row

    for i in range(len(genotypes_row)):
        if genotypes_row.count(genotypes_row[i]) < min_individuals_in_group:
            genotypes_row[i] = None
    found_genotypes = sorted(list(set([e for e in genotypes_row if e is not None])))
    print('found_genotypes:', found_genotypes)
    if len(found_genotypes) < 2:
        # only one genotype for this VNTR
        return

    to_drop_cols = []
    for i, col in enumerate(gene_df.columns[1:]):
        if genotypes_row[i] == None:
            to_drop_cols.append(col)
    gene_df = gene_df.drop(columns=to_drop_cols)

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

    if (np.median(temp[gene_name]) == 0):
        # gene is not expressed
        return

    # normalize expression values by fitting them to a normal distribution
    from sklearn.preprocessing import quantile_transform
    normalied_expr = quantile_transform(np.array([[e] for e in temp[gene_name]]), random_state=0, copy=True)
    temp[gene_name] = list(pd.DataFrame(normalied_expr)[0])

    groups = [list(temp.loc[temp['%s' % vntr_genotype_title] == repeat_count]['%s' % gene_name]) for repeat_count in found_genotypes]

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
    print('summary printed for ', vntr_mod.fvalue, vntr_mod.f_pvalue, tissue_name)
#    aov_table = sm.stats.anova_lm(vntr_mod)
#    print(aov_table)

    print('DONE with OLS, starting ols for SNPs')

    p_value_file = 'all_vntr_pvalues/%s/%s/pvalues.txt' % (tissue_name, vntr_id)
    if not os.path.exists(os.path.dirname(p_value_file)):
        os.makedirs(os.path.dirname(p_value_file))
    with open(p_value_file, 'w') as outfile:
        outfile.write('%s\n' % vntr_mod.f_pvalue)

    effect_size = get_caviar_zscore(vntr_mod, vntr_genotype_title)
    effect_size_file = 'all_vntrs_zscore/%s/%s/Z.txt' % (tissue_name, vntr_id)
    if not os.path.exists(os.path.dirname(effect_size_file)):
        os.makedirs(os.path.dirname(effect_size_file))
    with open(effect_size_file, 'w') as outfile:
        outfile.write('%s\n' % effect_size)

    regression_results = 'regression_results/%s/%s.txt' % (tissue_name.replace(' ', '-'), vntr_id)
    if not os.path.exists(os.path.dirname(regression_results)):
        os.makedirs(os.path.dirname(regression_results))
    with open(regression_results, 'w') as outfile:
        outfile.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (gene_name, reference_vntrs[vntr_id].chromosome,
                                    reference_vntrs[vntr_id].start_point, vntr_mod.params[vntr_genotype_title],
                                    vntr_mod.f_pvalue, vntr_mod.bse[vntr_genotype_title]))

    return
    if vntr_mod.f_pvalue > significance_threshold:
        print('not a significant VNTR for sure')
        return

    global low_pvalue_vntrs
    add_tissue(low_pvalue_vntrs, vntr_id, tissue_name)
    best_pvalues[vntr_id] = min(best_pvalues[vntr_id], vntr_mod.f_pvalue)

    variant_pvalues = [(vntr_genotype_title, vntr_mod.f_pvalue)]

    # add rows for each SNP genotype
    snp_map, snps = load_snp_file(vntr_id)
    if snps == []:
        print('PLINK hasnt been run')
        return
    print('snps are loaded')
    snp_titles = []
    for i in range(len(snps)):
        snp_column = [snp_map[individual_id][snps[i]] for individual_id in temp.index]
        snp_titles.append('SNP_%s' % snps[i].split('_')[1])
        temp[snp_titles[-1]] = snp_column

    best_snp_p = 1e9
    best_snp_f = 0
    best_snp = None

    snp_linear_models = []
    for snp_id in snp_titles:
        snp_mod = ols('%s ~ %s' % (anova_target, snp_id), data=temp).fit()
        snp_mod = sm.OLS(temp[anova_target], temp[[snp_id, 'const']]).fit()
        snp_linear_models.append((snp_mod.f_pvalue, snp_id, snp_mod))
        variant_pvalues.append((snp_id, snp_mod.f_pvalue))

    p_value_file = 'pvalues/%s/%s/pvalues.txt' % (tissue_name, vntr_id)
    if not os.path.exists(os.path.dirname(p_value_file)):
        os.makedirs(os.path.dirname(p_value_file))
    with open(p_value_file, 'w') as outfile:
        for vid, v_pvalue in variant_pvalues:
            outfile.write('%s\t%s\n' % (vid, v_pvalue))

    beat_10 = True
    beat_20 = True
    beat_100 = True

    caviar_variants = [vntr_genotype_title]
    caviar_zscores = [get_caviar_zscore(vntr_mod, vntr_genotype_title)]
    snp_linear_models = sorted(snp_linear_models)
    counter = 0
    for snp_pvalue, snp_id, snp_mod in snp_linear_models:
        snp_vntr_lm = ols('%s ~ %s + %s' % (anova_target, snp_id, vntr_genotype_title), data=temp).fit()
        snp_vntr_lm = sm.OLS(temp[anova_target], temp[[snp_id, vntr_genotype_title, 'const']]).fit()

        # f-value and p-value of adding VNTR to SNP linear model
        vntr_contribution_result = snp_vntr_lm.compare_f_test(snp_mod)[0:2]
#        print(snp_mod.fvalue, snp_mod.f_pvalue, snp_vntr_lm.fvalue, snp_vntr_lm.f_pvalue, vntr_contribution_result, snp_vntr_lm.compare_f_test(vntr_mod)[0:2])
        if snp_mod.fvalue > best_snp_f:
            best_snp_f = snp_mod.fvalue
            best_snp_p = snp_mod.f_pvalue
            best_snp = snp_id

        if vntr_contribution_result[1] > 0.05 and counter < 10:
            beat_10 = False
        if vntr_contribution_result[1] > 0.05 and counter < 20:
            beat_20 = False
        if vntr_contribution_result[1] > 0.05 and counter < 100:
            beat_100 = False

        if counter < 100:
            caviar_variants.append(snp_id)
            caviar_zscores.append(get_caviar_zscore(snp_mod, snp_id))

#        anova_results = sm.stats.anova_lm(snp_mod, snp_vntr_lm)
        counter += 1

    causality_rank = run_caviar(caviar_variants, caviar_zscores, temp, tissue_name, vntr_id)
    global caviar_top_5
    global caviar_top_1
    if causality_rank <= 5:
        add_tissue(caviar_top_5, vntr_id, tissue_name)
    if causality_rank == 1:
        add_tissue(caviar_top_1, vntr_id, tissue_name)

    significant_vntr = vntr_mod.f_pvalue < significance_threshold
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
    if significant_vntr and smallest_group >= min_individuals_in_group:
        global significant_vntrs
        add_tissue(significant_vntrs, vntr_id, tissue_name)

        global top_p_value_vntrs
        if pv < best_snp_p:
            add_tissue(top_p_value_vntrs, vntr_id, tissue_name)
        global beat_top_10_snps_vntrs
        global beat_top_20_snps_vntrs
        global beat_top_100_snps_vntrs
        if beat_10:
            add_tissue(beat_top_10_snps_vntrs, vntr_id, tissue_name)
        if beat_20:
            add_tissue(beat_top_20_snps_vntrs, vntr_id, tissue_name)
        if beat_100:
            add_tissue(beat_top_100_snps_vntrs, vntr_id, tissue_name)

        expression_correlation_file = 'genotype_expression_correlation/%s/%s/correlation.txt' % (tissue_name, vntr_id)
        if not os.path.exists(os.path.dirname(expression_correlation_file)):
            os.makedirs(os.path.dirname(expression_correlation_file))
        with open(expression_correlation_file, 'w') as outfile:
            for i, g in enumerate(groups):
                outfile.write('%s %s\n' % (found_genotypes[i], ','.join([str(e) for e in g])))
                print(g)
                if len(g) > 0:
                    print(get_average(g))
        print(vntr_id, gene_name)
#        exit(0)


if __name__ == '__main__':
    run_anova_for_bjarni()
    exit(0)
    run_anova()
    print(highest_fs)
    print(lowest_p)
    print('low_pvalue_vntrs: %s:' % len(low_pvalue_vntrs.keys()), low_pvalue_vntrs)
    print('significant vntrs: %s:' % len(significant_vntrs.keys()), significant_vntrs)
    print('top_pvalue_vntrs: %s:' % len(top_p_value_vntrs.keys()), top_p_value_vntrs)
    print('beat_10: %s:' % len(beat_top_10_snps_vntrs.keys()), beat_top_10_snps_vntrs)
    print('beat_20: %s:' % len(beat_top_20_snps_vntrs.keys()), beat_top_20_snps_vntrs)
    print('beat_100: %s:' % len(beat_top_100_snps_vntrs.keys()), beat_top_100_snps_vntrs)
    print('caviar_top_1: %s: ' % len(caviar_top_1.keys()))
    print('caviar_top_5: %s: ' % len(caviar_top_5.keys()))

    all_vntrs = sorted(list(set(top_p_value_vntrs.keys() + beat_top_10_snps_vntrs.keys() + beat_top_20_snps_vntrs.keys() + beat_top_100_snps_vntrs.keys())))
    for vntr_id in all_vntrs:
        print('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (vntr_id, reference_vntrs[vntr_id].gene_name, reference_vntrs[vntr_id].annotation,
              reference_vntrs[vntr_id].pattern, best_pvalues[vntr_id], ','.join(top_p_value_vntrs[vntr_id]) if vntr_id in top_p_value_vntrs.keys() else '-',
              ','.join(beat_top_100_snps_vntrs[vntr_id]) if vntr_id in beat_top_100_snps_vntrs.keys() else '-',
              ','.join(beat_top_20_snps_vntrs[vntr_id]) if vntr_id in beat_top_20_snps_vntrs.keys() else '-',
              ','.join(beat_top_10_snps_vntrs[vntr_id]) if vntr_id in beat_top_10_snps_vntrs.keys() else '-'))

