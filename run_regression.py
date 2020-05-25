import glob
import os
import operator
import pandas as pd
import sys
import numpy as np
import math

import statsmodels.api as sm
from statsmodels.formula.api import ols

from advntr.models import load_unique_vntrs_data

from load_confounders import load_peer_factors, load_population_pcs, load_genders
from gene_locations import GeneLocations
from hwe_test import get_loci_in_HWE

VNTR_genotypes_dir = '../gtex_genotypes/'
genotypes_dir_1kg = '../1kg_advntr_genotype/'
#VNTR_genotypes_dir = '/pedigree2/projects/adVNTR/gtex_genotypes/'
wgs_id_gtex_id_file = '../GTEX_sample_id_conversion.txt'
# vntr_models_dir = '/pedigree2/projects/adVNTR/vntr_data/hg38_selected_VNTRs_Illumina.db'
vntr_models_dir = '/home/mehrdad/workspace/adVNTR/vntr_data/hg38_selected_VNTRs_Illumina.db.bck'

snp_directory = '../files/phg001219.v1.GTEx_v8_WGS.genotype-calls-vcf.c1/'
original_snp_file = snp_directory + 'plink_866Ind_filtered_biallelic_snps_0.05maf'

rpkm_directory = '../Expression_by_Subtissue/'
caviar_result_dir = 'caviar_inputs/'

if len(sys.argv) > 2:
    VNTR_genotypes_dir = sys.argv[1] + '/'
    rpkm_directory = sys.argv[2] + '/'

low_pvalue_vntrs = {}
significant_vntrs = {}
top_p_value_vntrs = {}
beat_top_10_snps_vntrs = {}
beat_top_20_snps_vntrs = {}
beat_top_100_snps_vntrs = {}
best_pvalues = {}
best_pvalue_ranks = {}
caviar_top_5 = {}
caviar_top_1 = {}
best_caviar_probs = {}
best_caviar_ranks = {}

min_individuals_in_group = 4
significance_threshold = 0.0005
thresholds = {}
try:
    with open('thresholds.txt') as infile:
        lines = infile.readlines()
        for l in lines:
            thresholds[l.split('\t')[0].replace(' ', '-')] = float(l.split('\t')[1])
except:
    # need to rerun after finding thresholds from initial results
    rpkms = glob.glob(rpkm_directory + '/*')
    for rpkm in rpkms:
        thresholds[rpkm.split('.')[0]] = 0.0005

run_permutation_test = False
bootstrapping = False

highest_fs = 0
lowest_p = 1e10

try:
    gene_locations_obj = GeneLocations()
except:
    pass

if __name__ == '__main__':
    ref_vntrs = load_unique_vntrs_data(vntr_models_dir)
    reference_vntrs = {}
    for ref_vntr in ref_vntrs:
        reference_vntrs[ref_vntr.id] = ref_vntr


def get_average(lst):
    return sum(lst) / len(lst)


def get_wgs_id_to_individual_id_map():
    try:
        with open(wgs_id_gtex_id_file) as infile:
            lines = infile.readlines()
    except:
        lines = []
    result = {}
    for line in lines:
        line = line.strip().split()
        result[line[1]] = line[0]
    return result


def load_1k_genotypes(reference_vntrs, average=True, limit=None, individual_ids=None):
    res = {}
    # res['HG03193'][527655] = 2.5
    genotype_files = glob.glob(genotypes_dir_1kg + '*.genotype')
    if limit:
        genotype_files = genotype_files[:limit]

    for genotype_file in genotype_files:
        individual_id = os.path.basename(genotype_file).split('.')[0]
        if individual_ids is not None and individual_id not in individual_ids:
            continue
        res[individual_id] = {}
        with open(genotype_file) as infile:
            lines = infile.readlines()
        if lines[0].startswith('This is an modifie'):
            lines = lines[1:]
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


def load_individual_genotypes(reference_vntrs, average=True, limit=None):
    res = {}
    # res['GTEX-QWERT'][527655] = 2.5
    wgs_id_to_gtex_id = get_wgs_id_to_individual_id_map()
    genotype_files = glob.glob(VNTR_genotypes_dir + '*.out')
    if limit:
        genotype_files = genotype_files[:limit]

    for genotype_file in genotype_files:
        wgs_id = os.path.basename(genotype_file).split('.')[0]
#        wgs_id = os.path.basename(genotype_file).split('_')[1]
        if wgs_id in wgs_id_to_gtex_id.keys():
            individual_id = wgs_id_to_gtex_id[wgs_id]
        else:
            individual_id = wgs_id
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
        try:
            start, end = gene_locations_obj.get_gene_coordinates(gene_name)
        except:
            start = reference_vntrs[vntr_id].start_point
            end = reference_vntrs[vntr_id].start_point + reference_vntrs[vntr_id].get_length()
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
            if X.corr(Y) is np.nan or math.isnan(X.corr(Y)):
                result[i].append(0.0)
            else:
                result[i].append(X.corr(Y))
    return pd.DataFrame(result, columns=variant_titles, index=variant_titles)


def lookfor(x, p):
    """look for occurence of x in frame column p and output now numb, ID and score"""
    for row in range(0,len(p.index)):
        if x in p.values[row][0]:
            top = p.values[row][0]
            score = p.values[row][2]
            return(row, top, score)


def run_caviar(caviar_variants, caviar_zscores, gene_df, tissue_name, vntr_id):
    rank = 1e10
    top_vntr_score = 0

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
        caviar_cmd = 'CAVIAR -l "%s" -z "%s" -o "%s" -c 1 -f 1 > "%s"' % (ld_file, z_file, caviar_output, temp_dir+"log")
        os.system(caviar_cmd)

    if not os.path.exists(caviar_post_file):
        print('caviar didnt produce posterior probability file')
    else:
        post = pd.read_csv(caviar_post_file, sep="\t", header=0)
        post = post.sort_values(post.columns[2], ascending=False)
        post = post.reset_index(drop=True)
        I, topvntr, top_vntr_score = lookfor('%s' % vntr_id, post)
        rank = I + 1

    return rank, top_vntr_score


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
    global best_pvalue_ranks
    global best_caviar_ranks
    global best_caviar_probs
    best_pvalues = {_id: 1 for _id in range(1000000)}
    best_pvalue_ranks = {_id: 1e10 for _id in range(1000000)}
    best_caviar_ranks = {_id: 1e10 for _id in range(1000000)}
    best_caviar_probs = {_id: 0 for _id in range(1000000)}
#    global beat_top_100_snps_vntrs
#    global top_p_value_vntrs
#    beat_top_100_snps_vntrs = {_id: set() for _id in range(1000000)}
#    top_p_value_vntrs = {_id: set() for _id in range(1000000)}

    vntr_genotypes = get_vntr_alleles(genotypes)
    print(vntr_genotypes[0:10], vntr_genotypes[-1])

    peer_files = glob.glob('PEER_results/*')
    tissue_names = [pf.split('/')[-1][13:-3] for pf in peer_files]
    print(tissue_names)

#    with open('important_vntr_ids.txt') as infile:
#        lines = infile.readlines()
#    important_vntr_ids = set([int(line.strip()) for line in lines if line.strip() != ''])
#    print(important_vntr_ids)
    eqtl_targets = get_loci_in_HWE(reference_vntrs, load_individual_genotypes(reference_vntrs, False))

    computed_tissues = [e[19:] for e in glob.glob('regression_results/*')]
    print computed_tissues

    for tissue_name in tissue_names:
        tissue_rpkm_file = rpkm_directory + tissue_name + '.rpkm'
        df = pd.read_csv(tissue_rpkm_file, delimiter='\t', header=1)
        if len(df.columns) < 100:
            print('skip %s as it has too few individuals' % tissue_name)
            continue
        for vntr_id, number_of_genotypes in vntr_genotypes:
            if vntr_id not in eqtl_targets:
                continue
            if reference_vntrs[vntr_id].chromosome[3:] == 'Y':
                continue
            if number_of_genotypes <= 1:
                continue
#            if vntr_id not in important_vntr_ids:
#                continue
            if bootstrapping:
                res = []
                for k in range(100):
                    res.append(run_anova_for_vntr(df, genotypes, vntr_id, tissue_name))
                print(sum([1 if (e < significance_threshold * 2) else 0 for e in res]), vntr_id,
                '%s: %s' % (reference_vntrs[vntr_id].chromosome, reference_vntrs[vntr_id].start_point))
            else:
                run_anova_for_vntr(df, genotypes, vntr_id, tissue_name)


def load_bjarni_genotypes():
#    # res['BJARNI_0'][527655] = 2.5
    genotype_file = 'Bjarni/toMehrdad/VNTR_genotypes'
    genotype_file = 'Bjarni/toMehrdad/GRCh38_Blood_VNTR_genotypes.txt'
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
    expression_file = 'Bjarni/toMehrdad/Blood_VNTR-eQTLs_ENSEMBL87.table'
    df = pd.read_csv(expression_file, delimiter=' ', header=None)
    df = df.drop(columns=[1, 2, 3])
    df.columns = ['Gene'] + ['Bjarni_%s' % i for i in range(len(df.columns[1:]))]

    # add a row for each peer factor
    start_row = len(df.index)
    peer_factor_map, peer_factors = load_peer_factors(file_name='peer_factors_Bjarni_3')
    for i in range(0, len(peer_factors)):
        peer_row = []
        for j in range(1, len(df.columns)):
            individual_id = str(int(df.columns[j][7:]) + 4)
            peer_row.append(peer_factor_map[individual_id][peer_factors[i]])
        df.loc[i + start_row] = [peer_factors[i]] + peer_row

    df = df.set_index([df.columns[0]])
    df = df.transpose()
    pvalues = {}

    gene_name_vntr_map = {}
    with open('Blood_VNTR-eQTLs.txt') as infile:
        blood_vntrs = [int(line.strip()) for line in infile.readlines()]
    for r in reference_vntrs.values():
        if r.id in blood_vntrs:
            gene_name_vntr_map[r.gene_name] = r.id
    df.columns = [gene_name_vntr_map[g] if not g.startswith('peer') else g for g in df.columns]

    het_test = {}
    with open('heterozygosity_test.txt') as infile:
        lines = infile.readlines()
        for l in lines:
            vid, pval = l.split()
            vid = int(vid)
            pval = float(pval)
            het_test[vid] = pval

    positives = 0
    negatives = 0
    for vntr_id in df.columns:
        if str(vntr_id).startswith('peer_'):
            continue
        if het_test[vntr_id] < 0.05:
            continue
        if reference_vntrs[vntr_id].chromosome[3:] == 'Y':
            continue
        number_of_alleles = len(set([e for e in genotypes_df[vntr_id] if e not in ['None', None]]))
        if number_of_alleles <= 1:
            continue
        gene_df = df.drop(columns=list([e for e in df.columns if e != vntr_id and not str(e).startswith('peer_')]))
        anova_target = 'expression_%s' % vntr_id
        gene_df.columns = [anova_target] + list(gene_df.columns[1:])
        gene_df[anova_target] = get_normalized_gene_expression(gene_df[anova_target])
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
        gene_df.drop(gene_df.index[to_remove], inplace=True)
        gene_df[vntr_genotype_title] = np.array(gene_df[vntr_genotype_title], dtype='float')

        confoudners_str = ' + '.join(peer_factors)#+pop_pcs
        confounders_lm = ols('%s ~ %s' % (anova_target, confoudners_str), data=gene_df).fit()
        gene_df['residuals'] = np.array(confounders_lm.resid, dtype='float')

        gene_df['const'] = 1
        vntr_mod = ols('%s ~ %s' % ('residuals', vntr_genotype_title), data=gene_df).fit()
        pvalues[vntr_id] = vntr_mod.f_pvalue

        with open('regression_results/Whole-Blood/%s.txt' % vntr_id) as infile:
            original_effect_size = float(infile.readlines()[0].split('\t')[3])
        effect_size = vntr_mod.params[vntr_genotype_title]
        print(vntr_mod.f_pvalue < significance_threshold * 2, vntr_id,
              '%s: %s' % (reference_vntrs[vntr_id].chromosome, reference_vntrs[vntr_id].start_point))
        if vntr_mod.f_pvalue < significance_threshold * 2:
            print(vntr_id)
            if original_effect_size * effect_size < 0:
                print('not matching effect size direction for %s', vntr_id)
            positives += 1
        else:
            negatives += 1
        continue
        run_anova_for_vntr(df, genotypes, vntr_id)
    print('------')
    print(positives, negatives)


def run_anova_for_geuvadis():
    with open('ensemble87_to_transcript_ens.txt') as infile:
        lines = infile.readlines()[1:]
    gene_ens_to_trans_ens = {line.strip().split()[0]: [] for line in lines}
    for line in lines:
        g, t = line.strip().split()
        gene_ens_to_trans_ens[g].append(t)

    with open('ensemblToGeneName.txt') as infile:
        lines = infile.readlines()
    trans_ens_to_gene_name = {line.strip().split()[0]: line.strip().split()[1] for line in lines}
    expression_df = pd.read_csv('Geuvadis/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt', delimiter='\t', header=0)

    gene_ensembles = expression_df['Gene_Symbol']
    gene_names = []
    trans_ens_keys = set(trans_ens_to_gene_name.keys())
    gene_ens_keys = set(gene_ens_to_trans_ens.keys())
    for g_ens in gene_ensembles:
        found = False
        g_ens = g_ens.split('.')[0]
        if g_ens in gene_ens_keys:
            for t_ens in gene_ens_to_trans_ens[g_ens]:
                if t_ens in trans_ens_keys:
                    gene_names.append(trans_ens_to_gene_name[t_ens])
                    found = True
                    break
        if not found:
            gene_names.append(g_ens)

    ensembl_ids = {gene_names[i]: expression_df['Gene_Symbol'][i] for i in range(len(gene_names))}
    expression_df['Gene_Symbol'] = gene_names
    individual_ids = set(list(expression_df.columns[4:]))

    genotypes = load_1k_genotypes(reference_vntrs, individual_ids=individual_ids)
    vntr_genotypes = get_vntr_alleles(genotypes)

    with open('Blood_VNTR-eQTLs.txt') as infile:
        blood_vntrs = [int(line.strip()) for line in infile.readlines()]
    print(len(blood_vntrs), len(vntr_genotypes))

    het_test = {}
    with open('heterozygosity_test.txt') as infile:
        lines = infile.readlines()
        for l in lines:
            vid, pval = l.split()
            vid = int(vid)
            pval = float(pval)
            het_test[vid] = pval

    positive, negative = 0, 0
    for vntr_id, number_of_genotypes in vntr_genotypes:
        if reference_vntrs[vntr_id].chromosome[3:] == 'Y':
            continue
        if vntr_id not in blood_vntrs: # TODO: temporary
            continue
        if het_test[vntr_id] < 0.05:
            continue
        if number_of_genotypes <= 1:
            continue
        tissue_name = 'Whole Blood'
        p_value_file = 'geuvadis_all_vntr_pvalues/%s/%s/pvalues.txt' % (tissue_name, vntr_id)
#        if os.path.exists(p_value_file) and os.path.getsize(p_value_file) > 0: #TODO: Temporary
#            continue
        gene_name = reference_vntrs[vntr_id].gene_name.replace('-', '__')

        gene_df = expression_df.loc[expression_df['Gene_Symbol'] == gene_name]
        if gene_df.shape[0] == 0:
            print('no expression')
            # don't have expression for this gene
            continue
        while gene_df.shape[0] != 1:
            gene_df.drop(gene_df.index[[1]], inplace=True)

        gene_df = gene_df.reset_index()
        gene_df.drop(columns=['TargetID', 'Chr', 'Coord', 'index'], inplace=True)

        genotypes_row = get_genotypes_row_for_df(gene_df, genotypes, vntr_id)
        vntr_genotype_title = '%s_%s_Genotype' % (gene_name, vntr_id)
        gene_df.loc[len(gene_df.index)] = [vntr_genotype_title] + genotypes_row

        to_drop_cols = []
        for i, col in enumerate(gene_df.columns[1:]):
            if genotypes_row[i] == None:
                to_drop_cols.append(col)
        gene_df = gene_df.drop(columns=to_drop_cols)

        # add a row for each peer factor
        start_row = 2
        peer_factor_map, peer_factors = load_peer_factors(file_name='peer_factors_Geuvadis_15')
        for i in range(0, len(peer_factors)):
            peer_row = []
            for j in range(1, len(gene_df.columns)):
                individual_id = gene_df.columns[j]
                peer_row.append(peer_factor_map[individual_id][peer_factors[i]])
            gene_df.loc[i + start_row] = [peer_factors[i]] + peer_row

        # add a row for each population pc
        start_row = 2 + len(peer_factors)
        pop_pc_map, pop_pcs = load_population_pcs('1KG_PCA_results/pca.gtex.pca.evec')
        for i in range(0, len(pop_pcs)):
            pop_row = []
            for j in range(1, len(gene_df.columns)):
                individual_id = gene_df.columns[j]
                if individual_id not in pop_pc_map.keys():
                    pop_row.append(0)
                else:
                    pop_row.append(pop_pc_map[individual_id][pop_pcs[i]])
            gene_df.loc[i + start_row] = [pop_pcs[i]] + pop_row

        gene_df = gene_df.set_index('Gene_Symbol').transpose()

        if (np.median(gene_df[gene_name]) == 0):
            print(ensembl_ids[gene_name], 'gene not expression')
            # gene is not expressed
            continue
        gene_df[gene_name] = get_normalized_gene_expression(gene_df[gene_name])

        confoudners_str = ' + '.join(peer_factors + pop_pcs)
        confounders_lm = ols('%s ~ %s' % (gene_name, confoudners_str), data=gene_df).fit()
        gene_df['residuals'] = confounders_lm.resid


        tissue_name = 'Whole Blood'
        if run_permutation_test:
            for i in range(100):
                gene_df['shuffled_genotypes'] = np.random.permutation(list(gene_df[vntr_genotype_title]))
                permutated_mod = ols('%s ~ %s' % ('residuals', 'shuffled_genotypes'), data=gene_df).fit()
                permutated_pvalue_file = '/mnt/geuvadis_permutated_pvalues_multiple_tests/permutated_pvalues_%s/%s/%s/permutated_pvalues.txt' % (i, tissue_name, vntr_id)
                if not os.path.exists(os.path.dirname(permutated_pvalue_file)):
                    os.makedirs(os.path.dirname(permutated_pvalue_file))
                with open(permutated_pvalue_file, 'w') as outfile:
                    outfile.write('%s\n' % permutated_mod.f_pvalue)

        vntr_mod = ols('%s ~ %s' % ('residuals', vntr_genotype_title), data=gene_df).fit()
        print(vntr_mod.f_pvalue < 0.0023370968265774856, vntr_mod.f_pvalue, ensembl_ids[gene_name], vntr_id, '%s: %s' % (reference_vntrs[vntr_id].chromosome, reference_vntrs[vntr_id].start_point))

        with open('regression_results/Whole-Blood/%s.txt' % vntr_id) as infile:
            original_effect_size = float(infile.readlines()[0].split('\t')[3])
        effect_size = vntr_mod.params[vntr_genotype_title]
        if effect_size * original_effect_size < 0:
            print('vntr %s have different association direction' % vntr_id)
        if vntr_mod.f_pvalue < 0.0023370968265774856:
            positive += 1
        else:
            negative += 1

        p_value_file = 'geuvadis_all_vntr_pvalues/%s/%s/pvalues.txt' % (tissue_name, vntr_id)
        if not os.path.exists(os.path.dirname(p_value_file)):
            os.makedirs(os.path.dirname(p_value_file))
        with open(p_value_file, 'w') as outfile:
            outfile.write('%s\n' % vntr_mod.f_pvalue)

        regression_results = 'geuvadis_regression_results/%s/%s.txt' % ('Whole Blood'.replace(' ', '-'), vntr_id)
        if not os.path.exists(os.path.dirname(regression_results)):
            os.makedirs(os.path.dirname(regression_results))
        with open(regression_results, 'w') as outfile:
            outfile.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (gene_name, reference_vntrs[vntr_id].chromosome,
                                                        reference_vntrs[vntr_id].start_point,
                                                        vntr_mod.params[vntr_genotype_title],
                                                        vntr_mod.f_pvalue, vntr_mod.bse[vntr_genotype_title]))

        # run_anova_for_vntr(df, genotypes, vntr_id, tissue_name)
    print(positive, negative)


# normalize expression values by fitting them to a normal distribution
def get_normalized_gene_expression(gene_expression):
    from sklearn.preprocessing import quantile_transform
    normalied_expr = quantile_transform(np.array([[e] for e in gene_expression]), random_state=0, copy=True)
    return list(pd.DataFrame(normalied_expr)[0])


def get_genotypes_row_for_df(gene_df, genotypes, vntr_id):
    genotypes_row = []
    for i in range(1, len(gene_df.columns)):
        individual_id = gene_df.columns[i]
        if individual_id in genotypes.keys() and vntr_id in genotypes[individual_id].keys() and genotypes[individual_id][vntr_id] is not None:
            genotypes_row.append(genotypes[individual_id][vntr_id])
        else:
            genotypes_row.append(None)
    return genotypes_row


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
    seq = [i for i in range(1, len(gene_df.columns))]
    if bootstrapping:
        from random import sample
        to_drop_cols = sample(seq, len(seq) / 5)
        gene_df = gene_df.drop(columns=gene_df.columns[to_drop_cols])

    genotypes_row = get_genotypes_row_for_df(gene_df, genotypes, vntr_id)
    vntr_genotype_title = '%s_%s_Genotype' % (gene_name, vntr_id)
    gene_df.loc[1] = [vntr_genotype_title] + genotypes_row

    normalized_expressions = get_normalized_gene_expression(gene_df.loc[0][1:])
    all_genotypes = sorted(list(set([e for e in genotypes_row if e is not None])))
    plot_data = [[normalized_expressions[i] for i in range(len(normalized_expressions)) if genotypes_row[i] == repeat_count] for repeat_count in all_genotypes]

    for i in range(len(genotypes_row)):
        if genotypes_row.count(genotypes_row[i]) < min_individuals_in_group:
            genotypes_row[i] = None
    found_genotypes = sorted(list(set([e for e in genotypes_row if e is not None])))
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

    # add a row for gender
    start_row = 2 + len(peer_factors) + len(pop_pcs)
    gender_map = load_genders()
    gender_row = []
    for j in range(1, len(gene_df.columns)):
        individual_id = gene_df.columns[j]
        if individual_id not in gender_map.keys():
            gender_row.append(0)
        else:
            gender_row.append(0 if gender_map[individual_id] == 'M' else 1)
    gene_df.loc[start_row] = ['Sex'] + gender_row

    temp = gene_df.set_index('Description').transpose()
    # temp: (except SNPs that will be added later)
    # Description   gene_name   vntr_genotype_title SNP1    SNP2    ... SNPn    Peer1   Peer2   PopPC1  PopPC2  Sex
    # GTEX-QWETY    0.6         2.5                 2       1       ... 1       0.5     0.6     0.4     0.8     1

    if (np.median(temp[gene_name]) == 0):
        # gene is not expressed
        return

    temp[gene_name] = get_normalized_gene_expression(temp[gene_name])
    groups = [list(temp.loc[temp['%s' % vntr_genotype_title] == repeat_count]['%s' % gene_name]) for repeat_count in found_genotypes]
    anova_target = 'residuals'

    # find residuals
    confoudners_str = ' + '.join(pop_pcs + peer_factors + ['Sex'])
    confounders_lm = ols('%s ~ %s' % (gene_name, confoudners_str), data=temp).fit()
    temp['residuals'] = confounders_lm.resid

    temp['const'] = 1
    vntr_mod = ols('%s ~ %s' % (anova_target, vntr_genotype_title), data=temp).fit()

    if bootstrapping:
        return vntr_mod.f_pvalue

    if run_permutation_test:
        for i in range(100):
            temp['shuffled_genotypes'] = np.random.permutation(list(temp[vntr_genotype_title]))

            permutated_mod = ols('%s ~ %s' % (anova_target, 'shuffled_genotypes'), data=temp).fit()

            permutated_pvalue_file = '/mnt/permutated_pvalues_multiple_tests/permutated_pvalues_%s/%s/%s/permutated_pvalues.txt' % (i, tissue_name, vntr_id)
            if not os.path.exists(os.path.dirname(permutated_pvalue_file)):
                os.makedirs(os.path.dirname(permutated_pvalue_file))
            with open(permutated_pvalue_file, 'w') as outfile:
                outfile.write('%s\n' % permutated_mod.f_pvalue)

    print('summary printed for ', vntr_mod.fvalue, vntr_mod.f_pvalue, tissue_name, vntr_id)
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

    if vntr_mod.f_pvalue > thresholds[tissue_name.replace(' ', '-')]:
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
    drop_list = ["GTEX-11DXY", "GTEX-12BJ1", "GTEX-13NYS", "GTEX-13O1R", "GTEX-14A5I", "GTEX-14ICL", "GTEX-14PHW", "GTEX-16Z82", "GTEX-OHPK", "GTEX-QLQW", "GTEX-YFCO", "GTEX-ZVTK"]
    drop_list = [e for e in drop_list if e in temp.index]
    temp = temp.drop(drop_list)
    for i in range(len(snps)):
        snp_column = [snp_map[individual_id][snps[i]] for individual_id in temp.index]
        snp_titles.append('SNP_%s' % snps[i].split('_')[1])
        temp[snp_titles[-1]] = snp_column

    best_snp_p = 1e9
    best_snp_f = 0

    snp_linear_models = []
    vntr_pvalue_rank = 1
    for snp_id in snp_titles:
        snp_mod = ols('%s ~ %s' % (anova_target, snp_id), data=temp).fit()
        snp_mod = sm.OLS(temp[anova_target], temp[[snp_id, 'const']]).fit()
        snp_linear_models.append((snp_mod.f_pvalue, snp_id, snp_mod))
        if snp_mod.f_pvalue < vntr_mod.f_pvalue:
            vntr_pvalue_rank += 1
        variant_pvalues.append((snp_id, snp_mod.f_pvalue))
    best_pvalue_ranks[vntr_id] = min(vntr_pvalue_rank, best_pvalue_ranks[vntr_id])

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

        import math
        if counter < 100 and not math.isnan(get_caviar_zscore(snp_mod, snp_id)):
            caviar_variants.append(snp_id)
            caviar_zscores.append(get_caviar_zscore(snp_mod, snp_id))

#        anova_results = sm.stats.anova_lm(snp_mod, snp_vntr_lm)
        counter += 1

    causality_rank, causality_prob = run_caviar(caviar_variants, caviar_zscores, temp, tissue_name, vntr_id)
    best_caviar_ranks[vntr_id] = min(causality_rank, best_caviar_ranks[vntr_id])
    best_caviar_probs[vntr_id] = max(causality_prob, best_caviar_probs[vntr_id])
    global caviar_top_5
    global caviar_top_1
    if causality_rank <= 5:
        add_tissue(caviar_top_5, vntr_id, tissue_name)
    if causality_rank == 1:
        add_tissue(caviar_top_1, vntr_id, tissue_name)

    significant_vntr = vntr_mod.f_pvalue < thresholds[tissue_name.replace(' ', '-')]
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
            for i, g in enumerate(plot_data):
                outfile.write('%s %s\n' % (all_genotypes[i], ','.join([str(e) for e in g])))
                print(g)
                if len(g) > 0:
                    print(get_average(g))
        print(vntr_id, gene_name)
#        exit(0)


if __name__ == '__main__':
    # run_anova_for_geuvadis()
    # run_anova_for_bjarni()
    # exit(0)
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

    all_vntrs = sorted(list(set(significant_vntrs.keys() + top_p_value_vntrs.keys() + beat_top_10_snps_vntrs.keys() + beat_top_20_snps_vntrs.keys() + beat_top_100_snps_vntrs.keys())))

    # with open('causal_vntrs.txt', 'w') as outfile:
    #     for vntr_id in all_vntrs:
    #         outfile.write('%s\t%s\t%s\n' % (vntr_id, best_caviar_ranks[vntr_id], best_caviar_probs[vntr_id]))

    for vntr_id in all_vntrs:
        print('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t' % (vntr_id, reference_vntrs[vntr_id].gene_name, reference_vntrs[vntr_id].annotation,
              reference_vntrs[vntr_id].pattern, best_pvalues[vntr_id], ','.join(top_p_value_vntrs[vntr_id]) if vntr_id in top_p_value_vntrs.keys() else '-',
              ','.join(beat_top_100_snps_vntrs[vntr_id]) if vntr_id in beat_top_100_snps_vntrs.keys() else '-',
              ','.join(beat_top_20_snps_vntrs[vntr_id]) if vntr_id in beat_top_20_snps_vntrs.keys() else '-',
              ','.join(beat_top_10_snps_vntrs[vntr_id]) if vntr_id in beat_top_10_snps_vntrs.keys() else '-',
              ','.join(caviar_top_1[vntr_id]) if vntr_id in caviar_top_1.keys() else '-',
              best_pvalue_ranks[vntr_id], best_caviar_ranks[vntr_id], best_caviar_probs[vntr_id]))

