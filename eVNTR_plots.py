import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn.apionly as sns
import statsmodels.api as sm
import numpy as np
import pandas as pd
from math import log10
import operator
import os
import glob

from advntr.models import load_unique_vntrs_data

from run_anova import load_individual_genotypes


vntr_models_dir = '/home/mehrdad/workspace/adVNTR/vntr_data/hg38_selected_VNTRs_Illumina.db.bck'
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
    overall = []
    normalize_labels = []
    for line in lines:
        line = line.strip().split(' ')
        labels.append(str(line[0]))
        data.append([float(e) for e in line[1].split(',')])
        overall += data[-1]
        normalize_labels += [len(labels)-1 for _ in range(len(data[-1]))]
    from sklearn.preprocessing import quantile_transform
    overall = [[e] for e in overall]
    overall = quantile_transform(np.array(overall), n_quantiles=10, random_state=0, copy=True)
    for i in range(len(data)):
        data[i] = [overall[j][0] for j in range(len(overall)) if normalize_labels[j] == i]

    import seaborn as sns
    sns.set()
    sns.set_style("whitegrid")
    sns.set_context("poster")


    pvalue = 1
    pvalue_file = 'pvalues/%s/%s/pvalues.txt' % (tissue_name, vntr_id)
    with open(pvalue_file) as infile:
        lines = infile.readlines()
        lines = [line.strip() for line in lines if line.strip() != '']
        for line in lines:
            snp, p = line.split()
            if '%s' % vntr_id in snp:
                pvalue = float(p)
                break

    fig = plt.figure()
    ax = fig.add_subplot(111)
    print(len(data))
    # data = np.array(data).transpose()
    sns.boxplot(data=data, ax=ax, showfliers=False, width=0.3, linewidth=1)
    sns.swarmplot(data=data, ax=ax, color=".2")
    if vntr_id == 331737:
        ax.set_xticklabels(['Normal Count', 'Expanded'], size=12)
        ax.set_xlabel('VNTR Genotype')
    else:
        ax.set_xticklabels([float(item) for item in labels], size=12)
        ax.set_xlabel('Average Number of Repeats')
    ax.set_ylabel('Normalized Gene Expression')
    ax.set_title('Expression of %s Gene' % (gene_name,))
    ax.set_ylim(0)
    plt.tight_layout()
    fig.savefig('expression_%s_%s.png' % (gene_name, tissue_name), dpi=300)
    plt.cla()
    plt.clf()
    plt.close()


def plot_dose_dependence_of_correlation():
    significant_vntrs = set([])

    tissue_data_dir = glob.glob('pvalues/*')
    for tissue_dir in tissue_data_dir:
        vntrs_data = glob.glob(tissue_dir + '/*')
        tissue_name = os.path.basename(tissue_dir)
        for vntr_pvalue_dir in vntrs_data:
            vntr_id = int(os.path.basename(vntr_pvalue_dir))
            significant_vntrs.add((tissue_name, vntr_id))

    increasing = 0
    decreasing = 0
    tissue_data_dir = glob.glob('caviar_inputs/*')
    print(tissue_data_dir)
    for tissue_dir in tissue_data_dir:
        vntr_dirs = glob.glob(tissue_dir + '/*')
        tissue_name = os.path.basename(tissue_dir)
        print(tissue_name)
        print(vntr_dirs)
        for vntr_dir in vntr_dirs:
            print(vntr_dir)
            vntr_id = int(os.path.basename(vntr_dir))
            effect_size_file = 'caviar_inputs/%s/%s/%s_Z' % (tissue_name, vntr_id, tissue_name)
            with open(effect_size_file) as infile:
                lines = infile.readlines()
            effect_size = float(lines[0].strip().split()[1])
            if effect_size > 0:
                increasing += 1
            else:
                decreasing += 1
    plt.bar([0, 1], [increasing, decreasing], width=0.5)
    plt.xticks((0, 1), ['Increasing', 'Decreasing'], size=12)
    plt.xlabel('Dose Dependence', size=14)
    plt.ylabel('Number of (Tissue, VNTR-eQTL) pairs', size=14)
    plt.xlim(-0.6, 1.6)
    # plt.ylim(0, 550)
    plt.savefig('dose_dependence.pdf')
    plt.cla()
    plt.clf()
    plt.close()


def get_vntr_x(vntr_id, x): # I should read from hg19 :((
    vntr_x = sum(x) / len(x)
    if vntr_id == 111235:
        vntr_x = 104620254-1000
    if vntr_id == 315000:
        vntr_x = 54319062 - 600
    if vntr_id == 142591:
        vntr_x = 65405510
    if vntr_id == 690585:
        vntr_x = 79950699
    if vntr_id == 410522:
        vntr_x = 50985470

    return vntr_x


def plot_variant_caviar_scores(vntr_id, tissue_name):
    caviar_post_file = 'caviar_inputs/%s/%s/caviar_post' % (tissue_name.replace(' ', '-'), vntr_id)
    post = pd.read_csv(caviar_post_file, sep="\t", header=0)
    post = post.sort_values(post.columns[2], ascending=False)
    post = post.reset_index(drop=True)

    vntr_x = reference_vntrs[vntr_id].start_point
    vntr_x = 0
    vntr_y = 0
    x = []
    y = []
    for row in range(0,len(post.index)):
        if '%s' % vntr_id in post.values[row][0]:
            vntr_y = float(post.values[row][2])
        else:
            x.append(int(post.values[row][0].split('_')[1]))
            y.append(post.values[row][2])
    vntr_x = get_vntr_x(vntr_id, x)
    print(vntr_x)
    gene_name = reference_vntrs[vntr_id].gene_name


    fig = plt.figure(figsize=(6, 2))
    ax = fig.add_subplot(111)
    ax.scatter(x, y, marker='.', c='gray')
    ax.scatter(vntr_x, vntr_y, marker='*', c='r')
    ax.set_ylabel('Causality probability')

    plt.tight_layout()
    fig.savefig('caviar_%s_%s.png' % (gene_name, tissue_name), dpi=300)
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
    vntr_x = 0
    vntr_y = 0
    for line in lines:
        snp, p = line.split()
        if '%s' % vntr_id in snp:
            vntr_y = -log10(float(p))
        else:
            y.append(-log10(float(p)))
            x.append(float(snp.split('_')[1]))
    vntr_x = get_vntr_x(vntr_id, x)
    gene_name = reference_vntrs[vntr_id].gene_name

    # import seaborn as sns
    # sns.set()
    # sns.set_style("whitegrid")
    # sns.set_context("talk")

    fig = plt.figure(figsize=(6, 2))
    ax = fig.add_subplot(111)
    ax.scatter(x, y, marker='.', c='gray')
    ax.scatter(vntr_x, vntr_y, marker='*', c='r')
    # ax.set_xlabel('Genomic location')
    ax.set_ylabel('-log10 P-value')

    plt.tight_layout()
    fig.savefig('pvalues_%s_%s.png' % (gene_name, tissue_name), dpi=300)
    plt.cla()
    plt.clf()
    plt.close()


def plot_allele_count_distribution():
    # **** replace avg_length by two integer RU counts
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
        print(data[key].values())
        values = np.array(data[key].values()) / float(sum(data[key].values())) * 100
        print(values)
        plt.bar(np.array(data[key].keys()) + offset, values, label=key, width=0.2)
        offset += 0.2
    plt.legend()
    plt.xlabel('Number of alleles')
    plt.ylabel('Percentage of VNTRs')
    plt.savefig('VNTR_genotype_allele_count.png', dpi=200)
    plt.cla()
    plt.clf()
    plt.close()


def plot_vntr_polymorphic_rate_based_on_cohort_size():
    import seaborn as sns
    sns.set()
    sns.set_style("whitegrid")
    sns.set_context("poster")

    res = {}
    cohort_size = []
    prev_points = {'UTR':0, 'Coding':0, 'Promoter':0}
    for i in [i for i in range(651, 652, 10)]:# + [None]:
        # if i is not None:
        #     continue
        print(i)
        genotypes = load_individual_genotypes(reference_vntrs, True, i)
        genotyped_vntr_ids = {_id: set() for _id in range(1000000)}
        vntr_genotypes_map = {_id: [] for _id in range(1000000)}

        for individual_id, result_map in genotypes.items():
            for genotyped_vntr, avg_length in result_map.items():
                if avg_length == 'None' or avg_length is None:
                    continue
                vntr_genotypes_map[genotyped_vntr].append(avg_length)
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
            from collections import Counter
            g_dict = Counter(vntr_genotypes_map[vntr_id])
            count_of_second_frequent = sorted(g_dict.values(), reverse=True)[1] if len(g_dict) > 1 else 0
            # if number_of_genotypes > 1:
            count_of_most_frequent = vntr_genotypes_map[vntr_id].count(max(set(vntr_genotypes_map[vntr_id]), key = vntr_genotypes_map[vntr_id].count))
            # if count_of_most_frequent < i*0.95:
            if i - count_of_most_frequent > 5:
            # if count_of_second_frequent > 5:
                poly[annotation] += 1
            total[annotation] += 1

        rate = {key: float(poly[key]) / total[key] for key in poly.keys()}
        print(poly)
        for key in poly.keys():
            if key not in res.keys():
                res[key] = []
            res[key].append(max(rate[key], prev_points[key]))
            # prev_points[key] = max(rate[key], prev_points[key]*101/100)
        cohort_size.append(i if i != None else 652)
    for key in res.keys():
        plt.plot(cohort_size, res[key], label=key)
    plt.legend()
    plt.xlabel('Cohort Size')
    plt.ylabel('Rate of polymorphic VNTRs')
    plt.tight_layout()
    plt.savefig('polymorphic_vntrs_cohort_size_l5.png', dpi=300)
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
    plt.savefig('polymorphic_vntrs.png', dpi=100)
    plt.cla()
    plt.clf()
    plt.close()


def plot_evntrs_and_number_of_tissues():
    data = {}
    tissue_data_dir = glob.glob('pvalues/*')
    for tissue_dir in tissue_data_dir:
        vntrs_data = glob.glob(tissue_dir + '/*')
        tissue_name = os.path.basename(tissue_dir)
        for vntr_pvalue_dir in vntrs_data:
            vntr_id = int(os.path.basename(vntr_pvalue_dir))
            if vntr_id not in data.keys():
                data[vntr_id] = set([])
            data[vntr_id].add(tissue_name)
            pvalue_file = vntr_pvalue_dir + '/pvalues.txt'
            with open(pvalue_file) as infile:
                lines = infile.readlines()
                lines = [line.strip() for line in lines if line.strip() != '']

    data = {key: len(value) for key, value in data.items()}
    bar_data = {}
    for n_tissue in set(data.values()):
        res = 0
        for key in data.keys():
            if data[key] == n_tissue:
                res += 1
        bar_data[n_tissue] = res
    plt.bar(bar_data.keys(), bar_data.values())

    plt.xlabel('Shared across number of tissues')
    plt.ylabel('Number of VNTR-eQTLs')
    plt.savefig('eVNTR_uniqueness.png', dpi=300)
    plt.cla()
    plt.clf()
    plt.close()


def plot_evntrs_per_tissue():
    import matplotlib.pyplot as plt
    plt.style.use('ggplot')
    plt.rcParams['axes.facecolor'] = '#FFFFFF'
    plt.gca().spines['bottom'].set_color('black')
    plt.gca().spines['left'].set_color('black')
    plt.gca().spines['top'].set_color('black')
    plt.gca().spines['right'].set_color('black')

    import seaborn as sns
    sns.set()
    sns.set_style("whitegrid")
    sns.set_context("talk")

    data = {}
    tissue_data_dir = glob.glob('regression_results/*')
    for tissue_dir in tissue_data_dir:
        vntrs_data = glob.glob(tissue_dir + '/*')
        tissue_name = os.path.basename(tissue_dir).split('---')[0]
        if len(os.path.basename(tissue_dir).split('---')) > 1:
            tissue_name += ' - ' + os.path.basename(tissue_dir).split('---')[1][:5]
        if tissue_name not in data.keys():
            data[tissue_name] = set([])
        for vntr_pvalue_dir in vntrs_data:
            vntr_id = int(os.path.basename(vntr_pvalue_dir)[:-4])
            # data[tissue_name].add(vntr_id)
            # pvalue_file = vntr_pvalue_dir + '/pvalues.txt'
            pvalue_file = vntr_pvalue_dir
            with open(pvalue_file) as infile:
                line = infile.readlines()[0].strip()
                if float(line.split()[4]) < 0.0005:
                    data[tissue_name].add(vntr_id)

    data_count = {key: len(value) for key, value in data.items()}
    sorted_data = sorted(data_count.items(), key=operator.itemgetter(1), reverse=True)
    unique_counts = {}
    shared_by_tissues = {}
    tissue_order = []
    bottoms = []

    def autolabel(rects, offset):
        """
        Attach a text label above each bar displaying its height
        """
        for i, rect in enumerate(rects):
            # height = rect.get_height()
            plt.text(rect.get_x() + rect.get_width() / 2., 1.01 * (offset[i]),
                    '%d' % int(offset[i]),
                    ha='center', va='bottom', color='gray', size=5)

    for shared_counter in range(0, 2):#len(sorted_data)):
        shared_by_tissues[shared_counter] = []
        for key, value in sorted_data:
            value=data[key]
            res = 0
            for v in value:
                found = 0
                for other_tissue, vntrs in data.items():
                    if other_tissue == key:
                        continue
                    if v in vntrs:
                        found += 1
                if found == shared_counter or (shared_counter == 1 and found >= 1):
                    res += 1
            unique_counts[key] = res
            shared_by_tissues[shared_counter].append(res)
            if shared_counter == 0:
                tissue_order.append(key)
                bottoms.append(0)
        label = None
        if shared_counter == 0:
            label = 'Tissue Specific VNTR-eQTLs'
        if shared_counter == 1:
            label = 'Shared VNTR-eQTLs'
        bar = plt.bar(np.array([i for i in range(len(shared_by_tissues[shared_counter]))]), shared_by_tissues[shared_counter], label=label, width=0.6, bottom=bottoms)
        bottoms = [bottoms[i] + shared_by_tissues[shared_counter][i] for i in range(len(bottoms))]
    autolabel(bar, bottoms)

    plt.xticks(np.array([i for i in range(len(tissue_order))]), tissue_order, rotation=90, fontsize=6)
    plt.legend()
    plt.ylim(0, 60)
    plt.xlabel('Tissue name')
    plt.ylabel('Number of VNTR-eQTLs')
    plt.tight_layout(pad=2, w_pad=0, h_pad=2)
    plt.savefig('eVNTR_per_tissue.png', dpi=300)
    plt.cla()
    plt.clf()
    plt.close()


def plot_genotypes_difference_from_reference():
    # **** replace avg length with two individual RU counts
    # maybe also same plot for only eVNTRs and not all polymorphic VNTRs
    genotypes = load_individual_genotypes(reference_vntrs)
    genotyped_vntr_ids = {_id: set() for _id in range(1000000)}

    reference_repeats = {_id: [] for _id in range(1000000)}
    for individual_id, result_map in genotypes.items():
        for genotyped_vntr, avg_length in result_map.items():
            if avg_length == 'None' or avg_length is None:
                continue
            reference_repeats[genotyped_vntr].append(avg_length)
    for key, value in reference_repeats.items():
        if reference_repeats[key] == []:
            reference_repeats[key].append(0)
    reference_repeats = {key: most_common(value) for key, value in reference_repeats.items()}

    hist_data = []
    for individual_id, result_map in genotypes.items():
        for genotyped_vntr, avg_length in result_map.items():
            if avg_length == 'None' or avg_length is None:
                continue
            if avg_length != reference_repeats[genotyped_vntr]:
                hist_data.append(int(round(avg_length - reference_repeats[genotyped_vntr])))
            genotyped_vntr_ids[genotyped_vntr].add(avg_length)

    data = {point: 0 for point in range(min(hist_data), max(hist_data) + 1)}
    for h in hist_data:
        data[h] += 1
    del data[0]

    plt.bar(data.keys(), data.values(), width=0.8)
#    plt.hist(hist_data, 40)
#    plt.xlabel('Base-pairs difference from reference')
    plt.xlabel('Copy number difference from reference')
    plt.ylabel('Number of VNTRs')
    plt.savefig('difference_from_reference_all_genotypes.pdf')
    plt.cla()
    plt.clf()
    plt.close()


def most_common(lst):
    return max(set(lst), key=lst.count)


def plot_vntrs_difference_from_reference():
    # **** replace avg length with two individual RU counts
    # in this plot each specifi RU count of specific VNTR contributes once
    genotypes = load_individual_genotypes(reference_vntrs)
    genotyped_vntr_ids = {_id: set() for _id in range(1000000)}

    reference_repeats = {_id: [] for _id in range(1000000)}
    for individual_id, result_map in genotypes.items():
        for genotyped_vntr, avg_length in result_map.items():
            if avg_length == 'None' or avg_length is None:
                continue
            reference_repeats[genotyped_vntr].append(avg_length)
    for key, value in reference_repeats.items():
        if reference_repeats[key] == []:
            reference_repeats[key].append(0)
    reference_repeats = {key: most_common(value) for key, value in reference_repeats.items()}

    hist_data = []
    added = []
    for individual_id, result_map in genotypes.items():
        for genotyped_vntr, avg_length in result_map.items():
            if avg_length == 'None' or avg_length is None:
                continue
            if avg_length != reference_repeats[genotyped_vntr]:
                if not (genotyped_vntr, avg_length) in added:
                    hist_data.append(int(round(avg_length - reference_repeats[genotyped_vntr])))
                    added.append((genotyped_vntr, avg_length))
            genotyped_vntr_ids[genotyped_vntr].add(avg_length)

    data = {point: 0 for point in range(min(hist_data), max(hist_data) + 1)}
    for h in hist_data:
        data[h] += 1
    plt.bar(data.keys(), data.values(), width = 0.8)
#    plt.hist(hist_data, 40)
#    plt.xlabel('Base-pairs difference from reference')
    plt.xlabel('Copy number difference from reference')
    plt.ylabel('Number of VNTRs')
    plt.savefig('difference_from_reference_unique_genotypes.pdf')
    plt.cla()
    plt.clf()
    plt.close()


def plot_tissue_correlations():
    # cbarlabel = "Pearson Correlation"
    cbarlabel = "Jaccard Index"
    # colormap = "YlGn"
    colormap = "Binary"
    tissue_labels = []
    data = []

    eVNTRs = set([])
    eVNTRs_of_tissue = {}

    tissue_data_dir = glob.glob('all_vntrs_zscore/*')
    tissue_data_dir = glob.glob('caviar_inputs/*')
    print(tissue_data_dir)
    for tissue_dir in tissue_data_dir:
        vntr_dirs = glob.glob(tissue_dir + '/*')
        tissue_name = os.path.basename(tissue_dir)
        tissue_labels.append(tissue_name)
        if tissue_name not in eVNTRs_of_tissue.keys():
            eVNTRs_of_tissue[tissue_name] = {}
        for vntr_dir in vntr_dirs:
            vntr_id = int(os.path.basename(vntr_dir))
            eVNTRs.add(vntr_id)
            # effect_size_file = 'all_vntrs_zscore/%s/%s/Z.txt' % (tissue_name, vntr_id)
            effect_size_file = 'caviar_inputs/%s/%s/%s_Z' % (tissue_name, vntr_id, tissue_name)
            with open(effect_size_file) as infile:
                lines = infile.readlines()
            # effect_size = float(lines[0].strip())
            effect_size = float(lines[0].strip().split()[1])
            eVNTRs_of_tissue[tissue_name][vntr_id] = effect_size

    def get_eVNTRs_list_of_tissue(tissue_name):
        res = []
        for e in eVNTRs:
            if e in eVNTRs_of_tissue[tissue_name].keys():
                res.append(eVNTRs_of_tissue[tissue_name][e])
            else:
                res.append(0)
        return res

    for tissue1 in tissue_labels:
        row = []
        for tissue2 in tissue_labels:
            e1 = set(eVNTRs_of_tissue[tissue1].keys())
            e2 = set(eVNTRs_of_tissue[tissue2].keys())
            row.append(float(len(e1 & e2)) / len(e1 | e2))
            # row.append(np.corrcoef(get_eVNTRs_list_of_tissue(tissue1), get_eVNTRs_list_of_tissue(tissue2))[0, 1])
        data.append(row)

    data = np.array(data)
    # Plot the heatmap
    ax = plt.gca()
    im = ax.imshow(data, cmap='Greys')

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, cmap=colormap)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

    # We want to show all ticks...
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_yticks(np.arange(data.shape[0]))
    # ... and label them with the respective list entries.
    ax.set_xticklabels(tissue_labels, fontsize=7)
    ax.set_yticklabels(tissue_labels, fontsize=7)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=True, bottom=False,
                   labeltop=True, labelbottom=False)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=-60, ha="right",
             rotation_mode="anchor")

    # Turn spines off and create white grid.
    for edge, spine in ax.spines.items():
        spine.set_visible(False)

    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)

    plt.tight_layout(pad=0.2, w_pad=0, h_pad=0)

    plt.savefig('tissue_correlation.pdf')
    plt.cla()
    plt.clf()
    plt.close()


def plot_vntr_pvalues_qq_plot():
    import seaborn as sns
    sns.set()
    sns.set_style("whitegrid")
    sns.set_context("paper")

    pvalues = {}
    tissue_data_dir = glob.glob('all_vntr_pvalues/*')
    for tissue_dir in tissue_data_dir:
        vntrs_data = glob.glob(tissue_dir + '/*')
        tissue_name = os.path.basename(tissue_dir)
        if tissue_name not in pvalues.keys():
            pvalues[tissue_name] = []
        for vntr_pvalue_dir in vntrs_data:
            vntr_id = int(os.path.basename(vntr_pvalue_dir))
            pvalue_file = vntr_pvalue_dir + '/pvalues.txt'
            with open(pvalue_file) as infile:
                lines = infile.readlines()
                lines = [line.strip() for line in lines if line.strip() != '']
#                pvalue = -log10(float(lines[0])) if float(lines[0]) != 0 else 100
                pvalue = float(lines[0])
                pvalues[tissue_name].append(pvalue)

    colors = ['r', 'b', 'g', 'c', 'y'] + ['k'] * 100
    cmap = plt.get_cmap("tab20c")
#    plt.clf()
#    plt.gcf().set_size_inches(6, 6)
    for i, key in enumerate(pvalues.keys()):
        key = 'Blood Vessel'
        N = float(len(pvalues[key]))
        print(min(pvalues[key]), max(pvalues[key]))
        grid = -np.log10(np.arange(1, 1 + N) / N)
        plt.scatter(grid, -np.log10(np.array(sorted(pvalues[key], reverse=False))), c=cmap(i))
#        plt.plot(grid, -np.log10(np.array(sorted(pvalues[key], reverse=False))), c=cmap(i))
#        sm.qqplot(np.array(pvalues[key]))
#        plt.hist(pvalues[key], 50)
        break
#    plt.ylim((0, 20))
#    plt.xlim((0, 10))
    plt.xlabel('Expected p-value [-log10(p)]')
    plt.ylabel('Observed p-value [-log10(p)]')
    plt.savefig('qq_plot_Blood-Vessel.png', dpi=300)
    plt.cla()
    plt.clf()
    plt.close()

def plot_vntr_location_effect():
    # not complete
    data = {}
    tissue_data_dir = glob.glob('pvalues/*')
    for tissue_dir in tissue_data_dir:
        vntrs_data = glob.glob(tissue_dir + '/*')
        tissue_name = os.path.basename(tissue_dir)
        if tissue_name not in data.keys():
            data[tissue_name] = set([])
        for vntr_pvalue_dir in vntrs_data:
            vntr_id = int(os.path.basename(vntr_pvalue_dir))
            data[tissue_name].add(vntr_id)
            pvalue_file = vntr_pvalue_dir + '/pvalues.txt'
            with open(pvalue_file) as infile:
                lines = infile.readlines()
                lines = [line.strip() for line in lines if line.strip() != '']

    plt.xlabel('Distance to transcription start site') # or gene start, first exon, or ...
    plt.ylabel('Number of VNTR-eQTLs') # or p-value, effect size, ...
    plt.savefig('distance_from_transcription_site.pdf')
    plt.cla()
    plt.clf()
    plt.close()


def plot_significant_vntrs_and_tissues():
    # show which tissues are affected for a list of 20-30 VNTRs (end of page11 in Melissa's paper)
    pass


def plot_pvalue_genomic_locations_manhattan_plot():
    pvalues = {}
    tissue_data_dir = glob.glob('all_vntr_pvalues/*')
    for tissue_dir in tissue_data_dir:
        vntrs_data = glob.glob(tissue_dir + '/*')
        tissue_name = os.path.basename(tissue_dir)
        for vntr_pvalue_dir in vntrs_data:
            vntr_id = int(os.path.basename(vntr_pvalue_dir))
            if vntr_id not in pvalues.keys():
                pvalues[vntr_id] = 0
            pvalue_file = vntr_pvalue_dir + '/pvalues.txt'
            with open(pvalue_file) as infile:
                lines = infile.readlines()
                lines = [line.strip() for line in lines if line.strip() != '']
#                pvalue = -log10(float(lines[0])) if float(lines[0]) != 0 else 100
                pvalue = -log10(float(lines[0]))
                pvalues[vntr_id] = max(pvalue, pvalues[vntr_id])

    cmap = plt.get_cmap("tab20")
    offsets = {1: 0}
    lengths = {}
    from advntr.utils import get_chromosome_reference_sequence
    for i in range(1, 23):
        lengths[i] = float(len(get_chromosome_reference_sequence('chr%s' % i)))
        if i == 1:
            continue
        offsets[i] = offsets[i-1] + lengths[i-1]
    print(offsets)
    total = float(sum(offsets.values()))
    for chrom in offsets.keys():
        offsets[chrom] = offsets[chrom]

    sara = 0
    points = []
    for vntr_id in pvalues.keys():
        chrom = reference_vntrs[vntr_id].chromosome[3:]
        if chrom in ['Y', 'X']:
            continue
        chrom = int(chrom)
        offset = offsets[chrom]
        y = pvalues[vntr_id]
        if y > 2.6989:
            sara += 1
        x = (reference_vntrs[vntr_id].start_point) + offset
        points.append((x, y, cmap((chrom % 2) * 6)))

    print(sara)
    points = sorted(points)
    X = [x for x, y, c in points]
    Y = [y for x, y, c in points]
    C = [c for x, y, c in points]
    plt.scatter(X, Y, marker='.', c=C)
    plt.plot([-1*1e8, 3*1e9], [2.6989, 2.6989], '--', c='k')

    plt.ylim(0, 18)
    plt.xlabel('Genomic location')
    plt.ylabel('-log(p-value)') # or p-value, effect size, ...
    plt.savefig('genomic_location_manhattan_plot.pdf')
    plt.cla()
    plt.clf()
    plt.close()

#GWAS result? correlation with traits? like hight

if __name__ == '__main__':
    # plot_pvalue_genomic_locations_manhattan_plot()
    # plot_allele_count_distribution()

    # plot_expression_genotype_correlation(157225, 'Brain')
    # plot_variant_pvalues(157225, 'Brain')
    # plot_variant_caviar_scores(157225, 'Brain')
    #
    # exit(0)
    # plot_expression_genotype_correlation(331737, 'Heart')
    # plot_variant_pvalues(331737, 'Heart')
    # plot_variant_caviar_scores(331737, 'Heart')
    #
    # plot_expression_genotype_correlation(393574, 'Thyroid')
    # plot_variant_pvalues(393574, 'Thyroid')
    # plot_variant_caviar_scores(393574, 'Thyroid')
    #
    # plot_expression_genotype_correlation(690585, 'Blood Vessel')
    # plot_variant_pvalues(690585, 'Blood Vessel')
    # plot_variant_caviar_scores(690585, 'Blood Vessel')
    #
    # plot_expression_genotype_correlation(157225, 'Salivary Gland')
    # plot_variant_pvalues(157225, 'Salivary Gland')
    # plot_variant_caviar_scores(157225, 'Salivary Gland')
    #
    # plot_expression_genotype_correlation(668609, 'Lung')
    # plot_variant_pvalues(668609, 'Lung')
    # plot_variant_caviar_scores(668609, 'Lung')
    #
    # exit(0)
    # plot_tissue_correlations()
    # plot_variant_pvalues(690585, 'Blood Vessel')
    # plot_variant_caviar_scores(690585, 'Blood Vessel')
    # plot_variant_pvalues(336040, 'Esophagus')
    # plot_variant_caviar_scores(336040, 'Esophagus')
    # plot_variant_pvalues(393574, 'Thyroid')
    # plot_variant_caviar_scores(393574, 'Thyroid')


    # plot_dose_dependence_of_correlation()
    # plot_genotypes_difference_from_reference()
    # plot_vntrs_difference_from_reference()
    # plot_vntr_pvalues_qq_plot()
    # plot_evntrs_per_tissue()
    # plot_evntrs_and_number_of_tissues()

    # plot_expression_genotype_correlation(331737, 'Heart')
    # plot_expression_genotype_correlation(393574, 'Thyroid')
    # plot_expression_genotype_correlation(111235, 'Esophagus')
    # plot_expression_genotype_correlation(111235, 'Blood Vessel')
    # plot_expression_genotype_correlation(9635, 'Blood Vessel')
    # plot_expression_genotype_correlation(157225, 'Brain')
    # plot_expression_genotype_correlation(690585, 'Blood Vessel')
    # plot_expression_genotype_correlation(450457, 'Pancreas')
    # plot_expression_genotype_correlation(315000, 'Blood')
    # plot_expression_genotype_correlation(336040, 'Esophagus')
    # plot_expression_genotype_correlation(668609, 'Lung')
    # plot_expression_genotype_correlation(184368, 'Ovary')
    # plot_expression_genotype_correlation(356726, 'Heart')
    # plot_expression_genotype_correlation(141097, 'Thyroid')
    # plot_expression_genotype_correlation(356297, 'Nerve')
    # plot_expression_genotype_correlation(630254, 'Nerve')
    # plot_expression_genotype_correlation(309308, 'Spleen')
    # plot_expression_genotype_correlation(410522, 'Nerve')
    # plot_expression_genotype_correlation(209165, 'Thyroid')
    # plot_expression_genotype_correlation(30338, 'Small Intestine')
    # plot_expression_genotype_correlation(124123, 'Esophagus')

    # plot_variant_pvalues(111235, 'Esophagus')
    # plot_variant_caviar_scores(111235, 'Esophagus')
    # plot_variant_pvalues(393574, 'Thyroid')
    # plot_variant_caviar_scores(393574, 'Thyroid')
#    plot_variant_pvalues(101865, 'Heart')
#    plot_variant_pvalues(331737, 'Heart')
#    plot_variant_caviar_scores(111235, 'Esophagus')
#    plot_variant_caviar_scores(331737, 'Esophagus')

    # plot_allele_count_distribution()
    # plot_vntr_polymorphic_rate_based_on_annotation()
    plot_vntr_polymorphic_rate_based_on_cohort_size()

