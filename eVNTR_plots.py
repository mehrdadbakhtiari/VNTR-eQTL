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
    if vntr_id == 450457:
        gene_name = 'RPA2'
    data_file = 'genotype_expression_correlation/%s/%s/correlation.txt' % (tissue_name, vntr_id)
    if not os.path.exists(data_file):
        return
    with open(data_file) as infile:
        lines = infile.readlines()
    data = []
    labels = []
    overall = []
    normalize_labels = []
    for line in lines:
        print(line.strip())
        line = line.strip().split(' ')
        if len(line[1].split(',')) < 5:
            print('skip')
            continue
        labels.append(str(line[0]))
        data.append([float(e) for e in line[1].split(',')])
        overall += data[-1]
        normalize_labels += [len(labels)-1 for _ in range(len(data[-1]))]
    from sklearn.preprocessing import quantile_transform
    overall = [[e] for e in overall]
    overall = quantile_transform(np.array(overall), n_quantiles=10, random_state=0, copy=True)
    for i in range(len(data)):
        data[i] = [overall[j][0] for j in range(len(overall)) if normalize_labels[j] == i]

    print('points', sum([len(data[i]) for i in range(len(data))]))
    import seaborn as sns
    sns.set()
    sns.set_style("whitegrid")
    sns.set_context("paper")

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

    print(tissue_name, pvalue)
    # return
    fig = plt.figure(figsize=(6, 4))
    ax = fig.add_subplot(111)
    # data = np.array(data).transpose()
    if len(labels) < 4:
        sns.set_palette(sns.color_palette([sns.xkcd_rgb["denim blue"], sns.xkcd_rgb['medium green'], sns.xkcd_rgb['pale red']]))
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
    #ax.set_ylim(0)
    plt.tight_layout()

    plt.gca().spines['bottom'].set_color('black')
    plt.gca().spines['left'].set_color('black')
    plt.gca().spines['top'].set_color('black')
    plt.gca().spines['right'].set_color('black')

    fig.savefig('expression_%s_%s.pdf' % (gene_name, tissue_name))
    plt.cla()
    plt.clf()
    plt.close()


def plot_dose_dependence_of_correlation(use_beta=True):
    import seaborn as sns
    sns.set()
    sns.set_style("whitegrid")
    sns.set_context("paper")

    significant_vntrs = set([])

    tissue_data_dir = glob.glob('pvalues/*')
    for tissue_dir in tissue_data_dir:
        vntrs_data = glob.glob(tissue_dir + '/*')
        tissue_name = os.path.basename(tissue_dir)
        for vntr_pvalue_dir in vntrs_data:
            vntr_id = int(os.path.basename(vntr_pvalue_dir))
            significant_vntrs.add((tissue_name.replace(' ', '-'), vntr_id))
            # significant_vntrs.add(vntr_id)
    print(len(significant_vntrs))

    increasing = 0
    decreasing = 0
    hist_data = []
    tissue_data_dir = glob.glob('caviar_inputs/*')
    print(tissue_data_dir)
    for tissue_dir in tissue_data_dir:
        vntr_dirs = glob.glob(tissue_dir + '/*')
        tissue_name = os.path.basename(tissue_dir)
        for vntr_dir in vntr_dirs:
            vntr_id = int(os.path.basename(vntr_dir))
            effect_size_file = 'caviar_inputs/%s/%s/%s_Z' % (tissue_name, vntr_id, tissue_name)
            with open(effect_size_file) as infile:
                lines = infile.readlines()
            effect_size = float(lines[0].strip().split()[1])
            hist_data.append(effect_size)
            if effect_size > 0:
                increasing += 1
            else:
                decreasing += 1

    if use_beta:
        hist_data = []
        increasing = 0
        decreasing = 0
        tissue_data_dir = glob.glob('regression_results/*')
        print(tissue_data_dir)
        for tissue_dir in tissue_data_dir:
            print(tissue_dir)
            tissue_name = os.path.basename(tissue_dir)
            vntr_files = glob.glob(tissue_dir + '/*')
            for vntr_file in vntr_files:
                vntr_id = int(os.path.basename(vntr_file[:-4]))
                if (tissue_name, vntr_id) not in significant_vntrs:
                    continue
                with open(vntr_file) as infile:
                    lines = infile.readlines()
                effect_size = float(lines[0].strip().split('\t')[3])
                hist_data.append(effect_size)
                if effect_size > 0:
                    increasing += 1
                elif effect_size < -0:
                    decreasing += 1

    print(increasing, decreasing)
    from pylab import rcParams
    # rcParams['figure.figsize'] = 8, 5

    # plt.bar([0, 1], [increasing, decreasing], width=0.5)
    # plt.xticks((0, 1), ['Increasing', 'Decreasing'], size=12)
    data = {}
    print(data)
    for h in hist_data:
        h = round(10 * h) / 10.0
        if h not in data.keys():
            data[h] = 0
        data[h] += 1
    print(data)
    # plt.bar(data.keys(), data.values(), width=0.10)
    # plt.hist(hist_data, bins=100, histtype='step')
    sns.set_palette(sns.hls_palette(8, l=.3, s=.8))
    sns.set_palette(sns.color_palette("dark", 40))
    sns.distplot(hist_data, bins=100)
    plt.xlabel('Dose Dependence', size=14)
    plt.ylabel('Significant associations', size=14)
    # plt.xlim(-0.6, 1.6)
    plt.tight_layout()
    plt.gca().spines['bottom'].set_color('black')
    plt.gca().spines['left'].set_color('black')
    plt.gca().spines['top'].set_color('black')
    plt.gca().spines['right'].set_color('black')

    plt.savefig('dose_dependence.png', dpi=300)
    plt.cla()
    plt.clf()
    plt.close()


def plot_variant_caviar_scores(vntr_id, tissue_name, xlim=None):
    caviar_post_file = 'caviar_inputs/%s/%s/caviar_post' % (tissue_name.replace(' ', '-'), vntr_id)
    post = pd.read_csv(caviar_post_file, sep="\t", header=0)
    post = post.sort_values(post.columns[2], ascending=False)
    post = post.reset_index(drop=True)

    vntr_x = reference_vntrs[vntr_id].start_point
    # vntr_x = 0
    vntr_y = 0
    x = []
    y = []
    for row in range(0,len(post.index)):
        if '%s' % vntr_id in post.values[row][0]:
            vntr_y = float(post.values[row][2])
        else:
            x.append(int(post.values[row][0].split('_')[1]))
            y.append(post.values[row][2])
    # vntr_x = get_vntr_x(vntr_id, x)
    print(vntr_x)
    gene_name = reference_vntrs[vntr_id].gene_name


    fig = plt.figure(figsize=(6, 1.8))
    ax = fig.add_subplot(111)
    ax.scatter(x, y, marker='.', c='gray')
    ax.scatter(vntr_x, vntr_y, marker='*', c='r')
    ax.set_ylabel('Causality probability')

    if xlim:
        ax.set_xlim(xlim[0], xlim[1])
    plt.tight_layout()

    plt.gca().spines['bottom'].set_color('black')
    plt.gca().spines['left'].set_color('black')
    plt.gca().spines['top'].set_color('black')
    plt.gca().spines['right'].set_color('black')

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
    vntr_x = reference_vntrs[vntr_id].start_point
    # vntr_x = 0
    vntr_y = 0
    for line in lines:
        snp, p = line.split()
        if '%s' % vntr_id in snp:
            vntr_y = -log10(float(p))
        else:
            y.append(-log10(float(p)))
            x.append(float(snp.split('_')[1]))
    # vntr_x = get_vntr_x(vntr_id, x)
    gene_name = reference_vntrs[vntr_id].gene_name

    # import seaborn as sns
    # sns.set()
    # sns.set_style("whitegrid")
    # sns.set_context("talk")

    fig = plt.figure(figsize=(6, 1.8))
    ax = fig.add_subplot(111)
    ax.scatter(x, y, marker='.', c='gray')
    ax.scatter(vntr_x, vntr_y, marker='*', c='r')
    # ax.set_xlabel('Genomic location')
    ax.set_ylabel('-log10 P-value')

    plt.tight_layout()
    xlim = plt.xlim()

    plt.gca().spines['bottom'].set_color('black')
    plt.gca().spines['left'].set_color('black')
    plt.gca().spines['top'].set_color('black')
    plt.gca().spines['right'].set_color('black')

    fig.savefig('pvalues_%s_%s.png' % (gene_name, tissue_name), dpi=300)
    plt.cla()
    plt.clf()
    plt.close()
    return xlim


def plot_pvalue_and_caviar(vntr_id, tissue_name):
    xlim = plot_variant_pvalues(vntr_id, tissue_name)
    plot_variant_caviar_scores(vntr_id, tissue_name, xlim=xlim)


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


def plot_effect_per_allele_frequency():
    import matplotlib.pyplot as plt
    # plt.style.use('ggplot')
    # plt.rcParams['axes.facecolor'] = '#FFFFFF'

    import seaborn as sns
    sns.set()
    sns.set_style("whitegrid")
    sns.set_context("paper")

    plt.gca().spines['bottom'].set_color('black')
    plt.gca().spines['left'].set_color('black')
    plt.gca().spines['top'].set_color('black')
    plt.gca().spines['right'].set_color('black')

    genotypes = load_individual_genotypes(reference_vntrs, False)
    genotyped_vntr_ids = {_id: set() for _id in range(1000000)}
    vntr_genotypes_map = {_id: [] for _id in range(1000000)}

    for individual_id, result_map in genotypes.items():
        for genotyped_vntr, alleles in result_map.items():
            if alleles == 'None' or alleles is None or None in alleles:
                continue
            vntr_genotypes_map[genotyped_vntr] += alleles
            for allele in alleles:
                genotyped_vntr_ids[genotyped_vntr].add(allele)
    vntr_genotypes = {_vntr_id: len(lengths) for _vntr_id, lengths in genotyped_vntr_ids.items() if len(lengths) != 0}

    maf = {}
    for vntr_id, number_of_genotypes in vntr_genotypes.items():
        from collections import Counter
        g_dict = Counter(vntr_genotypes_map[vntr_id])
        count_of_second_frequent = sorted(g_dict.values(), reverse=True)[1] if len(g_dict) > 1 else 0
        count_of_most_frequent = vntr_genotypes_map[vntr_id].count(
            max(set(vntr_genotypes_map[vntr_id]), key=vntr_genotypes_map[vntr_id].count))
        maf[vntr_id] = count_of_second_frequent / float(len(vntr_genotypes_map[vntr_id]))

    import math
    print('computed minor allele frequencies')

    thresholds = {}
    with open('thresholds.txt') as infile:
        lines = infile.readlines()
        for l in lines:
            thresholds[l.split('\t')[0].replace(' ', '-')] = float(l.split('\t')[1])

    x = []
    y = []
    xpvalue = []
    ypvalue = []
    gene_names = []
    vntr_bestpvalue = {}
    vntr_bestbeta = {}
    tissue_data_dir = glob.glob('regression_results/*')
    for tissue_dir in tissue_data_dir:
        vntrs_data = glob.glob(tissue_dir + '/*')
        tissue_name = os.path.basename(tissue_dir).split('---')[0]
        if len(os.path.basename(tissue_dir).split('---')) > 1:
            tissue_name += ' - ' + os.path.basename(tissue_dir).split('---')[1][:5]
        for vntr_pvalue_dir in vntrs_data:
            vntr_id = int(os.path.basename(vntr_pvalue_dir)[:-4])
            if maf[vntr_id] not in vntr_bestpvalue.keys():
                vntr_bestpvalue[maf[vntr_id]] = 0
            if maf[vntr_id] not in vntr_bestbeta.keys():
                vntr_bestbeta[maf[vntr_id]] = 0
            pvalue_file = vntr_pvalue_dir
            with open(pvalue_file) as infile:
                line = infile.readlines()[0].strip()
                effect_size = float(line.split()[3])
                pvalue = float(line.split()[4])

                gene_names.append(reference_vntrs[vntr_id].gene_name)
                x.append(maf[vntr_id])
                y.append(abs(effect_size))
                vntr_bestbeta[maf[vntr_id]] = max(vntr_bestbeta[maf[vntr_id]], abs(effect_size))
                xpvalue.append(maf[vntr_id])
                ypvalue.append(-math.log(pvalue, 10))
                vntr_bestpvalue[maf[vntr_id]] = max(vntr_bestpvalue[maf[vntr_id]], -math.log(pvalue, 10))
    plt.ylabel('Effect Size')
    plt.xlabel('MAF: minor allele frequency (second most common allele)')
    # x = []
    # y = []
    # for key, value in vntr_bestbeta.items():
    #     x.append(maf[key])
    #     y.append(value)
    def is_peak(maf, beta):
        # beta = vntr_bestbeta[maf]
        peaks = []
        for _maf, _beta in vntr_bestbeta.items():
            if _maf == maf:
                continue
            if maf - 0.05 < _maf < maf + 0.05:
                peaks.append(_beta)
        return beta > 2.0 * sum(peaks) / float(len(peaks))
    # colors = []
    # for i in range(len(x)):
    #     if x[i] > 0.05 and is_peak(x[i], y[i]):
    #         colors.append('tab:red')
    #     else:
    #         colors.append('tab:blue')
    plt.scatter(x, y, marker='.', alpha=0.4)
    for i in range(len(x)):
        if x[i] > 0.05 and y[i] == vntr_bestbeta[x[i]] and (y[i] > 0.45 or (y[i] > 0.29 and x[i] > 0.1)):
            if x[i] > 0.17 or y[i] > 0.48 or gene_names[i] in ['FAM71E2', 'EPS8L2', 'SRRD', 'RAET1G', 'ANKRD30BL']:
                plt.text(x[i], y[i], gene_names[i], {'ha': 'left', 'va': 'bottom'}, rotation=45, fontsize=6)
            else:
                print(x[i], y[i], gene_names[i])
                #if gene_names[i] == 'FAM71E2':
                #    text_pos = (x[i] + 0.05, y[i] + 0.07)
                #    arrow_end = (x[i] - text_pos[0] + 0.003, y[i] - text_pos[1] + 0.005)
                if gene_names[i] == 'ZNF232':
                    text_pos = (x[i], y[i] + 0.1)
                    arrow_end = (x[i] - text_pos[0], y[i] - text_pos[1] + 0.01)
                elif gene_names[i] == 'CRYBB2P1':
                    text_pos = (x[i] + 0.03, y[i]+0.08)
                    arrow_end = (x[i] - text_pos[0] + 0.001, y[i] - text_pos[1] + 0.008)
                # plt.annotate(gene_names[i], xy=(x[i], y[i]), xytext=text_pos, fontsize=6, arrowprops=dict(facecolor='black', shrink=0.05))
                arrow_color = 'k'
                plt.arrow(text_pos[0], text_pos[1], arrow_end[0], arrow_end[1], fc=arrow_color, ec=arrow_color, width=0.001, lw=0.5)
                plt.text(text_pos[0], text_pos[1], gene_names[i], {'ha': 'left', 'va': 'bottom'}, rotation=45, fontsize=6)
    # sns.distplot(x)
    # sns.kdeplot(x)
    plt.tight_layout(pad=2, w_pad=0, h_pad=2)
    plt.savefig('maf_effect_size.png', dpi=300)

    plt.cla()
    plt.clf()
    plt.close()

    plt.ylabel('$\\mathdefault{-log_{10}(Pvalue)}$', fontsize=12)
    plt.xlabel('MAF: minor allele frequency (second most common allele)', fontsize=12)
    # xpvalue = []
    # ypvalue = []
    # for key, value in vntr_bestpvalue.items():
    #     xpvalue.append(maf[key])
    #     ypvalue.append(value)
    plt.scatter(xpvalue, ypvalue, marker='.', alpha=0.4)

    for i in range(len(xpvalue)):
        if xpvalue[i] > 0.05 and ypvalue[i] == vntr_bestpvalue[xpvalue[i]] and ypvalue[i] > 10:
            print(gene_names[i])
            if gene_names[i] not in ['RMDN1', 'CADM1', 'COLQ', 'SNHG16', 'ANKEF1', 'ZNF736', 'CCDC57', 'PTTG1', 'EPDR1', 'ALKBH5', 'BAIAP2L2']:
                plt.text(xpvalue[i], ypvalue[i], gene_names[i], {'ha': 'left', 'va': 'bottom'}, rotation=45, fontsize=6)
            else:
                print(xpvalue[i], ypvalue[i], gene_names[i])
                if gene_names[i] == 'RMDN1':
                    text_pos = (x[i], ypvalue[i] + 3)
                    arrow_end = (x[i] - text_pos[0], ypvalue[i] - text_pos[1] + 1)
                    plt.text(xpvalue[i]-0.005, ypvalue[i]+0.04, gene_names[i], {'ha': 'left', 'va': 'bottom'}, rotation=45,
                             fontsize=6)
                    continue
                elif gene_names[i] == 'CADM1':
                    plt.text(xpvalue[i] + 0.0003, ypvalue[i] - 1, gene_names[i], {'ha': 'left', 'va': 'bottom'},
                             rotation=45, fontsize=6)
                    continue
                elif gene_names[i] == 'COLQ':
                    plt.text(xpvalue[i] - 0.003, ypvalue[i] - 0.5, gene_names[i], {'ha': 'left', 'va': 'bottom'},
                             rotation=45, fontsize=6)
                    continue
                elif gene_names[i] == 'BAIAP2L2' and xpvalue[i] < 0.37:
                    plt.text(xpvalue[i], ypvalue[i], gene_names[i], {'ha': 'left', 'va': 'bottom'}, rotation=45,
                             fontsize=6)
                    continue
                elif gene_names[i] == 'BAIAP2L2':
                    text_pos = (x[i] + 0.020, ypvalue[i]+8)
                    arrow_end = (x[i] + 0.0020, ypvalue[i]+0.6)
                elif gene_names[i] == 'SNHG16':
                    text_pos = (x[i], ypvalue[i] + 8)
                    arrow_end = (x[i], ypvalue[i] + 0.7)
                elif gene_names[i] == 'ANKEF1':
                    text_pos = (x[i] + 0.022, ypvalue[i]+4)
                    arrow_end = (x[i] + 0.0025, ypvalue[i]+0.03)
                elif gene_names[i] == 'ZNF736':
                    text_pos = (x[i] + 0.02, ypvalue[i] - 1.5)
                    arrow_end = (x[i] + 0.0025, ypvalue[i]-0.02)
                elif gene_names[i] == 'CCDC57':
                    text_pos = (x[i] + 0.013, ypvalue[i] + 2)
                    arrow_end = (x[i] + 0.0025, ypvalue[i]+0.02)
                elif gene_names[i] == 'PTTG1':
                    text_pos = (x[i], ypvalue[i] + 4)
                    arrow_end = (x[i], ypvalue[i] + 0.7)
                elif gene_names[i] == 'EPDR1':
                    text_pos = (x[i], ypvalue[i] + 12)
                    arrow_end = (x[i], ypvalue[i] + 0.7)
                elif gene_names[i] == 'ALKBH5':
                    text_pos = (x[i] + 0.02, ypvalue[i])
                    arrow_end = (x[i] + 0.0025, ypvalue[i])
                else:
                    continue
                # plt.annotate(gene_names[i], xy=(x[i], y[i]), xytext=text_pos, fontsize=6, arrowprops=dict(facecolor='black', shrink=0.05))
                arrow_color = 'k'
                plt.annotate(gene_names[i], xy=arrow_end, xytext=text_pos, va='bottom', bbox=dict(fc='none', ec='none', pad=0), rotation=45,
                             arrowprops=dict(facecolor='black', shrink=0, width=1.5, headlength=5, headwidth=4.5), fontsize=6)
                #plt.arrow(text_pos[0], text_pos[1], arrow_end[0], arrow_end[1], fc=arrow_color, ec=arrow_color, width=0.002)
                #plt.text(text_pos[0], text_pos[1], gene_names[i], {'ha': 'left', 'va': 'bottom'}, rotation=45, fontsize=6)

    #plt.plot([0, 0.5], [-math.log(0.0005, 10), -math.log(0.0005, 10)], '--', color='tab:red')
    lower = -math.log(0.00101531415963, 10)
    upper = -math.log(3.79592740513e-05, 10)
#    plt.axhspan(lower, upper, facecolor='tab:red', alpha=0.4, xmin=0, xmax=0.5)
    plt.fill([0,0.5,0.5,0], [lower,lower,upper,upper], 'tab:red', alpha=0.3, edgecolor='red')
    plt.tight_layout(pad=2, w_pad=0, h_pad=2)
    plt.savefig('maf_pvalue.pdf')


def plot_evntrs_per_tissue():
    import matplotlib.pyplot as plt
#    plt.style.use('ggplot')
    plt.rcParams['axes.facecolor'] = '#FFFFFF'
    plt.gca().spines['bottom'].set_color('black')
    plt.gca().spines['left'].set_color('black')
    plt.gca().spines['top'].set_color('black')
    plt.gca().spines['right'].set_color('black')

    import seaborn as sns
    sns.set()
    sns.set_style("whitegrid")
    sns.set_context("paper")

    thresholds = {}
    with open('thresholds.txt') as infile:
        lines = infile.readlines()
        for l in lines:
            thresholds[l.split('\t')[0].replace(' ', '-')] = float(l.split('\t')[1])

    data = {}
    tissue_data_dir = glob.glob('regression_results/*')
    for tissue_dir in tissue_data_dir:
        vntrs_data = glob.glob(tissue_dir + '/*')
        threshold = thresholds[os.path.basename(tissue_dir)]
        tissue_name = os.path.basename(tissue_dir).split('---')[0].replace('-', ' ')
        if len(os.path.basename(tissue_dir).split('---')) > 1:
            tissue_name += ' - ' + os.path.basename(tissue_dir).split('---')[1].replace('-', ' ')
        if tissue_name not in data.keys():
            data[tissue_name] = set([])
        for vntr_pvalue_dir in vntrs_data:
            vntr_id = int(os.path.basename(vntr_pvalue_dir)[:-4])
            # data[tissue_name].add(vntr_id)
            # pvalue_file = vntr_pvalue_dir + '/pvalues.txt'
            pvalue_file = vntr_pvalue_dir
            with open(pvalue_file) as infile:
                line = infile.readlines()[0].strip()
                if float(line.split()[4]) <= threshold:
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
        red_color = '#d62728'
        blue_color = "#1f77b4"
        if shared_counter == 0:
            label = 'Tissue Specific VNTR-eQTLs'
            color = red_color
        if shared_counter == 1:
            label = 'Shared VNTR-eQTLs'
            color = blue_color
        bar = plt.bar(np.array([i for i in range(len(shared_by_tissues[shared_counter]))]), shared_by_tissues[shared_counter], label=label, width=0.6, bottom=bottoms, color=color)
        bottoms = [bottoms[i] + shared_by_tissues[shared_counter][i] for i in range(len(bottoms))]
        print(bottoms)
        print(sum(bottoms))
    autolabel(bar, bottoms)
    print(shared_by_tissues)

    for i in range(len(tissue_order)):
        tissue_order[i] = tissue_order[i].replace('Gastroesophageal Junction', 'GE Junction')
        tissue_order[i] = tissue_order[i].replace('Cells - ', '')
        p = tissue_order[i].find('(')
        if p != -1:
            tissue_order[i] = tissue_order[i][:p].strip()

    plt.xticks(np.array([i for i in range(len(tissue_order))]), tissue_order, rotation=90, fontsize=6)
    plt.legend()
    plt.ylim(0, 45)
    plt.xlabel('Tissue name')
    plt.ylabel('Number of VNTR-eQTLs')
    plt.tight_layout(pad=2, w_pad=0, h_pad=2)
    plt.savefig('eVNTR_per_tissue.png', dpi=300)
    plt.cla()
    plt.clf()
    plt.close()


def plot_affected_bp_per_individual():
#    import matplotlib.pyplot as plt

    import seaborn as sns
    sns.set()
    sns.set_style("whitegrid")
    sns.set_context("paper")

    plt.gca().spines['bottom'].set_color('black')
    plt.gca().spines['left'].set_color('black')
    plt.gca().spines['top'].set_color('black')
    plt.gca().spines['right'].set_color('black')

    genotypes = load_individual_genotypes(reference_vntrs, False)
    genotyped_vntr_ids = {_id: set() for _id in range(1000000)}

    reference_repeats = {_id: [] for _id in range(1000000)}
    hist_data = []
    cp_diff = []
    for individual_id, result_map in genotypes.items():
        affected = 0
        different_sites = 0
        for genotyped_vntr, genotype in result_map.items():
            if genotype == 'None' or genotype is None or None in genotype or 'None' in genotype:
                continue
            ref_vntr_length = reference_vntrs[genotyped_vntr].get_length()
            vntr_pattern_length = len(reference_vntrs[genotyped_vntr].pattern)
            ref_copy_number = len(reference_vntrs[genotyped_vntr].get_repeat_segments())
            affected += sum([abs(ref_vntr_length - allele * vntr_pattern_length) for allele in genotype if allele != ref_copy_number])
            different_sites += sum([1 for allele in genotype if allele != ref_copy_number] + [0])
            reference_repeats[genotyped_vntr] += genotype
        hist_data.append(affected)
        cp_diff.append(different_sites)

    # data = {point: 0 for point in range(min(hist_data), max(hist_data) + 1)}
    # for h in hist_data:
    #     data[h] += 1
    # del data[0]
    #
    # plt.bar(data.keys(), data.values(), width=0.8)
    print('different sites')
    print(sum(cp_diff) / float(len(cp_diff)))
    print(min(cp_diff))
    print(cp_diff)
    hist_mean = sum(hist_data) / float(len(hist_data))
    print(hist_mean)
    plt.hist(hist_data, 100)
    plt.axvline(hist_mean, color='k', linestyle='dashed', linewidth=1)
    #    plt.xlabel('Base-pairs difference from reference')
    plt.xlabel('Total BP affected by VNTR variation')
    plt.ylabel('Number of individuals')
    plt.xscale('log')
    ax = plt.gca()
    ax.xaxis.grid(False)
    plt.savefig('affected_bp_per_individual.pdf')
    plt.cla()
    plt.clf()
    plt.close()


def plot_most_common_allele_difference_from_reference():
    from brokenaxes import brokenaxes

    genotypes = load_individual_genotypes(reference_vntrs, False)
    genotyped_vntrs = set([])

    most_common_allele = {_id: [] for _id in range(1000000)}
    for individual_id, result_map in genotypes.items():
        for genotyped_vntr, alleles in result_map.items():
            if alleles == 'None' or alleles is None or alleles[0] is None:
                continue
            genotyped_vntrs.add(genotyped_vntr)
            most_common_allele[genotyped_vntr] += alleles
    for key, value in most_common_allele.items():
        if most_common_allele[key] == []:
            most_common_allele[key].append(0)
    most_common_allele = {key: most_common(value) for key, value in most_common_allele.items() if key in genotyped_vntrs}

    hist_data = []
    for vntr_id, major_allele in most_common_allele.items():
        hist_data.append(major_allele - len(reference_vntrs[vntr_id].get_repeat_segments()))

    bar_ranges = 6
    print('avg:', sum([abs(e) for e in hist_data if abs(e) > 0]) / float(len([abs(e) for e in hist_data if abs(e) > 0])))
    print('avg:', sum([abs(e) for e in hist_data]) / float(len(hist_data)))
    print(len([1 for e in hist_data if abs(e) > 3]), len([1 for e in hist_data if abs(e) > 0]), len(hist_data))
    data = {point: 0 for point in range(-bar_ranges, bar_ranges + 1)}
    for h in hist_data:
        h = min(bar_ranges, h)
        h = max(h, -bar_ranges)
        data[h] += 1
    # del data[0]

    print(data)
    second_peak = max([data[i] for i in data.keys() if i !=0 ])
    bax = brokenaxes(ylims=((-1, second_peak + second_peak * 0.4), (data[0]-second_peak*0.3, data[0]+second_peak*0.2)), hspace=.09)

    bax.bar(data.keys(), data.values(), width=0.8)
#    plt.hist(hist_data, 40)
    bax.set_ylabel('Number of VNTR Loci', labelpad=40)

    # bax.canvas.draw()
    # print(bax.get_xticklabels())
    # print(bax.get_xticklabels()[0])

    # labels = [item.get_text() for item in bax.get_xticklabels()[1]]
    # labels[0] = '<=6'
    # labels[-1] = '>=6'
    bax.set_xticklabels(['', u'<=6', u'-4', u'-2', u'0', u'2', u'4', u'>=6', ''])

    bax.set_xlabel('Most common allele repeat difference from GRCh38', labelpad=20)
    plt.savefig('common_allele_difference_from_reference.pdf')
    plt.cla()
    plt.clf()
    plt.close()


def plot_bp_difference_from_reference():
    from brokenaxes import brokenaxes

    genotypes = load_individual_genotypes(reference_vntrs, False)

    hist_data = []
    for individual_id, result_map in genotypes.items():
        for genotyped_vntr, alleles in result_map.items():
            if alleles == 'None' or alleles is None or alleles[0] is None:
                continue
            ref_repeats = len(reference_vntrs[genotyped_vntr].get_repeat_segments())
            pattern_len = len(reference_vntrs[genotyped_vntr].pattern)
            for allele in alleles:
                # if allele != ref_repeats:
                hist_data.append(int(round(allele * pattern_len - ref_repeats * pattern_len))) # bp difference
                # hist_data.append(int(round(allele - ref_repeats))) # copy number difference

    bar_ranges = 10
    print('avg:', sum([abs(e) for e in hist_data if abs(e) > 0]) / float(len([abs(e) for e in hist_data if abs(e) > 0])))
    print('avg:', sum([abs(e) for e in hist_data]) / float(len(hist_data)))
    print(len([1 for e in hist_data if abs(e) > 3]), len([1 for e in hist_data if abs(e) > 0]))
    data = {point: 0 for point in range(-bar_ranges, bar_ranges + 1)}
    for h in hist_data:
        h /= 10
        h = min(bar_ranges, h)
        h = max(h, -bar_ranges)
        data[h] += 0.001
    # del data[0]

    print(data)
    second_peak = max([data[i] for i in data.keys() if i !=0 ])
    bax = brokenaxes(ylims=((-1, second_peak + second_peak * 0.4), (data[0]-second_peak*0.3, data[0]+second_peak*0.2)), hspace=.09)

    bax.bar(data.keys(), data.values(), width=0.8)
#    plt.hist(hist_data, 40)
    bax.set_ylabel('Number of VNTRs (x 1000)', labelpad=40)

    # bax.canvas.draw()
    # print(bax.get_xticklabels())
    # print(bax.get_xticklabels()[0])
    labels = [item.get_text() for item in bax.get_xticklabels()[1]]
    labels[0] = '<=6'
    labels[-1] = '>=6'
    bax.set_xticklabels(['<=60', u'<=100', u'-50', u'0', '50', u'>=100', '>=6'])

    bax.set_xlabel('Base-pairs difference from reference', labelpad=20)
    plt.savefig('bp_difference_from_reference_all_genotypes.pdf')
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

    import scipy
    corrMatrix = np.full((len(tissue_labels), len(tissue_labels)), np.nan)
    TISSUES = list(tissue_labels)
    all_vntrs = list(eVNTRs)
    for tissue1 in TISSUES:
        for tissue2 in TISSUES:
            t1ind = TISSUES.index(tissue1)
            t2ind = TISSUES.index(tissue2)
            e1 = set(eVNTRs_of_tissue[tissue1].keys())
            e2 = set(eVNTRs_of_tissue[tissue2].keys())
            eff = float(len(e1 & e2)) / len(e1 | e2)
            eff1 = [1 if e in e1 else 0 for e in all_vntrs]
            eff2 = [1 if e in e2 else 0 for e in all_vntrs]
            eff1 = [eVNTRs_of_tissue[tissue1][e] if e in e1 else 0 for e in all_vntrs]
            eff2 = [eVNTRs_of_tissue[tissue2][e] if e in e2 else 0 for e in all_vntrs]
            eff = scipy.stats.spearmanr(eff1, eff2)[0]
#            eff = len(e1 & e2) / 163.0
            corrMatrix[t1ind, t2ind] = eff
            corrMatrix[t2ind, t1ind] = eff
    TISSUES = [t.replace('-', ' ').replace('   ', ' - ') for t in TISSUES]
    for i in range(len(TISSUES)):
        TISSUES[i] = TISSUES[i].replace('Gastroesophageal Junction', 'GE Junction')
        p = TISSUES[i].find('(')
        if p != -1:
            TISSUES[i] = TISSUES[i][:p].strip()

    columns = [e.replace('Cells - ', '') for e in TISSUES]
    corrMatrix = pd.DataFrame(corrMatrix, columns=columns, index=TISSUES)
    cg = sns.clustermap(corrMatrix, cmap="Blues", yticklabels=1, xticklabels=2)#mako, PuBu
    cg.ax_col_dendrogram.set_visible(False)
    cg.savefig("tissue_correlation2.pdf")#, dpi=300)
    return

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
    ax.set_xticklabels(tissue_labels, fontsize=5)
    ax.set_yticklabels(tissue_labels, fontsize=5)

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
    genes_info = {}
    with open('/home/mehrdad/workspace/adVNTR/results/hg38_annotation/refseq_genes.bed') as infile:
        genes_lines = infile.readlines()
        for line in genes_lines:
            line = line.strip().split()[:6]
            chromosome, start, end, identifier, _, direction = line
            start = int(start)
            end = int(end)
            if chromosome not in genes_info.keys():
                genes_info[chromosome] = []
            genes_info[chromosome].append((start, end, identifier, direction))
    genes_coordinates = {}
    for chromosome, coordinates in genes_info.items():
        genes_coordinates[chromosome] = sorted(coordinates)

    tss_distances = {}
    computed = set()

    def get_tss_distance(ref_vntr, coordinates):
        if ref_vntr.id in computed:
            return tss_distances[ref_vntr.id]
        res = 1e9
        for start, end, identifier, direction in coordinates:
            if direction == '+':
                TSS = start
            else:
                TSS = end
            new_res = abs(ref_vntr.start_point - TSS)
            if ref_vntr.start_point < TSS < ref_vntr.start_point + ref_vntr.get_length():
                new_res = 0
            new_res2 = abs(ref_vntr.start_point + ref_vntr.get_length() - end)
            new_res = min(new_res, new_res2)
            if new_res < res:
                res = new_res
            # elif res != 1e9 and new_res > 20 * res:
            #     break
        computed.add(ref_vntr.id)
        tss_distances[ref_vntr.id] = res
        return res

    thresholds = {}
    with open('thresholds.txt') as infile:
        lines = infile.readlines()
        for l in lines:
            thresholds[l.split('\t')[0].replace(' ', '-')] = float(l.split('\t')[1])

    het_test = {}
    with open('heterozygosity_test.txt') as infile:
        lines = infile.readlines()
        for l in lines:
            vid, pval = l.split()
            vid = int(vid)
            pval = float(pval)
            het_test[vid] = pval

    data = {}
    all_data = {}
    hist_data = []
    tissue_data_dir = glob.glob('regression_results/*')
    both_annotation_vntrs = set()
    tss_significant_pvalue = set()
    close_vntrs = set()
    print(tissue_data_dir)
    other_vntrs = set()
    other_sig_vntrs = set()
    for tissue_dir in tissue_data_dir:
#        print(tissue_dir)
        tissue_name = os.path.basename(tissue_dir)
        vntr_files = glob.glob(tissue_dir + '/*')
        for vntr_file in vntr_files:
            vntr_id = int(os.path.basename(vntr_file[:-4]))
            if het_test[vntr_id] < 0.05:
                continue
#            if reference_vntrs[vntr_id].annotation not in ['UTR']:
#                continue
            with open(vntr_file) as infile:
                lines = infile.readlines()
            effect_size = float(lines[0].strip().split('\t')[3])
            pvalue = float(lines[0].strip().split('\t')[4])
            dist_to_tss = get_tss_distance(reference_vntrs[vntr_id],
                                           genes_coordinates[reference_vntrs[vntr_id].chromosome])
            if dist_to_tss < 1:
                close_vntrs.add(vntr_id)
            if dist_to_tss < 50:
                both_annotation_vntrs.add(vntr_id)
                if pvalue < thresholds[tissue_name]:
                    tss_significant_pvalue.add(vntr_id)
            else:
                other_vntrs.add(vntr_id)
                if pvalue < thresholds[tissue_name]:
                    other_sig_vntrs.add(vntr_id)
            denom = 100
            if dist_to_tss / denom not in all_data.keys():
                all_data[dist_to_tss / denom] = 0
            all_data[dist_to_tss / denom] += 1
            if pvalue < 0.0005:
                if dist_to_tss < 1:
                    print(vntr_id, tissue_name, pvalue, reference_vntrs[vntr_id].gene_name)
                if dist_to_tss > 12000:
                    print(vntr_id, reference_vntrs[vntr_id].gene_name, reference_vntrs[vntr_id].annotation)
                if dist_to_tss/denom not in data.keys():
                    data[dist_to_tss/denom] = 0
                data[dist_to_tss/denom] += 1
                hist_data.append(effect_size)

    print(data)
    print('close vntrs:', len(close_vntrs))
    print('both annotation:', len(both_annotation_vntrs))
    print('significant TSS vntrs', len(tss_significant_pvalue))
    print('other vntrs and other sig vntrs', len(other_vntrs), len(other_sig_vntrs))
    ratio = {}
    for key in all_data.keys():
        sig = data[key] if key in data.keys() else 0
        ratio[key] = sig / float(all_data[key])
    plt.bar(ratio.keys(), ratio.values())
    # plt.bar(data.keys(), data.values())
    # plt.bar(all_data.keys(), all_data.values())
    plt.xlabel('Distance to transcription start site') # or gene start, first exon, or ...
    plt.ylabel('Number of VNTR-eQTLs') # or p-value, effect size, ...
    plt.savefig('distance_from_transcription_site.pdf')
    plt.cla()
    plt.clf()
    plt.close()


def plot_significant_vntrs_and_tissues():
    # show which tissues are affected for a list of 20-30 VNTRs (end of page11 in Melissa's paper)
    pass


def plot_eqtl_significance_per_annotation():
    import matplotlib.pyplot as plt

    # plt.style.use('ggplot')
    # plt.rcParams['axes.facecolor'] = '#FFFFFF'
    # plt.gca().spines['bottom'].set_color('black')
    # plt.gca().spines['left'].set_color('black')
    # plt.gca().spines['top'].set_color('black')
    # plt.gca().spines['right'].set_color('black')

    import seaborn as sns
    sns.set()
    sns.set_style("whitegrid")
    sns.set_context("paper")

    data = {}
    tissue_data_dir = glob.glob('regression_results/*')
    for tissue_dir in tissue_data_dir:
        vntrs_data = glob.glob(tissue_dir + '/*')
        tissue_name = os.path.basename(tissue_dir).split('---')[0]
        if len(os.path.basename(tissue_dir).split('---')) > 1:
            tissue_name += ' - ' + os.path.basename(tissue_dir).split('---')[1][:5]
        for vntr_pvalue_dir in vntrs_data:
            vntr_id = int(os.path.basename(vntr_pvalue_dir)[:-4])
            pvalue_file = vntr_pvalue_dir
            with open(pvalue_file) as infile:
                line = infile.readlines()[0].strip()
                pvalue = float(line.split()[4])
                if pvalue > 0.0005:
                    continue
                pvalue = max(pvalue, 1e-45)
                if vntr_id not in data.keys():
                    data[vntr_id] = []
                data[vntr_id].append(pvalue)

    genes_info = {}
    with open('/home/mehrdad/workspace/adVNTR/results/hg38_annotation/refseq_genes.bed') as infile:
        genes_lines = infile.readlines()
        for line in genes_lines:
            line = line.strip().split()[:6]
            chromosome, start, end, identifier, _, direction = line
            start = int(start)
            end = int(end)
            if chromosome not in genes_info.keys():
                genes_info[chromosome] = []
            genes_info[chromosome].append((start, end, identifier, direction))
    genes_coordinates = {}
    for chromosome, coordinates in genes_info.items():
        genes_coordinates[chromosome] = sorted(coordinates)

    tss_distances = {}
    computed = set()

    def get_tss_distance(ref_vntr, coordinates):
        if ref_vntr.id in computed:
            return tss_distances[ref_vntr.id]
        res = 1e9
        for start, end, identifier, direction in coordinates:
            if direction == '+':
                TSS = start
            else:
                TSS = end
            new_res = abs(ref_vntr.start_point - TSS)
            if ref_vntr.start_point < TSS < ref_vntr.start_point + ref_vntr.get_length():
                new_res = 0
            new_res2 = abs(ref_vntr.start_point + ref_vntr.get_length() - end)
            new_res = min(new_res, new_res2)
            if new_res < res:
                res = new_res
            # elif res != 1e9 and new_res > 20 * res:
            #     break
        computed.add(ref_vntr.id)
        tss_distances[ref_vntr.id] = res
        return res



    hist_data = {}
    for vntr_id in data.keys():
        annotation = reference_vntrs[vntr_id].annotation
        if get_tss_distance(reference_vntrs[vntr_id], genes_coordinates[reference_vntrs[vntr_id].chromosome]) < 10:
            annotation = 'Promoter'
        if annotation not in hist_data.keys():
            hist_data[annotation] = []
        hist_data[annotation].append(min(data[vntr_id]))

    for index, annotation in enumerate(hist_data.keys()):
        bins = np.logspace(np.log10(1e-50),np.log10(1.0), 20)
        values, base = np.histogram(hist_data[annotation], bins=bins)
        i = index * index * 2
        print(i, annotation)
        print(sum([1 if el < (1e-11) else 0 for el in hist_data[annotation]]) / float(len(hist_data[annotation])))
        # plt.bar(bins[:-1]+i*bins[:-1], values, width=bins[:-1]+i*bins[:-1], label=annotation)

        #option2:
        # plt.hist(hist_data[annotation], bins=np.logspace(np.log10(1e-50),np.log10(1.0), 200), label=annotation, histtype='step', cumulative=True)
        pass
    #option3
    red_color = '#d62728'
    blue_color = "#1f77b4"
    green_color = '#2ca02c'
    plt.hist(hist_data.values(), bins=np.logspace(np.log10(1e-50),np.log10(1.0), 200), label=hist_data.keys(), cumulative=True, histtype='step', normed=True, color=[red_color, blue_color, green_color], linewidth=1.2)

    plt.xscale('log')
    plt.xlabel('P-value', fontsize=12)

    plt.draw()
    labels = [item.get_text() for item in plt.gca().get_yticklabels()]
    plt.gca().set_yticklabels([str(float(e)*100) + '%' if e != '' else e for e in labels])

    labels = [item.get_text() for item in plt.gca().get_xticklabels()]
    print(labels)
    plt.gca().set_xticklabels([e.replace('10^{-46}', '<10^{-46}') for e in labels])

    plt.gca().tick_params(labelsize=10)

    plt.ylabel('Percent of VNTR', fontsize=12)
    plt.legend(loc=2, fontsize=10)
    plt.tight_layout()
    plt.savefig('significance_for_annotation.pdf')
    plt.savefig('significance_for_annotation.pdf')
    plt.cla()
    plt.clf()
    plt.close()


def plot_causality_ranks():
    pvalue_ranks = '''1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n2\n2\n2\n2\n2\n2\n2\n3\n3\n3\n3\n3\n3\n4\n4\n4\n4\n5\n5\n6\n6\n7\n7\n8\n8\n8\n9\n9\n10\n10\n11\n12\n12\n13\n15\n15\n17\n19\n23\n23\n24\n25\n25\n25\n27\n29\n29\n30\n31\n31\n32\n32\n33\n37\n41\n44\n46\n47\n47\n47\n54\n55\n58\n61\n61\n90\n90\n93\n97\n100\n111\n116\n118\n134\n136\n141\n142\n151\n151\n190\n346\n994'''
    caviar_ranks = '''1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1\n2\n2\n2\n3\n2\n2\n2\n3\n3\n3\n3\n3\n3\n4\n4\n4\n4\n5\n5\n6\n6\n7\n7\n8\n8\n8\n9\n9\n10\n10\n11\n12\n12\n14\n15\n15\n17\n19\n23\n23\n24\n25\n25\n25\n27\n29\n31\n30\n31\n31\n32\n32\n33\n37\n41\n44\n46\n47\n47\n47\n54\n55\n58\n61\n61\n90\n90\n93\n97\n100\n111\n116\n118\n134\n136\n141\n142\n151\n151\n190\n346\n994'''
    harmomnic_ranks = ''''''
    pvalue_ranks = [float(e) for e in pvalue_ranks.split('\n') if e != '']
    caviar_ranks = [float(e) for e in caviar_ranks.split('\n') if e != '']
    harmomnic_ranks = [(1.0/(0.5*(1.0/pvalue_ranks[i]+1.0/caviar_ranks[i])), caviar_ranks[i], pvalue_ranks[i]) for i in range(len(caviar_ranks))]
    harmomnic_ranks = sorted(harmomnic_ranks)
    pvalue_ranks = [e[2] for e in harmomnic_ranks]
    caviar_ranks = [e[1] for e in harmomnic_ranks]
    harmomnic_ranks = [e[0] for e in harmomnic_ranks]

    discrepancies = [2*abs(pvalue_ranks[i] - caviar_ranks[i]) / float((pvalue_ranks[i] + caviar_ranks[i])) for i in range(len(pvalue_ranks))]
    print('mean dicrepancy:', sum(discrepancies) / len(discrepancies))

    print('T <= 10:', sum([1 for e in harmomnic_ranks if e <= 10]))
    print('harmonic <= 1:', sum([1 for e in harmomnic_ranks if e <= 1]))

    plt.xlabel('VNTR Number (ordered by rank)', fontsize=12)
    plt.ylabel('Rank', fontsize=12)
    plt.yscale('log')
    plt.plot([i + 1 for i in range(len(pvalue_ranks))], pvalue_ranks, label='P-value Ranks', linewidth=1.2)
    plt.plot([i + 1 for i in range(len(caviar_ranks))], caviar_ranks, label='CAVIAR Ranks', linewidth=1.2)
    plt.plot([i + 1 for i in range(len(harmomnic_ranks))], harmomnic_ranks, label='Harmonic Ranks', linewidth=1.2)
    plt.legend(fontsize=10)
    plt.draw()
    plt.gca().tick_params(labelsize=10)
    plt.savefig('VNTR_causality_ranks.png')
    plt.cla()
    plt.clf()
    plt.close()


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

    # plot_tissue_correlations()

    # plot_affected_bp_per_individual()

    # plot_eqtl_significance_per_annotation()
    # plot_bp_difference_from_reference()
    # plot_most_common_allele_difference_from_reference()

    # plot_dose_dependence_of_correlation()
    # plot_genotypes_difference_from_reference()
    # plot_vntrs_difference_from_reference()
    # plot_vntr_pvalues_qq_plot()
    # plot_evntrs_per_tissue()

    # plot_effect_per_allele_frequency()
    # plot_evntrs_and_number_of_tissues()

    peer_files = glob.glob('PEER_results/*')
    tissue_names = [pf.split('/')[-1][13:-3] for pf in peer_files]
    # for tissue in tissue_names:
    #     plot_expression_genotype_correlation(423956, tissue)

    #plot_expression_genotype_correlation(450457, 'Pancreas')
    #plot_pvalue_and_caviar(450457, 'Pancreas')

    # plot_expression_genotype_correlation(450457, 'Colon - Transverse')
    # plot_pvalue_and_caviar(450457, 'Colon - Transverse')

    #plot_expression_genotype_correlation(331737, 'Brain - Hippocampus')
    #plot_pvalue_and_caviar(331737, 'Brain - Hippocampus')

    # plot_expression_genotype_correlation(690585, 'Nerve - Tibial')
    # plot_pvalue_and_caviar(690585, 'Nerve - Tibial')

    # plot_allele_count_distribution()
    # plot_vntr_polymorphic_rate_based_on_annotation()
    # plot_vntr_polymorphic_rate_based_on_cohort_size()
    plot_causality_ranks()
#    plot_vntr_location_effect()
