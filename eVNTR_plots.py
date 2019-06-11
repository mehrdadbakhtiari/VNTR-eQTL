import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn.apionly as sns
import numpy as np
import pandas as pd
from math import log10

from advntr.models import load_unique_vntrs_data

vntr_models_dir = '/home/mehrdad/workspace/adVNTR/vntr_data/hg38_selected_VNTRs_Illumina.db'
ref_vntrs = load_unique_vntrs_data(vntr_models_dir)
reference_vntrs = {}
for ref_vntr in ref_vntrs:
        reference_vntrs[ref_vntr.id] = ref_vntr


def plot_expression_genotype_correlation():
    gene_name = 'AS3MT'
    tissue_name = 'Blood Vessel'
    data = [[2.03838086128, 2.7370822429699997, 2.05135184526, 1.5914614200600001, 0.958451598883, 1.1516859531399999, 1.3514041900600002, 1.74997264147, 2.07785964012, 0.278437435627],
    [2.7890409231200004, 2.66182696819, 2.8365943431900003, 1.85346615314, 3.58080442746, 4.039058208469999, 2.99621518453, 2.16904437542, 1.97231662273, 1.4924068749000001, 2.8461830616, 3.74974346161, 4.124842882159999, 3.01081669331, 6.583733081819999, 2.89373207092, 2.00292903185, 1.64222568274, 1.8419418335, 3.3163380622900003, 2.60120218992, 1.2146719694099999, 1.8719456791900002, 2.87023258209, 4.13667408625, 4.193371772769999, 2.12403798103, 3.3078730106400003, 1.65636360645, 4.9908556938199995, 2.46760606766, 1.40581941605, 1.6427392959600002, 2.5626500845, 2.07668371995, 3.36892787615],
    [3.32425562541, 3.5104527473400005, 3.08298095067, 4.673421621319999, 4.17712974548, 3.19682518641, 4.23930692673, 3.14745998383, 3.64102172852, 2.1144301891299997, 4.476975917819999, 4.137458801269999, 3.69442009926, 3.72947589556, 3.40179014206, 4.59486802419, 5.2096545696300005, 4.51752734184, 2.93920099735, 4.5821495056199995, 4.823669274649999, 3.88390767574, 1.60192894936, 3.4485193491, 3.98949337006, 1.9181122779799997, 2.9278344710699997],]
    labels = ['2', '2.5', '3']

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
    fig.savefig('expression_%s.pdf' % gene_name)
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


def plot_variant_pvalues():
    fig = plt.figure(figsize=(6, 2))
    ax = fig.add_subplot(111)

    gene_name = 'AS3MT'

    data = [0.0001, 0.000005, 0.001, 0.01, 0.001, 0.0001, 0.001, 0.0001]
    data = [-log10(e) for e in data]
    x = [i for i in range(len(data))]
    ax.scatter(x, data, marker='.', c='gray')
    vntr_x = 4.5
    vntr_p_value = -log10(0.000000000001)
    ax.scatter(vntr_x, vntr_p_value, marker='*', c='r')
    # ax.set_xlabel('Genomic location')
    ax.set_ylabel('-log10 P-value')
    fig.savefig('pvalues_%s.pdf' % gene_name)
    plt.cla()
    plt.clf()
    plt.close()


if __name__ == '__main__':
    # plot_expression_genotype_correlation()
    # plot_variant_pvalues()
    plot_variant_caviar_scores(111235, 'Esophagus')
    plot_variant_caviar_scores(331737, 'Esophagus')

