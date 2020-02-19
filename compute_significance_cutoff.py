
import glob
import numpy as np
from bisect import bisect


FDR_RATE = 0.05

def load_pvalues_in_directory(dir, file_name_suffix, tissue_name):
    res = []
    tissues_dir = glob.glob(dir + '/*')
    for tissue_dir in tissues_dir:
        if tissue_dir.split('/')[-1] != tissue_name:
            continue
        vntr_dirs = glob.glob(tissue_dir + '/*')
        for vntr_dir in vntr_dirs:
            file_name = vntr_dir + '/' + file_name_suffix
            with open(file_name) as input:
                pvalue = float(input.readlines()[0].strip())
                res.append(pvalue)
    return res


def plot_qq_plot(pvalues, permuted_pvalues_list, significance_cutoff, tissue_name=None, filling='full'):
    import matplotlib
    matplotlib.use('Agg')
    import seaborn as sns
    sns.set()
    import matplotlib.pyplot as plt
    sns.set_style("whitegrid")
    sns.set_context("paper")
    sns.set_palette(sns.color_palette("deep", 40))
    sns.set_palette(sns.color_palette("Paired", 46))


    from pandas import DataFrame
    permuted_x = []
    permuted_y = []
    permutation_averages = None
    for i, permuted_pvalues in enumerate(permuted_pvalues_list):
        label = 'Permuted Genotypes' if i == 0 else None
        N = float(len(permuted_pvalues))
        # print(min(permuted_pvalues), max(permuted_pvalues))
        grid = -np.log10(np.arange(1, 1 + N) / N)
        permuted_x += list(grid)
        permuted_y += list(-np.log10(np.array(sorted(permuted_pvalues, reverse=False))))
        if permutation_averages is None:
            permutation_averages = np.array(-np.log10(np.array(sorted(permuted_pvalues, reverse=False))))
        else:
            permutation_averages += np.array(-np.log10(np.array(sorted(permuted_pvalues, reverse=False))))
        plt.scatter(grid, -np.log10(np.array(sorted(permuted_pvalues, reverse=False))), color='gray', marker='.', alpha=0.15, zorder=2)
    permutation_averages /= len(permuted_pvalues_list)

    df = DataFrame(permuted_y, permuted_x, columns=['observed_pvalues'])
    df['expected_pvalues'] = df.index

    permuted_label = 'Permuted Genotypes' if tissue_name == 'Minor Salivary Gland' else None
    # ax = sns.lineplot(x="expected_pvalues", y="observed_pvalues", data=df, n_boot=30, sort=True, color='gray')
    # plt.plot(grid, permutation_averages, label=permuted_label, color='gray', zorder=1)
    limit = grid[[0]]
    plt.plot([0, limit], [0, limit], zorder=1, color='gray')

    N = float(len(pvalues))
    print(min(pvalues), max(pvalues))
    grid = -np.log10(np.arange(1, 1 + N) / N)
    identified_label = tissue_name #'identified genotypes'
    marker_style = dict(marker='.', fillstyles='left')
    # marker = matplotlib.markers.MarkerStyle(marker='.', fillstyle='right')
    # plt.scatter(grid, -np.log10(np.array(sorted(pvalues, reverse=False))), label=identified_label, zorder=3, marker=marker, color='tab:red')#,marker='.')
    marker = matplotlib.markers.MarkerStyle(marker='.', fillstyle=filling)
    plt.scatter(grid, -np.log10(np.array(sorted(pvalues, reverse=False))), label=identified_label, zorder=3, marker=marker)
    plt.scatter(grid, permutation_averages, label=permuted_label, marker='.', zorder=1, color='gray')

    # ax = sns.scatterplot(grid, -np.log10(np.array(sorted(pvalues, reverse=False))), label='identified genotypes', marker='.')

    plt.legend(ncol=2, fontsize=4)
    plt.xlabel('Expected p-value [-log10(p)]')
    plt.ylabel('Observed p-value [-log10(p)]')
    plt.savefig('qq_plots_aggregate/qq_plot_p_values_%s.png' % tissue_name.replace(' ', '-'), dpi=300)
    # plt.cla()
    # plt.clf()
    # plt.close()


def find_cutoff(tissue_name, tissue_number=None, geuvadis=False):
    geuvadis_dir = 'geuvadis_' if geuvadis else ''
    pvalues = load_pvalues_in_directory('%sall_vntr_pvalues/' % geuvadis_dir, 'pvalues.txt', tissue_name)
    # print('loaded p-values')
    permuted_pvalues = []
    for i in range(0, 50):
        permuted_pvalues.append(load_pvalues_in_directory('/mnt/%spermutated_pvalues_multiple_tests/permutated_pvalues_%s/' % (geuvadis_dir, i), 'permutated_pvalues.txt', tissue_name))
    # print('loaded permuted p-values')
    # print(np.median(np.array(pvalues)), min(pvalues))
    # print(np.median(np.array(permuted_pvalues)), min(permuted_pvalues))

    sorted_pvalues = sorted(pvalues)
    sorted_permuted = sorted(permuted_pvalues[0])
    print('all tests:', len(sorted_pvalues))

    significance_cutoff = 0
    for theta in sorted(pvalues + permuted_pvalues[0], reverse=True):
        fdr = float(bisect(sorted_permuted, theta)) / bisect(sorted_pvalues, theta)
        if fdr <= FDR_RATE:
            significance_cutoff = theta
            print('found fdr', fdr, theta)
            # print('n1 ', bisect(sorted_permuted, theta))
            # print('n2 ', bisect(sorted_pvalues, theta))
            break

    # Benjamini-Hochberg
    res = 0
    for rank, theta in enumerate(sorted_pvalues):
        rank += 1
        im_q = float(rank) / len(pvalues) * FDR_RATE
        if theta < im_q:
            res = theta
    print('benjamini: ', res)
    print('benjamini significants ', bisect(sorted_pvalues, res))
    if geuvadis:
        return res

    from math import log
    fillings = ['left', 'right', 'top', 'full']
    filling = fillings[tissue_number / 12] if tissue_number is not None else 'full'
    filling = 'full'
    plot_qq_plot(pvalues, permuted_pvalues, -log(res, 10), tissue_name, filling=filling)
    return res

tissue_dirs = glob.glob('all_vntr_pvalues/*')
tissue_names = [e.split('/')[-1] for e in tissue_dirs]

find_cutoff('Whole Blood', 0, geuvadis=True)

for i, tname in enumerate(tissue_names):
    print(tname)
    cutoff = find_cutoff(tname, i)
    with open('thresholds.txt', 'a') as outfile:
        outfile.write('%s\t%s\n' % (tname, cutoff))
