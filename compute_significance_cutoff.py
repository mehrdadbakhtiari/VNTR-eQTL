
import glob
import numpy as np
from bisect import bisect


FDR_RATE = 0.05

def load_pvalues_in_directory(dir, file_name_suffix):
    res = []
    tissues_dir = glob.glob(dir + '/*')
    for tissue_dir in tissues_dir:
        vntr_dirs = glob.glob(tissue_dir + '/*')
        for vntr_dir in vntr_dirs:
            file_name = vntr_dir + '/' + file_name_suffix
            with open(file_name) as input:
                pvalue = float(input.readlines()[0].strip())
                res.append(pvalue)
    return res


def plot_qq_plot(pvalues, permuted_pvalues, significance_cutoff):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import seaborn as sns
    sns.set()
    sns.set_style("whitegrid")
    sns.set_context("paper")

    N = float(len(pvalues))
    print(min(pvalues), max(pvalues))
    grid = -np.log10(np.arange(1, 1 + N) / N)
    plt.scatter(grid, -np.log10(np.array(sorted(pvalues, reverse=False))), label='identified genotypes')

    N = float(len(permuted_pvalues))
    print(min(permuted_pvalues), max(permuted_pvalues))
    grid = -np.log10(np.arange(1, 1 + N) / N)
    plt.scatter(grid, -np.log10(np.array(sorted(permuted_pvalues, reverse=False))), label='permuted genotypes')

    plt.plot([0, 5], [0, 5], label='x = y', color='red')
    # plt.plot([0, 5], [significance_cutoff, significance_cutoff], label='Cutoff', color='black')

    #    plt.ylim((0, 20))
    #    plt.xlim((0, 10))
    plt.legend()
    plt.xlabel('Expected p-value [-log10(p)]')
    plt.ylabel('Observed p-value [-log10(p)]')
    plt.savefig('qq_plot_p_values.png', dpi=300)
    plt.cla()
    plt.clf()
    plt.close()


def find_cutoff():
    pvalues = load_pvalues_in_directory('all_vntr_pvalues/', 'pvalues.txt')
    permuted_pvalues = load_pvalues_in_directory('permutated_pvalues/', 'permutated_pvalues.txt')
    print(np.median(np.array(pvalues)), min(pvalues))
    print(np.median(np.array(permuted_pvalues)), min(permuted_pvalues))

    sorted_pvalues = sorted(pvalues)
    sorted_permuted = sorted(permuted_pvalues)
    print('all tests:', len(sorted_pvalues))

    significance_cutoff = 0
    for theta in sorted(pvalues + permuted_pvalues, reverse=True):
        fdr = float(bisect(sorted_permuted, theta)) / bisect(sorted_pvalues, theta)
        if fdr <= FDR_RATE:
            significance_cutoff = theta
            print('found fdr', fdr, theta)
            print('n1 ', bisect(sorted_permuted, theta))
            print('n2 ', bisect(sorted_pvalues, theta))
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

    from math import log
    plot_qq_plot(pvalues, permuted_pvalues, -log(significance_cutoff, 10))


find_cutoff()
