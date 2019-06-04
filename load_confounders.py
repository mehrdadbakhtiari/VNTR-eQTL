

def load_peer_factors(tissue_name='Blood Vessel'):
    K = 15
    peer_factor_file = 'PEER_results/peer_factors_%s_%s' % (tissue_name, K)
    with open(peer_factor_file) as infile:
        lines = infile.readlines()

    peer_factors = {}
    factors = []
    for j, line in enumerate(lines):
        line = line.strip().split('\t')
        peer_factors[line[0]] = {}
        for i in range(1, len(line)):
            peer_factor = 'peer_%s' % (i-1)
            peer_factors[line[0]][peer_factor] = float(line[i])
            if j == 0:
                factors.append(peer_factor)

    return peer_factors, factors

def load_population_pcs():
    population_pc_file = 'PCA_results/pca.gtex.pca.evec'
    with open(population_pc_file) as infile:
        lines = infile.readlines()

    pop_structure_factors = {}
    pcs = []
    for j in range(len(lines)):
        if j == 0 or len(lines[j].strip()) == 0:
            continue
        line = lines[j].strip().split()
        pop_structure_factors[line[0]] = {}
        for i in range(1, len(line)-1):
            pc = 'pop_pc_%s' % (i-1)
            pop_structure_factors[line[0]][pc] = float(line[i])
            if j == 1:
                pcs.append(pc)

    return pop_structure_factors, pcs

