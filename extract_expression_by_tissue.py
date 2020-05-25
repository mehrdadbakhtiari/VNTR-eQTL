rna_table_file = '../Sra_table_RNA-Seq_only'
rpkm_file = '../files/phe000020.v1.GTEx_RNAseq.expression-data-matrixfmt.c1/GTEx_Data_20160115_v7_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct'
expression_by_tissue_dir = '../Expression_by_Subtissue/'

import sys
if (sys.argv) >= 4:
    rna_table_file = sys.argv[1]
    rpkm_file = sys.argv[2]
    expression_by_tissue_dir = sys.argv[3]

protein_id_to_gene = {}


def get_average(lst):
    return sum(lst) / len(lst)


def write_expression_for_tissue(tissue_name, data):
    output_file = expression_by_tissue_dir + tissue_name + '.rpkm'
    rows = len(data.keys())
    columns = len(data[data.keys()[0]][tissue_name].keys())
    dirname = os.path.dirname(output_file)
    if not os.path.exists(dirname):
        os.makedirs(dirname)
    with open(output_file, 'w') as outfile:
        outfile.write('%s\t%s\n' % (rows, columns))
        outfile.write('Name\tDescription\t')
        outfile.write('\t'.join(data[data.keys()[0]][tissue_name].keys()))
        outfile.write('\n')
        for protein in data.keys():
            outfile.write('%s\t%s' % (protein, protein_id_to_gene[protein]))
            for individual in data[protein][tissue_name].keys():
                expression = get_average(data[protein][tissue_name][individual])
                outfile.write('\t%s' % expression)
            outfile.write('\n')


experiment_tissue = {}
main_tissue = {}


with open(rna_table_file) as infile:
    lines = infile.readlines()
for line in lines:
    experiment_id = line.split('\t')[24]
    exact_tissue = line.split('\t')[25]
    tissue = line.split('\t')[28]
    experiment_tissue[experiment_id] = exact_tissue
    main_tissue[exact_tissue] = tissue

expression_by_tissue = {}
# expression_by_tissue['protein1']['brain']['individual_id'] = [0.5]

with open(rpkm_file) as infile:
    lines = infile.readlines()[2:]

rpkm_header = lines[0].strip().split('\t')[2:]
print(rpkm_header)
print(len(lines))
for count, line in enumerate(lines[1:]):
    print(count)
    line = line.strip().split('\t')
    protein_id, description = line[0:2]
    protein_id_to_gene[protein_id] = description
    expression_by_tissue[protein_id] = {}
    line = line[2:]
    for i in range(len(line)):
        expression = float(line[i])
        tissue = experiment_tissue[rpkm_header[i]]
#        tissue = main_tissue[tissue]
        individual_id = "-".join(rpkm_header[i].split("-")[:2])
        if tissue not in expression_by_tissue[protein_id].keys():
            expression_by_tissue[protein_id][tissue] = {}
#        if protein_id not in expression_by_tissue[tissue].keys():
#            expression_by_tissue[tissue][protein_id] = {}
        if individual_id not in expression_by_tissue[protein_id][tissue].keys():
            expression_by_tissue[protein_id][tissue][individual_id] = []
        expression_by_tissue[protein_id][tissue][individual_id].append(expression)


for tissue in expression_by_tissue[expression_by_tissue.keys()[0]].keys():
    write_expression_for_tissue(tissue, expression_by_tissue)

