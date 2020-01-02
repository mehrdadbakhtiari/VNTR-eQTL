

class GeneLocations:
    def __init__(self, directory='/pedigree2/projects/adVNTR/hg38_annotation/', knownToEnsembl_file='knownToEnsembl.txt', ensemblToGeneName_file='ensemblToGeneName.txt', genes_coordinates_file='refseq_genes.bed', refseq_to_gene_file='Refseq2Gene.txt'):
        self.ensembl_to_gene_name = {}
        self.ucsc_to_gene_name = {}
        self.gene_coordinates = {}
        self.refseq_to_gene = {}

        knownToEnsembl_file = directory + knownToEnsembl_file
        ensemblToGeneName_file = directory + ensemblToGeneName_file
        genes_coordinates_file = directory + genes_coordinates_file
        refseq_to_gene_file = directory + refseq_to_gene_file

#        with open(ensemblToGeneName_file) as infile:
#            lines = infile.readlines()
#        for line in lines:
#            line = line.strip().split()
#            self.ensembl_to_gene_name[line[0]] = line[1]
#
#        with open(knownToEnsembl_file) as infile:
#            lines = infile.readlines()
#        for line in lines:
#            line = line.strip().split()
#            self.ucsc_to_gene_name[line[0]] = self.ensembl_to_gene_name[line[1]]

        with open(refseq_to_gene_file) as infile:
            lines = infile.readlines()
        for line in lines:
            line = line.strip().split()
            self.refseq_to_gene[line[0]] = line[1]

        with open(genes_coordinates_file) as infile:
            lines = infile.readlines()
#        known_ucsc_ids = set(self.ucsc_to_gene_name.keys())
        known_refseq_ids = set(self.refseq_to_gene.keys())
        for line in lines:
            line = line.strip().split()
            if len(line[0]) > len('chrXX'):
                continue
            if line[3].split('.')[0] not in known_refseq_ids:
                continue
#            self.gene_coordinates[self.ucsc_to_gene_name[line[3]]] = (int(line[1]), int(line[2]))
            self.gene_coordinates[self.refseq_to_gene[line[3].split('.')[0]]] = (int(line[1]), int(line[2]))

    def get_gene_coordinates(self, gene_name):
        if gene_name not in self.gene_coordinates.keys():
            return None, None
        return self.gene_coordinates[gene_name]

