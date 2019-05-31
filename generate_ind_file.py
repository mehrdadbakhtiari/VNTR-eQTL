"""
Generate ind file with population labels

Usage: python generate_ind_file.py sample_file > ind_file
"""

import sys

try:
    indfile = sys.argv[1]
    sampfile = sys.argv[2]
except:
    sys.stderr.write(__doc__)
    sys.exit(1)

sample_to_gender = {}
sample_to_pop = {}
labels = ["Asian", "AfricanAmerican", "European", "Amerindian", "NotReported", "Unknown"]
nums = [1, 2, 3, 4, 98, 99]
num_to_pop = dict(zip(nums, labels))

samples = []

with open(sampfile, "r") as f:
    for line in f:
        if line.startswith("#"): continue
        if line.strip() == "": continue
        if line.startswith("dbGaP_Subject_ID"): continue
        items = line.strip().split("\t")
        try:
            sample = items[1]
            samples.append(sample)
            pop = int(items[5])
            sex = 'M' if int(items[3]) == 1 else 'F'
            sample_to_pop[sample] = num_to_pop[int(pop)]
            sample_to_gender[sample] = sex
        except IndexError: continue

with open(indfile, "r") as f:
    for line in f:
        items = line.strip().split()
        sampid = "-".join(items[0].split("-")[0:2])
        pop = sample_to_pop.get(sampid, "NA")
        items[2] = pop
        sex = sample_to_gender.get(sampid, "NA")
        items[1] = sex
        sys.stdout.write("\t".join(items)+"\n")

