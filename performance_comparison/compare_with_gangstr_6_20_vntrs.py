# Script for comparison between adVNTR-NN and GangSTR.
# Reference VNTRs for adVNTR is required.
# Input: adVNTR genotype result and GangSTR's output and target VNTR file.
# Output: Accuracy for each heterozygous scenario.

from advntr import models
from collections import defaultdict

advntr_db = '/home/jonghun/advntr_vcf/adVNTR/vntr_data/hg38_selected_VNTRs_Illumina.db'

####################################################################
# Get similar VNTR IDs for filtering
similar_vntrs_hg38 = set()
with open("similar_vntrs_hg38_maxlen1000.txt", "r") as f:
    for line in f:
        similar_vntrs_hg38.add(int(line))

# Read genotype results for those and the target VNTRs
reference_vntrs = models.load_unique_vntrs_data(advntr_db)
ref_vntrs = {ref_vntr.id: ref_vntr for ref_vntr in reference_vntrs}
target_vntrs = [ref_vntr for ref_vntr in reference_vntrs if 6 <= len(ref_vntr.pattern) <= 20 and ref_vntr.id not in similar_vntrs_hg38]
target_vids = set([ref_vntr.id for ref_vntr in target_vntrs])

####################################################################
# Calculate accuracy for each secnario
print("This script calculates accuracies for each scenarios for adVNTR-NN and GangSTR")

# chromosome + start_point to vid map for GangSTR
chr_start_to_vid = defaultdict(int)
for t in target_vntrs:
    chr_start_to_vid[t.chromosome + str(t.start_point)] = t.id

def convert_vid_to_chr_position(ref_vntr):
    return ref_vntr.chromosome+str(ref_vntr.start_point)

def get_ref_repeat(ref_vntr):
    return ref_vntr.estimated_repeats

gangstr_genotypes = defaultdict(lambda: defaultdict(str))
advntr_genotypes = defaultdict(lambda: defaultdict(str))

gangstr_acc = defaultdict(lambda: defaultdict(bool))
advntr_acc = defaultdict(lambda: defaultdict(bool))

for mode in ['more', 'less']:
    for edited_copy in range(1, 4):
        scenario = str(edited_copy) + mode

        gangstr_result = "{}_{}_gangstr.genotype".format(edited_copy, mode)
        advntr_result = "{}_{}_advntr.genotype".format(edited_copy, mode)
        # Load adVNTR genotype
        with open(advntr_result, "r") as f:
            for line in f:
                if "None" in line or "/" in line:
                    if vid not in target_vids:
                        continue
                    genotype = line.strip()
                    advntr_genotypes[scenario][vid] = genotype

                    ref_repeat = get_ref_repeat(ref_vntrs[vid])
                    correct_genotype = ""
                    if mode == "more":
                        correct_genotype = str(ref_repeat) + '/' + str(ref_repeat + edited_copy)
                    else:
                        if ref_repeat - edited_copy < 0:
                            continue
                        correct_genotype = str(ref_repeat - edited_copy) + '/' + str(ref_repeat)

                    if correct_genotype == genotype:
                        advntr_acc[scenario][vid] = True

                else:
                    vid = int(line.strip())

        # Load GangSTR genotype
        with open(gangstr_result, "r") as f:
            for line in f:
                if "Processing" in line:
                    vntr_locus = line.strip().split(" ")[-1].split(":")
                    chromosome = vntr_locus[0]
                    start_point = vntr_locus[1]
                if "WARNING:" in line:
                    formated_genotype = "None"
                    vid = chr_start_to_vid[chromosome+start_point]
                    if vid not in target_vids:
                        continue
                if "Not enough reads extracted" in line:
                    formated_genotype = "None"
                    vid = chr_start_to_vid[chromosome+start_point]
                    if vid not in target_vids:
                        continue
                if "Genotyper" in line:
                    split = line.strip().split("\t")
                    genotype = split[1].split(" ")[3:]
                    formated_genotype = genotype[0][:-1] + "/" + genotype[1]
                    vid = chr_start_to_vid[chromosome+start_point]
                    if vid not in target_vids:
                        continue
                    gangstr_genotypes[scenario][chromosome+start_point] = formated_genotype

                    ref_repeat = get_ref_repeat(ref_vntrs[vid])
                    correct_genotype = ""
                    if mode == "more":
                        correct_genotype = str(ref_repeat) + '/' + str(ref_repeat + edited_copy)
                    else:
                        if ref_repeat - edited_copy < 0:
                            continue
                        correct_genotype = str(ref_repeat - edited_copy) + '/' + str(ref_repeat)

                    if correct_genotype == formated_genotype:
                        gangstr_acc[scenario][vid] = True

ru2less_target_vids = set()
ru3less_target_vids = set()
for vid in target_vids:
    if ref_vntrs[vid].estimated_repeats >= 3:
        ru2less_target_vids.add(vid)
    if ref_vntrs[vid].estimated_repeats >= 4:
        ru3less_target_vids.add(vid)

for scenario in advntr_genotypes.keys():
    print("Scenario: {}".format(scenario))

    edited_copy = int(scenario[0])
    mode = scenario[1:]
    target_count = len(target_vntrs)
    if scenario == "2less":
        target_count = len(ru2less_target_vids)
    elif scenario == "3less":
        target_count = len(ru3less_target_vids)
    else:
        target_count = len(target_vntrs)

    advntr_correct_count = 0
    for vid in advntr_acc[scenario].keys():
        if scenario == "2less":
            if vid in ru2less_target_vids and advntr_acc[scenario][vid]:
                advntr_correct_count += 1
        elif scenario == "3less" and advntr_acc[scenario][vid]:
            if vid in ru3less_target_vids:
                advntr_correct_count += 1
        else:
            if advntr_acc[scenario][vid]:
                advntr_correct_count += 1

    gangstr_correct_count = 0
    for vid in gangstr_acc[scenario].keys():
        if scenario == "2less":
            if vid in ru2less_target_vids and gangstr_acc[scenario][vid]:
                gangstr_correct_count += 1
        elif scenario == "3less" and gangstr_acc[scenario][vid]:
            if vid in ru3less_target_vids:
                gangstr_correct_count += 1
        else:
            if gangstr_acc[scenario][vid]:
                gangstr_correct_count += 1

    print("adVNTR-NN accuracy: {}".format(float(advntr_correct_count)/(target_count)))
    print("GangSTR accuracy: {}".format(float(gangstr_correct_count)/(target_count)))