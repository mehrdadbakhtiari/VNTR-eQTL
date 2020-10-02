# Script for comparison between adVNTR-NN and VNTRseek.
# Reference VNTRs for adVNTR and VNTRseek are required to run this module.
# Input: adVNTR genotype result and VNTRseek output and target VNTR files.
# Output: Accuracies for each heterozygous scenario.

from advntr.models import *
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns

advntr_db = "/Users/jonghun/git/vntr/adVNTR/vntr_data/hg19_genic_VNTRs.db"


class VNTRSeekTarget:
    def __init__(self, repeat_id, first_idx, last_idx, copy_number, chr, left_f, pattern, arrayseq, right_f):
        self.repeat_id = repeat_id
        self.first_index = int(first_idx)
        self.last_index = int(last_idx)
        self.copy_number = copy_number
        self.chromosome = chr
        self.left_flanking = left_f
        self.pattern = pattern
        self.array_seq = arrayseq
        self.right_flanking = right_f


# Find common VNTRs in DB (hg19)
def find_common_VNTRs_in_VNTRseek_and_adVNTR(outfile, advntr_db):

    out = open(outfile, "w")
    ignoring_ids = check_same_starting_position_vntrs()

    # Read VNTRSeek database
    # Repeatid, FirstIndex, LastIndex, CopyNumber, FastaHeader, FlankingLeft1000, Pattern, ArraySequence, FlankingRight1000, Conserved
    vntrseek_db = "./vntrseek_reference/t26__.seq"
    vntrseek_reference = defaultdict()
    with open(vntrseek_db, "r") as f:
        header = f.readline()
        for line in f:
            repeat_id, first_index, last_index, copy_number, chromosome, left_flanking, pattern, arrayseq, right_flanking, _ = line.split(",")
            db_element = VNTRSeekTarget(repeat_id, first_index, last_index, copy_number, chromosome, left_flanking, pattern, arrayseq, right_flanking)
            vntrseek_reference[repeat_id] = db_element

    # advntr_db = "/Users/jonghun/git/vntr/adVNTR/vntr_data/hg19_genic_VNTRs.db"
    reference_vntrs = load_unique_vntrs_data(db_file=advntr_db)

    vntrseek_vntr_ids_in_common = []
    advntr_vntr_ids_in_common = []
    for repeat_id, vntrseek_ref_vntr in vntrseek_reference.items():
        if repeat_id not in ignoring_ids:
            for ref_vntr in reference_vntrs:
                if ref_vntr.pattern == vntrseek_ref_vntr.pattern and ref_vntr.chromosome == vntrseek_ref_vntr.chromosome:
                    if vntrseek_ref_vntr.first_index == ref_vntr.start_point + 1:
                        # if abs(vntrseek_ref_vntr.last_index - ref_vntr.start_point + ref_vntr.get_length()-1) < len(ref_vntr.pattern):
                        vntrseek_vntr_ids_in_common.append(vntrseek_ref_vntr.repeat_id)
                        advntr_vntr_ids_in_common.append(ref_vntr.id)
                        out.write(vntrseek_ref_vntr.repeat_id)
                        out.write("\n")
                        out.write(str(ref_vntr.id))
                        out.write("\n")
                        out.write(vntrseek_ref_vntr.chromosome)
                        out.write("\n")
                        out.write(vntrseek_ref_vntr.pattern)
                        out.write("\n")
                        out.write(vntrseek_ref_vntr.array_seq)
                        out.write("\n")
                        out.write(ref_vntr.pattern)
                        out.write("\n")
                        out.write(''.join(ref_vntr.get_repeat_segments()))
                        out.write("\n")
                        out.write(str(vntrseek_ref_vntr.first_index) + " ")
                        out.write(str(vntrseek_ref_vntr.last_index))
                        out.write("\n")
                        out.write(str(ref_vntr.start_point + 1) + " ")
                        out.write(str(ref_vntr.start_point + ref_vntr.get_length()-1))
                        out.write("\n")
                        break

    return vntrseek_vntr_ids_in_common, advntr_vntr_ids_in_common


def create_vntrseek_database(target_ids, vntrseek_ref_name, new_ref_name):
    # We don't change indist file
    # # .indist file
    # with open(new_ref_name + ".indist", "w") as outfile:
    #     with open(vntrseek_ref_name + ".indist", "r") as infile:
    #         line_count = 0
    #         for line in infile:
    #             line_count += 1
    #             if -int(line.strip()) in target_ids:
    #                 outfile.write(line)
    # print("total line", line_count)

    # .leb36 file
    with open(new_ref_name + ".leb36", "w") as outfile:
        with open(vntrseek_ref_name + ".leb36", "r") as infile:
            line_count = 0
            for line in infile:
                line_count += 1
                if int(line.strip().split(" ")[0]) in target_ids:
                    outfile.write(line)

    # .seq file
    with open(new_ref_name + ".seq", "w") as outfile:
        with open(vntrseek_ref_name + ".seq", "r") as infile:
            line_count = 0
            header = infile.readline()
            outfile.write(header)
            for line in infile:
                line_count += 1
                if int(line.split(",")[0]) in target_ids:
                    outfile.write(line)


def get_target_ids(infile):
    with open(infile, "r") as f:
        line_count = 0
        for line in f:
            line_count += 1
            if line.startswith("["):
                ids = [int(s[1:-1]) for s in line.strip()[1:-1].split(", ")]
    return ids


def check_same_starting_position_vntrs():
    vntrseek_db = "./vntrseek_reference/t26__.seq"
    vntrseek_reference = defaultdict()
    first_index_dict = defaultdict(list)
    repeat_id_list = []
    ref_count = 0
    with open(vntrseek_db, "r") as f:
        header = f.readline()
        for line in f:
            ref_count += 1
            repeat_id, first_index, last_index, copy_number, chromosome, left_flanking, pattern, arrayseq, right_flanking, _ = line.split(
                ",")
            db_element = VNTRSeekTarget(repeat_id, first_index, last_index, copy_number, chromosome, left_flanking,
                                        pattern, arrayseq, right_flanking)
            vntrseek_reference[repeat_id] = db_element
            first_index_dict[first_index+chromosome].append(repeat_id)

    print("Check the list size and set size, if it is diff, there is a duplicates")
    print("Total reference count", ref_count)
    print(len(first_index_dict.keys()))
    # Total reference count is 230306
    # Considering the chromosome and starting index, 230228.
    # There are duplicates on starting points
    # That is, 78 are redundant (230306 - 230228)

    # Checked: Repeat ids are unique
    dulpitcate_count = 0
    duplicated_ids = []
    for first_index, repeat_id_list in first_index_dict.items():
        if len(repeat_id_list) > 1:
            dulpitcate_count += 1
            duplicated_ids.append(repeat_id_list)
    # ignore the second ids (only take the first id)
    ignoring_ids = [ids[0] for ids in duplicated_ids]
    return ignoring_ids


def write_vntr_ids_in_common_hg19(advntr_db):
    outfile = "vntrs_in_common_hg19genic.txt"
    vntrseek_vntr_ids_in_common, advntr_ids_in_common = find_common_VNTRs_in_VNTRseek_and_adVNTR(outfile, advntr_db)
    with open("common_vntr_ids_vntrseek_advntr_hg19genic.txt", "w") as f:
        f.write("VNTRseek IDs\n")
        for id in vntrseek_vntr_ids_in_common:
            f.write(str(id) + "\n")

        f.write("adVNTR IDs\n")
        for id in advntr_ids_in_common:
            f.write(str(id) + "\n")


def get_target_vntr_ids(infile):
    vntrseek_ids = []
    advntr_ids = []
    with open(infile, "r") as f:
        is_vntr_seek_ids = False
        for line in f:
            if "VNTRseek IDs" in line:
                is_vntr_seek_ids = True
                continue
            if "adVNTR IDs" in line:
                is_vntr_seek_ids = False
                continue
            if is_vntr_seek_ids:
                vntrseek_ids.append(int(line.strip()))
            else:
                advntr_ids.append(int(line.strip()))

    assert len(vntrseek_ids) == len(advntr_ids), "The size of targets shold be same"
    return vntrseek_ids, advntr_ids


def create_databases_for_runningtime_comp(vntrseek_targets, advntr_targets, target_vntr_counts):
    vntrseek_reference = "./vntrseek_reference/t26__"
    for target_vntr_count in target_vntr_counts:
        new_ref_name = "./vntrseek_reference/t26_" + str(target_vntr_count)
        # create database for vntrseek
        create_vntrseek_database(vntrseek_targets[:target_vntr_count], vntrseek_reference, new_ref_name)

        # write advntr targets as a file
        with open("advntr_vntrseek_common_target_{}.txt".format(target_vntr_count), "w") as f:
            targets = advntr_targets[:target_vntr_count]
            for target in targets:
                f.write(str(target) + "\n")


def create_databases_for_simulation_comp(vntrseek_targets, advntr_targets):
    vntrseek_reference = "./vntrseek_reference/t26__"

    new_ref_name = "./vntrseek_reference/t26_simulation"
    # create database for vntrseek
    create_vntrseek_database(vntrseek_targets, vntrseek_reference, new_ref_name)

    # write advntr targets as a file
    with open("advntr_vntrseek_common_target_simulation.txt", "w") as f:
        targets = advntr_targets
        for target in targets:
            f.write(str(target) + "\n")


def get_targets_for_simulations(advntr_targets):
    advntr_target_list_less = [-1, [], [], []]
    # Sort target vntrs for each chromosome
    advntr_targets.sort(key=lambda x: x.chromosome)

    # modify each chromosome
    previous_chromosome = advntr_targets[0].chromosome
    vntr_in_same_chromosome = []

    total_vntr_count = 0
    for vntr in advntr_targets:
        current_chromosome = vntr.chromosome
        if current_chromosome == previous_chromosome:
            vntr_in_same_chromosome.append(vntr)
        else:
            # print("Processing " + previous_chromosome, len(vntr_in_same_chromosome))
            vntr_in_same_chromosome.sort(key=lambda x: x.start_point)
            for copy_number_change in range(1, 4):
                # Base repeat unit count is 2 (to be called as a tandem repeat)
                target_ids = [ref_vntr.id for ref_vntr in vntr_in_same_chromosome if ref_vntr.estimated_repeats > copy_number_change]
                advntr_target_list_less[copy_number_change].extend(target_ids)

            # Reset for the next
            total_vntr_count += len(vntr_in_same_chromosome)
            vntr_in_same_chromosome = [vntr]
        previous_chromosome = current_chromosome

    # Last chromosome
    # print("Processing " + previous_chromosome, len(vntr_in_same_chromosome))
    vntr_in_same_chromosome.sort(key=lambda x: x.start_point)
    for copy_number_change in range(1, 4):
        # Base repeat unit count is 2 (to be called it is tandem repeat)
        target_ids = [ref_vntr.id for ref_vntr in vntr_in_same_chromosome if ref_vntr.estimated_repeats > copy_number_change]
        advntr_target_list_less[copy_number_change].extend(target_ids)
    total_vntr_count += len(vntr_in_same_chromosome)

    # print("Total processed VNTR count", total_vntr_count)

    return advntr_target_list_less

# 1. Find common VNTR IDs form adVNTR and VNTRseek
write_vntr_ids_in_common_hg19(advntr_db)

# 2. load target VNTR Ids
vntrseek_targets, advntr_targets = get_target_vntr_ids("common_vntr_ids_vntrseek_advntr_hg19genic.txt")
print("Loaded target VNTRs of VNTRseek and adVNTR")
print("Total vntrseek targets", len(vntrseek_targets))

# Load id mapping
advntr_to_vntrseek = dict()
vntrseek_to_advntr = dict()
with open("vntrseek_advntr_id_map.tsv", "r") as f:
    for line in f:
        split = line.strip().split("\t")
        vntrseek_id = int(split[0])
        advntr_id = int(split[1])
        advntr_to_vntrseek[advntr_id] = vntrseek_id
        vntrseek_to_advntr[vntrseek_id] = advntr_id

# 4. Crate target VNTRs for simulation tests (Only difference with the running time comparison is the total length)
# Since both tools are using short-reads for the accuracy comparison, we filter out VNTRs having > 140 bp
reference_vntrs = load_unique_vntrs_data(db_file=advntr_db)
ref_vntr_dict = dict()
for ref_vntr in reference_vntrs:
    ref_vntr_dict[ref_vntr.id] = ref_vntr

vntrseek_simulation_target_ids = []
advntr_simulation_target_ids = []
target_advntr_ref_vntrs = []
for i, target_id in enumerate(advntr_targets):
    ref_vntr = ref_vntr_dict[target_id]
    if ref_vntr.get_length() <= 140:
        vntrseek_simulation_target_ids.append(vntrseek_targets[i])
        advntr_simulation_target_ids.append(advntr_targets[i])
        target_advntr_ref_vntrs.append(ref_vntr)

vntrseek_simulation_target_ids = vntrseek_simulation_target_ids[:10000]
advntr_simulation_target_ids = advntr_simulation_target_ids[:10000]
target_advntr_ref_vntrs = target_advntr_ref_vntrs[:10000]

# print("list sizes")
# print(len(vntrseek_simulation_target_ids))
# print(len(advntr_simulation_target_ids))
# print(len(target_advntr_ref_vntrs))

# 5. Check how many are in indist DB among 100000
indist_ids = set()
with open("vntrseek_reference/t26__.indist", "r") as f:
    for line in f:
        indist_ids.add(-int(line.strip()))

# 6. Add indistinguishable VNTRs
advntr_indist = [555011, 138884, 518405, 282503, 400652, 298261, 200860, 164381, 308514, 248885, 206646, 468922,\
                 46164, 535133, 112479, 186978, 132967, 422253, 355696, 232564, 476537, 330362, 92668]

distinguishable_vntrseek_ids = set()
distinguishable_advntr_ids = set()
distinguishable_advntr_ref_vntrs = []

# Filter IDs  they are in indist
for i, seek_target_id in enumerate(vntrseek_simulation_target_ids):
    if seek_target_id in indist_ids:
        # print(i, seek_target_id)
        continue
    if vntrseek_to_advntr[seek_target_id] in advntr_indist:
        continue
    else:
        distinguishable_vntrseek_ids.add(seek_target_id)
        distinguishable_advntr_ids.add(vntrseek_to_advntr[seek_target_id])
        distinguishable_advntr_ref_vntrs.append(ref_vntr_dict[vntrseek_to_advntr[seek_target_id]])

# print("final distinguishable advntr target count", len(distinguishable_advntr_ids))
# print("final distinguishable vntrseek target count", len(distinguishable_vntrseek_ids))
# print("final distinguishable advntr ref target count", len(distinguishable_advntr_ref_vntrs))

advntr_target_list_less = get_targets_for_simulations(distinguishable_advntr_ref_vntrs)
# print(distinguishable_advntr_ids - set([vntr.id for vntr in distinguishable_advntr_ref_vntrs]))
# print(distinguishable_advntr_ids - set(advntr_target_list_less[1]))
# print(set(advntr_target_list_less[1]) - distinguishable_advntr_ids)

print("number of target VNTRs for contraction scenarios")
print(len(advntr_target_list_less[1]))
print(len(advntr_target_list_less[2]))
print(len(advntr_target_list_less[3]))

# Calculate accuracy for adVNTR-NN
modes = ['more', 'less']
max_edited_unit = 3
# Usage: result[1more][vntr_id] = True
advntr_result_scenario_total_length_map = defaultdict(lambda: defaultdict(list))
advntr_result_scenario_repeat_length_map = defaultdict(lambda: defaultdict(list))

advntr_result_scenario_id_map = defaultdict(lambda: defaultdict(int))
advntr_result_id_scenario_map = defaultdict(lambda: defaultdict(int))
for mode in modes:
    for edited_unit in range(1, max_edited_unit + 1):
        result_file_path = "simul_genotype_result/"
        result_file = result_file_path + "{}{}.result".format(edited_unit, mode)
        target_advntr_set = distinguishable_advntr_ids
        if mode == "less":
            target_advntr_set = set(advntr_target_list_less[edited_unit])

        not_in_target_count = 0
        with open(result_file, "r") as f:
            for line in f:
                if "/" in line or "None" in line:
                    # Check if this vntr ID is in the targets
                    if vntr_id not in target_advntr_set:
                        # print("not in the target", vntr_id, advntr_to_vntrseek[vntr_id])
                        not_in_target_count += 1
                        continue

                    if "None" in line:
                        advntr_result_scenario_id_map[str(edited_unit) + mode][vntr_id] = False
                        advntr_result_id_scenario_map[vntr_id][str(edited_unit) + mode] = False

                        # Counts twice because it is haplotype-based result
                        advntr_result_scenario_repeat_length_map[str(edited_unit) + mode][len(ref_vntr_dict[vntr_id].pattern)].append(False)
                        advntr_result_scenario_repeat_length_map[str(edited_unit) + mode][len(ref_vntr_dict[vntr_id].pattern)].append(False)

                        hap1_length = ref_vntr_dict[vntr_id].get_length()
                        hap2_length = ref_vntr_dict[vntr_id].get_length()
                        if mode == "more":
                            hap2_length += len(ref_vntr_dict[vntr_id].pattern) * edited_unit
                        else:
                            hap2_length -= len(ref_vntr_dict[vntr_id].pattern) * edited_unit

                        advntr_result_scenario_total_length_map[str(edited_unit) + mode][hap1_length].append(False)
                        advntr_result_scenario_total_length_map[str(edited_unit) + mode][hap2_length].append(False)
                        continue

                    genotype = line.strip().split("/")
                    estimated_repeat_count = ref_vntr_dict[vntr_id].estimated_repeats
                    hap1_length = ref_vntr_dict[vntr_id].get_length()
                    hap2_length = ref_vntr_dict[vntr_id].get_length()
                    if mode == "more":
                        hap2_length += len(ref_vntr_dict[vntr_id].pattern) * edited_unit
                    else:
                        hap2_length -= len(ref_vntr_dict[vntr_id].pattern) * edited_unit

                    hap1 = int(genotype[0])
                    hap2 = int(genotype[1])

                    if mode == "more":
                        if hap1 == estimated_repeat_count:
                            advntr_result_scenario_total_length_map[str(edited_unit) + mode][hap1_length].append(True)
                            advntr_result_scenario_repeat_length_map[str(edited_unit) + mode][len(ref_vntr_dict[vntr_id].pattern)].append(True)
                        else:
                            advntr_result_scenario_total_length_map[str(edited_unit) + mode][hap1_length].append(False)
                            advntr_result_scenario_repeat_length_map[str(edited_unit) + mode][
                                len(ref_vntr_dict[vntr_id].pattern)].append(False)
                        if hap2 == estimated_repeat_count + edited_unit:
                            advntr_result_scenario_total_length_map[str(edited_unit) + mode][hap2_length].append(True)
                            advntr_result_scenario_repeat_length_map[str(edited_unit) + mode][
                                len(ref_vntr_dict[vntr_id].pattern)].append(True)
                        else:
                            advntr_result_scenario_total_length_map[str(edited_unit) + mode][hap2_length].append(False)
                            advntr_result_scenario_repeat_length_map[str(edited_unit) + mode][
                                len(ref_vntr_dict[vntr_id].pattern)].append(False)

                        if int(genotype[0]) == estimated_repeat_count and int(genotype[1]) == estimated_repeat_count + edited_unit:
                            advntr_result_scenario_id_map[str(edited_unit) + mode][vntr_id] = True
                            advntr_result_id_scenario_map[vntr_id][str(edited_unit) + mode] = True
                        else:
                            advntr_result_scenario_id_map[str(edited_unit) + mode][vntr_id] = False
                            advntr_result_id_scenario_map[vntr_id][str(edited_unit) + mode] = False
                    else:
                        if hap1 == estimated_repeat_count - edited_unit:
                            advntr_result_scenario_total_length_map[str(edited_unit) + mode][hap1_length].append(True)
                            advntr_result_scenario_repeat_length_map[str(edited_unit) + mode][
                                len(ref_vntr_dict[vntr_id].pattern)].append(True)
                        else:
                            advntr_result_scenario_total_length_map[str(edited_unit) + mode][hap1_length].append(False)
                            advntr_result_scenario_repeat_length_map[str(edited_unit) + mode][
                                len(ref_vntr_dict[vntr_id].pattern)].append(False)
                        if hap2 == estimated_repeat_count:
                            advntr_result_scenario_total_length_map[str(edited_unit) + mode][hap2_length].append(True)
                            advntr_result_scenario_repeat_length_map[str(edited_unit) + mode][
                                len(ref_vntr_dict[vntr_id].pattern)].append(True)
                        else:
                            advntr_result_scenario_total_length_map[str(edited_unit) + mode][hap2_length].append(False)
                            advntr_result_scenario_repeat_length_map[str(edited_unit) + mode][
                                len(ref_vntr_dict[vntr_id].pattern)].append(False)

                        if int(genotype[0]) == estimated_repeat_count - edited_unit and int(genotype[1]) == estimated_repeat_count:
                            advntr_result_scenario_id_map[str(edited_unit) + mode][vntr_id] = True
                            advntr_result_id_scenario_map[vntr_id][str(edited_unit) + mode] = True
                        else:
                            advntr_result_scenario_id_map[str(edited_unit) + mode][vntr_id] = False
                            advntr_result_id_scenario_map[vntr_id][str(edited_unit) + mode] = False
                else:
                    vntr_id = int(line.strip())

        # print(not_in_target_count, mode, edited_unit)

advntr_accuracy_per_vntr = open("advntr_accuracy_per_vntr.csv", "w")
advntr_accuracy_list_per_vntr = []
for advntr_id in distinguishable_advntr_ids:
    advntr_true_count = advntr_result_id_scenario_map[advntr_id].values().count(True)
    possible_scenario_count = len(advntr_result_id_scenario_map[advntr_id].keys())

    accuracy = str(float(advntr_true_count)/possible_scenario_count*100)
    advntr_accuracy_list_per_vntr.append(accuracy)

    # Check when adVNTR's performance is poor
    # if accuracy == "0.0":
    #     print("Accuracy is zero!")
    #     print(advntr_id, advntr_to_vntrseek[advntr_id])

    advntr_accuracy_per_vntr.write(str(advntr_id) + "," + accuracy + "\n")
advntr_accuracy_per_vntr.close()

for mode in modes:
    for edited_unit in range(1, max_edited_unit + 1):
        correct_count = advntr_result_scenario_id_map[str(edited_unit) + mode].values().count(True)
        total_target_count = len(advntr_result_scenario_id_map[str(edited_unit) + mode].keys())
        accuracy = float(correct_count) / total_target_count * 100

        print("{} correct count".format(correct_count))
        print("{} total VNTR count".format(total_target_count))
        print("advntr {}{} accuracy: {}".format(edited_unit, mode, accuracy))

def plot_advntr_accuracy():
    global scenario, vntr_id, ref_vntr, acc, length, result_bool_list, red_color, blue_color
    # For revision (eVNTR paper)
    # Accuracy stratified by Repeat unit size and total length (accuracy)
    repeat_unit_bin_size = 6
    total_length_bin_size = 30
    accuracy_stratified_by_repeat_unit = defaultdict(lambda: defaultdict(list))
    accuracy_stratified_by_total_length = defaultdict(lambda: defaultdict(list))
    for scenario in advntr_result_scenario_id_map.keys():
        for vntr_id in advntr_result_scenario_id_map[scenario].keys():
            ref_vntr = ref_vntr_dict[vntr_id]
            pattern_length = len(ref_vntr.pattern)
            total_length = ref_vntr.get_length()

            accuracy_stratified_by_repeat_unit[scenario][pattern_length // repeat_unit_bin_size].append(
                advntr_result_scenario_id_map[scenario][vntr_id])

            # Bin - modified length?
            how_many = int(scenario[0])
            if "more" in scenario:
                total_length += pattern_length * how_many
            elif "less" in scenario:
                total_length -= pattern_length * how_many
            else:
                print("ERROR in scenario name", scenario)
                exit(-1)
            accuracy_stratified_by_total_length[scenario][total_length // total_length_bin_size].append(
                advntr_result_scenario_id_map[scenario][vntr_id])
    print("repeat unit accuracy")
    for scenario in accuracy_stratified_by_repeat_unit.keys():
        print("\nscenario {}".format(scenario))
        for bin in accuracy_stratified_by_repeat_unit[scenario].keys():
            # print("bin: {}".format(bin * repeat_unit_bin_size))
            acc = float(accuracy_stratified_by_repeat_unit[scenario][bin].count(True)) / len(
                accuracy_stratified_by_repeat_unit[scenario][bin])
            print("{}".format(acc))
    # Averaging for all VNTRs in each bin
    print("Average of repeat unit accuracy")
    for bin in accuracy_stratified_by_repeat_unit[scenario].keys():
        print("bin {}".format(bin * repeat_unit_bin_size))
        total_bin_count = 0.0
        weighted_average = .0
        for scenario in accuracy_stratified_by_repeat_unit.keys():
            total_bin_count += len(accuracy_stratified_by_repeat_unit[scenario][bin])

        for scenario in accuracy_stratified_by_repeat_unit.keys():
            # print("scenario {}".format(scenario))
            # print("bin: {}".format(bin * repeat_unit_bin_size))
            current_bin_count = float(len(accuracy_stratified_by_repeat_unit[scenario][bin]))
            if current_bin_count == 0:
                continue
            weighted_average += float(
                accuracy_stratified_by_repeat_unit[scenario][bin].count(True)) / current_bin_count * (
                                            current_bin_count / total_bin_count)

        print("Averaged bin accuracy\t{}".format(weighted_average))
    print("Total length accuracy")
    for scenario in accuracy_stratified_by_total_length.keys():
        print("\nscenario {}".format(scenario))
        print("vntr id num {}".format(len(accuracy_stratified_by_total_length[scenario].keys())))
        for bin in accuracy_stratified_by_total_length[scenario].keys():
            # print("bin: {}".format(bin * total_length_bin_size))
            acc = float(accuracy_stratified_by_total_length[scenario][bin].count(True)) / len(
                accuracy_stratified_by_total_length[scenario][bin])
            print("{}".format(acc))
    print("Accuracy by repeat length - Hayplotype-based")
    repeat_length_count = defaultdict(int)
    repeat_length_acc = defaultdict(float)
    for scenario in advntr_result_scenario_repeat_length_map.keys():
        for length in sorted(advntr_result_scenario_repeat_length_map[scenario].keys()):
            repeat_length_count[length] += len(advntr_result_scenario_repeat_length_map[scenario][length])
    for scenario in advntr_result_scenario_repeat_length_map.keys():
        for length in sorted(advntr_result_scenario_repeat_length_map[scenario].keys()):
            result_bool_list = advntr_result_scenario_repeat_length_map[scenario][length]
            acc = float(result_bool_list.count(True)) / float(len(result_bool_list))
            repeat_length_acc[length] += float(result_bool_list.count(True)) / repeat_length_count[length] * 100
    sns.set()
    sns.set_style("ticks")
    sns.set_context("paper")
    red_color = '#d62728'
    blue_color = "#1f77b4"
    plt.bar(repeat_length_acc.keys(), repeat_length_acc.values())
    # plt.scatter(repeat_length_acc.keys(), repeat_length_acc.values(), s=4, color='k')
    # plt.xlim((0, 200))
    plt.xlabel("Repeat unit length", fontsize=12)
    plt.ylabel("Accuracy (%)", fontsize=12)
    plt.xticks(fontsize=11)
    plt.yticks(fontsize=11)
    # Black border
    plt.gca().spines['bottom'].set_color('black')
    plt.gca().spines['left'].set_color('black')
    plt.gca().spines['top'].set_color('black')
    plt.gca().spines['right'].set_color('black')
    plt.savefig("revision_adVNTR_accuracy_by_repeat_length_bar.png", dpi=300)

    print("Accuracy by total length - Separated by haplotype")
    total_length_count = defaultdict(int)
    total_length_acc = defaultdict(float)

    for scenario in advntr_result_scenario_total_length_map.keys():
        for length in sorted(advntr_result_scenario_total_length_map[scenario].keys()):
            total_length_count[length] += len(advntr_result_scenario_total_length_map[scenario][length])

    for scenario in advntr_result_scenario_total_length_map.keys():
        hap_lengths = []
        hap_accuracies = []
        for length in sorted(advntr_result_scenario_total_length_map[scenario].keys()):
            result_bool_list = advntr_result_scenario_total_length_map[scenario][length]
            acc = float(result_bool_list.count(True)) / float(len(result_bool_list))
            total_length_acc[length] += float(result_bool_list.count(True)) / total_length_count[length] * 100

            hap_lengths.append(length)
            hap_accuracies.append(acc)

    sns.set()
    sns.set_style("ticks")
    sns.set_context("paper")
    red_color = '#d62728'
    blue_color = "#1f77b4"

    plt.figure()
    plt.scatter(total_length_acc.keys(), total_length_acc.values(), s=4, color='k')
    plt.xlim((0, 200))
    plt.xlabel("VNTR Length", fontsize=12)
    plt.ylabel("Accuracy (%)", fontsize=12)
    plt.xticks(fontsize=11)
    plt.yticks(fontsize=11)

    # Black border
    plt.gca().spines['bottom'].set_color('black')
    plt.gca().spines['left'].set_color('black')
    plt.gca().spines['top'].set_color('black')
    plt.gca().spines['right'].set_color('black')

    # plt.show()
    plt.savefig("revision_adVNTR_accuracy_by_total_length_dotted.png", dpi=300)

# plot_advntr_accuracy()

# Calculate accuracy for VNTRseek
modes = ['more', 'less']
max_edited_unit = 3
# Usage: result[1more][vntr_id] = True
vntrseek_result_scenario_id_map = defaultdict(lambda: defaultdict(int))
vntrseek_result_id_scenario_map = defaultdict(lambda: defaultdict(int))
for mode in modes:
    for edited_unit in range(1, max_edited_unit + 1):
        # vcf_file_path = "simul_genotype_result/"
        vcf_file_path = "./"
        vcf = vcf_file_path + "VNTRPIPE_simul_{}{}.span2.vcf".format(edited_unit, mode)
        vcf = vcf_file_path + "VNTRPIPE_simul_{}{}.allwithsupport.span2.vcf".format(edited_unit, mode)

        target_advntr_set = distinguishable_advntr_ids
        if mode == "less":
            target_advntr_set = set(advntr_target_list_less[edited_unit])
        not_in_target_count = 0

        with open(vcf, "r") as f:
            for line in f:
                if line.startswith("#"):
                    if "##numVNTRs=" in line:
                        num_of_VNTRs = line.strip().split("=")[-1]
                        # print(num_of_VNTRs, edited_unit, mode)
                    continue

                # #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  simul_1more
                split = line.strip().split("\t")
                vid = int(split[2][2:])
                advntr_id = vntrseek_to_advntr[int(vid)]
                if advntr_id not in target_advntr_set:
                    # print("not in the target", vntr_id, advntr_to_vntrseek[vntr_id])
                    not_in_target_count += 1
                    continue

                ref_vntr = ref_vntr_dict[int(vntrseek_to_advntr[int(vid)])]
                ref = split[3]
                alt = split[4]
                info = split[7]
                pattern = info.split(";")[3].split("=")[-1]
                # GT:SP:CGL	0/1:7,4:0,3
                format_col = split[9]
                genotype = format_col.split(":")[0]
                copy_changes = format_col.split(":")[-1].split(",")

                alt_list = alt.strip().split(",")
                if len(alt_list) >= 2:
                    vntrseek_result_scenario_id_map[str(edited_unit) + mode][vid] = False
                    vntrseek_result_id_scenario_map[vid][str(edited_unit) + mode] = False
                    continue
                if '.' in copy_changes:
                    vntrseek_result_scenario_id_map[str(edited_unit) + mode][vid] = False
                    vntrseek_result_id_scenario_map[vid][str(edited_unit) + mode] = False
                    continue
                if len(copy_changes) == 1:
                    vntrseek_result_scenario_id_map[str(edited_unit) + mode][vid] = False
                    vntrseek_result_id_scenario_map[vid][str(edited_unit) + mode] = False
                    continue

                if mode == 'more':
                    if int(copy_changes[0]) == 0 and int(copy_changes[1]) == edited_unit:
                        vntrseek_result_scenario_id_map[str(edited_unit) + mode][vid] = True
                        vntrseek_result_id_scenario_map[vid][str(edited_unit) + mode] = True
                    # if abs(len(alt) - len(ref)) == len(ref_vntr.pattern) * edited_unit:
                    #     vntrseek_result[str(edited_unit) + mode][vid] = True
                    else:
                        vntrseek_result_scenario_id_map[str(edited_unit) + mode][vid] = False
                        vntrseek_result_id_scenario_map[vid][str(edited_unit) + mode] = False
                else:  # mode == 'less'
                    if int(copy_changes[0]) == 0 and int(copy_changes[1]) == -edited_unit:
                        vntrseek_result_scenario_id_map[str(edited_unit) + mode][vid] = True
                        vntrseek_result_id_scenario_map[vid][str(edited_unit) + mode] = True
                    else:
                        vntrseek_result_scenario_id_map[str(edited_unit) + mode][vid] = False
                        vntrseek_result_id_scenario_map[vid][str(edited_unit) + mode] = False

            # For VNTRs not appeared in result vcf
            for advntr_id in target_advntr_set:
                vntrseek_id = advntr_to_vntrseek[advntr_id]
                if vntrseek_id not in vntrseek_result_scenario_id_map[str(edited_unit) + mode]:
                    vntrseek_result_scenario_id_map[str(edited_unit) + mode][vntrseek_id] = False
                    vntrseek_result_id_scenario_map[vntrseek_id][str(edited_unit) + mode] = False

vntrseek_accuracy_per_vntr = open("vntrseek_accuracy_per_vntr.csv", "w")
vntrseek_accuracy_list_per_vntr = []
for advntr_id in distinguishable_advntr_ids:
    vntrseek_id = advntr_to_vntrseek[advntr_id]
    advntr_true_count = vntrseek_result_id_scenario_map[vntrseek_id].values().count(True)
    possible_scenario_count = len(vntrseek_result_id_scenario_map[vntrseek_id].keys())

    accuracy = str(float(advntr_true_count)/possible_scenario_count*100)
    vntrseek_accuracy_list_per_vntr.append(accuracy)

    vntrseek_accuracy_per_vntr.write(str(vntrseek_id) + "," + accuracy + "\n")
vntrseek_accuracy_per_vntr.close()

for mode in modes:
    for edited_unit in range(1, max_edited_unit + 1):
        correct_count = vntrseek_result_scenario_id_map[str(edited_unit) + mode].values().count(True)
        total_target_count = len(vntrseek_result_scenario_id_map[str(edited_unit) + mode].keys())
        accuracy = float(correct_count) / total_target_count * 100

        print("{} correct count".format(correct_count))
        print("{} total VNTR count".format(total_target_count))
        print("vntrseek {}{} accuracy: {}".format(edited_unit, mode, accuracy))

