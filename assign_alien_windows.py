import sys

prefix = sys.argv[1]
snp_threshold = int(sys.argv[2])
hom_het_ratio = int(sys.argv[3])
cov_dev_threshold = float(sys.argv[4])
distance_for_cov_dev_block_merging = int(sys.argv[5])
block_proportion = float(sys.argv[6])
cov_dev_block_proportion = float(sys.argv[7])

###load chromosome lengths
chromosome_lengths = {}
with open(sys.argv[8], 'r') as ref_fai:
    for line in ref_fai:
        temp = line.rstrip( ).split( )
        chromosome_lengths[temp[0]] = int(temp[1])

cov_dev_dict = {}
cov_dev_blocks = []
current_block = []

with open(prefix + "_cov_dev_1Mb.tsv", 'r') as cov_dev_file:
    for line in cov_dev_file:
        if not "chrUn" in line:
            temp = line.rstrip( ).split( )
            cov_dev_dict[temp[0]+":"+temp[1]] = float(temp[2])

            if float(temp[2]) < cov_dev_threshold:
                if current_block:
                    previous_cov_dev_chromosome = current_block[-1].split(":")[0]
                    previous_cov_dev_window = current_block[-1].split(":")[1]
                    if previous_cov_dev_chromosome == temp[0] and (int(temp[1]) - int(previous_cov_dev_window)) <= distance_for_cov_dev_block_merging:
                        current_block.append(temp[0]+":"+temp[1])
                    else:
                        cov_dev_blocks.append(current_block)
                        current_block = []
                        current_block.append(temp[0]+":"+temp[1])
                else:
                    current_block.append(temp[0]+":"+temp[1])

    cov_dev_blocks.append(current_block)

merged_cov_dev_blocks = [[i[0], i[-1]] for i in cov_dev_blocks if int(i[-1].split(":")[1]) - int(i[0].split(":")[1]) >= 1000000]

window_ids_within_merged_cov_dev_blocks = []

filtered_merged_cov_dev_blocks = []

current_block_windows = []

for i in merged_cov_dev_blocks:
    chromosome = i[0].split(":")[0]
    start = int(i[0].split(":")[1])
    end = int(i[1].split(":")[1])

    count = 0
    for window in range(start, end+1000000, 1000000):
        if cov_dev_dict.get(chromosome+":"+str(window)) <= cov_dev_threshold:
            count += 1

        current_block_windows.append(chromosome+":"+str(window))

    if count >= (cov_dev_block_proportion * ((end+1000000 - start) / 1000000)):
        for window in current_block_windows:
            window_ids_within_merged_cov_dev_blocks.append(window)
        filtered_merged_cov_dev_blocks.append(i)

muticum_windows = []

for line in sys.stdin:
    if not "chrUn" in line:
        temp = line.rstrip( ).split( )

        if int(temp[2]) >= snp_threshold:
            if int(temp[3]) == 0:
                if temp[0]+":"+temp[1] in window_ids_within_merged_cov_dev_blocks:
                    muticum_windows.append(temp[0]+":"+temp[1])
            elif int(temp[2]) / int(temp[3]) >= hom_het_ratio:
                if temp[0]+":"+temp[1] in window_ids_within_merged_cov_dev_blocks:
                    muticum_windows.append(temp[0]+":"+temp[1])

final_segments_coarse = []
count = 0
for i in filtered_merged_cov_dev_blocks:
    chromosome = i[0].split(":")[0]
    start = int(i[0].split(":")[1])
    end = int(i[1].split(":")[1])
    running_mean = []
    for window in range(start, end+1000000, 1000000):
        if chromosome+":"+str(window) in muticum_windows:
            count += 1
        running_mean.append([window, count])
    if (end - start) / 1000000 <= 5:
        if count > (0.50 * ((end - start) / 1000000)):
            final_segments_coarse.append([chromosome, start, end])
    elif count > (block_proportion * ((end - start) / 1000000)):
        final_segments_coarse.append([chromosome, start, end])
    else:
        last_increase = ""
        block_from_snps_coarse = []
        for n, i in enumerate(running_mean):
            if i[0] - running_mean[0][0] == 0:
                if i[1] == 1:
                    last_increase = i[0]
            else:
                if i[1] - running_mean[n-1][1] == 1:
                    last_increase = i[0]

                if i[1] / ((i[0] - running_mean[0][0]) / 1000000) >= block_proportion:
                    block_from_snps_coarse.append(i)
                else:
                    if last_increase:
                        final_segments_coarse.append([chromosome, start, last_increase])
                    break

        ###same but from the end of chromosome walking backwards
        running_mean = []
        for window in range(start, end+1000000, 1000000):
            if chromosome+":"+str(window) in muticum_windows:
                count += 1
        running_mean.append([window, count])
        last_increase = ""
        block_from_snps_coarse = []
        for n, i in enumerate(running_mean[::-1]):
            if running_mean[-1][0] - i[0] == 0:
                if i[1] == 1:
                    last_increase = i[0]
            else:
                if i[1] - running_mean[n+1][1] == 1:
                    last_increase = i[0]

                if i[1] / ((running_mean[-1][0] - i[0]) / 1000000) >= block_proportion:
                    block_from_snps_coarse.append(i)
                else:
                    if last_increase:
                        final_segments_coarse.append([chromosome, start, last_increase])
                    break

    count = 0

final_segments_fine = []

###make cov dev dict with 100Kbp windows
cov_dev_dict_100Kb = {}

cov_dev_blocks = []
current_block = []

with open(prefix + "_cov_dev_100Kb.tsv", 'r') as cov_dev_file:
    for line in cov_dev_file:
        if not "chrUn" in line:
            temp = line.rstrip( ).split( )
            cov_dev_dict_100Kb[temp[0]+":"+temp[1]] = float(temp[2])

with open(prefix + "_segments.txt", 'w') as out_file:
    for segment in final_segments_coarse:
        chromosome = segment[0]
        start = int(segment[1])
        end = int(segment[2])

        chromosome_length = chromosome_lengths.get(chromosome)

        out_file.write("\t".join(["Coarse segment identification:", chromosome, str(start), str(end)]) + "\n")
        left_junction_candidates = []
        right_junction_candidates = []

        if start == 0:
            for i in range(start, start+1000000, 100000):
                if cov_dev_dict_100Kb.get(chromosome +":"+str(i)):
                    if cov_dev_dict_100Kb.get(chromosome +":"+str(i)) < 0.5:
                        left_junction_candidates.append([chromosome, str(i)])
        else:
            for i in range(start-2000000, start+1000000, 100000):
                if cov_dev_dict_100Kb.get(chromosome +":"+str(i)):
                    if cov_dev_dict_100Kb.get(chromosome +":"+str(i)) < 0.5:
                        left_junction_candidates.append([chromosome, str(i)])

        for i in range(end, end+2000000, 100000):
            if cov_dev_dict_100Kb.get(chromosome +":"+str(i)):
                if cov_dev_dict_100Kb.get(chromosome +":"+str(i)) < 0.5:
                    right_junction_candidates.append([chromosome, str(i)])

        if left_junction_candidates:
            if int(left_junction_candidates[0][1]) - 200000 < 0:
                out_file.write("\t".join(["left junction within:", left_junction_candidates[0][0], "0", str(int(left_junction_candidates[0][1]))]) + "\n")    
            else:
                out_file.write("\t".join(["left junction within:", left_junction_candidates[0][0], str(int(left_junction_candidates[0][1]) - 200000), str(int(left_junction_candidates[0][1]))]) + "\n")
            left = [left_junction_candidates[0][0], str(int(left_junction_candidates[0][1]) - 100000)]
        else:
            out_file.write("\t".join(["left junction within:", chromosome, str(start), str(start+100000)]) + "\n")
            left = [chromosome, str(start)]

        if right_junction_candidates:
            if int(right_junction_candidates[-1][1]) + 200000 > chromosome_length:
                out_file.write("\t".join(["right junction within:", right_junction_candidates[-1][0], str(int(right_junction_candidates[-1][1])), str(chromosome_length)]) + "\n")    
            else: 
                out_file.write("\t".join(["right junction within:", right_junction_candidates[-1][0], str(int(right_junction_candidates[-1][1])), str(int(right_junction_candidates[-1][1]) + 200000)]) + "\n")
            right = [chromosome, str(int(right_junction_candidates[-1][1]))]
        else:
            out_file.write("\t".join(["right junction within:", chromosome, str(end-100000), str(end)]) + "\n")
            right = [chromosome, str(end)]

        final_segments_fine.append([left[0], left[1], right[1]])
        out_file.write("\n")

windows_in_final_segments = []
for segment in final_segments_fine:
    for i in range(int(segment[1]), int(segment[2])+100000, 100000):
        windows_in_final_segments.append(segment[0]+":"+str(i))

with open(prefix + "_cov_dev_1Mb.tsv", 'r') as cov_dev_file:
    with open(prefix + "_muticum_assignments_1Mb.tsv", 'w') as out_file:
        for line in cov_dev_file:
            if not "chrUn" in line:
                temp = line.rstrip( ).split( )
                score = 0

                for i in range(int(temp[1]), int(temp[1])+1000000, 100000):
                    if temp[0]+":"+str(i) in windows_in_final_segments:
                        score += 0.1

                if score >= 0.1:
                    out_file.write(temp[0] + "\t" + temp[1] + "\t" + temp[2] + "\t" + "muticum" + "\n")
                else:
                    out_file.write(temp[0] + "\t" + temp[1] + "\t" + temp[2] + "\t" + "wheat" + "\n")
