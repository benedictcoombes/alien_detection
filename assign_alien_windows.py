import sys

prefix = sys.argv[1]
snp_threshold = int(sys.argv[2])
hom_het_ratio = int(sys.argv[3])

cov_dev_dict = {}

cov_dev_blocks = []
current_block = []

with open(prefix + "_cov_dev_1Mb.tsv", 'r') as cov_dev_file:
    for line in cov_dev_file:
        if not "chrUn" in line:
            temp = line.rstrip( ).split( )
            # print(temp)
            cov_dev_dict[temp[0]+":"+temp[1]] = float(temp[2])

            if float(temp[2]) < 0.7:
                # print(temp)
                if current_block:
                    previous_cov_dev_chromosome = current_block[-1].split(":")[0]
                    previous_cov_dev_window = current_block[-1].split(":")[1]
                    if previous_cov_dev_chromosome == temp[0] and (int(temp[1]) - int(previous_cov_dev_window)) <= 5000000:
                        current_block.append(temp[0]+":"+temp[1])
                    else:
                        cov_dev_blocks.append(current_block)
                        current_block = []
                        current_block.append(temp[0]+":"+temp[1])
                else:
                    current_block.append(temp[0]+":"+temp[1])

    cov_dev_blocks.append(current_block)

merged_cov_dev_blocks = [[i[0], i[-1]] for i in cov_dev_blocks if int(i[-1].split(":")[1]) - int(i[0].split(":")[1]) >= 1000000]

# print(cov_dev_blocks)

window_ids_within_merged_cov_dev_blocks = []

filtered_merged_cov_dev_blocks = []

current_block_windows = []

for i in merged_cov_dev_blocks:
    chromosome = i[0].split(":")[0]
    start = int(i[0].split(":")[1])
    end = int(i[1].split(":")[1])

    print(chromosome, start, end)

    count = 0
    for window in range(start, end+1000000, 1000000):
        # print(chromosome, start, end, window, cov_dev_dict.get(chromosome+":"+str(window)))
        if cov_dev_dict.get(chromosome+":"+str(window)) <= 0.7:
            count += 1

        current_block_windows.append(chromosome+":"+str(window))

    # print(chromosome, start, end, count, 0.8 * ((end+1000000 - start) / 1000000))
    if count >= (0.8 * ((end+1000000 - start) / 1000000)):
        for window in current_block_windows:
            window_ids_within_merged_cov_dev_blocks.append(window)
        filtered_merged_cov_dev_blocks.append(i)

# print(window_ids_within_merged_cov_dev_blocks)

muticum_windows = []

# for line in sys.stdin:
#     temp = line.rstrip( ).split( )

#     if int(temp[2]) >= snp_threshold:
#         if int(temp[3]) == 0:
#             if cov_dev_dict[temp[0]+":"+temp[1]] < 0.8:
#                 #print("\t".join(temp))
#                 print(temp[0] + "\t" + temp[1] + "\t" + str(cov_dev_dict[temp[0]+":"+temp[1]]) + "\t" + "muticum")
#                 muticum_windows.append(temp[0]+":"+temp[1])
#             else:
#                 print(temp[0] + "\t" + temp[1] + "\t" + str(cov_dev_dict[temp[0]+":"+temp[1]]) + "\t" + "wheat")
#         elif int(temp[2]) / int(temp[3]) >= hom_het_ratio:
#             if cov_dev_dict[temp[0]+":"+temp[1]] < 0.8:
#                 #print("\t".join(temp))     
#                 print(temp[0] +	"\t" + temp[1] + "\t" +	str(cov_dev_dict[temp[0]+":"+temp[1]]) + "\t" + "muticum")
#                 muticum_windows.append(temp[0]+":"+temp[1])
#             else:
#                 print(temp[0] + "\t" + temp[1] + "\t" + str(cov_dev_dict[temp[0]+":"+temp[1]]) + "\t" + "wheat")
#         else:
#             print(temp[0] + "\t" + temp[1] + "\t" + str(cov_dev_dict[temp[0]+":"+temp[1]]) + "\t" + "wheat")
#     else:
#         print(temp[0] +	"\t" + temp[1] + "\t" +	str(cov_dev_dict[temp[0]+":"+temp[1]]) + "\t" + "wheat")

for line in sys.stdin:
    if not "chrUn" in line:
        temp = line.rstrip( ).split( )

        if int(temp[2]) >= snp_threshold:
            if int(temp[3]) == 0:
                if temp[0]+":"+temp[1] in window_ids_within_merged_cov_dev_blocks:
                    #print("\t".join(temp))
                    # print(temp[0] + "\t" + temp[1] + "\t" + str(cov_dev_dict[temp[0]+":"+temp[1]]) + "\t" + "muticum")
                    muticum_windows.append(temp[0]+":"+temp[1])
                # else:
                #     print(temp[0] + "\t" + temp[1] + "\t" + str(cov_dev_dict[temp[0]+":"+temp[1]]) + "\t" + "wheat")
            elif int(temp[2]) / int(temp[3]) >= hom_het_ratio:
                if temp[0]+":"+temp[1] in window_ids_within_merged_cov_dev_blocks:
                    #print("\t".join(temp))     
                    # print(temp[0] + "\t" + temp[1] + "\t" + str(cov_dev_dict[temp[0]+":"+temp[1]]) + "\t" + "muticum")
                    muticum_windows.append(temp[0]+":"+temp[1])
        #         else:
        #             print(temp[0] + "\t" + temp[1] + "\t" + str(cov_dev_dict[temp[0]+":"+temp[1]]) + "\t" + "wheat")
        #     else:
        #         print(temp[0] + "\t" + temp[1] + "\t" + str(cov_dev_dict[temp[0]+":"+temp[1]]) + "\t" + "wheat")
        # else:
        #     print(temp[0] + "\t" + temp[1] + "\t" + str(cov_dev_dict[temp[0]+":"+temp[1]]) + "\t" + "wheat")

# for i in muticum_windows:
#     print(i,cov_dev_dict.get(i))



# for i in merged_cov_dev_blocks:
#     print(i)

# for i in filtered_merged_cov_dev_blocks:
#     print(i)



# for i in window_ids_within_merged_cov_dev_blocks:
#     print(i)

final_segments_coarse = []
# print("FINAL SEGMENTS COARSE")
count = 0
for i in filtered_merged_cov_dev_blocks:
    chromosome = i[0].split(":")[0]
    start = int(i[0].split(":")[1])
    end = int(i[1].split(":")[1])
    running_mean = []
    for window in range(start, end+1000000, 1000000):
        if chromosome+":"+str(window) in muticum_windows:
            count += 1
        # print(window, count)
        running_mean.append([window, count])
    # print(chromosome, start, end, count)
    if (end - start) / 1000000 <= 5:
        # print(i, count)
        if count > (0.50 * ((end - start) / 1000000)):
            final_segments_coarse.append([chromosome, start, end])
    elif count > (0.14 * ((end - start) / 1000000)):
        # print(i[0], i[-1])
        final_segments_coarse.append([chromosome, start, end])
    else:
        last_increase = ""
        block_from_snps_coarse = []
        for n, i in enumerate(running_mean):
            # print(i)
            if i[0] - running_mean[0][0] == 0:
                if i[1] == 1:
                    # print("yes")
                    last_increase = i[0]
            else:
                if i[1] - running_mean[n-1][1] == 1:
                    last_increase = i[0]

                # print(i[1] / ((i[0] - running_mean[0][0]) / 1000000))
                if i[1] / ((i[0] - running_mean[0][0]) / 1000000) >= 0.14:
                    # print("yes")
                    block_from_snps_coarse.append(i)
                else:
                    if last_increase:
                        print([chromosome, start, last_increase])
                        final_segments_coarse.append([chromosome, start, last_increase])
                    break

        ###same but from the end of chromosome walking backwards
        running_mean = []
        for window in range(start, end+1000000, 1000000):
            if chromosome+":"+str(window) in muticum_windows:
                count += 1
            # print(window, count)
        running_mean.append([window, count])
        last_increase = ""
        block_from_snps_coarse = []
        for n, i in enumerate(running_mean[::-1]):
            print(i)
            if running_mean[-1][0] - i[0] == 0:
                if i[1] == 1:
                    # print("yes")
                    last_increase = i[0]
            else:
                if i[1] - running_mean[n+1][1] == 1:
                    last_increase = i[0]

                # print(i[1] / ((i[0] - running_mean[0][0]) / 1000000))
                if i[1] / ((running_mean[-1][0] - i[0]) / 1000000) >= 0.14:
                    # print("yes")
                    block_from_snps_coarse.append(i)
                else:
                    if last_increase:
                        print([chromosome, start, last_increase])
                        final_segments_coarse.append([chromosome, start, last_increase])
                    break

    count = 0

# print( )


final_segments_fine = []

###make cov dev dict with 100Kbp windows

cov_dev_dict_100Kb = {}

cov_dev_blocks = []
current_block = []

with open(prefix + "_cov_dev_100Kb.tsv", 'r') as cov_dev_file:
    for line in cov_dev_file:
        if not "chrUn" in line:
            temp = line.rstrip( ).split( )
            # print(temp)
            cov_dev_dict_100Kb[temp[0]+":"+temp[1]] = float(temp[2])

with open(prefix + "_segments.txt", 'w') as out_file:
    for segment in final_segments_coarse:
        chromosome = segment[0]
        start = int(segment[1])
        end = int(segment[2])

        # print(chromosome, start, end)
        out_file.write("\t".join(["Coarse segment identification:", chromosome, str(start), str(end)]) + "\n")
        left_junction_candidates = []
        right_junction_candidates = []

        # print("LEFT JUNCTION")

        if start == 0:
            for i in range(start, start+1000000, 100000):
            # print(i)
                if cov_dev_dict_100Kb.get(chromosome +":"+str(i)):
                    if cov_dev_dict_100Kb.get(chromosome +":"+str(i)) < 0.5:
                        # print(i)
                        left_junction_candidates.append([chromosome, str(i)])
                        # print(chromosome, str(i), cov_dev_dict_100Kb.get(chromosome +":"+str(i)))
        else:
            for i in range(start-2000000, start+1000000, 100000):
                # print(i)
                if cov_dev_dict_100Kb.get(chromosome +":"+str(i)):
                    if cov_dev_dict_100Kb.get(chromosome +":"+str(i)) < 0.5:
                        # print(i)
                        left_junction_candidates.append([chromosome, str(i)])
                        # print(chromosome, str(i), cov_dev_dict_100Kb.get(chromosome +":"+str(i)))

        # if end == 

        # print("RIGHT JUNCTION")
        for i in range(end, end+2000000, 100000):
            if cov_dev_dict_100Kb.get(chromosome +":"+str(i)):
                if cov_dev_dict_100Kb.get(chromosome +":"+str(i)) < 0.5:
                    right_junction_candidates.append([chromosome, str(i)])
                    # print(chromosome, str(i), cov_dev_dict_100Kb.get(chromosome +":"+str(i)))

        # print(left_junction_candidates)
        # print(right_junction_candidates)

        if left_junction_candidates:
            # print("left junction within: ", "\t".join([left_junction_candidates[0][0], str(int(left_junction_candidates[0][1]) - 100000), str(int(left_junction_candidates[0][1]))]))
            out_file.write("\t".join(["left junction within:", left_junction_candidates[0][0], str(int(left_junction_candidates[0][1]) - 200000), str(int(left_junction_candidates[0][1]))]) + "\n")
            left = [left_junction_candidates[0][0], str(int(left_junction_candidates[0][1]) - 100000)]
        else:
            # print("left junction within: ", "\t".join([chromosome, str(start+100000)]))
            out_file.write("\t".join(["left junction within:", chromosome, str(start), str(start+100000)]) + "\n")
            left = [chromosome, str(start)]

        if right_junction_candidates:
            # print("right junction within: ", "\t".join([right_junction_candidates[-1][0], str(int(right_junction_candidates[-1][1])), str(int(right_junction_candidates[-1][1]) + 100000)]))
            out_file.write("\t".join(["right junction within:", right_junction_candidates[-1][0], str(int(right_junction_candidates[-1][1])), str(int(right_junction_candidates[-1][1]) + 200000)]) + "\n")
            right = [chromosome, str(int(right_junction_candidates[-1][1]))]
        else:
            # print("right junction within: ", "\t".join([chromosome, str(end-100000)]))
            out_file.write("\t".join(["right junction within:", chromosome, str(end-100000), str(end)]) + "\n")
            right = [chromosome, str(end)]

        final_segments_fine.append([left[0], left[1], right[1]])
        # print(left[0], left[1], right[1])
        # print( )
        out_file.write("\n")

windows_in_final_segments = []
for segment in final_segments_fine:
    for i in range(int(segment[1]), int(segment[2])+100000, 100000):
        windows_in_final_segments.append(segment[0]+":"+str(i))




# with open(prefix + "_cov_dev_100Kb.tsv", 'r') as cov_dev_file:
#     with open(prefix + "_muticum_assignments_100Kb.tsv", 'w') as out_file:
#         for line in cov_dev_file:
#             if not "chrUn" in line:
#                 temp = line.rstrip( ).split( )
#                 if temp[0]+":"+temp[1] in windows_in_final_segments:
#                     out_file.write(temp[0] + "\t" + temp[1] + "\t" + temp[2] + "\t" + "muticum" + "\n")
#                 else:
#                     out_file.write(temp[0] + "\t" + temp[1] + "\t" + temp[2] + "\t" + "wheat" + "\n")


with open(prefix + "_cov_dev_1Mb.tsv", 'r') as cov_dev_file:
    with open(prefix + "_muticum_assignments_1Mb.tsv", 'w') as out_file:
        for line in cov_dev_file:
            if not "chrUn" in line:
                temp = line.rstrip( ).split( )
                # print(temp)
                score = 0

                for i in range(int(temp[1]), int(temp[1])+1000000, 100000):
                    if temp[0]+":"+str(i) in windows_in_final_segments:
                        score += 0.1

                if score >= 0.2:
                    out_file.write(temp[0] + "\t" + temp[1] + "\t" + temp[2] + "\t" + "muticum" + "\n")
                else:
                    out_file.write(temp[0] + "\t" + temp[1] + "\t" + temp[2] + "\t" + "wheat" + "\n")

                # if temp[0]+":"+temp[1] in windows_in_final_segments:
                #     out_file.write(temp[0] + "\t" + temp[1] + "\t" + temp[2] + "\t" + "muticum" + "\n")
                # else:
                #     out_file.write(temp[0] + "\t" + temp[1] + "\t" + temp[2] + "\t" + "wheat" + "\n")





