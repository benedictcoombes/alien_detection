import sys

window_size=sys.argv[4]
prefix=sys.argv[5]

with open(sys.argv[1], 'r') as int_line_cov_file:
	with open(sys.argv[2], 'r') as wheat_parent_1_cov_file:
		with open(sys.argv[3], 'r') as wheat_parent_2_cov_file:
			wheat_parent_1 = [[line.split( )[0], line.split( )[1], line.split( )[4]] for line in wheat_parent_1_cov_file]
			wheat_parent_2 = [[line.split( )[0], line.split( )[1], line.split( )[4]] for line in wheat_parent_2_cov_file]
			int_line = [[line.split( )[0], line.split( )[1], line.split( )[3]] for line in int_line_cov_file]


with open(prefix + "_cov_dev_" + window_size + ".tsv", 'w') as out_file:
	for n, i in enumerate(int_line):
		cov = float(i[2]) * 100
		wheat_parent_1_cov = float(wheat_parent_1[n][2])
		wheat_parent_2_cov = float(wheat_parent_2[n][2])
		if wheat_parent_1_cov == 0:
			vs_wheat_parent_1 = 0
		else:
			vs_wheat_parent_1 = cov / wheat_parent_1_cov
		if wheat_parent_2_cov == 0:
			vs_wheat_parent_2 = 0
		else:
			vs_wheat_parent_2 = cov / wheat_parent_2_cov

		if abs(1-vs_wheat_parent_1) <= abs(1-vs_wheat_parent_2):
			out_file.write("\t".join([str(i[0]), str(i[1]), str(vs_wheat_parent_1)]) + "\n")

		else:
			out_file.write("\t".join([str(i[0]), str(i[1]), str(vs_wheat_parent_2)]) + "\n")
