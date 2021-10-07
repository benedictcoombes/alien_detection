import sys

prefix=sys.argv[2]
window_size=sys.argv[3]

with open('/ei/projects/c/c38b7c13-1c82-4b59-ae07-b631556928a3/scratch/introgression_work/introgression_detection_25_lines/introgression_lines/muticum_lines/parent_coverage/paragon/paragon_cov_norm_' + window_size + '.tsv', 'r') as paragon_cov_file:
	with open('/ei/projects/c/c38b7c13-1c82-4b59-ae07-b631556928a3/scratch/introgression_work/introgression_detection_25_lines/introgression_lines/muticum_lines/parent_coverage/pavon/pavon_cov_norm_' + window_size + '.tsv', 'r') as pavon_cov_file:
	# with open('/ei/projects/c/c38b7c13-1c82-4b59-ae07-b631556928a3/scratch/introgression_work/introgression_detection_25_lines/introgression_lines/muticum_lines/parent_coverage/CS/CS_cov_norm_' + window_size + '.tsv', 'r') as pavon_cov_file:
		with open(sys.argv[1], 'r') as int_line_cov_file:
			paragon = [[line.split( )[0], line.split( )[1], line.split( )[4]] for line in paragon_cov_file]
			pavon = [[line.split( )[0], line.split( )[1], line.split( )[4]] for line in pavon_cov_file]
			int_line = [[line.split( )[0], line.split( )[1], line.split( )[3]] for line in int_line_cov_file]


with open(prefix + "_cov_dev_" + window_size + ".tsv", 'w') as out_file:
	for n, i in enumerate(int_line):
		cov = float(i[2]) * 100
		paragon_cov = float(paragon[n][2])
		pavon_cov = float(pavon[n][2])
		if paragon_cov == 0:
			vs_paragon = 0
		else:
			vs_paragon = cov / paragon_cov
		if pavon_cov == 0:
			vs_pavon = 0
		else:
			vs_pavon = cov / pavon_cov

		if abs(1-vs_paragon) <= abs(1-vs_pavon):
			out_file.write("\t".join([str(i[0]), str(i[1]), str(vs_paragon)]) + "\n")

		else:
			out_file.write("\t".join([str(i[0]), str(i[1]), str(vs_pavon)]) + "\n")
