import sys
from collections import defaultdict

line_prefix = sys.argv[3]
alien = sys.argv[4]
homo_threshold = int(sys.argv[5])
het_threshold = int(sys.argv[6])

with open(sys.argv[1], 'r') as alien_specific_snp_file:
    alien_dict = defaultdict()
    for line in alien_specific_snp_file:
        if not line.startswith('#'):
            temp = line.split( )
            alien_dict[temp[0]+':'+temp[1]] = temp[2]
    print('Alien-specific SNPs loaded')

previous_chr = "chr1A"
with open(sys.argv[2], 'r') as introgression_line_vcf:
    with open(line_prefix + "_" + alien + "_specific_SNP_assignments_homo.bed", 'w') as homo_out:
        with open(line_prefix + "_" + alien + "_specific_SNP_assignments_het.bed", 'w') as het_out:
            for line in introgression_line_vcf:
                if line and not line.startswith("#") and (line.startswith("chr") or line.startswith("chr")) and "PASS" in line:
                    temp = line.split( )
                    if temp[0] != previous_chr:
                        print("Processed " + previous_chr)
                        previous_chr = temp[0]
                    if "BQB" in line:
                        DP = temp[7].split(";")[3]
                    else:
                        DP = temp[7].split(";")[2]

                    if "BQB" in line:
                        DP4 = temp[7].split(";")[4].split("=")[1].split(",")
                    else:
                        DP4 = temp[7].split(";")[3].split("=")[1].split(",")

                    alt_allele_count = int(DP4[2]) + int(DP4[3])
                    total_allele_count = int(DP4[0]) + int(DP4[1]) + int(DP4[2]) + int(DP4[3])

                    ref_allele = temp[3]
                    alt_allele = temp[4]

                    x = temp[0]+":"+temp[1]

                    if alien_dict.get(x):
                        if "," in alien_dict.get(x):
                            alien_alleles = alien_dict.get(x).split(",")
                        else:
                            alien_alleles = [alien_dict.get(x)]
                    else:
                        alien_alleles = ""

                    if alien_alleles:
                        ###site in int line vcf present in one or more of the genome specific vcfs
                        if "Homo" in line and alt_allele_count == total_allele_count and alt_allele_count >= homo_threshold and int(DP.split("=")[1]) <= (homo_threshold * 4):
                            if len(alien_alleles) == 1:
                                if temp[4] == alien_alleles[0]:
                                    out = [temp[0], temp[1], str(int(temp[1]) + 1)]
                                    homo_out.write("\t".join(out) + "\n")
                            elif len(alien_alleles) == 2:
                                if temp[4] == alien_alleles[0] or temp[4] == alien_alleles[1]:
                                    out = [temp[0], temp[1], str(int(temp[1]) + 1)]
                                    homo_out.write("\t".join(out) + "\n")

                        elif "Het" in line and int(DP.split("=")[1]) <= (homo_threshold * 4):
                            if "," in temp[4]:
                                if len(temp[4].split(",")) == 2 and total_allele_count >= (het_threshold * 2):
                                    het_alleles = temp[4].split(",")
                                    if len(alien_alleles) == 1:
                                        if het_alleles[0] == alien_alleles[0] or het_alleles[1] == alien_alleles[0]:
                                            out = [temp[0], temp[1], str(int(temp[1]) + 1)]
                                            het_out.write("\t".join(out) + "\n")
                                    elif len(alien_alleles) == 2:
                                        if het_alleles[0] == alien_alleles[0] or het_alleles[0] == alien_alleles[1] or het_alleles[1] == alien_alleles[0] or het_alleles[1] == alien_alleles[1]:
                                            out = [temp[0], temp[1], str(int(temp[1]) + 1)]
                                            het_out.write("\t".join(out) + "\n")
                            else:
                                het_allele = temp[4]
                                if len(alien_alleles) == 1:
                                    if het_allele == alien_alleles[0] and alt_allele_count >= het_threshold and total_allele_count - alt_allele_count >= het_threshold:
                                        out = [temp[0], temp[1], str(int(temp[1]) + 1)]
                                        het_out.write("\t".join(out) + "\n")
                                elif len(alien_alleles) == 2:
                                    if het_allele == alien_alleles[0] or het_allele == alien_alleles[1] and alt_allele_count >= het_threshold and total_allele_count - alt_allele_count >= het_threshold:
                                        out = [temp[0], temp[1], str(int(temp[1]) + 1)]
                                        het_out.write("\t".join(out) + "\n")
