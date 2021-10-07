import sys
from collections import OrderedDict

alien = sys.argv[4]

with open(sys.argv[1], 'r') as muticum_vcf:
    with open(sys.argv[2], 'r') as paragon_vcf:
        with open(sys.argv[3], 'r') as pavon_vcf:

            paragon_dict = OrderedDict()
            muticum_dict = OrderedDict()
            pavon_dict = OrderedDict()

            for line in paragon_vcf:
                if not line.startswith("#"):
                    temp = line.split( )
                    if "," in temp[4]:
                        alleles = temp[4].split(",")
                        if len(alleles) == 2:
                            paragon_dict[temp[0]+":"+temp[1]] = [alleles[0], alleles[1], "Het"]
                    else:
                        if "Homo" in line:
                            paragon_dict[temp[0]+":"+temp[1]] = [temp[4], "", "Homo"]
                        else:    
                            paragon_dict[temp[0]+":"+temp[1]] = [temp[4], "", "Het"]

            for line in pavon_vcf:
                if not line.startswith("#"):
                    temp = line.split( )
                    if "," in temp[4]:
                        alleles = temp[4].split(",")
                        if len(alleles) == 2:
                            pavon_dict[temp[0]+":"+temp[1]] = [alleles[0], alleles[1], "Het"]
                    else:
                        if "Homo" in line:
                            pavon_dict[temp[0]+":"+temp[1]] = [temp[4], "", "Homo"]
                        else:    
                            pavon_dict[temp[0]+":"+temp[1]] = [temp[4], "", "Het"]

            for line in muticum_vcf:
                if not line.startswith("#"):
                    temp = line.split( )
                    if "," in temp[4]:
                        alleles_tmp = temp[4]+",Het"
                        alleles = alleles_tmp.split(",")
                        muticum_dict[temp[0]+":"+temp[1]] = alleles
                    else:
                        if "Homo" in line:
                            muticum_dict[temp[0]+":"+temp[1]] = [temp[4], "", "Homo"]
                        else:    
                            muticum_dict[temp[0]+":"+temp[1]] = [temp[4], "", "Het"]

print("vcfs loaded") 

###get muticum specific SNPs
with open(alien + "_specific_snps.vcf", 'w') as muticum_out:
    for x in muticum_dict.keys():
        chromosome = x.split(":")[0]
        base = x.split(":")[1]
        alleles = muticum_dict.get(x)[0:2]
        zygosity = muticum_dict.get(x)[2]

        if paragon_dict.get(x):
            paragon_alleles = paragon_dict.get(x)
        else:
            paragon_alleles = ["", ""]
        if pavon_dict.get(x):
            pavon_alleles = pavon_dict.get(x)
        else:
            pavon_alleles = ["", ""]

        if paragon_alleles[0] and not pavon_alleles[0]:
        ##SNP at site in paragon but not in pavon
            if alleles[1]:
            #muticum has 2 alt alleles
                if paragon_alleles[1]:
                    #paragon has 2 alt alleles
                    if alleles[0] != paragon_alleles[0] and alleles[0] != paragon_alleles[1] and alleles[1] != paragon_alleles[0] and alleles[1] != paragon_alleles[1]:
                        muticum_out.write("\t".join([chromosome, base, ".", ".", alleles[0]+","+alleles[1], zygosity]) + "\n")
                    if alleles[0] != paragon_alleles[0] and alleles[0] != paragon_alleles[1] and (alleles[1] == paragon_alleles[0] or alleles[1] == paragon_alleles[1]):
                        muticum_out.write("\t".join([chromosome, base, ".", ".", alleles[0], zygosity]) + "\n")
                    if (alleles[0] == paragon_alleles[0] or alleles[0] == paragon_alleles[1]) and alleles[1] != paragon_alleles[0] and alleles[1] == paragon_alleles[1]:
                        muticum_out.write("\t".join([chromosome, base, ".", ".", alleles[1], zygosity]) + "\n")
                else:
                    #paragon has 1 alt allele
                    if alleles[0] != paragon_alleles[0] and alleles[1] != paragon_alleles[0]:
                        muticum_out.write("\t".join([chromosome, base, ".", ".", alleles[0]+","+alleles[1], zygosity]) + "\n")
                    elif alleles[0] != paragon_alleles[0] and alleles[1] == paragon_alleles[0]:
                        muticum_out.write("\t".join([chromosome, base, ".", ".", alleles[0], zygosity]) + "\n")
                    elif alleles[0] == paragon_alleles[0] and alleles[1] != paragon_alleles[0]:
                        muticum_out.write("\t".join([chromosome, base, ".", ".", alleles[1], zygosity]) + "\n")
            else:
                #muticum has 1 alt allele
                if paragon_alleles[1]:
                    #paragon has 2 alt alleles
                    if alleles[0] != paragon_alleles[0] and alleles[0] != paragon_alleles[1]:
                        muticum_out.write("\t".join([chromosome, base, ".", ".", alleles[0], zygosity]) + "\n")



        elif not paragon_alleles[0] and pavon_alleles[0]:
            ##SNP at site in pavon but not in paragon
            if alleles[1]:
            #muticum has 2 alt alleles
                if pavon_alleles[1]:
                    #pavon has 2 alt alleles
                    if alleles[0] != pavon_alleles[0] and alleles[0] != pavon_alleles[1] and alleles[1] != pavon_alleles[0] and alleles[1] != pavon_alleles[1]:
                        muticum_out.write("\t".join([chromosome, base, ".", ".", alleles[0]+","+alleles[1], zygosity]) + "\n")
                    if alleles[0] != pavon_alleles[0] and alleles[0] != pavon_alleles[1] and (alleles[1] == pavon_alleles[0] or alleles[1] == pavon_alleles[1]):
                        muticum_out.write("\t".join([chromosome, base, ".", ".", alleles[0], zygosity]) + "\n")
                    if (alleles[0] == pavon_alleles[0] or alleles[0] == pavon_alleles[1]) and alleles[1] != pavon_alleles[0] and alleles[1] == pavon_alleles[1]:
                        muticum_out.write("\t".join([chromosome, base, ".", ".", alleles[1], zygosity]) + "\n")
                else:
                    #pavon has 1 alt allele
                    if alleles[0] != pavon_alleles[0] and alleles[1] != pavon_alleles[0]:
                        muticum_out.write("\t".join([chromosome, base, ".", ".", alleles[0]+","+alleles[1], zygosity]) + "\n")
                    elif alleles[0] != pavon_alleles[0] and alleles[1] == pavon_alleles[0]:
                        muticum_out.write("\t".join([chromosome, base, ".", ".", alleles[0], zygosity]) + "\n")
                    elif alleles[0] == pavon_alleles[0] and alleles[1] != pavon_alleles[0]:
                        muticum_out.write("\t".join([chromosome, base, ".", ".", alleles[1], zygosity]) + "\n")
            else:
                #muticum has 1 alt allele
                if pavon_alleles[1]:
                    #pavon has 2 alt alleles
                    if alleles[0] != pavon_alleles[0] and alleles[0] != pavon_alleles[1]:
                        muticum_out.write("\t".join([chromosome, base, ".", ".", alleles[0], zygosity]) + "\n")

        elif paragon_alleles[0] and pavon_alleles[0]:
            ##SNP at site in both paragon and in pavon
            if alleles[1]:
            #muticum has 2 alt alleles
                if not alleles[0] in pavon_alleles and not alleles[0] in paragon_alleles:
                    if alleles[1] in pavon_alleles or alleles[1] in paragon_alleles:
                        muticum_out.write("\t".join([chromosome, base, ".", ".", alleles[0], zygosity]) + "\n")
                    else:
                        muticum_out.write("\t".join([chromosome, base, ".", ".", alleles[0]+","+alleles[1], zygosity]) + "\n")
                elif alleles[1] in pavon_alleles and not alleles[1] in paragon_alleles:
                    muticum_out.write("\t".join([chromosome, base, ".", ".", alleles[1], zygosity]) + "\n")

            else:
                #muticum has 1 alt allele
                if not alleles[0] in pavon_alleles and not alleles[0] in paragon_alleles:
                    muticum_out.write("\t".join([chromosome, base, ".", ".", alleles[0], zygosity]) + "\n")

        else:
            if alleles[1]:
                muticum_out.write("\t".join([chromosome, base, ".", ".", alleles[0]+","+alleles[1], zygosity]) + "\n")
            else:
                muticum_out.write("\t".join([chromosome, base, ".", ".", alleles[0], zygosity]) + "\n")
