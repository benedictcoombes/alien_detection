alien_species=$1
parent_1=$2
parent_2=$3
int_line_name=$4

bam_suffix=$5
int_line_vcf_suffix=$6
parent_lines_vcf_suffix=$7

#Index the wheat reference genome
samtools faidx 161010_Chinese_Spring_v1.0_pseudomolecules.fasta

###Produce genome file for the reference genome to be used as input to bedtools
awk '{OFS="\t"; print $1,$2}' 161010_Chinese_Spring_v1.0_pseudomolecules.fasta.fai > RefSeqv1.0_genome_file.txt

###Produce 1Mbp and 100Kbp window bed files
bedtools makewindows -w 1000000 -g RefSeqv1.0_genome_file.txt > RefSeqv1.0_1Mb_windows.bed 
bedtools makewindows -w 100000 -g RefSeqv1.0_genome_file.txt > RefSeqv1.0_100Kb_windows.bed

###Count number of reads in each 1Mbp window for the introgression line and the parent lines
hts_nim_tools count-reads RefSeqv1.0_1Mb_windows.bed ${int_line_name}${bam_suffix} | sort -k1,1 -k2,2n | awk '{OFS="\t"; print}' > ${int_line_name}_cov_1Mb_windows.tsv 
hts_nim_tools count-reads RefSeqv1.0_1Mb_windows.bed ${parent_1}${bam_suffix} | sort -k1,1 -k2,2n | awk '{OFS="\t"; print}' > ${parent_1}_cov_1Mb_windows.tsv 
hts_nim_tools count-reads RefSeqv1.0_1Mb_windows.bed ${parent_2}${bam_suffix} | sort -k1,1 -k2,2n | awk '{OFS="\t"; print}' > ${parent_2}_cov_1Mb_windows.tsv

###Count number of reads in each 100Kbp window for the introgression line and the parent lines
hts_nim_tools count-reads RefSeqv1.0_100Kb_windows.bed ${int_line_name}${bam_suffix} | sort -k1,1 -k2,2n | awk '{OFS="\t"; print}' > ${int_line_name}_cov_100Kb_windows.tsv 
hts_nim_tools count-reads RefSeqv1.0_100Kb_windows.bed ${parent_1}${bam_suffix} | sort -k1,1 -k2,2n | awk '{OFS="\t"; print}' > ${parent_1}_cov_100Kb_windows.tsv 
hts_nim_tools count-reads RefSeqv1.0_100Kb_windows.bed ${parent_2}${bam_suffix} | sort -k1,1 -k2,2n | awk '{OFS="\t"; print}' > ${parent_2}_cov_100Kb_windows.tsv

###Computing coverage deviation in each genomic window between the introgression line and the parent lines
python3 cov_deviation.py ${int_line_name}_cov_1Mb_windows.tsv ${parent_1}_cov_1Mb_windows.tsv ${parent_2}_cov_1Mb_windows.tsv 1Mb ${int_line_name}
python3 cov_deviation.py ${int_line_name}_cov_100Kb_windows.tsv ${parent_1}_cov_100Kb_windows.tsv ${parent_2}_cov_100Kb_windows.tsv 100Kb ${int_line_name}

###Producing species_specific_SNPs. This takes the vcfs of the introgressed species and the 2 wheat parents and returns 
###a file of species-specfic SNPs for the intrgoressed species
python3 alien_specific_snps.py ${int_line_name}${int_line_vcf_suffix} ${parent_1}${parent_lines_vcf_suffix} ${parent_2}${parent_lines_vcf_suffix} ${alien_species}

###matching alien specific SNPs with introgression line SNPs
python3 match_alien_specific_snps.py muticum_specific_snps.tsv ${int_line_name}${int_line_vcf_suffix} ${int_line_name} ${alien_species} 5 3
bedtools coverage -a RefSeqv1.0_1Mb_windows.bed -b ${int_line_name}_${alien_species}_specific_SNP_assignments_homo.bed | cut -f 1,2,4 > ${int_line_name}_homo_mut_1Mb_windows.tsv 
bedtools coverage -a RefSeqv1.0_1Mb_windows.bed -b ${int_line_name}_${alien_species}_specific_SNP_assignments_het.bed | cut -f 1,2,4 > ${int_line_name}_het_mut_1Mb_windows.tsv

###Produce coverage deviation plot and homozygous/heterozygous SNP density plot assigned introgression blocks plot. 
###The coverage deviation and SNP plots are usually sufficient to spot where the introgressions are.
###You are looking for a block with coverage deviation below 1 (depending on the genetic distance of the donor chromosome from the wheat chromosome
###at the introgression location, this will vary. Introgressions from non-homologous chromosomes will usually have deviation below 0.5. Primary genepool
###introgression, such as tauschii, urartu, or durum, if inserted within their homologous subgenome, will have a much more subtle drop in coverage).
Rscript plot_coverage_deviation_SNPs.R ${int_line_name}

###Assigning introgression blocks based on coverage deviation blocks that overlap with a region of high homozygous species-specific SNP density
###and low density of heterozygous species-specific SNP density
###The assigned blocks script is designed to automate what can usually be deduced by eye from the coverage and SNP plots. This works most of the time
###but the parameters sometimes need tweaking based on a number of variables, including the number of species-specific SNPs, the SNP density of the introgressed 
###region, and how contiguous the block is
paste ${int_line_name}_homo_mut_1Mb_windows.tsv ${int_line_name}_het_mut_1Mb_windows.tsv | cut -f 1,2,3,6 | python3 assign_alien_windows.py ${int_line_name} 55 4 0.7 5000000 0.14 0.6 RefSeqv1.0.fasta.fai

###Plotting coverage deviation colour coded based on assigned introgression blocks
Rscript plot_assigned_introgression_blocks.R ${int_line_name}