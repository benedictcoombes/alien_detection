# alien_detection

These scripts were used during the analysis for "Whole genome sequencing reveals the structural and functional impact of synthetic alien introgression in wheat" to identify _Am. muticum_ introgression regions within a hexaploid wheat background.

Input for these scripts are the following:

BAM alignment files and vcfs produced from mapping WGS illumina reads from the wheat parents, _Am. muticum_, and the introgression lines to the wheat reference genome (we used RefSeq v1.0). These are produced and processed using details included in the methods section of the paper.

samtools faidx 161010_Chinese_Spring_v1.0_pseudomolecules.fasta
awk '{OFS="\t"; print $1,$2}' 161010_Chinese_Spring_v1.0_pseudomolecules.fasta.fai > RefSeqv1.0_genome_file.txt
bedtools makewindows -w 1000000 -g RefSeqv1.0_genome_file.txt > RefSeqv1.0_1Mb_windows.bed
bedtools makewindows -w 100000 -g RefSeqv1.0_genome_file.txt > RefSeqv1.0_100Kb_windows.bed

num_reads_post_duplicate_removal = n

prefix = line_name

hts_nim_tools count-reads RefSeqv1.0_1Mb_windows.bed $bam | sort -k1,1 -k2,2n | awk -v num_reads=$num_reads_post_duplicate_removal '{OFS="\t"; $5=$4/num_reads; print}' > ${prefix}_cov_1Mb_windows.tsv

python3 cov_deviation.py int_line_cov_file wheat_parent_1_cov_file wheat_parent_2_cov_file window_size introgression_line_name


###making alien_specific_SNPs
python3 alien_specific_snps.py muticum.vcf paragon.vcf pavon.vcf muticum


###matching alien specific SNPs with introgression line SNPs
python3 match_alien_specific_snps.py muticum_specific_snps.tsv ${prefix}_vfallelecalls_dp5.vcf ${prefix} muticum 5 3

bedtools coverage -a RefSeqv1.0_1Mb_windows.bed -b ${prefix}_muticum_specific_SNP_assignments_homo.bed | cut -f 1,2,4 > ${prefix}_homo_mut_1Mb_windows.tsv
bedtools coverage -a RefSeqv1.0_1Mb_windows.bed -b ${prefix}_muticum_specific_SNP_assignments_het.bed | cut -f 1,2,4 > ${prefix}_het_mut_1Mb_windows.tsv

###making introgression blocks
paste ${prefix}_homo_mut_1Mb_windows.tsv ${prefix}_het_mut_1Mb_windows.tsv | cut -f 1,2,3,6 | python3 assign_alien_windows.py ${prefix} 55 4





