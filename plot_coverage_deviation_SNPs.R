library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

prefix=args[1]

#####plot coverage deviation from parents
dat_cov_dev <- read.table(paste(prefix, '_cov_dev_1Mb.tsv', sep = ''))

colnames(dat_cov_dev) <- c("chr", "position", "cov_dev")
dat_cov_dev <- dat_cov_dev[dat_cov_dev$chr != "chrUn",]
dat_cov_dev$position <- dat_cov_dev$position / 1000000

cov_dev_plot <- ggplot(data = dat_cov_dev) +
  geom_point(aes(x=position, y=cov_dev), fill = "#386cb0", colour = "#386cb0", size = 0.5) +
  # scale_fill_manual(values = c("#ef3b2c", "#386cb0")) +
  # scale_colour_manual(values= c("#ef3b2c", "#386cb0")) +
  scale_x_continuous(breaks=c(0,100,200,300,400,500,600,700,800,900,1000)) +
  # geom_rect(data=centromere_df, aes(xmin=xmin,ymin=ymin,ymax=ymax,xmax=xmax),inherit.aes=FALSE) +
  facet_wrap(~ chr, ncol = 3, strip.position="left") + 
  labs(y="Mapping coverage deviation", x="Chromosomal Position (Mbp)") +
  theme_bw(base_size = 13) +
  theme(strip.background =element_rect(fill="lightgrey"), legend.position = "None") +
  coord_cartesian(ylim = c(0,2))

ggsave(paste(prefix, '_cov_dev_1Mbp.png', sep=''), cov_dev_plot, width=16, height=8, device='png')

#####plot homo and het muticum-matched SNPs
dat_homo <- read.table(paste(prefix, '_homo_mut_1Mb_windows.tsv', sep = ''))
colnames(dat_homo) <- c("chr", "position", "SNPs")
dat_homo <- dat_homo[dat_homo$chr != "chrUn",]
dat_homo$position <- dat_homo$position / 1000000

dat_het <- read.table(paste(prefix, '_het_mut_1Mb_windows.tsv', sep = ''))
colnames(dat_het) <- c("chr", "position", "SNPs")
dat_het <- dat_het[dat_het$chr != "chrUn",]
dat_het$position <- dat_het$position / 1000000

snp_df <- data.frame(dat_homo$chr, dat_homo$position, dat_homo$SNPs, dat_het$SNPs)
colnames(snp_df) <- c("chr", "position", "homo", "het")

SNP_plot <- ggplot(data = snp_df, aes(x=position)) +
  geom_col(aes(y=homo), fill = "#386cb0", alpha=0.5, width = 1) +
  geom_col(aes(y=het), fill = "#ef3b2c", alpha=0.5, width = 1) +
  # scale_fill_manual(values = c("#ef3b2c", "#386cb0")) +
  # scale_colour_manual(values= c("#ef3b2c", "#386cb0")) +
  scale_x_continuous(breaks=c(0,100,200,300,400,500,600,700,800,900,1000)) +
  # geom_rect(data=centromere_df, aes(xmin=xmin,ymin=ymin,ymax=ymax,xmax=xmax),inherit.aes=FALSE) +
  facet_wrap(~ chr, ncol = 3, strip.position="left") + 
  labs(y="Number of SNPs", x="Chromosomal Position (Mbp)") +
  theme_bw(base_size = 13) +
  theme(strip.background =element_rect(fill="lightgrey"), legend.position = "None")
# coord_cartesian(ylim = c(0,2))

ggsave(paste(prefix, '_SNP_plot_1Mbp.png', sep = ''), SNP_plot, width=16, height=8, device='png')






