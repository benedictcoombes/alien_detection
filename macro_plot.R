prefix="DH65"

dat <- read.table(paste(prefix, '_muticum_assignments_1Mb.tsv',sep = ''))
colnames(dat) <- c("chr", "position", "cov", "genotype")
dat$position <- dat$position / 1000000

plot <- ggplot(data = dat, aes(x=position, y=cov/median(cov), fill = genotype, colour = genotype)) +
  geom_point(size=0.5) +
  scale_fill_manual(values = c("#ef3b2c", "#386cb0")) +
  scale_colour_manual(values= c("#ef3b2c", "#386cb0")) +
  scale_x_continuous(breaks=c(0,100,200,300,400,500,600,700,800,900,1000)) +
  geom_rect(data=centromere_df, aes(xmin=xmin,ymin=ymin,ymax=ymax,xmax=xmax),inherit.aes=FALSE) +
  facet_wrap(. ~ chr, ncol = 3, strip.position="left") + 
  labs(y="Mapping coverage deviation", x="Chromosomal Position (Mbp)", colour ="Genotype", fill = "Genotype") +
  theme_bw(base_size = 13) +
  theme(strip.background =element_rect(fill="lightgrey"), legend.position = "None") +
  coord_cartesian(ylim = c(0,2))

ggsave(paste(prefix, '_macro_structure_1Mbp.png', sep=''), plot, width=16, height=8, device='png')
ggsave(paste(prefix, '_macro_structure_1Mbp.svg', sep=''), plot, width=16, height=8, device='svg')
