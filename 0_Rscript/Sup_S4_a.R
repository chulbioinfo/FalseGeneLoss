library(dplyr)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(scales)
library(grid)
library(gtable)
library(gridExtra)
library(vroom)
library(ggnewscale)

zf_path <- "/disk1/juwan_d1/cat_rerun/zebra_finch/tm/CAT2MISSING_pre_taegut/chromosomeMap/out_circos/gene_density.200kb.rename.nuc"
an_path <- "/disk1/juwan_d1/cat_rerun/anna/tm/CAT2MISSING_pre_calann/chromosomeMap/out_circos/gene_density.200kb.rename.nuc"
pl_path <- "/disk1/juwan_d1/cat_rerun/platypus/tm/CAT2MISSING_pre_ornana/chromosomeMap/out_circos/gene_density.200kb.rename.nuc"
cp_path <- "/disk1/juwan_d1/cat_rerun/climbing_perch/tm/CAT2MISSING_pre_anates/chromosomeMap/out_circos/gene_density.200kb.rename.nuc"

zf_df <- read.csv(zf_path, sep="\t")
an_df <- read.csv(an_path, sep="\t")
pl_df <- read.csv(pl_path, sep="\t")
cp_df <- read.csv(cp_path, sep="\t")

zf_df <- zf_df %>%
  filter(grepl('^chr',X.1_usercol))

an_df <- an_df %>%
  filter(grepl('^chr',X.1_usercol))

pl_df <- pl_df %>%
  filter(grepl('^chr',X.1_usercol))

cp_df <- cp_df %>%
  filter(grepl('^chr',X.1_usercol))

zf_plot <- ggplot(zf_df) +
  facet_wrap(.~X.1_usercol, ncol=7) + 
  geom_point(aes(x=X4_usercol, y=X5_usercol), color="red", size=0.1, alpha=0.3) + # gene_density, missing ratio
  geom_point(aes(x=X4_usercol, y=X7_pct_gc*100), color="blue", size=0.1, alpha=0.3) + # gene_density, GC content
  stat_cor(data=zf_df, inherit.aes = FALSE, aes(x=X4_usercol, y=X5_usercol, group=X.1_usercol), label.x=c(15), label.y = c(95), 
           method="pearson", color="red", cex=2) + 
  stat_cor(data=zf_df, inherit.aes = FALSE, aes(x=X4_usercol, y=X7_pct_gc*100, group=X.1_usercol), label.x=c(15), label.y = c(90), 
           method="pearson", color="blue",cex=2) + 
  theme(text=element_text(size=6))
zf_plot

an_plot <- ggplot(an_df) +
  facet_wrap(.~X.1_usercol, ncol=7) + 
  geom_point(aes(x=X4_usercol, y=X5_usercol), color="red", size=0.1, alpha=0.3) + # gene_density, missing ratio
  geom_point(aes(x=X4_usercol, y=X7_pct_gc*100), color="blue", size=0.1, alpha=0.3) + # gene_density, GC content
  stat_cor(data=an_df, inherit.aes = FALSE, aes(x=X4_usercol, y=X5_usercol, group=X.1_usercol), label.x=c(15), label.y = c(95), 
           method="pearson", color="red", cex=2) + 
  stat_cor(data=an_df, inherit.aes = FALSE, aes(x=X4_usercol, y=X7_pct_gc*100, group=X.1_usercol), label.x=c(15), label.y = c(90), 
           method="pearson", color="blue",cex=2) + 
  theme(text=element_text(size=6))
an_plot

pl_plot <- ggplot(pl_df) +
  facet_wrap(.~X.1_usercol, ncol=7) + 
  geom_point(aes(x=X4_usercol, y=X5_usercol), color="red", size=0.1, alpha=0.3) + # gene_density, missing ratio
  geom_point(aes(x=X4_usercol, y=X7_pct_gc*100), color="blue", size=0.1, alpha=0.3) + # gene_density, GC content
  stat_cor(data=pl_df, inherit.aes = FALSE, aes(x=X4_usercol, y=X5_usercol, group=X.1_usercol), label.x=c(15), label.y = c(95), 
           method="spearman", color="red", cex=2) + 
  stat_cor(data=pl_df, inherit.aes = FALSE, aes(x=X4_usercol, y=X7_pct_gc*100, group=X.1_usercol), label.x=c(15), label.y = c(90), 
           method="spearman", color="blue",cex=2) + 
  theme(text=element_text(size=6))
pl_plot

cp_plot <- ggplot(cp_df) +
  facet_wrap(.~X.1_usercol, ncol=7) + 
  geom_point(aes(x=X4_usercol, y=X5_usercol), color="red", size=0.1, alpha=0.3) + # gene_density, missing ratio
  geom_point(aes(x=X4_usercol, y=X7_pct_gc*100), color="blue", size=0.1, alpha=0.3) + # gene_density, GC content
  stat_cor(data=cp_df, inherit.aes = FALSE, aes(x=X4_usercol, y=X5_usercol, group=X.1_usercol), label.x=c(15), label.y = c(95), 
           method="pearson", color="red", cex=2) + 
  stat_cor(data=cp_df, inherit.aes = FALSE, aes(x=X4_usercol, y=X7_pct_gc*100, group=X.1_usercol), label.x=c(15), label.y = c(90), 
           method="pearson", color="blue",cex=2) + 
  theme(text=element_text(size=6))
cp_plot

zf_chroms_c <- rev(c("chr.1", "chr.1A", "chr.2", "chr.3", "chr.4", "chr.4A", "chr.5", "chr.6", "chr.7", "chr.8", "chr.9", "chr.10", "chr.11", "chr.12", "chr.13", "chr.14", "chr.15", "chr.16", "chr.17", 
                 "chr.18", "chr.19", "chr.20", "chr.21", "chr.22", "chr.23", "chr.24", "chr.25", "chr.26", "chr.27", "chr.28", "chr.Z", "chr.29", "chr.30", "chr.5,u1", "chr.31", "chr.35,u1", "chr.11,u1", 
                 "chr.2,u1", "chr.26,u1", "chr.29,u1", "chr.35,u2", "chr.2,u2", "chr.1A,u1", "chr.32", "chr.16,u1", "chr.34", "chr.33", "chr.36", "chr.34,u1", "chr.35", "chr.31,u1", "chr.3,u1"))

an_chrom_c <- rev(c("chr.1", "chr.2", "chr.3", "chr.4", "chr.4A", "chr.4B", "chr.5", "chr.5A", "chr.6", "chr.7", "chr.8", "chr.9", "chr.10", "chr.11", "chr.12", "chr.13", "chr.14",
                    "chr.15", "chr.17", "chr.18", "chr.19", "chr.20", "chr.21", "chr.22", "chr.23", "chr.24", "chr.25", "chr.26", "chr.27", "chr.28", "chr.Z", "chr.33", "chr.W"))

pl_chrom_c <- rev(c("chr.1", "chr.2", "chr.3", "chr.4", "chr.5", "chr.6", "chr.7", "chr.8", "chr.9", "chr.10", "chr.11", "chr.12", "chr.13", "chr.14", "chr.15", 
                    "chr.16", "chr.17", "chr.18", "chr.19", "chr.20", "chr.21", "chr.X1", "chr.X2", "chr.X3", "chr.X4", "chr.X5", "chr.Y1", "chr.Y2", "chr.Y3", "chr.Y4", "chr.Y5"))

cp_chrom_c <- rev(c("chr.1", "chr.2", "chr.3", "chr.4", "chr.5", "chr.6", "chr.7", "chr.8", "chr.9", "chr.10", "chr.11", "chr.12", 
                    "chr.13", "chr.14", "chr.15", "chr.16", "chr.17", "chr.18", "chr.19", "chr.21", "chr.22", "chr.23", "chr.24"))

zf_df$X.1_usercol <- factor(zf_df$X.1_usercol, levels=zf_chroms_c)
an_df$X.1_usercol <- factor(an_df$X.1_usercol, levels=an_chrom_c)
pl_df$X.1_usercol <- factor(pl_df$X.1_usercol, levels=pl_chrom_c)
cp_df$X.1_usercol <- factor(cp_df$X.1_usercol, levels=cp_chrom_c)

zf_plot <- ggplot(zf_df) +
  geom_bar(aes(x=X.1_usercol, y=X4_usercol, fill=X.1_usercol), width=0.7, alpha=0.7, stat = "summary", fun = "mean") + # gene_density
  geom_hline(aes(yintercept=3.29), lwd=0.2, color="black", linetype="dashed")+
  labs(title="Zebra finch") + 
  theme(plot.title=element_text(size=8), axis.text.x=element_text(size=6, angle=0), axis.text.y=element_text(size=6, angle=0), legend.position = "none", axis.title=element_blank()) + 
  coord_flip(expand=FALSE)
zf_plot

an_plot <- ggplot(an_df) +
  geom_bar(aes(x=X.1_usercol, y=X4_usercol, fill=X.1_usercol), width=0.7, alpha=0.7, stat = "summary", fun = "mean") + # gene_density
  geom_hline(aes(yintercept=2.77), lwd=0.2, color="black", linetype="dashed")+
  labs(title="Anna's hummingbird") + 
  theme(plot.title=element_text(size=8), axis.text.x=element_text(size=6, angle=0), axis.text.y=element_text(size=6, angle=0), legend.position = "none", axis.title=element_blank()) + 
  coord_flip(expand=FALSE)
an_plot

pl_plot <- ggplot(pl_df) +
  geom_bar(aes(x=X.1_usercol, y=X4_usercol, fill=X.1_usercol), width=0.7, alpha=0.7, stat = "summary", fun = "mean") + # gene_density
  geom_hline(aes(yintercept=1.95), lwd=0.2, color="black", linetype="dashed")+
  labs(title="Platypus") + 
  theme(plot.title=element_text(size=8), axis.text.x=element_text(size=6, angle=0), axis.text.y=element_text(size=6, angle=0), legend.position = "none", axis.title=element_blank()) + 
  coord_flip(expand=FALSE)
pl_plot

cp_plot <- ggplot(cp_df) +
  geom_bar(aes(x=X.1_usercol, y=X4_usercol, fill=X.1_usercol), width=0.7, alpha=0.7, stat = "summary", fun = "mean") + # gene_density
  geom_hline(aes(yintercept=8.63), lwd=0.2, color="black", linetype="dashed")+
  labs(title="Climbing perch") + 
  theme(plot.title=element_text(size=8), axis.text.x=element_text(size=6, angle=0), axis.text.y=element_text(size=6, angle=0), legend.position = "none", axis.title=element_blank()) + 
  coord_flip(expand=FALSE)
cp_plot

my_all_plot <- c(list(zf_plot), list(an_plot), list(pl_plot), list(cp_plot))
my_layout <- rbind(c(1,2,3,4))
sup.fig.S3.plot <- arrangeGrob(grid.arrange(grobs = my_all_plot, layout_matrix = my_layout))
ggsave("/disk1/juwan_d1/cat_rerun/1_Main_figures/Sup/sup.fig.S4.a-d.png", sup.fig.S3.plot, width=6, height=4.5,units = "in")