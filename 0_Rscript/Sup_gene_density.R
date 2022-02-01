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

zf_path <- "/disk1/juwan_d1/cat_rerun/zebra_finch/tm/CAT2MISSING_pre_taegut/chromosomeMap/out_correlation/vgp_taegut.stat_summary.tabs"
an_path <- "/disk1/juwan_d1/cat_rerun/anna/tm/CAT2MISSING_pre_calann/chromosomeMap/out_correlation/vgp_calann.stat_summary.tabs"
pl_path <- "/disk1/juwan_d1/cat_rerun/platypus/tm/CAT2MISSING_pre_ornana/chromosomeMap/out_correlation/vgp_ornana.stat_summary.tabs"
cp_path <- "/disk1/juwan_d1/cat_rerun/climbing_perch/tm/CAT2MISSING_pre_anates/chromosomeMap/out_correlation/vgp_anates.stat_summary.tabs"

zf_df <- read.csv(zf_path, sep="\t")
an_df <- read.csv(an_path, sep="\t")
pl_df <- read.csv(pl_path, sep="\t")
cp_df <- read.csv(cp_path, sep="\t")

zf_df <- zf_df %>%
  filter(grepl('^chr',scaffold))
an_df <- an_df %>%
  filter(grepl('^chr',scaffold))
pl_df <- pl_df %>%
  filter(grepl('^chr',scaffold))
cp_df <- cp_df %>%
  filter(grepl('^chr',scaffold))

zf_df$scaffold <- gsub("chr.", "", zf_df$scaffold)
an_df$scaffold <- gsub("chr.", "", an_df$scaffold)
pl_df$scaffold <- gsub("chr.", "", pl_df$scaffold)
cp_df$scaffold <- gsub("chr.", "", cp_df$scaffold)

zf_chroms_c <- rev(c("1", "1A", "2", "3", "4", "4A", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", 
                 "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "Z", "29", "30", "5,u1", "31", "35,u1", "11,u1", 
                 "2,u1", "26,u1", "29,u1", "35,u2", "2,u2", "1A,u1", "32", "16,u1", "34", "33", "36", "34,u1", "35", "31,u1", "3,u1"))

an_chrom_c <- rev(c("1", "2", "3", "4", "4A", "4B", "5", "5A", "6", "7", "8", "9", "10", "11", "12", "13", "14",
                    "15", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "Z", "33", "W"))

pl_chrom_c <- rev(c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", 
                    "16", "17", "18", "19", "20", "21", "X1", "X2", "X3", "X4", "X5", "Y1", "Y2", "Y3", "Y4", "Y5"))

cp_chrom_c <- rev(c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", 
                    "13", "14", "15", "16", "17", "18", "19", "21", "22", "23", "24"))

zf_df$scaffold <- factor(zf_df$scaffold, levels=zf_chroms_c)
an_df$scaffold <- factor(an_df$scaffold, levels=an_chrom_c)
pl_df$scaffold <- factor(pl_df$scaffold, levels=pl_chrom_c)
cp_df$scaffold <- factor(cp_df$scaffold, levels=cp_chrom_c)

zf_df$gene_density <- zf_df$num_of_genes/zf_df$size * 1000000
an_df$gene_density <- an_df$num_of_genes/an_df$size * 1000000
pl_df$gene_density <- pl_df$num_of_genes/pl_df$size * 1000000
cp_df$gene_density <- cp_df$num_of_genes/cp_df$size * 1000000

zf_average_gene_density <- sum(zf_df$num_of_genes) / sum(zf_df$size) * 1000000
an_average_gene_density <- sum(an_df$num_of_genes) / sum(an_df$size) * 1000000
pl_average_gene_density <- sum(pl_df$num_of_genes) / sum(pl_df$size) * 1000000
cp_average_gene_density <- sum(cp_df$num_of_genes) / sum(cp_df$size) * 1000000

zf_plot <- ggplot(zf_df) +
  geom_bar(aes(x=scaffold, y=gene_density, fill=scaffold), width=0.7, alpha=0.7, stat = "summary", fun = "mean") + # gene_density
  geom_hline(aes(yintercept=zf_average_gene_density), lwd=0.2, color="black", linetype="dashed")+
  labs(title="Zebra finch") + 
  theme(plot.title=element_text(size=8), axis.text.x=element_text(size=6, angle=0), 
        axis.text.y=element_text(size=6, angle=0), legend.position = "none", axis.title=element_blank()) + 
  coord_flip(expand=FALSE)
zf_plot

an_plot <- ggplot(an_df) +
  geom_bar(aes(x=scaffold, y=gene_density, fill=scaffold), width=0.7, alpha=0.7, stat = "summary", fun = "mean") + # gene_density
  geom_hline(aes(yintercept=an_average_gene_density), lwd=0.2, color="black", linetype="dashed")+
  labs(title="Anna's hummingbird") + 
  theme(plot.title=element_text(size=8), axis.text.x=element_text(size=6, angle=0), 
        axis.text.y=element_text(size=6, angle=0), legend.position = "none", axis.title=element_blank()) + 
  coord_flip(expand=FALSE)
an_plot

pl_plot <- ggplot(pl_df) +
  geom_bar(aes(x=scaffold, y=gene_density, fill=scaffold), width=0.7, alpha=0.7, stat = "summary", fun = "mean") + # gene_density
  geom_hline(aes(yintercept=pl_average_gene_density), lwd=0.2, color="black", linetype="dashed")+
  labs(title="Platypus") + 
  theme(plot.title=element_text(size=8), axis.text.x=element_text(size=6, angle=0), 
        axis.text.y=element_text(size=6, angle=0), legend.position = "none", axis.title=element_blank()) + 
  coord_flip(expand=FALSE)
pl_plot

cp_plot <- ggplot(cp_df) +
  geom_bar(aes(x=scaffold, y=gene_density, fill=scaffold), width=0.7, alpha=0.7, stat = "summary", fun = "mean") + # gene_density
  geom_hline(aes(yintercept=cp_average_gene_density), lwd=0.2, color="black", linetype="dashed")+
  labs(title="Climbing perch") + 
  theme(plot.title=element_text(size=8), axis.text.x=element_text(size=6, angle=0), 
        axis.text.y=element_text(size=6, angle=0), legend.position = "none", axis.title=element_blank()) + 
  coord_flip(expand=FALSE)
cp_plot

my_all_plot <- c(list(zf_plot), list(an_plot), list(pl_plot), list(cp_plot))
my_layout <- rbind(c(1,2,3,4))
sup.fig.S3.plot <- arrangeGrob(grid.arrange(grobs = my_all_plot, layout_matrix = my_layout))
ggsave("/disk1/juwan_d1/cat_rerun/1_Main_figures/Sup/sup.fig.S3.a-d.png", sup.fig.S3.plot, width=3.5, height=4,units = "in")