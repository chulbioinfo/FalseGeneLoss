library(dplyr)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(scales)
library(grid)
library(gtable)
library(gridExtra)

#########################
### Shared variables ####
#########################

species_list <- c("Zebra finch", "Anna's hummingbird", "Platypus", "Climbing perch")

#########################
######## Fig 3.a ########
#########################
# load input data
zf_gene_df <- read.csv("./input/zebra_finch.gc_and_repeat_of_genes.stats", sep="\t")
an_gene_df <- read.csv("./input/anna.gc_and_repeat_of_genes.stats", sep="\t")
pl_gene_df <- read.csv("./input/platypus.gc_and_repeat_of_genes.stats", sep="\t")
cl_gene_df <- read.csv("./input/climbing_perch.gc_and_repeat_of_genes.stats", sep="\t")

# set the boundary
all_gene_df <- rbind(zf_gene_df, an_gene_df, pl_gene_df, cl_gene_df)
gc_min <- max(min(all_gene_df$GCcontents)-5,0)
gc_max <- min(max(all_gene_df$GCcontents)+5,100)
repeat_min <- max(min(all_gene_df$repeat_percent)-5,0)
repeat_max <- min(max(all_gene_df$repeat_percent)+5,100)

# generate plot
fig.3.a.maker <- function(species_s, gene_df) {
  fgl_df <- gene_df %>%
    filter(type=="fgl")
  gene_df$type <- factor(gene_df$type, levels=c("fgl","control"))
  scatter <- ggplot(gene_df,aes(x=GCcontents, y=repeat_percent,color=type))+
    geom_point(size=0.01, alpha=0.05) +
    geom_point(data=fgl_df,aes(x=GCcontents, y=repeat_percent),color="brown",size=0.01) + 
    geom_point(data=subset(gene_df, type=="fgl"),aes(x=mean(GCcontents), y=mean(repeat_percent)),color="black", fill="brown",lwd=2, size=8, shape=21) + 
    geom_point(data=subset(gene_df, type=="control"),aes(x=mean(GCcontents), y=mean(repeat_percent)),color="black", fill="grey",lwd=2, size=8, shape=21) + 
    scale_color_manual(values=c("brown","grey")) +
    labs(title = species_s,x="GC-content",y="Repeat content") +
    coord_cartesian(xlim=c(gc_min,gc_max),ylim=c(repeat_min,repeat_max),expand=FALSE) + 
    theme(legend.position="none",
          plot.title = element_text(size=7),
          axis.text=element_text(size=6),
          axis.line = element_line(size=0.01),
          axis.ticks = element_line(size=0.05),
          axis.title=element_blank())
  p2 <- ggMarginal(scatter, type="density", size=6, groupFill = TRUE, lwd=0.05, outlier.shape=NA)
  my_list <- list(p2)
  return(my_list)}

zf_list <- fig.3.a.maker(species_s="Zebra finch",gene_df=zf_gene_df)
an_list <- fig.3.a.maker(species_s="Anna's hummingbird",gene_df=an_gene_df)
pl_list <- fig.3.a.maker(species_s="Platypus",gene_df=pl_gene_df)
cl_list <- fig.3.a.maker(species_s="Climbing perch",gene_df=cl_gene_df)
gene_list <- c(zf_list, an_list, pl_list, cl_list)
my_layout <- rbind(c(1,2,3,4))
fig.3.a.plot <- arrangeGrob(grobs=gene_list, layout_matrix=my_layout)
fig.3.a.plot_path <- "./Fig3.a.scatter_plot.png"
ggsave(fig.3.a.plot_path, fig.3.a.plot, width=5.8, height=1.8) 

################################
######## Fig 3.a (done) ########
################################



#########################
######## Fig 3.b ########
#########################

# load input data + bind them as one single dataframe
input_df <- read.csv("./input/compare.missing.tsv", sep="\t")
input_df$species <- factor(input_df$species, levels=c("Zebra finch", "Platypus", "Anna's hummingbird", "Climbing perch"))
input_df$platform <- factor(input_df$platform, levels=c("Sanger-based", "Illumina-based"))
# only leave cds, intron, and intergenic regions 
input_df <- input_df %>%
  filter(type %in% c("cds", "intron","intergenic", "noncoding_genes"))
# rename CDS
input_df$type <- ifelse(input_df$type == "cds", "CDS (exon)", input_df$type)
input_df$type <- ifelse(input_df$type == "intron", "intronic", input_df$type)
input_df$type <- ifelse(input_df$type == "noncoding_genes", "noncoding genes", input_df$type)
input_df$type <- factor(input_df$type, 
                                      levels=c("intergenic", "noncoding genes", "intronic","CDS (exon)"))
# draw a plot
fig3.b.plot <- ggplot(input_df, aes(x=type, y=missing_percent, fill=species))+
  facet_wrap(species~.,nrow=4)+
  geom_bar(stat="identity", width=0.6) +
  labs(y="Missing (%)", fill="Types") + 
  scale_fill_manual(values=c("#592800", "#b35600",  "#00695c", "#4db8ab"))+
  theme_pubr()+
  theme(axis.title.y=element_blank(),
        text=element_text(size=6),
        strip.text = element_blank(),
        strip.background = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size=6),
        legend.key.size = unit(0.2, "cm"),
        legend.key.width = unit(0.2,"cm")) +
  coord_flip(ylim=c(0,20), xlim=c(0, 4.5), expand=FALSE)

fig3.b.plot_filepath <- "./fig3.b.missing_ratio.png"
fig3.b.legend_filepath <- "./fig3.b.legend.png"
ggsave(fig3.b.plot_filepath, fig3.b.plot + theme(legend.position = "none"), width=2.5, height=3.2)
fig3.b.legend = gtable_filter(ggplotGrob(fig3.b.plot), "guide-box")
ggsave(fig3.b.legend_filepath, fig3.b.legend, width=3, height=0.5,units = "in")

################################
######## Fig 3.b (done) ########
################################


#########################
######## Fig 3.c ########
#########################

# load input data + bind them as one single dataframe
missing_gene_colnames <- c("chrom","start","end","gene","missing_length","missing_ratio")
zf_missing_gene_df <- read.csv("./input/zebra_finch.ref.genes.not_aligned.bed", sep="\t", col.names = missing_gene_colnames)
an_missing_gene_df <- read.csv("./input/anna.ref.genes.not_aligned.bed", sep="\t", col.names = missing_gene_colnames)
pl_missing_gene_df <- read.csv("./input/platypus.ref.genes.not_aligned.bed", sep="\t", col.names = missing_gene_colnames)
cl_missing_gene_df <- read.csv("./input/climbing_perch.ref.genes.not_aligned.bed", sep="\t", col.names = missing_gene_colnames)
zf_missing_gene_df$species <- "Zebra finch"
an_missing_gene_df$species <- "Anna's hummingbird"
pl_missing_gene_df$species <- "Platypus"
cl_missing_gene_df$species <- "Climbing perch"
all_missing_gene_df <- rbind(zf_missing_gene_df, an_missing_gene_df, pl_missing_gene_df, cl_missing_gene_df)
all_missing_gene_df$species <- factor(all_missing_gene_df$species, levels=species_list)
all_missing_gene_df$platform <- ifelse(all_missing_gene_df$species %in% c("Zebra finch", "Platypus"), "Sanger-based", "Illumina-based")
all_missing_gene_df$platform <- factor(all_missing_gene_df$platform, levels=c("Sanger-based", "Illumina-based"))
tail(all_missing_gene_df)
# generate plot
fig.3.c.plot <- ggplot(all_missing_gene_df, aes(x=missing_ratio)) + 
  geom_vline(aes(xintercept=10),color="grey",linetype="dotted")+
  scale_x_continuous(breaks=c(0,10,25,50,75,100))+
  scale_y_continuous(breaks=c(0,0.20,0.40,0.60,0.80,1), labels=c(0,20,40,60,80,100))+
  scale_color_manual(values=c("#592800", "#b35600",  "#00695c", "#4db8ab"))+
  stat_ecdf(aes(color=species, linetype=platform), geom="step") + 
  scale_linetype_manual(values=c("dashed","solid"))+
  coord_cartesian(xlim=c(0,100), ylim=c(0,1),expand=FALSE)+
  theme_pubr() +  
  labs(x="Missing (%)", y="Cumulative density of genes (%)", color="Species", linetype="Platform")+
  theme(axis.text=element_text(size=6),
        axis.title=element_text(size=6),
        strip.text = element_text(size=6),
        strip.background = element_blank(),
        legend.position = "bottom",
        legend.text=element_text(size=6),
        legend.direction = "horizontal",
        legend.box="vertical",
        legend.title = element_text(size=7)) + 
  guides(color = guide_legend(ncol = 2,  keyheight = unit(1.5, "mm")))

# save output
fig.3.c.plot_filepath = "./fig3.c.cumulative_density.png"
fig.3.c.legend_filepath = "./fig3.c.legend.png"
ggsave(fig.3.c.plot_filepath, fig.3.c.plot + theme(legend.position = "none"), width=3, height=3) 
fig.3.c.legend = gtable_filter(ggplotGrob(fig.3.c.plot), "guide-box")
ggsave(fig.3.c.legend_filepath, fig.3.c.legend, width=4, height=3,units = "in")

################################
######## Fig 3.c (done) ########
################################
