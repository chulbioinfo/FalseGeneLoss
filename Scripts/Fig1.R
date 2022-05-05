library(ggpubr)
library(ggrepel)
library(ggplot2)
library(gridExtra)
library(vroom)
library(reshape2)
library(ggnewscale)
library(ggExtra)
##################
#### Figure 1 ####
##################

##################
### Fig.1. a-d ###
##################

rename_newly_identified <- function(scaffold_df, chrs_c) {
  scaffold_df$scaffold <- ifelse(scaffold_df$scaffold %in% chrs_c, gsub("chr.", "chr.*", scaffold_df$scaffold), scaffold_df$scaffold)
  scaffold_df$scaffold <- ifelse(scaffold_df$scaffold %in% chrs_c, gsub("up", "*up", scaffold_df$scaffold), scaffold_df$scaffold)
  return(scaffold_df)
}

fig.1.a_d.maker <- function(species_s) {
  # load input file
  # Fig. 1. a-d: scaffold > 100kbp + remove Y chromosome result for platypus
  scaffold_stat_df <- read.csv(paste0("./input/", species_s,".stat_summary.tsv"), sep="\t", stringsAsFactors = FALSE)
  scaffold_stat_df <- scaffold_stat_df %>%
    filter(scaffold_stat_df[["size"]] >=100000)
  if (species_s == "vgp_ornana") {scaffold_stat_df <- scaffold_stat_df %>% filter(!scaffold_stat_df[["scaffold"]] %like% "^chr.Y")} else {scaffold_stat_df <- scaffold_stat_df}
  # rename newly identified chromosomes
  newly_identified_chrs_c <- c()
  if (species_s %in% c("vgp_anates", "vgp_calann")) {
    newly_identified_chrs_c <- c() #levels(factor(scaffold_stat_df$scaffold))
  } else if (species_s == "vgp_taegut") {
    newly_identified_chrs_c <- c("chr.29", "chr.30", "chr.31", "chr.32", "chr.33", "chr.34", "chr.35", "chr.36",
                                 "chr.29,u1", "chr.30,u1", "chr.31,u1", "chr.32,u1", "chr.33,u1", "chr.34,u1", "chr.35,u1", "chr.36,u1",
                                 "chr.29,u2", "chr.30,u2", "chr.31,u2", "chr.32,u2", "chr.33,u2", "chr.34,u2", "chr.35,u2", "chr.36,u2")
  } else if (species_s == "vgp_ornana") {
    newly_identified_chrs_c <- c() #c("chr.8", "chr.9", "chr.16", "chr.19", "chr.21", "chr.X4")
  }
  newly_identified_chrs_c <- subset(scaffold_stat_df$scaffold, scaffold_stat_df$missing_ratio >= 30)
  print(newly_identified_chrs_c)
  scaffold_stat_df <- rename_newly_identified(scaffold_df=scaffold_stat_df, chrs_c=newly_identified_chrs_c)
  # Calculation of missing ratio of each scaffold
  all_genome_length = sum(scaffold_stat_df$size)
  all_unaligned_length = sum(scaffold_stat_df$missing_length)
  all_aligned_length = sum(scaffold_stat_df$aligned_length)
  all_unaligned_pct = all_unaligned_length/all_genome_length*100
  all_aligned_pct =  all_aligned_length/all_genome_length*100
  print(head(scaffold_stat_df))
  # Transform the size and missing ratio as log-scaled
  scaffold_stat_df$log10_size <- log10(scaffold_stat_df$size)
  scaffold_stat_df$final_missing <- (scaffold_stat_df$missing_ratio)/100*scaffold_stat_df$log10_size
  scaffold_stat_df$final_aligned <- (scaffold_stat_df$aligned_ratio)/100*scaffold_stat_df$log10_size
  scaffold_stat_df$final_gap <- (scaffold_stat_df$gap_ratio)/100*scaffold_stat_df$log10_size
  # Reorder input dataframe (descending order of its size)
  scaffold_stat_df <- scaffold_stat_df[with(scaffold_stat_df, order(-size)), ]
  # Rename unplaced scaffolds (u1, u2, ...)
  unplaced_scaffolds <- scaffold_stat_df%>%
    filter(type=="unplaced_scaffold")
  unplaced_scaffolds <- unplaced_scaffolds$scaffold
  unplaced_scaffolds_name_dict <- list()
  int_start=0
  for (i in unplaced_scaffolds) {
    int_start <- int_start + 1
    unplaced_scaffolds_name_dict[i] = int_start}
  scaffold_stat_df$scaffold <- ifelse(scaffold_stat_df$scaffold  %like% "^chr.", 
                                      gsub("chr.","", scaffold_stat_df$scaffold),
                                      paste0("up",unplaced_scaffolds_name_dict[scaffold_stat_df$scaffold]))
  newly_identified_chrs_c <- subset(scaffold_stat_df$scaffold, scaffold_stat_df$missing_ratio >= 30)
  scaffold_stat_df <- rename_newly_identified(scaffold_df=scaffold_stat_df, chrs_c=newly_identified_chrs_c)
  scaffold_stat_df$gene_density <- scaffold_stat_df$num_of_genes/scaffold_stat_df$size*1000000
  # type of scaffolds: chromosome, unplaced, and unlocalized scaffold
  scaffold_stat_df.m <- melt(scaffold_stat_df,
                             id.vars = c("scaffold", "type", "size", "log10_size", "scaffold_gc","scaffold_repeat", "missing_ratio"), 
                             measure.vars=c("final_missing", "final_aligned", "final_gap"))
  print(head(scaffold_stat_df))
  scaffold_stat_df$type <- factor(scaffold_stat_df$type, levels=c("chromosome", "unplaced_scaffold", "unlocalized_scaffold"))
  scaffold_stat_df <- scaffold_stat_df[order(-scaffold_stat_df$size),]
  scaffold_stat_df$scaffold <- gsub("chr.","", scaffold_stat_df$scaffold)
  scaffold_order <- scaffold_stat_df$scaffold
  print(scaffold_order)
  scaffold_stat_df$scaffold <- factor(scaffold_stat_df$scaffold, levels=scaffold_order)
  print(head(scaffold_stat_df))
  # generate plot
  fig.1.a_d <- ggplot(data=scaffold_stat_df) + 
    geom_bar(mapping=aes(x=scaffold, y=size/10000000), width=.7, stat="identity", lwd=0.2) + 
    scale_fill_manual(values=c("brown", "grey", "grey")) + 
    scale_y_continuous(limits=c(-9,20), breaks=c(-7.5,-5.5,-3.5,-1.5,0,10,20), labels=c("Gene density","Missing","Repeat","GC","0","10", "20"), expand=c(0,0)) + 
    new_scale_fill()+
    geom_tile(mapping=aes(x=scaffold, y=-1.5, group=1, fill=scaffold_gc, height=2)) + 
    scale_fill_gradientn(colors = brewer.pal(9, "Reds"), limits=c(30,70), name="GC content   ")+
    new_scale_fill()+
    geom_tile(mapping=aes(x=scaffold, y=-3.5, group=1, fill=scaffold_repeat, height=2)) + 
    scale_fill_gradientn(colors = brewer.pal(9, "Blues"), limits=c(10,100),  name="Repeat content   ")+
    new_scale_fill()+
    geom_tile(mapping=aes(x=scaffold, y=-5.5, group=1, fill=missing_ratio, height=2)) + 
    scale_fill_gradientn(colors = brewer.pal(9, "Greens"), limits=c(0,100),  name="Missing ratio   ")+
    new_scale_fill()+
    geom_tile(mapping=aes(x=scaffold, y=-7.5, group=1, fill=gene_density, height=2)) + 
    scale_fill_gradientn(colors = brewer.pal(9, "Oranges"), limits=c(0,120),  name="Gene density (1 Mbp)")+
    labs(x="", y="Size (10 Mbp)") + 
    theme_classic()+ 
    theme(legend.position="bottom",
          legend.key.width = unit(0.6,"cm"),
          legend.key.size = unit(0.2, "cm"),
          legend.key.height = unit(0.2,"cm"),
          legend.box="vertical",
          legend.box.just = "left",
          text=element_text(size=6),
          axis.text.x=element_text(size=6, angle=90, hjust=1),
          panel.grid.major.y = element_line())
  fig.1.a_d_filepath <-paste0("./fig1.a_d.scaffold_aligned_plot_missing.",species_s,".png")
  ggsave(fig.1.a_d_filepath, fig.1.a_d + theme(legend.position="none"), width=6, height=2.2)
  
  # save legend as file as well
  if (species_s == "vgp_taegut") {
    fig.1.a_d_legendpath <-paste0("./fig1.a_d.scaffold_aligned_plot_missing.legend.png")
    fig.1.a_d_legend = gtable_filter(ggplotGrob(fig.1.a_d), "guide-box")
    ggsave(fig.1.a_d_legendpath, fig.1.a_d_legend, width=3, height=6)}
}

fig.1.a_d.maker(species_s="vgp_taegut")
fig.1.a_d.maker(species_s="vgp_calann")
fig.1.a_d.maker(species_s="vgp_ornana")
fig.1.a_d.maker(species_s="vgp_anates")


#########################
### Fig.1. a-d (done) ###
#########################


##################
#### Fig.1. e ####
##################

species_list <- c("Zebra finch", "Anna's hummingbird", "Platypus", "Climbing perch")

# generate plot
fig.1.e.maker <- function(species_s) {
  ref_df <- read.csv(paste0("PATH_TO_REF_DATA"),sep="\t", header=FALSE, col.names=c("chrom", "type", "repeat", "gc", "size"))
  aligned_df <- read.csv(paste0("PATH_TO_ALN_DATA"), sep="\t", header=FALSE, col.names=c("chrom", "type", "repeat", "gc", "size"))
  missing_df <- read.csv(paste0("PATH_TO_MISSING_DATA"), sep="\t", header=FALSE, col.names=c("chrom", "type", "repeat", "gc", "size"))
  aligned_n_missing_df <- rbind(missing_df, aligned_df)
  aligned_n_missing_df$type <- factor(aligned_n_missing_df$type, levels=c("aligned", "missing"))
  gc_wilcox_result <- wilcox.test(aligned_df$gc, missing_df$gc)
  repeat_wilcox_result <- wilcox.test(aligned_df$repeat., missing_df$repeat.)
  scatter <- ggplot(aligned_n_missing_df, aes(x=gc, y=repeat., color=type))+
    geom_point(alpha=0.05, size=0.01)+ 
    geom_point(data=subset(aligned_n_missing_df, type=="missing"),aes(x=mean(gc), y=mean(repeat.)),color="black", fill="brown",lwd=2, size=8, shape=21) + 
    geom_point(data=subset(aligned_n_missing_df, type=="aligned"),aes(x=mean(gc), y=mean(repeat.)),color="black", fill="grey",lwd=2, size=8, shape=21) + 
    scale_color_manual(values=c("grey", "brown")) +
    labs(title = species_s,x="GC-content",y="Repeat content") +
    coord_cartesian(xlim=c(0,100),ylim=c(0,100),expand=FALSE) + 
    theme(legend.position="none",
          plot.title = element_text(size=7),
          axis.text=element_text(size=6),
          axis.line = element_line(size=0.01),
          axis.ticks = element_line(size=0.05),
          axis.title=element_blank())
  p2 <- ggMarginal(scatter, type="density", size=6, groupFill = TRUE, lwd=0.05)
  my_list <- list(p2)
  return(my_list)}

zf_list <- fig.1.e.maker(species_s="Zebra finch")
an_list <- fig.1.e.maker(species_s="Anna's hummingbird")
pl_list <- fig.1.e.maker(species_s="Platypus")
cl_list <- fig.1.e.maker(species_s="Climbing perch")
gene_list <- c(zf_list, an_list, pl_list, cl_list)
my_layout <- rbind(c(1,2,3,4))
fig.1.e.plot <- arrangeGrob(grobs=gene_list, layout_matrix=my_layout)
fig.1.e.plot_path <- "./Fig1.e.scatter_plot.png"
ggsave(fig.1.e.plot_path, fig.1.e.plot, width=5.8, height=1.8) 


#######################
### Fig.1. e (done) ###
#######################
