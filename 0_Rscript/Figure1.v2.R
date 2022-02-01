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


fig.1.a_d.maker <- function(my_path,species_s) {
  # load input file
  scaffold_stat_df <- read.csv(paste0(my_path, species_s,".stat_summary.tabs"), sep="\t", stringsAsFactors = FALSE)
  # Fig. 1. a-d: scaffold > 100kbp + remove Y chromosome result for platypus
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
    # manual replacement of chr 29 of zebra finch
    # Size
    scaffold_stat_df[32,3] <- 1995839
    # GC content
    scaffold_stat_df[32,4] <- 58.6618434516287
    # repeat content
    scaffold_stat_df[32,5] <- 31.93894898
    # missing ratio
    scaffold_stat_df[32,11] <- 39.4  
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
  # type of scaffolds: chromosome, unplaced, and unlocalized scaffold
  scaffold_stat_df.m <- melt(scaffold_stat_df,
                             id.vars = c("scaffold", "type", "size", "log10_size", "scaffold_gc","scaffold_repeat"), 
                             measure.vars=c("final_missing", "final_aligned", "final_gap"))
  print(head(scaffold_stat_df.m))
  scaffold_stat_df.m$variable <- factor(scaffold_stat_df.m$variable, levels=c("final_missing", "final_aligned", "final_gap"))
  scaffold_stat_df.m$type <- factor(scaffold_stat_df.m$type, levels=c("chromosome", "unplaced_scaffold", "unlocalized_scaffold"))
  # generate plot
  fig.1.a_d <- ggplot(data=scaffold_stat_df.m) + 
    geom_bar(mapping=aes(x=reorder(gsub("chr.","", scaffold), -size), y=value, fill=variable), width=.7, stat="identity", lwd=0.2) + 
    scale_fill_manual(values=c("brown", "grey", "grey")) + 
    scale_y_continuous(limits=c(-2.5,8.5), breaks=c(-2,-1,0,5), labels=c("Repeat","GC","0","5")) + 
    new_scale_fill()+
    geom_tile(mapping=aes(x=reorder(gsub("chr.","", scaffold), -size), y=-1, group=1, fill=scaffold_gc)) + 
    scale_fill_gradientn(colors = brewer.pal(9, "Reds"), limits=c(30,70), name="GC content   ")+
    new_scale_fill()+
    geom_tile(mapping=aes(x=reorder(gsub("chr.","", scaffold), -size), y=-2, group=1, fill=scaffold_repeat)) + 
    scale_fill_gradientn(colors = brewer.pal(9, "Blues"), limits=c(10,100),  name="Repeat content   ")+
    labs(x="", y="Log10(size)") + 
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
  fig.1.a_d_filepath <-paste0("/disk1/juwan_d1/cat_rerun/1_Main_figures/Fig1/fig1.a_d.scaffold_aligned_plot_missing.",species_s,".png")
  ggsave(fig.1.a_d_filepath, fig.1.a_d + theme(legend.position="none"), width=6, height=1.2)
  # save legend as file as well
  if (species_s == "vgp_taegut") {
    fig.1.a_d_legendpath <-paste0("/disk1/juwan_d1/cat_rerun/1_Main_figures/Fig1/fig1.a_d.scaffold_aligned_plot_missing.legend.png")
    fig.1.a_d_legend = gtable_filter(ggplotGrob(fig.1.a_d), "guide-box")
    ggsave(fig.1.a_d_legendpath, fig.1.a_d_legend, width=3, height=3)}
}

# input file paths
zebra_finch_path_s="/disk1/juwan_d1/cat_rerun/zebra_finch/tm/CAT2MISSING_pre_taegut/chromosomeMap/out_correlation/"
anna_hummingbird_path_s="/disk1/juwan_d1/cat_rerun/anna/tm/CAT2MISSING_pre_calann/chromosomeMap/out_correlation/"
platypus_path_s="/disk1/juwan_d1/cat_rerun/platypus/tm/CAT2MISSING_pre_ornana/chromosomeMap/out_correlation/"
climbing_perch_path_s="/disk1/juwan_d1/cat_rerun/climbing_perch/tm/CAT2MISSING_pre_anates/chromosomeMap/out_correlation/"

mylist <- list('vgp_taegut'=zebra_finch_path_s,
               "vgp_calann"=anna_hummingbird_path_s,
               'vgp_ornana'=platypus_path_s,
               'vgp_anates'=climbing_perch_path_s)
for (i in seq_along(mylist)){
  species_s = names(mylist)[i]
  scaffold_alignment_path_s = mylist[[i]]
  fig.1.a_d.maker(my_path=scaffold_alignment_path_s,species_s=species_s)
}

#########################
### Fig.1. a-d (done) ###
#########################


##################
#### Fig.1. e ####
##################

species_list <- c("Zebra finch", "Anna's hummingbird", "Platypus", "Climbing perch")

zf_path <- "/disk1/juwan_d1/cat_rerun/zebra_finch/tm/CAT2MISSING_pre_taegut/"
an_path <- "/disk1/juwan_d1/cat_rerun/anna/tm/CAT2MISSING_pre_calann/"
pl_path <- "/disk1/juwan_d1/cat_rerun/platypus/tm/CAT2MISSING_pre_ornana/"
cl_path <- "/disk1/juwan_d1/cat_rerun/climbing_perch/tm/CAT2MISSING_pre_anates/"

output_path <- "/disk1/juwan_d1/cat_rerun/1_Main_figures/Fig1/"

# generate plot
fig.1.f.maker <- function(species_s, species_path) {
  ref_df <- read.csv(paste0(species_path, "window_size_gc_and_repeat/ref.summary.tsv"),sep="\t", header=FALSE, col.names=c("chrom", "type", "repeat", "gc", "size"))
  aligned_df <- read.csv(paste0(species_path, "window_size_gc_and_repeat/aligned.summary.tsv"), sep="\t", header=FALSE, col.names=c("chrom", "type", "repeat", "gc", "size"))
  missing_df <- read.csv(paste0(species_path, "window_size_gc_and_repeat/missing.summary.tsv"), sep="\t", header=FALSE, col.names=c("chrom", "type", "repeat", "gc", "size"))
  aligned_n_missing_df <- rbind(missing_df, aligned_df)
  aligned_n_missing_df$type <- factor(aligned_n_missing_df$type, levels=c("aligned", "missing"))
  gc_wilcox_result <- wilcox.test(aligned_df$gc, missing_df$gc)
  repeat_wilcox_result <- wilcox.test(aligned_df$repeat., missing_df$repeat.)
  print(gc_wilcox_result$p.value)
  print(paste0("Aligned (GC): ", mean(aligned_df$gc)))
  print(paste0("Missing (GC): ", mean(missing_df$gc)))
  print(repeat_wilcox_result$p.value)
  print(paste0("Aligned (repeat): ", mean(aligned_df$repeat.)))
  print(paste0("Missing (repeat): ", mean(missing_df$repeat.)))
  scatter <- ggplot(aligned_n_missing_df, aes(x=gc, y=repeat., color=type))+
    geom_point(alpha=0.03, size=0.01)+ 
    geom_point(data=subset(aligned_n_missing_df, type=="aligned"),aes(x=mean(gc), y=mean(repeat.)),color="black", fill="grey",lwd=2, size=8, shape=21) + 
    geom_point(data=subset(aligned_n_missing_df, type=="missing"),aes(x=mean(gc), y=mean(repeat.)),color="black", fill="brown",lwd=2, size=8, shape=21) + 
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

zf_list <- fig.1.f.maker(species_s="Zebra finch",species_path=zf_path)
an_list <- fig.1.f.maker(species_s="Anna's hummingbird",species_path=an_path)
pl_list <- fig.1.f.maker(species_s="Platypus",species_path=pl_path)
cl_list <- fig.1.f.maker(species_s="Climbing perch",species_path=cl_path)
gene_list <- c(zf_list, an_list, pl_list, cl_list)
my_layout <- rbind(c(1,2,3,4))
fig.1.f.plot <- arrangeGrob(grobs=gene_list, layout_matrix=my_layout)
fig.1.f.plot_path <- paste0(output_path, "Fig1.f.scatter_plot.png")
ggsave(fig.1.f.plot_path, fig.1.f.plot, width=5.8, height=1.8) 

#######################
### Fig.1. e (done) ###
#######################

