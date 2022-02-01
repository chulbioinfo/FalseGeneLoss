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
  #scaffold_stat_df$variable <- factor(scaffold_stat_df$variable, levels=c("final_missing", "final_aligned", "final_gap"))
  scaffold_stat_df$type <- factor(scaffold_stat_df$type, levels=c("chromosome", "unplaced_scaffold", "unlocalized_scaffold"))
  print(1)
  scaffold_stat_df <- scaffold_stat_df[order(-scaffold_stat_df$size),]
  print(2)
  scaffold_stat_df$scaffold <- gsub("chr.","", scaffold_stat_df$scaffold)
  scaffold_order <- scaffold_stat_df$scaffold
  print(scaffold_order)
  scaffold_stat_df$scaffold <- factor(scaffold_stat_df$scaffold, levels=scaffold_order)
  print(head(scaffold_stat_df))
  # generate plot
  fig.1.a_d <- ggplot(data=scaffold_stat_df) + 
    geom_bar(mapping=aes(x=scaffold, y=size/10000000), width=.7, stat="identity", lwd=0.2) + 
    #geom_bar(mapping=aes(x=reorder(gsub("chr.","", scaffold), -size), y=value, fill=variable), width=.7, stat="identity", lwd=0.2) + 
    scale_fill_manual(values=c("brown", "grey", "grey")) + 
    #scale_y_continuous(limits=c(-2.5,16), breaks=c(-2,-1,0,5), labels=c("Repeat","GC","0","5")) + 
    scale_y_continuous(limits=c(-9,20), breaks=c(-7.5,-5.5,-3.5,-1.5,0,10,20), labels=c("Gene density","Missing","Repeat","GC","0","10", "20"), expand=c(0,0)) + 
    new_scale_fill()+
    geom_tile(mapping=aes(x=scaffold, y=-1.5, group=1, fill=scaffold_gc, height=2)) + #, height=1)) + #
    scale_fill_gradientn(colors = brewer.pal(9, "Reds"), limits=c(30,70), name="GC content   ")+
    new_scale_fill()+
    geom_tile(mapping=aes(x=scaffold, y=-3.5, group=1, fill=scaffold_repeat, height=2)) + #, height=1)) + 
    scale_fill_gradientn(colors = brewer.pal(9, "Blues"), limits=c(10,100),  name="Repeat content   ")+
    #
    new_scale_fill()+
    geom_tile(mapping=aes(x=scaffold, y=-5.5, group=1, fill=missing_ratio, height=2)) + #, height=1)) + 
    scale_fill_gradientn(colors = brewer.pal(9, "Greens"), limits=c(0,100),  name="Missing ratio   ")+
    #
    new_scale_fill()+
    geom_tile(mapping=aes(x=scaffold, y=-7.5, group=1, fill=gene_density, height=2)) + #, height=1)) + 
    scale_fill_gradientn(colors = brewer.pal(9, "Oranges"), limits=c(0,120),  name="Gene density (1 Mbp)")+
    #
    #labs(x="", y="Log10(size)") + 
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
  fig.1.a_d_filepath <-paste0("/disk1/juwan_d1/cat_rerun/1_Main_figures/Fig1/fig1.a_d.scaffold_aligned_plot_missing.",species_s,".png")
  ggsave(fig.1.a_d_filepath, fig.1.a_d + theme(legend.position="none"), width=6, height=2.2)
  # save legend as file as well
  if (species_s == "vgp_taegut") {
    fig.1.a_d_legendpath <-paste0("/disk1/juwan_d1/cat_rerun/1_Main_figures/Fig1/fig1.a_d.scaffold_aligned_plot_missing.legend.png")
    fig.1.a_d_legend = gtable_filter(ggplotGrob(fig.1.a_d), "guide-box")
    ggsave(fig.1.a_d_legendpath, fig.1.a_d_legend, width=3, height=6)}
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

# ASSEMBLY COMPARISON
# (2-1) VGP vs prior (1kbp)
vgp_vs_prev_df <- vroom("/disk1/juwan_d1/cat_rerun/zebra_finch/tm/CAT2MISSING_pre_taegut/window_size_gc_and_repeat/summary_gc_n_repeat.all.tsv", delim="\t")
vgp_vs_prev_df = as.data.frame(vgp_vs_prev_df)
head(vgp_vs_prev_df)
colnames(vgp_vs_prev_df) <- c("chrom", "assembly_id", "Repeat", "GC","Size")
vgp_vs_prev_df$species <- ifelse(vgp_vs_prev_df$assembly_id %in% c("vgp_taegut", "pre_taegut"), "Zebra finch",
                                 ifelse(vgp_vs_prev_df$assembly_id %in% c("vgp_calann", "pre_calann"), "Anna's\nhummingbird",
                                        ifelse(vgp_vs_prev_df$assembly_id %in% c("vgp_ornana", "pre_ornana"), "Platypus", "Climbing\nperch")))
vgp_vs_prev_df$type <- ifelse(grepl("^vgp_", vgp_vs_prev_df$assembly_id), "VGP","Prior")
vgp_vs_prev_df <- vgp_vs_prev_df %>%
  select(chrom,type,Repeat,GC,species)

# (2-2) Missing vs. present (1Kbp)
# input file paths
zebra_finch_path_s="/disk1/juwan_d1/cat_rerun/zebra_finch/tm/CAT2MISSING_pre_taegut/missing_chrom/3_consensus/"
anna_hummingbird_path_s="/disk1/juwan_d1/cat_rerun/anna/tm/CAT2MISSING_pre_calann_backup/missing_chrom/3_consensus/"
platypus_path_s="/disk1/juwan_d1/cat_rerun/platypus/tm/CAT2MISSING_pre_ornana/missing_chrom/3_consensus/"
climbing_perch_path_s="/disk1/juwan_d1/cat_rerun/climbing_perch/tm/CAT2MISSING_pre_anates/missing_chrom/3_consensus/"

zebra_finch_missing_vs_present_df <- read.csv(paste0(zebra_finch_path_s, "GC_n_repeat_summary.tsv"), sep="\t", stringsAsFactors = FALSE)
anna_hummingbird_missing_vs_present_df <- read.csv(paste0(anna_hummingbird_path_s, "GC_n_repeat_summary.tsv"), sep="\t", stringsAsFactors = FALSE)
platypus_missing_vs_present_df <- read.csv(paste0(platypus_path_s, "GC_n_repeat_summary.tsv"), sep="\t", stringsAsFactors = FALSE)
climbing_perch_missing_vs_present_df <- read.csv(paste0(climbing_perch_path_s, "GC_n_repeat_summary.tsv"), sep="\t", stringsAsFactors = FALSE)

zebra_finch_missing_vs_present_df$species <- "Zebra finch"
anna_hummingbird_missing_vs_present_df$species <- "Anna's\nhummingbird"
platypus_missing_vs_present_df$species <- "Platypus"
climbing_perch_missing_vs_present_df$species <- "Climbing\nperch"

missing_vs_present_df <- rbind(zebra_finch_missing_vs_present_df,
                               anna_hummingbird_missing_vs_present_df,
                               platypus_missing_vs_present_df,
                               climbing_perch_missing_vs_present_df)

missing_vs_present_df$type <- factor(missing_vs_present_df$type, levels=c("aligned", "missing"))

# (2-3) merge (VGP vs prior) + (Missing vs. present)
merged_window_df <-rbind(vgp_vs_prev_df, missing_vs_present_df)
# melt and remove prior assembly result (prior: not in use)
gc_n_repeat_window_df <- melt(merged_window_df, id.vars=c("chrom", "type","species"), measure.vars = c("GC","Repeat"))
subset_gc_n_repeat_window_df <- subset(gc_n_repeat_window_df,type %in% c("VGP", "missing","aligned"))

subset_gc_n_repeat_window_df$type <- factor(subset_gc_n_repeat_window_df$type, 
                                            levels=c("VGP", "missing", "aligned"))
subset_gc_n_repeat_window_df$species <- factor(subset_gc_n_repeat_window_df$species, 
                                               levels=c("Zebra finch","Anna's\nhummingbird","Platypus", "Climbing\nperch"))

fig.1.e.plot <- ggplot(subset_gc_n_repeat_window_df,aes(x=value, group=type))+
  facet_grid(variable~species, scale="free")+
  geom_density(data=subset(subset_gc_n_repeat_window_df,type!="VGP"), aes(fill=type), alpha = 0.5, position = "identity", lwd=0)+
  geom_density(data=subset(subset_gc_n_repeat_window_df,type=="VGP"), alpha = 0.5, position = "identity", color="#00226b", linetype="dotted", lwd=0.5)+
  stat_central_tendency(type = "mean", aes(x=value, color=type), geom="line", lwd=0.5)+
  scale_fill_manual(values=c("brown","grey"))+
  scale_color_manual(values=c("#00226b","brown","grey"))+
  scale_y_continuous(breaks=c(0, 0.05,0.1))+
  theme_classic2()+
  labs(x="",y = "Density")+
  theme(text = element_text(size=6),
        strip.text=element_text(size=6),
        strip.background = element_blank(),
        axis.line=element_line(size=0.1),
        legend.position = "none",
        legend.text= element_text(size=6)) + 
  coord_cartesian(expand=FALSE)
# save output
fig.1.e.plot_filepath <-"/disk1/juwan_d1/cat_rerun/1_Main_figures/Fig1/fig1.e.gc_n_repeat.png"
ggsave(fig.1.e.plot_filepath, fig.1.e.plot, width=6, height=2.5)

#######################
### Fig.1. e (done) ###
#######################

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

zf_list <- fig.1.f.maker(species_s="Zebra finch",species_path=zf_path)
an_list <- fig.1.f.maker(species_s="Anna's hummingbird",species_path=an_path)
pl_list <- fig.1.f.maker(species_s="Platypus",species_path=pl_path)
cl_list <- fig.1.f.maker(species_s="Climbing perch",species_path=cl_path)
gene_list <- c(zf_list, an_list, pl_list, cl_list)
my_layout <- rbind(c(1,2,3,4))
fig.1.f.plot <- arrangeGrob(grobs=gene_list, layout_matrix=my_layout)
fig.1.f.plot_path <- paste0(output_path, "Fig1.f.scatter_plot.png")
ggsave(fig.1.f.plot_path, fig.1.f.plot, width=5.8, height=1.8) 
