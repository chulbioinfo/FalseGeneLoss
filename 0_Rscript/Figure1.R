library(ggpubr)
library(ggrepel)
library(ggplot2)
library(gridExtra)
library(vroom)
library(reshape2)
library(ggnewscale)

##################
#### Figure 1 ####
##################

##################
### Fig.1. a-d ###
##################

fig.1.a_d.maker <- function(my_path,species_s) {
  # load input file
  scaffold_stat_df <- read.csv(paste0(my_path, species_s,".stat_summary.tabs"), sep="\t", stringsAsFactors = FALSE)
  
  # Fig. 1. a-d: scaffold > 100kbp + remove Y chromosome result for platypus
  scaffold_stat_df <- scaffold_stat_df %>%
    filter(scaffold_stat_df[["size"]] >=100000)
  if (species_s == "vgp_ornana") {scaffold_stat_df <- scaffold_stat_df %>% filter(!scaffold_stat_df[["scaffold"]] %like% "^chr.Y")} else {scaffold_stat_df <- scaffold_stat_df}
  
  # Calculation of missing ratio of each scaffold
  all_genome_length = sum(scaffold_stat_df$size)
  all_unaligned_length = sum(scaffold_stat_df$unaligned_ratio)
  all_aligned_length = all_genome_length-all_unaligned_length
  all_unaligned_pct = all_unaligned_length/all_genome_length*100
  all_aligned_pct =  all_aligned_length/all_genome_length*100
  
  # Transform the size and missing ratio as log-scaled
  scaffold_stat_df$log10_size <- log10(scaffold_stat_df$size)
  scaffold_stat_df$log10_percent <- (scaffold_stat_df$size * (100-scaffold_stat_df$unaligned_ratio))/scaffold_stat_df$size/100*scaffold_stat_df$log10_size
  scaffold_stat_df$colors <-"#6b6b6b"
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
  
  # type of scaffolds: chromosome, unplaced, and unlocalized scaffold
  scaffold_stat_df$type <- factor(scaffold_stat_df$type, levels=c("chromosome", "unplaced_scaffold", "unlocalized_scaffold"))

  # generate plot
  fig.1.a_d <- ggplot(data=scaffold_stat_df) + 
    geom_bar(mapping=aes(x=reorder(gsub("chr.","", scaffold), -size), y=log10_size), width=.7, 
             fill="brown", stat="identity", lwd=0.2) + 
    geom_bar(mapping=aes(x=reorder(gsub("chr.","", scaffold), -size), y=log10_percent), width=.7, 
             fill="grey", stat="identity", lwd=0) +
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
  ggsave(fig.1.a_d_filepath, fig.1.a_d + theme(legend.position="none"), width=6, height=1.5)
  # save legend as file as well
  if (species_s == "vgp_taegut") {
    fig.1.a_d_legendpath <-paste0("/disk1/juwan_d1/cat_rerun/1_Main_figures/Fig1/fig1.a_d.scaffold_aligned_plot_missing.legend.png")
    fig.1.a_d_legend = gtable_filter(ggplotGrob(fig.1.a_d), "guide-box")
    ggsave(fig.1.a_d_legendpath, fig.1.a_d_legend, width=3, height=3)}
}

# input file paths
zebra_finch_path_s="/disk1/juwan_d1/cat_rerun/zebra_finch/tm/CAT2MISSING_pre_taegut/chromosomeMap/out_correlation/"
anna_hummingbird_path_s="/disk1/juwan_d1/cat_rerun/anna/tm/CAT2MISSING_pre_calann_backup/chromosomeMap/out_correlation/"
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
