library(ggpubr)
library(ggrepel)
library(ggplot2)
library(gridExtra)
library(reshape2)

##########################################################
##### Shared variables for supplementary figures ######
##########################################################

zf_path <- "/disk1/juwan_d1/cat_rerun/zebra_finch/tm/CAT2MISSING_pre_taegut/chromosomeMap/out_correlation/vgp_taegut.stat_summary.tabs"
an_path <- "/disk1/juwan_d1/cat_rerun/anna/tm/CAT2MISSING_pre_calann/chromosomeMap/out_correlation/vgp_calann.stat_summary.tabs"
pl_path <- "/disk1/juwan_d1/cat_rerun/platypus/tm/CAT2MISSING_pre_ornana/chromosomeMap/out_correlation/vgp_ornana.stat_summary.tabs"
cl_path <- "/disk1/juwan_d1/cat_rerun/climbing_perch/tm/CAT2MISSING_pre_anates/chromosomeMap/out_correlation/vgp_anates.stat_summary.tabs"
species_order <- c("Zebra finch", "Anna's hummingbird", "Platypus", "Climbing perch")

# load input files
zf_df <- read.csv(zf_path, sep="\t")
an_df <- read.csv(an_path, sep="\t")
pl_df <- read.csv(pl_path, sep="\t")
cl_df <- read.csv(cl_path, sep="\t")
zf_df$species <- "Zebra finch"
an_df$species <- "Anna's hummingbird"
pl_df$species <- "Platypus"
cl_df$species <- "Climbing perch"

# remove chromosome Y of platypus
pl_df <- pl_df %>%
  filter(!grepl("^chr.Y", scaffold))

# manual replacement of chr 29 of zebra finch
# Size
zf_df[32,3] <- 1995839
# GC content
zf_df[32,4] <- 58.6618434516287
# repeat content
zf_df[32,5] <- 31.93894898
# missing ratio
zf_df[32,11] <- 39.4  


#########################################
##### Supplementary Fig.S1.panel a ######
#########################################

sup.fig.s1.a.maker <- function(input_df, species_s) {
  input_df$type <-ifelse(input_df$type == "chromosome", "chromosome",
                          ifelse(input_df$type == "unlocalized_scaffold", "unlocalized","unplaced"))
  input_df$type <- factor(input_df$type, levels=c("chromosome","unplaced", "unlocalized"))
  #####
  df_subset <- input_df %>%
    filter(size>100000)
  #####
  melt.df_subset <- melt(df_subset, id.vars=c("scaffold", "type", "missing_ratio"),
                         measure.vars=c("scaffold_gc","scaffold_repeat"))
  sup.fig.S1.plot <- ggscatter(melt.df_subset, x="missing_ratio", y="value", color="variable", palette="variable", size=0.5, alpha=0.2) + 
    stat_cor(aes(color = variable), label.x.npc=c("center"), label.y = c(2,12), method="spearman", cex=1.5, hjust=0.5) + 
    scale_shape_manual(values=c(6,4,13))+
    coord_cartesian(xlim=c(0,100), ylim=c(0,100)) + 
    theme_pubr()+
    theme(text=element_text(size=6), 
          axis.text=element_text(size=6),
          axis.title=element_blank(),
          strip.background = element_blank(),
          strip.text=element_text(size=6),
          legend.position = "bottom",
          legend.box = "vertical",
          legend.box.spacing = unit(0.01, "mm"),
          legend.text= element_text(size=6)) +
    labs(x="Missing ratio (%)", y="Percentage (%)", color="Type", shape="Scaffold type") 
  if (species_s =="Zebra finch") {
    sup.fig.S1.plot_legend = gtable_filter(ggplotGrob(sup.fig.S1.plot), "guide-box")
    sup.fig.S1.plot_legend <- arrangeGrob(sup.fig.S1.plot_legend)
    sup.fig.S1.plot_legend_path <- "/disk1/juwan_d1/cat_rerun/1_Main_figures/Sup/sup.fig.S1.legend.png"
    ggsave(sup.fig.S1.plot_legend_path, sup.fig.S1.plot_legend, width=3, height=1,units = "in")
  }
  my_plot <- list(sup.fig.S1.plot + theme(legend.position = "none"))
  return(my_plot)
}

finch_s1.a.plot <- sup.fig.s1.a.maker(input_df=zf_df, "Zebra finch")
anna_s1.a.plot <- sup.fig.s1.a.maker(input_df=an_df, "Anna's hummingbird")
platypus_s1.a.plot <- sup.fig.s1.a.maker(input_df=pl_df, "Platypus")
perch_s1.a.plot <- sup.fig.s1.a.maker(input_df=cl_df, "Climbing perch")

sup.fig.S1.a.list <- c(finch_s1.a.plot, anna_s1.a.plot, platypus_s1.a.plot, perch_s1.a.plot)
my_layout<-rbind(c(1,2,3,4))
sup.fig.S1.a.plots <- arrangeGrob(grobs = sup.fig.S1.a.list, layout_matrix = my_layout)
sup.fig.S1.a.path <-"/disk1/juwan_d1/cat_rerun/1_Main_figures/Sup/sup.fig.S1.a.png"
ggsave(sup.fig.S1.a.path, sup.fig.S1.a.plots, width=5.5, height=1.5)

#########################################
#####              done            ######
#########################################

#########################################
##### Supplementary Fig.S1.panel b ######
#########################################

rename_newly_identified <- function(scaffold_df, chrs_c) {
  scaffold_df$scaffold <- ifelse(scaffold_df$scaffold %in% chrs_c, paste0("*",scaffold_df$scaffold), scaffold_df$scaffold)
  #scaffold_df$scaffold <- ifelse(scaffold_df$scaffold %in% chrs_c, gsub("up", "*up", scaffold_df$scaffold), scaffold_df$scaffold)
  return(scaffold_df)
}

sup.fig.s1.b.maker <- function(input_df) {
  input_df$type <-ifelse(input_df$type == "chromosome", "chromosome",
                         ifelse(input_df$type == "unlocalized_scaffold", "unlocalized","unplaced"))
  input_df$type <- factor(input_df$type, levels=c("chromosome","unplaced", "unlocalized"))
  #####
  df_subset <- input_df %>%
    filter(size>100000)
  #####
  df_subset <- df_subset[order(-df_subset$size),]
  unplaced_scaffolds <- df_subset%>%
    filter(type=="unplaced")
  unplaced_scaffolds <- unplaced_scaffolds$scaffold
  unplaced_scaffolds_name_dict <- list()
  int_start=0
  
  for (i in unplaced_scaffolds) {
    int_start <- int_start + 1
    unplaced_scaffolds_name_dict[i] = int_start}
  
  df_subset$scaffold <- ifelse(df_subset$scaffold  %like% "^chr.", 
                                      gsub("chr.","", df_subset$scaffold),
                                      paste0("up",unplaced_scaffolds_name_dict[df_subset$scaffold]))

  newly_identified_chrs_c <- subset(df_subset$scaffold, df_subset$missing_ratio >= 30)
  print(newly_identified_chrs_c)
  df_subset <- rename_newly_identified(scaffold_df=df_subset, chrs_c=newly_identified_chrs_c)

  df_subset$gene_density <- df_subset$num_of_genes/df_subset$size * 1000000
  average_gene_density <- sum(df_subset$num_of_genes) / sum(df_subset$size) * 1000000

  scaffold_order <- df_subset$scaffold
  df_subset$scaffold <- factor(df_subset$scaffold, levels=scaffold_order)
  sup.fig.s1.b.plot <- ggplot(df_subset) +
    geom_bar(aes(x=scaffold, y=gene_density, fill=missing_ratio), width=0.7, color="black", lwd=0.2, stat = "summary", fun = "mean") + 
    scale_fill_gradientn(colors = brewer.pal(9, "Greens"), limits=c(0,100), name="Missing ratio   ")+
    ###
    geom_hline(aes(yintercept=average_gene_density), lwd=0.2, color="black", linetype="dashed")+
    theme(plot.title=element_text(size=6), axis.text.x=element_text(size=6, angle=90), 
          axis.text.y=element_text(size=6, angle=0), axis.title=element_blank(),
          legend.text = element_text(size=6),
          legend.direction = "horizontal", 
          legend.key.width = unit(0.6,"cm"),
          legend.key.size = unit(0.2, "cm"),
          legend.key.height = unit(0.2,"cm"),) + 
    coord_cartesian(expand=FALSE)
  sup.fig.s1.b.legend = gtable_filter(ggplotGrob(sup.fig.s1.b.plot), "guide-box")
  sup.fig.s1.b.legend <- arrangeGrob(sup.fig.s1.b.legend)
  sup.fig.s1.b.legend_path <- "/disk1/juwan_d1/cat_rerun/1_Main_figures/Sup/sup.fig.S1.b.legend.png"
  ggsave(sup.fig.s1.b.legend_path, sup.fig.s1.b.legend, width=5, height=1,units = "in")
  sup.fig.s1.b.plot_list <- list(sup.fig.s1.b.plot + theme(legend.position = "none"))
  return(sup.fig.s1.b.plot_list)
}

finch_s1.b.plot <- sup.fig.s1.b.maker(input_df=zf_df)
anna_s1.b.plot <- sup.fig.s1.b.maker(input_df=an_df)
platypus_s1.b.plot <- sup.fig.s1.b.maker(input_df=pl_df)
perch_s1.b.plot <- sup.fig.s1.b.maker(input_df=cl_df)

sup.fig.s1.b.list <- c(finch_s1.b.plot, anna_s1.b.plot, platypus_s1.b.plot, perch_s1.b.plot)
my_layout <- rbind(c(1),c(2),c(3),c(4))
sup.fig.s1.b.plots <- arrangeGrob(grid.arrange(grobs = sup.fig.s1.b.list, layout_matrix = my_layout))
sup.fig.S1.b.path <-"/disk1/juwan_d1/cat_rerun/1_Main_figures/Sup/sup.fig.S1.b.png"
ggsave(sup.fig.S1.b.path, sup.fig.s1.b.plots, width=5.5, height=3.5, units = "in")

#########################################
#####              done            ######
#########################################