library(ggpubr)
library(ggrepel)
library(ggplot2)
library(gridExtra)
library(reshape2)

##########################################################
##### Shared variables for supplementary fig.S1,2,3 ######
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

#################################
##### Supplementary Fig.S1 ######
#################################

sup.fig.S1.maker <- function(df_subset, species_s) {
  df_subset$type <-ifelse(df_subset$type == "chromosome", "chromosome",
                                  ifelse(df_subset$type == "unlocalized_scaffold", "unlocalized","unplaced"))
  df_subset$type <- factor(df_subset$type, levels=c("chromosome","unplaced", "unlocalized"))
  melt.df_subset <- melt(df_subset, id.vars=c("scaffold", "type", "unaligned_ratio"),
                                 measure.vars=c("scaffold_gc","scaffold_repeat"))
  sup.fig.S1.plot <- ggscatter(melt.df_subset, x="unaligned_ratio", y="value", color="variable", palette="variable", size=0.5, alpha=0.2) + 
    stat_cor(aes(color = variable), label.x.npc=c("center"), label.y = c(2,12), method="pearson", cex=1.5) + 
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
    labs(x="Missing ratio (%)", y="Percentage (%)", color="Type", shape="Scaffold type", title=species_s) 
  if (species_s =="Zebra finch") {
    sup.fig.S1.plot_legend = gtable_filter(ggplotGrob(sup.fig.S1.plot), "guide-box")
    sup.fig.S1.plot_legend <- arrangeGrob(sup.fig.S1.plot_legend)
    sup.fig.S1.plot_legend_path <- "/disk1/juwan_d1/cat_rerun/1_Main_figures/Sup/sup.fig.S1.legend.png"
    ggsave(sup.fig.S1.plot_legend_path, sup.fig.S1.plot_legend, width=3, height=1,units = "in")
  }
  my_plot <- list(sup.fig.S1.plot + theme(legend.position = "none"))
  return(my_plot)
}

finch_S1_plot <- sup.fig.S1.maker(df_subset=zf_df, "Zebra finch")
anna_S1_plot <- sup.fig.S1.maker(df_subset=an_df, "Anna's hummingbird")
platypus_S1_plot <- sup.fig.S1.maker(df_subset=pl_df, "Platypus")
perch_S1_plot <- sup.fig.S1.maker(df_subset=cl_df, "Climbing perch")

sup.fig.S1.list <- c(finch_S1_plot, anna_S1_plot, platypus_S1_plot, perch_S1_plot)
my_layout<-rbind(c(1,2),
                 c(3,4))
sup.fig.S1.plots <- arrangeGrob(grobs = sup.fig.S1.list, layout_matrix = my_layout)
sup.fig.S1.path <-"/disk1/juwan_d1/cat_rerun/1_Main_figures/Sup/sup.fig.S1.png"
ggsave(sup.fig.S1.path, sup.fig.S1.plots, width=3.7, height=4)

#########################################
###### Supplementary Fig.S1 (done) ######
#########################################

#################################
##### Supplementary Fig.S2 ######
#################################

sup.fig.S2.maker <- function(df_subset, species_s) {
  df_subset$type <-ifelse(df_subset$type == "chromosome", "chromosome",
                          ifelse(df_subset$type == "unlocalized_scaffold", "unlocalized","unplaced"))
  df_subset$type <- factor(df_subset$type, levels=c("chromosome","unplaced", "unlocalized"))
  melt.df_subset <- melt(df_subset, id.vars=c("scaffold", "type", "size"),
                         measure.vars=c("scaffold_gc","scaffold_repeat"))
  sup.fig.S2.plot <- ggscatter(melt.df_subset, x="size", y="value", color="variable", palette="variable", size=0.5, alpha=0.2) + 
    stat_cor(aes(color = variable), label.x.npc=c("center"), label.y = c(2,12), method="spearman", cex=1.5) + 
    scale_shape_manual(values=c(6,4,13))+
    coord_cartesian(ylim=c(0,100)) + 
    scale_x_log10()+
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
    labs(x="Scaffold size (bp)", y="Percentage (%)", color="Type", shape="Scaffold type", title=species_s) 
  if (species_s =="Zebra finch") {
    sup.fig.S2.plot_legend = gtable_filter(ggplotGrob(sup.fig.S2.plot), "guide-box")
    sup.fig.S2.plot_legend <- arrangeGrob(sup.fig.S2.plot_legend)
    sup.fig.S2.plot_legend_path <- "/disk1/juwan_d1/cat_rerun/1_Main_figures/Sup/sup.fig.S2.legend.png"
    ggsave(sup.fig.S2.plot_legend_path, sup.fig.S2.plot_legend, width=3, height=1,units = "in")
  }
  my_plot <- list(sup.fig.S2.plot + theme(legend.position = "none"))
  return(my_plot)
}

finch_S2_plot <- sup.fig.S2.maker(df_subset=zf_df, "Zebra finch")
anna_S2_plot <- sup.fig.S2.maker(df_subset=an_df, "Anna's hummingbird")
platypus_S2_plot <- sup.fig.S2.maker(df_subset=pl_df, "Platypus")
perch_S2_plot <- sup.fig.S2.maker(df_subset=cl_df, "Climbing perch")

sup.fig.S2.list <- c(finch_S2_plot, anna_S2_plot, platypus_S2_plot, perch_S2_plot)
my_layout<-rbind(c(1,2),
                 c(3,4))
sup.fig.S2.plots <- arrangeGrob(grobs = sup.fig.S2.list, layout_matrix = my_layout)
sup.fig.S2.path <-"/disk1/juwan_d1/cat_rerun/1_Main_figures/Sup/sup.fig.S2.png"
ggsave(sup.fig.S2.path, sup.fig.S2.plots, width=3.7, height=4)

#########################################
###### Supplementary Fig.S2 (done) ######
#########################################

#################################
##### Supplementary Fig.S3 ######
#################################

sup.fig.S3.maker <- function(df_subset, species_s) {
  df_subset$type <-ifelse(df_subset$type == "chromosome", "chromosome",
                          ifelse(df_subset$type == "unlocalized_scaffold", "unlocalized","unplaced"))
  df_subset$type <- factor(df_subset$type, levels=c("chromosome","unplaced", "unlocalized"))
  df_subset$missing_type <- ifelse(df_subset$unaligned_ratio >= 50, "missing", "normal")
  df_subset$missing_type <- factor(df_subset$missing_type, levels=c("normal", "missing"))

  sup.fig.S3.plot <- ggplot(df_subset, aes(x=scaffold_gc, y=scaffold_repeat,group=missing_type)) + 
    geom_point(aes(shape=type, color=missing_type),size=0.3,alpha=0.6) + 
    stat_cor(data=df_subset, inherit.aes = FALSE, aes(x=scaffold_gc, y=scaffold_repeat), label.x=c(30), label.y = c(2), method="pearson", cex=1.5) + 
    scale_color_manual(values=c("grey", "brown"))+
    scale_shape_manual(values=c(6,4,13))+
    coord_cartesian(xlim=c(30,75), ylim=c(0,100)) + 
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
    labs(x="GC content (%)", y="Repeat content (%)", color="Missing >50%", shape="Scaffold type", title=species_s) 
  tmp_plot <- sup.fig.S3.plot + theme(legend.position = "none")
  sup.fig.S3.plot_with_marginal<- ggMarginal(tmp_plot, type="density", size=6, groupFill = TRUE, lwd=0.05)
  my_plot <- list(sup.fig.S3.plot_with_marginal)
  if (species_s =="Zebra finch") {
    sup.fig.S3.plot_legend = gtable_filter(ggplotGrob(sup.fig.S3.plot), "guide-box")
    sup.fig.S3.plot_legend_plot <- arrangeGrob(sup.fig.S3.plot_legend)
    sup.fig.S3.plot_legend_path <- "/disk1/juwan_d1/cat_rerun/1_Main_figures/Sup/sup.fig.S3.legend.png"
    ggsave(sup.fig.S3.plot_legend_path, sup.fig.S3.plot_legend_plot, width=3, height=1,units = "in")
  }
  return(my_plot)
}

finch_S3_plot <- sup.fig.S3.maker(df_subset=zf_df, "Zebra finch")
anna_S3_plot <- sup.fig.S3.maker(df_subset=an_df, "Anna's hummingbird")
platypus_S3_plot <- sup.fig.S3.maker(df_subset=pl_df, "Platypus")
perch_S3_plot <- sup.fig.S3.maker(df_subset=cl_df, "Climbing perch")

sup.fig.S3.plots <- c(finch_S3_plot, anna_S3_plot, platypus_S3_plot, perch_S3_plot)
my_layout<-rbind(c(1,2),c(3,4))
sup.fig.S3.plots <- arrangeGrob(grobs = sup.fig.S3.plots, layout_matrix = my_layout)
sup.fig.S3.path <-"/disk1/juwan_d1/cat_rerun/1_Main_figures/Sup/sup.fig.S3.png"
ggsave(sup.fig.S3.path, sup.fig.S3.plots, width=3.7, height=4)

#########################################
###### Supplementary Fig.S3 (done) ######
#########################################

#################################
##### Supplementary Fig.S4 ######
#################################

zf_df$log10_size <- log10(zf_df$size)
zf_df$log10_percent <- (zf_df$size * (100-zf_df$unaligned_ratio))/zf_df$size/100*zf_df$log10_size
zf_df <- zf_df[with(zf_df, order(-size)), ]
unplaced_scaffolds <- zf_df%>%
  filter(type=="unplaced_scaffold")
unplaced_scaffolds <- unplaced_scaffolds$scaffold
unplaced_scaffolds_name_dict <- list()
int_start=0

for (i in unplaced_scaffolds) {
  int_start <- int_start + 1
  unplaced_scaffolds_name_dict[i] = int_start}
zf_df$scaffold <- ifelse(zf_df$scaffold  %like% "^chr.", gsub("chr.","", zf_df$scaffold), paste0("up",unplaced_scaffolds_name_dict[zf_df$scaffold]))
zf_df$type <- factor(zf_df$type, levels=c("chromosome", "unplaced_scaffold", "unlocalized_scaffold"))
zf_df$scaffold_type <- ifelse(zf_df$log10_size>5, "scaffold > 100kbp", "scaffold <= 100kbp")
zf_df$scaffold_type <- factor(zf_df$scaffold_type, levels=c("scaffold > 100kbp", "scaffold <= 100kbp"))
sup.fig.S4 <- ggplot() + 
  facet_wrap(~scaffold_type, nrow=2, scales="free") +
  geom_bar(data=zf_df, mapping=aes(x=reorder(gsub("chr.","", scaffold), -size),y=log10_size, width=.7), fill="grey", stat="identity",lwd=0.2) + 
  geom_line(data=zf_df,mapping=aes(x=reorder(gsub("chr.","", scaffold), -size), y=scaffold_gc/10,group=1), color="red") + 
  geom_line(data=zf_df,mapping=aes(x=reorder(gsub("chr.","", scaffold), -size), y=scaffold_repeat/10,group=1), color="blue") + 
  scale_y_continuous(limits=c(0,8.5), breaks=c(-2,-1,0,5), labels=c("Repeat","GC","0","5"), sec.axis = sec_axis(~.*10, name="GC or repeat content(%)")) + 
  labs(x="", y="Log(Size) (bp)") + 
  theme_classic()+ 
  theme(legend.position="bottom",
        legend.key.width = unit(0.6,"cm"),
        legend.key.size = unit(0.2, "cm"),
        legend.key.height = unit(0.2,"cm"),
        legend.box="vertical",
        legend.box.just = "left",
        text=element_text(size=6),
        axis.text.x=element_text(size=6, color="#6b6b6b", angle=90, hjust=1),
        panel.grid.major.y = element_line()) 

sup.fig.S4_legend_path <- "/disk1/juwan_d1/cat_rerun/1_Main_figures/Sup/sup.fig.S4.legend.png"
sup.fig.S4_legend = gtable_filter(ggplotGrob(sup.fig.S4), "guide-box")
ggsave(sup.fig.S4_legend_path, sup.fig.S4_legend, width=3, height=3)
sup.fig.S4_path <- "/disk1/juwan_d1/cat_rerun/1_Main_figures/Sup/sup.fig.S4.png"
ggsave(sup.fig.S4_path, sup.fig.S4 + theme(legend.position="none"), width=6, height=5)


########################################
##### Supplementary Fig.S4 (done) ######
########################################
