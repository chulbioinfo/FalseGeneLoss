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

############################
######## Fig.7.a-d #########
############################

# get the summary file of ref (VGP), tgt (CAT), and prev (Previous) annotations
# shared theme
shared_theme <- function(my_col, y_breaks) {
  list(theme_pubr(),
       scale_y_continuous(breaks = y_breaks),
       theme(plot.title = element_text(size=7,face='bold', hjust = 0.5, margin=margin(0,0,0,0)),
             strip.text = element_blank(),
             strip.background = element_blank(),
             axis.line = element_line(size=0.1),
             legend.position = "bottom",
             legend.key.size = unit(3, "mm"),
             legend.direction = "horizontal",
             legend.text = element_text(size=6),
             legend.title = element_text(size=6),
             legend.margin=margin(c(0,0,0,0)),
             axis.text.y = element_text(size=6),
             axis.text.x =  element_text(size=6,hjust=0.5,vjust=0.4),
             plot.margin = margin(t = 1, r = 4, b = 1, l = 1, unit = "pt"),
             axis.ticks = element_line(size=0.1),
             axis.title = element_text(size=6, face='bold')),
       guides(linetype=guide_legend(ncol = 2, direction ="horizontal", byrow = TRUE, keyheight = unit(3, "mm")),
              shape=guide_legend(ncol = 2, direction ="horizontal", byrow = TRUE,keyheight = unit(3, "mm")),
              color = guide_legend(ncol = 2, direction ="horizontal", byrow = TRUE, keyheight = unit(3, "mm")))
  )
}
# function for order introns
transform_intron_position <- function(column) {
  column =   ifelse(column == 1, 1.5, # 5'UTR Intron
                    ifelse(column == 2, 2.5, # First Intron
                           ifelse(column == 3, 3,   # Internal Intron
                                  ifelse(column == 4, 3.5, # Last Intron
                                         ifelse(column == 5, 4.5, 100))))) # 3'UTR Intron
  return(column)
}

# gene body lines
draw_gene_body_lines <- list(
  geom_vline(xintercept = 1.5, linetype="dashed", color = "grey",    size=0.1),
  geom_vline(xintercept = 2.5, linetype="dashed", color = "grey",    size=0.1),
  geom_vline(xintercept = 3.5, linetype="dashed", color = "grey",    size=0.1),
  geom_vline(xintercept = 4.5, linetype="dashed", color = "grey",    size=0.1),
  geom_vline(xintercept = 1,   linetype="solid",  color = "#5e5e5e", size=0.1),
  geom_vline(xintercept = 2,   linetype="solid",  color = "#5e5e5e", size=0.1),
  geom_vline(xintercept = 3,   linetype="solid",  color = "#5e5e5e", size=0.1),
  geom_vline(xintercept = 4,   linetype="solid",  color = "#5e5e5e", size=0.1),
  geom_vline(xintercept = 5,   linetype="solid",  color = "#5e5e5e", size=0.1),
  scale_x_continuous(name="Intron", breaks = c(1.5, 2.5, 3, 3.5, 4.5), labels=c("5'", "F", "I", "L", "3'"),
                     sec.axis=sec_axis(~. + 0, name="Exon", breaks = c(1,2,3,4,5), labels=c("5'", "F", "I", "L", "3'")))
)

species_order <- c("zebra finch", "anna's hummingbird", "platypus", "climbing perch")
utr_n_cds_order <- c("five-prime-UTR", "FirstCDS",  "InternalCDS", "LastCDS", "three-prime-UTR")
intron_order <- c("5'UTR-intron","First-intron","Internal-intron", "Last-intron","3'UTR-intron")
genic_order <- c("five-prime-UTR", "5'UTR-intron", "FirstCDS","First-intron",  "InternalCDS", "Internal-intron","LastCDS","Last-intron", "three-prime-UTR", "3'UTR-intron")
assembly_order <- c("VGP", "Projected", "Prior")

upstream_3kb_order <- c()
for (i in c(1:30)) {
  upstream_3kb_order <- c(upstream_3kb_order, paste0("five-prime-", 3000 - 100*(i-1), "bp"))}
downstream_3kb_order <- c()
for (i in c(1:30)) {
  downstream_3kb_order <- c(downstream_3kb_order, paste0("three-prime-", 100*(i), "bp"))}
upstream_3kb_inside_order <- c()
for (i in c(1:30)) {
  upstream_3kb_inside_order <- c(upstream_3kb_inside_order, paste0("five-inside-", 100*(i), "bp"))
}
downstream_3kb_inside_order <- c()
for (i in c(1:30)) {
  downstream_3kb_inside_order <- c(downstream_3kb_inside_order, paste0("three-inside-", 3100-100*(i), "bp"))
}

# (1) load input dataframes
zf_path <- "/disk1/juwan_d1/cat_rerun/zebra_finch/tm/CAT2MISSING_pre_taegut/gene_structure/"
an_path <- "/disk1/juwan_d1/cat_rerun/platypus/tm/CAT2MISSING_pre_ornana_backup/gene_structure/"
pl_path <- "/disk1/juwan_d1/cat_rerun/anna/tm/CAT2MISSING_pre_calann/gene_structure/"
cl_path <- "/disk1/juwan_d1/cat_rerun/climbing_perch/tm/CAT2MISSING_pre_anates/gene_structure/"

# (1-1) ref
mode_s = "ref"
zf_ref_path <- paste0(zf_path, mode_s, "/", mode_s, ".gene_structure.tabs")
an_ref_path <- paste0(an_path, mode_s, "/", mode_s, ".gene_structure.tabs")
pl_ref_path <- paste0(pl_path, mode_s, "/", mode_s, ".gene_structure.tabs")
cl_ref_path <- paste0(cl_path, mode_s, "/", mode_s, ".gene_structure.tabs")

ZF_ref_f = vroom(zf_ref_path)
AN_ref_f = vroom(an_ref_path)
PL_ref_f = vroom(pl_ref_path)
CL_ref_f = vroom(cl_ref_path)
ref_df <- rbind(ZF_ref_f,PL_ref_f,AN_ref_f,CL_ref_f)
ref_df <- as.data.frame(ref_df)
ref_df$mode <- "VGP"

# (1-2) tgt (projected by CAT)
mode_s = "tgt"
zf_tgt_path <- paste0(zf_path, mode_s, "/", mode_s, ".gene_structure.tabs")
an_tgt_path <- paste0(an_path, mode_s, "/", mode_s, ".gene_structure.tabs")
pl_tgt_path <- paste0(pl_path, mode_s, "/", mode_s, ".gene_structure.tabs")
cl_tgt_path <- paste0(cl_path, mode_s, "/", mode_s, ".gene_structure.tabs")

ZF_tgt_f = vroom(zf_tgt_path)
AN_tgt_f = vroom(an_tgt_path)
PL_tgt_f = vroom(pl_tgt_path)
CL_tgt_f = vroom(cl_tgt_path)
tgt_df <- rbind(ZF_tgt_f,PL_tgt_f,AN_tgt_f,CL_tgt_f)
tgt_df <- as.data.frame(tgt_df)
tgt_df$mode <- "Projected"

# (1-3) prior
mode_s = "prev"
zf_prev_path <- paste0(zf_path, mode_s, "/", mode_s, ".gene_structure.tabs")
an_prev_path <- paste0(an_path, mode_s, "/", mode_s, ".gene_structure.tabs")
pl_prev_path <- paste0(pl_path, mode_s, "/", mode_s, ".gene_structure.tabs")
ZF_prev_f = vroom(zf_prev_path)
AN_prev_f = vroom(an_prev_path)
PL_prev_f = vroom(pl_prev_path)
prev_df <- rbind(ZF_prev_f,PL_prev_f,AN_prev_f)
prev_df <- as.data.frame(prev_df)
prev_df$mode <- "Projected"

# (2) manually rename intron records
rename_introns <- function(df){
  df$position <- ifelse(df$position %in% c("five-prime-UTR_FirstCDS_intron",
                                           "five-prime-UTR_five-prime-UTR_intron"),"5'UTR-intron",
                        ifelse(df$position %in% c("LastCDS_three-prime-UTR_intron",
                                                  "FirstCDS_three-prime-UTR_intron",
                                                  "three-prime-UTR_three-prime-UTR_intron"), "3'UTR-intron", 
                               ifelse(df$position  %in% c("FirstCDS_InternalCDS_intron",
                                                          "FirstCDS_LastCDS_intron"),"First-intron", 
                                      ifelse(df$position  %in% c("InternalCDS_InternalCDS_intron"), "Internal-intron",
                                             ifelse(df$position  %in% c("InternalCDS_LastCDS_intron"),"Last-intron", 
                                                    df$position)))))
  return(df)
}
ref_df <-rename_introns(ref_df)
tgt_df <-rename_introns(tgt_df)
prev_df <-rename_introns(prev_df)

colnames(ref_df)
colnames(tgt_df)
colnames(prev_df)

# (3) merge all dataframes into a single dataframe: for previous annotation, manual parsing would be performed
my_c <-colnames(ref_df)
prev_df_parse <- prev_df[,1:29]
colnames(prev_df_parse) <- my_c
merge_df <- rbind(ref_df,ALL_tgt_df,prev_df_parse)
merge_df$real_length <- sum(merge_df$numA,merge_df$numC,merge_df$numG,merge_df$numT)
merge_df$real_gc <- ifelse(merge_df$numN != merge_df$length, as.numeric(merge_df$GCcontent)*100, 0) # for the N filled regions

# (4) filter out some records from merged dataframe
flt_merge_df <- merge_df %>%
  filter(excepted == "TRUE") %>% # should not have exception flags
  filter(numN != length) %>% # should not be filled by N
  group_by(position,species, general_type,mode) %>% 
  summarise_at(vars("real_gc"),funs(mean,sd,se=sd(.)/sqrt(n())))
flt_merge_df$GCcontent <-flt_merge_df$mean
flt_merge_df$species <- factor(flt_merge_df$species, levels=c("zebrafinch","anna","platypus","perch"))


# 1-3. GC content and aligned ratio of previous/cat/VGP annotations
gene_position_plot_maker <- function(summary_df, species_s) {
  species_df <- summary_df %>%
    filter(species == species_s)
  species_df$mode <- factor(species_df$mode, levels=assembly_order)
  
  # (1) UTR and CDS
  utr_n_cds <-species_df %>%
    filter(position %in% utr_n_cds_order)
  utr_n_cds$position_rename <- as.numeric(factor(utr_n_cds$position, level=utr_n_cds_order))
  
  # (2) Intron
  introns <-species_df %>%
    filter(position %in% intron_order)
  introns$position_rename <- as.numeric(factor(introns$position, level=intron_order))
  introns$position_rename <- transform_intron_position(introns$position_rename)
  introns$position_rename <- as.numeric(introns$position_rename)
  utr_n_cds$position_rename <- as.numeric(utr_n_cds$position_rename)
  
  # (3) (UTR & CDS) + (Intron) --> gene body
  utr_n_cds$type <- "exon"
  introns$type <- "intron"
  utr_cds_introns <- rbind(utr_n_cds, introns)
  utr_cds_introns$type <- factor(utr_cds_introns$type, level=c("exon","intron"))
  
  # (4-1) upstream 3kb region:  +-3kb region centered from TSS
  upstream_3kb_region <- c(upstream_3kb_order,upstream_3kb_inside_order)
  upstream_3kb_df <- species_df %>%
    filter(position %in% upstream_3kb_region)
  upstream_3kb_df$position_rename <- as.numeric(factor(upstream_3kb_df$position, level=upstream_3kb_region))
  upstream_3kb_labels <- c("-3","TSS","+3")
  
  # (4-2) upstream 3kb region:  +-3kb region centered from TTS
  downstream_3kb_region <- c(downstream_3kb_inside_order, downstream_3kb_order)
  downstream_3kb_df <- species_df %>%
    filter(position %in% downstream_3kb_region)
  downstream_3kb_df$position_rename <- as.numeric(factor(downstream_3kb_df$position, level=downstream_3kb_region))
  downstream_3kb_labels <- c("-3","TTS","+3")

  # (5) upstream 3kb region:  +-3kb region centered from TTS
  utr_cds_introns$mode <- factor(utr_cds_introns$mode, levels=assembly_order)
  upstream_3kb_df$mode <- factor(upstream_3kb_df$mode, levels=assembly_order)
  downstream_3kb_df$mode <- factor(downstream_3kb_df$mode, levels=assembly_order)
  
  gene_body_plot <- ggplot(utr_cds_introns, aes(x =position_rename, y = GCcontent, color=mode)) +
    draw_gene_body_lines+
    geom_line(aes(group=interaction(mode,type),linetype = type, color=mode), lwd=0.2) +
    geom_point(aes(shape = type, size=type),size=0.3) +
    scale_x_continuous(name="Intron", breaks = c(1.5, 2.5, 3, 3.5, 4.5), labels=c("5'", "F", "I", "L", "3'"),
                       sec.axis=sec_axis(~. + 0, name="Exon", breaks = c(1,2,3,4,5), labels=c("5'", "F", "I", "L", "3'"))) + 
    scale_shape_manual(values=c(19,1), labels=c("Exon or upstream/downstream sequence", "Intron"))+
    scale_size_manual(values=c(1.5,2.5), labels=c("Exon or upstream/downstream sequence", "Intron"))+
    scale_color_manual(values=c("#0043c9","#e68200","#00c95e"))+
    scale_linetype(labels=c("Exon or upstream/downstream sequence", "Intron"))+
    coord_cartesian(ylim=c(30,75), xlim=c(0.5,5.5), expand=FALSE)+#, expand=FALSE)+#coord_cartesian(ylim=c(30,75), expand=FALSE)+
    labs(color="Color", title = "Gene body", linetype="Type", shape="Type", size="Type") + 
    shared_theme
  
  upstream_3kb_plot <- ggplot(upstream_3kb_df, aes(x =position_rename, y = GCcontent, color=mode)) + 
    geom_line(aes(group=mode, color=mode), lwd=0.2) +
    geom_vline(xintercept = 30.5, linetype="dashed", color = "grey", size=0.5) +
    scale_x_continuous(breaks=c(1,30.5,60), labels=up_2kb_labels, name=" ", #trans=reverselog_trans(10),
                       sec.axis=sec_axis(~. + 0, name=" ", breaks = c(1), labels=c(" "))) + 
    scale_color_manual(values=c("#0043c9","#e68200","#00c95e"))+
    coord_cartesian(ylim=c(30,75),expand=FALSE)+#coord_cartesian(ylim=c(30,75), expand=FALSE)+
    labs(color="Color", title="Upstream", x=" ") + 
    shared_theme
  
  downstream_3kb_plot <- ggplot(downstream_3kb_df, aes(x =position_rename, y = GCcontent, color=mode)) + 
    geom_line(aes(group=mode, color=mode), lwd=0.2) +
    geom_vline(xintercept = 30.5, linetype="dashed", color = "grey", size=0.5) +
    scale_x_continuous(breaks=c(1,30.5,60), labels=down_2kb_labels, name=" ",
                       sec.axis=sec_axis(~. + 0, name=" ", breaks = c(1), labels=c(" "))) + 
    scale_color_manual(values=c("#0043c9","#e68200","#00c95e"))+
    coord_cartesian(ylim=c(30,75),expand=FALSE)+
    labs(color="Color", title="Downstream", x="(Kbp)") + 
    shared_theme
  
  my_list <- list(upstream_3kb_plot + theme(legend.position = "none",axis.title.y=element_blank()),
                  gene_body_plot + theme(legend.position = "none",axis.title.y=element_blank(),axis.text.y=element_blank()),
                  downstream_3kb_plot + theme(legend.position = "none",axis.title.y=element_blank(),axis.text.y=element_blank()))
  if (species_s == "zebrafinch") {
    Fig7.legend = gtable_filter(ggplotGrob(gene_body_plot), "guide-box")
    ggsave("/disk1/juwan_d1/cat_rerun/1_Main_figures/Fig7/Fig7.e.legend.png", Fig7.legend, width=6, height=2,units = "in")
  }
  return(my_list)
}

finch_plot <- gene_position_plot_maker(summary_df=all_summary_df, species_s="zebrafinch")
anna_plot <- gene_position_plot_maker(summary_df=all_summary_df, species_s="anna")
platypus_plot <- gene_position_plot_maker(summary_df=all_summary_df, species_s="platypus")
perch_plot <- gene_position_plot_maker(summary_df=all_summary_df, species_s="perch")

my_all_plot <- c(finch_plot, anna_plot, platypus_plot, perch_plot)
length(my_all_plot)
my_layout <- rbind(c(1,1,1,1,1,1,
                     2,2,2,2,2,
                     3,3,3,3,3,
                     4,4,4,4,4,4,
                     5,5,5,5,5,
                     6,6,6,6,6),
                   c(7,7,7,7,7,7,
                     8,8,8,8,8,
                     9,9,9,9,9,
                     10,10,10,10,10,10,
                     11,11,11,11,11,
                     12,12,12,12,12))

Fig7.e.annot_comparison_plot <- arrangeGrob(grid.arrange(grobs = my_all_plot, layout_matrix = my_layout))
ggsave("/disk1/juwan_d1/cat_rerun/1_Main_figures/Fig7/fig7.a-d.png", Fig7.e.annot_comparison_plot, width=5.7, height=3.5,units = "in")
