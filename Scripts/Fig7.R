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
shared_theme <- function(y_breaks) {
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
             axis.title = element_text(size=6, face='bold'))
  )
}

species_order <- c("Zebra finch", "Anna's hummingbird", "Platypus", "Climbing perch")
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
zf_path <- "./input/zebra_finch/"
pl_path <- "./input/platypus/"
an_path <- "./input/anna/"
cl_path <- "./input/climbing_perch/"

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
ZF_ref_f$species <- "Zebra finch"
AN_ref_f$species <- "Anna's hummingbird"
PL_ref_f$species <- "Platypus"
CL_ref_f$species <- "Climbing perch"

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
ZF_tgt_f$species <- "Zebra finch"
AN_tgt_f$species <- "Anna's hummingbird"
PL_tgt_f$species <- "Platypus"
CL_tgt_f$species <- "Climbing perch"
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
ZF_prev_f$species <- "Zebra finch"
AN_prev_f$species <- "Anna's hummingbird"
PL_prev_f$species <- "Platypus"
prev_df <- rbind(ZF_prev_f,PL_prev_f,AN_prev_f)
prev_df <- as.data.frame(prev_df)
prev_df$mode <- "Prior"

# (3) merge all dataframes into a single dataframe: for previous annotation, manual parsing would be performed
merge_df <- rbind(ref_df,tgt_df,prev_df)
merge_df$real_length <- sum(merge_df$numA, merge_df$numC, merge_df$numG, merge_df$numT)
merge_df$real_gc <- ifelse(merge_df$numN != merge_df$length, as.numeric(merge_df$GCcontent)*100, 0) # for the N filled regions

# (4) filter out some records from merged dataframe
flt_merge_df <- merge_df %>%
  filter(excepted == "TRUE") %>% # should not have exception flags
  filter(numN != length) %>% # should not be filled by N
  group_by(position,species, general_type,mode) %>% 
  summarise_at(vars("real_gc"),funs(mean,sd,se=sd(.)/sqrt(n())))
flt_merge_df$GCcontent <-flt_merge_df$mean
flt_merge_df$species <- factor(flt_merge_df$species, levels=c("Zebra finch","Anna's hummingbird","Platypus","Climbing perch"))
write.csv(flt_merge_df, "/disk1/juwan_d1/cat_rerun/summary.csv")

# 1-3. GC content and aligned ratio of previous/cat/VGP annotations
gene_position_plot_maker <- function(summary_df) {
  # (4-1) upstream 3kb region:  +-3kb region centered from TSS
  upstream_3kb_region <- c(upstream_3kb_order,upstream_3kb_inside_order)
  upstream_3kb_df <- summary_df %>%
    filter(position %in% upstream_3kb_region)
  upstream_3kb_df$position_rename <- as.numeric(factor(upstream_3kb_df$position, level=upstream_3kb_region))
  upstream_3kb_labels <- c("-3","TSS","+3")
  upstream_3kb_df$general_type <- "upstream"
  
  # (4-2) upstream 3kb region:  +-3kb region centered from TTS
  downstream_3kb_region <- c(downstream_3kb_inside_order, downstream_3kb_order)
  downstream_3kb_df <- summary_df %>%
    filter(position %in% downstream_3kb_region)
  downstream_3kb_df$position_rename <- as.numeric(factor(downstream_3kb_df$position, level=downstream_3kb_region)) + 90
  downstream_3kb_labels <- c("  -3","TTS","+3 (kbp)          ")
  downstream_3kb_df$general_type <- "downstream"
  
  # (5) upstream 3kb region:  +-3kb region centered from TTS
  up_n_down_3kb_df <- rbind(upstream_3kb_df, downstream_3kb_df)
  up_n_down_3kb_df$mode <- factor(up_n_down_3kb_df$mode, levels=assembly_order)
  up_n_down_3kb_df$general_type <- factor(up_n_down_3kb_df$general_type, levels=c("upstream", "downstream"))
  y_breaks = c(30,40,50,60,70)
  print(levels(factor(up_n_down_3kb_df$position_rename)))

  up_n_down_3kb_plot <- ggplot(up_n_down_3kb_df, aes(x=position_rename, y=GCcontent, color=mode)) + 
    facet_grid(species~general_type, scales="free") + 
    geom_line(data=subset(up_n_down_3kb_df, general_type=="upstream"), aes(group=mode, color=mode), lwd=0.2) +
    geom_line(data=subset(up_n_down_3kb_df, general_type=="downstream"), aes(group=mode, color=mode), lwd=0.2) +
    geom_vline(data=filter(up_n_down_3kb_df, general_type=="upstream"), aes(xintercept = 30.5), linetype="dashed", color = "grey", size=0.5) +
    geom_vline(data=filter(up_n_down_3kb_df, general_type=="downstream"), aes(xintercept = 120.5), linetype="dashed", color = "grey", size=0.5) +
    scale_x_continuous(breaks=c(1,30.5, 60, 91, 120.5, 150), labels=c(upstream_3kb_labels, downstream_3kb_labels), name=" ") + 
    scale_color_manual(values=c("#0043c9","#e68200","#00c95e"))+
    coord_cartesian(ylim=c(30,75),expand=FALSE)+
    shared_theme(y_breaks) + 
    theme(plot.title=element_text(size=8,hjust=0))+
    labs(color="Annotation", x=" ", y="GC-content(%)")
  print(up_n_down_3kb_plot)
    
  Fig7.legend = gtable_filter(ggplotGrob(up_n_down_3kb_plot), "guide-box")
  ggsave("/disk1/juwan_d1/cat_rerun/1_Main_figures/Fig7/fig7.a-d.new.png", up_n_down_3kb_plot, width=2.4, height=5,units = "in")
  ggsave("/disk1/juwan_d1/cat_rerun/1_Main_figures/Fig7/Fig7.a-d.legend.png", Fig7.legend, width=1, height=2,units = "in")
  }

all_plot <- gene_position_plot_maker(summary_df=flt_merge_df)


