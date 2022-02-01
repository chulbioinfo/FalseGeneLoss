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


##########################
####### Fig.5.a-b ########
##########################

# load input data
input_path <- "/disk1/juwan_d1/cat_rerun/"
average_gc_miss_path <- paste0(input_path, "all.forR.gc_and_indel.tsv")
average_gc_miss_df <- read.csv(average_gc_miss_path, sep="\t")
utr_n_cds_order <- c("five-prime-UTR", "FirstCDS",  "InternalCDS", "LastCDS", "three-prime-UTR")
intron_order <- c("5'UTR-intron","First-intron","Internal-intron", "Last-intron","3'UTR-intron")
genic_order <- c("five-prime-UTR", "5'UTR-intron", "FirstCDS","First-intron",  "InternalCDS", "Internal-intron","LastCDS","Last-intron", "three-prime-UTR", "3'UTR-intron")
species_order <- c("zebra finch", "anna's hummingbird", "platypus", "climbing perch")
output_path <- "/disk1/juwan_d1/cat_rerun/1_Main_figures/Fig5/"
gene_biotype_order <- c("protein_coding", "lncRNA","rRNA","snoRNA","snRNA","tRNA")

upstream_3kb_order <- c()
for (i in c(1:30)) {
  upstream_3kb_order <- c(upstream_3kb_order, paste0("five-prime-", 3000 - 100*(i-1), "bp"))}
downstream_3kb_order <- c()
for (i in c(1:30)) {
  downstream_3kb_order <- c(downstream_3kb_order, paste0("three-prime-", 100*(i), "bp"))}

fig5.a_b.order <- c(upstream_3kb_order, utr_n_cds, downstream_3kb_order)

# function for the shared theme of plots
shared_theme <- function(my_col, y_breaks, species_s, variant_type, coeff) {
  list(theme_pubr(),
       scale_color_manual(name="", values = my_col, labels = c("GC content", paste0("False ", variant_type))),
       scale_y_continuous(breaks = y_breaks, 
                          limits =c(0,100), name = "GC-content (%)", sec.axis = sec_axis(~./coeff, name=paste0("False ", variant_type, " frequency (%)"))),
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
       guides(linetype=guide_legend(ncol = 1, keyheight = unit(3, "mm")), #direction ="horizontal", byrow = TRUE, 
              shape=guide_legend(ncol = 1, keyheight = unit(3, "mm")), #direction ="horizontal", byrow = TRUE, 
              color = guide_legend(ncol = 1, keyheight = unit(3, "mm"))) #direction ="horizontal", byrow = TRUE, 
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

# load data and bind them to a single dataframe

all_summary_f <- average_gc_miss_df %>% 
  filter(species %in% species_order)
all_summary_f$species <- factor(all_summary_f$species, levels=species_order)
all_summary_f$general_type <- ifelse(all_summary_f$position %in% genic_order, "genic",
                                     ifelse(all_summary_f$position %in% upstream_3kb_order , "5'",
                                            ifelse(all_summary_f$position %in% downstream_3kb_order,"3'", "ELSE")))
all_summary_f <- all_summary_f %>%
  filter(general_type %in% c("genic","5'","3'"))

fig.5.a_b.maker <- function(species_s, summary_df, variant_type) {
  all_summary_f <- summary_df %>%
    filter(species==species_s) %>%
    filter(CpG_island=="protein_coding")
  
  if (variant_type=="SNP") {
    if (species_s == "climbing perch") {coeff=8} else if (species_s=="anna's hummingbird") {coeff=200} else {coeff=80}
  } else if (variant_type=="indel") {
    if (species_s == "climbing perch") {coeff=50} else if (species_s=="anna's hummingbird") {coeff=1600} else {coeff=500}
  }
  
  # (1) UTR and CDS
  utr_n_cds <-all_summary_f %>%
    filter(position %in% utr_n_cds_order)
  utr_n_cds$position_rename <- as.numeric(factor(utr_n_cds$position, level=utr_n_cds_order))
  
  # (2) Intron
  introns <-all_summary_f %>%
    filter(position %in% intron_order)
  introns$position_rename <- as.numeric(factor(introns$position, level=intron_order))
  introns$position_rename <- transform_intron_position(introns$position_rename)
  
  # (3) (UTR and CDS) + (Intron) --> gene body 
  introns$position_rename <- as.numeric(introns$position_rename)
  utr_n_cds$position_rename <- as.numeric(utr_n_cds$position_rename)
  utr_n_cds$type <- "exon"
  introns$type <- "intron"
  utr_cds_introns <- rbind(utr_n_cds, introns)
  utr_cds_introns$type <- factor(utr_cds_introns$type, level=c("exon","intron"))
  
  # (4) +-3kb (upstream/downstream) 
  upstream_3kb <- all_summary_f %>%
    filter(position %in% upstream_3kb_order)
  upstream_3kb$position_rename <- as.numeric(factor(upstream_3kb$position, level=upstream_3kb_order))
  upstream_3kb_labels <- c("-3", "-1.5", "TSS   ")
  downstream_3kb <- all_summary_f %>%
    filter(position %in% downstream_3kb_order)
  downstream_3kb$position_rename <- as.numeric(factor(downstream_3kb$position, level=downstream_3kb_order))
  downstream_3kb_labels <- c("   TTS","+1.5", "+3")

  y_breaks = c(0,20,40,60,80)
  color_palette <- c("GC"="#d42700", "SNP"="#1a73c7", "indel"="#1a73c7")
  # (7) generate plots
  gene_body_plot <- ggplot(utr_cds_introns) +
    geom_line(aes(x=position_rename, y=GC_mean, group=type, linetype = type, color="GC"),  lwd=0.3) +
    geom_point(aes(x=position_rename, y=GC_mean, group=type, color="GC"), size=0.4) +
    geom_line(aes(x=position_rename, y=get(variant_type)*coeff, group=type, linetype=type, color=variant_type), lwd=0.3) +
    geom_point(aes(x=position_rename, y=get(variant_type)*coeff, group=type, color=variant_type), size=0.4) +
    draw_gene_body_lines +
    #scale_shape_manual(values=c(19,1), labels=c("Exon or upstream/downstream sequence", "Intron"))+
    #scale_linetype(labels=c("Exon or upstream/downstream sequence", "Intron"))+
    scale_linetype_manual(name="", values = c("solid", "dotted"), labels = c("Exon or upstream/downstream sequence", "Intron"))+
    coord_cartesian(ylim=c(0,80.5), xlim=c(0.5,5.5), expand=FALSE)+
    labs(color="Color", title = "Gene body", linetype="Type", y="Percentage (%)") + 
    shared_theme(color_palette, y_breaks, species_s, variant_type, coeff)
    
  up3kb_plot <- ggplot(upstream_3kb) +
    geom_line(aes(x=position_rename, y=GC_mean, color="GC"), lwd=0.3) +
    geom_point(aes(x=position_rename, y=GC_mean, color="GC"), size=0.4) +
    geom_line(aes(x=position_rename, y=get(variant_type)*coeff, color=variant_type), lwd=0.3) +
    geom_point(aes(x=position_rename, y=get(variant_type)*coeff, color=variant_type), size=0.4) +
    shared_theme(color_palette, y_breaks, species_s, variant_type, coeff) + 
    scale_x_continuous(breaks=c(1,15,30), labels=upstream_3kb_labels, name=" ",
                       sec.axis=sec_axis(~. + 0, name=" ", breaks = c(1), labels=c(" "))) + 
    coord_cartesian(ylim=c(0,80.5), xlim=c(0.8,30.5), expand=FALSE)+
    labs(color="Color", title="Upstream", y="Percentage (%)")
  
  down3kb_plot <- ggplot(downstream_3kb) +
    geom_line(aes(x=position_rename, y=GC_mean, color="GC"), lwd=0.3) +
    geom_point(aes(x=position_rename, y=GC_mean, color="GC"), size=0.4) +
    geom_line(aes(x=position_rename, y=get(variant_type)*coeff, color=variant_type), lwd=0.3) +
    geom_point(aes(x=position_rename, y=get(variant_type)*coeff, color=variant_type), size=0.4) +
    shared_theme(color_palette, y_breaks, species_s, variant_type, coeff) + 
    scale_x_continuous(breaks=c(1,15,30), labels=downstream_3kb_labels, name=" ",
                       sec.axis=sec_axis(~. + 0, name=" ", breaks = c(1), labels=c(" "))) + 
    coord_cartesian(ylim=c(0,80.5), xlim=c(0.8,30.5), expand=FALSE)+
    labs(color="Color", title="Downstream") 
  my_legend <- gtable_filter(ggplotGrob(gene_body_plot), "guide-box")
  my_list <- list(up3kb_plot + theme(legend.position = "none", axis.title.y.right=element_blank(), axis.text.y.right=element_blank(), axis.ticks.length.y.right=unit(0,"cm")),
                  gene_body_plot + theme(legend.position = "none",axis.title.y=element_blank(), axis.text.y=element_blank()),
                  down3kb_plot + theme(legend.position = "none",axis.title.y.left=element_blank(), axis.text.y.left=element_blank(), axis.ticks.length.y.left=unit(0,"cm")))
  fig.5.legend = gtable_filter(ggplotGrob(gene_body_plot), "guide-box")
  fig.5.legend <- arrangeGrob(fig.5.legend)
  fig.5.legend_path <- paste0(output_path, "fig5.legend.png")
  ggsave(fig.5.legend_path, fig.5.legend, width=6, height=0.6,units = "in")
  print(6)
  return(my_list)
}
my_layout <- rbind(c(1,1,1,1,1,2,2,2,2,3,3,3,3,3))

# generate plots
head(all_summary_f)
zf_gene_SNP_plot <- fig.5.a_b.maker(species_s="zebra finch", summary_df=all_summary_f, variant_type="SNP")
ZF.SNP.plots <- c(zf_gene_SNP_plot)
fig5.a_b.zf_SNP_plot <- arrangeGrob(grobs=ZF.SNP.plots,layout_matrix = my_layout)
fig5.a_b.zf_SNP_filepath <- paste0(output_path, "fig5.a_d.ZF.SNP.png")
ggsave(fig5.a_b.zf_SNP_filepath, fig5.a_b.zf_SNP_plot, width=3, height=1.3) 

an_gene_SNP_plot <- fig.5.a_b.maker(species_s="anna's hummingbird", summary_df=all_summary_f, variant_type="SNP")
AN.SNP.plots <- c(an_gene_SNP_plot)
fig5.a_b.an_SNP_plot <- arrangeGrob(grobs=AN.SNP.plots,layout_matrix = my_layout)
fig5.a_b.an_SNP_filepath <- paste0(output_path, "fig5.a_d.AN.SNP.png")
ggsave(fig5.a_b.an_SNP_filepath, fig5.a_b.an_SNP_plot, width=3, height=1.3) 

zf_gene_indel_plot <- fig.5.a_b.maker(species_s="zebra finch", summary_df=all_summary_f, variant_type="indel")
ZF.indel.plots <- c(zf_gene_indel_plot)
fig5.a_b.zf_indel_plot <- arrangeGrob(grobs=ZF.indel.plots,layout_matrix = my_layout)
fig5.a_b.zf_indel_filepath <- paste0(output_path, "fig5.a_d.ZF.indel.png")
ggsave(fig5.a_b.zf_indel_filepath, fig5.a_b.zf_indel_plot, width=3, height=1.3) 

an_gene_indel_plot <- fig.5.a_b.maker(species_s="anna's hummingbird", summary_df=all_summary_f, variant_type="indel")
AN.indel.plots <- c(an_gene_indel_plot)
fig5.a_b.an_indel_plot <- arrangeGrob(grobs=AN.indel.plots,layout_matrix = my_layout)
fig5.a_b.an_indel_filepath <- paste0(output_path, "fig5.a_d.AN.indel.png")
ggsave(fig5.a_b.an_indel_filepath, fig5.a_b.an_indel_plot, width=3, height=1.3) 


#################################
######## Fig.5.a_b (done) #########
#################################
