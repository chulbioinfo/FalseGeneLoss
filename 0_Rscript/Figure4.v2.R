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
######## Fig.4.a #########
##########################

input_path = "/disk1/juwan_d1/cat_rerun/"
output_path <- "/disk1/juwan_d1/cat_rerun/1_Main_figures/Fig4/"
species_order <- c("zebra finch", "anna's hummingbird", "platypus", "climbing perch")
utr_n_cds_order <- c("five-prime-UTR", "FirstCDS",  "InternalCDS", "LastCDS", "three-prime-UTR")
intron_order <- c("5'UTR-intron","First-intron","Internal-intron", "Last-intron","3'UTR-intron")
genic_order <- c("five-prime-UTR", "5'UTR-intron", "FirstCDS","First-intron",  "InternalCDS", "Internal-intron","LastCDS","Last-intron", "three-prime-UTR", "3'UTR-intron")

exon_order <-c("FirstExon","InternalExon","LastExon")
gene_biotype_order <- c("protein_coding", "lncRNA","rRNA","snoRNA","snRNA","tRNA")

upstream_3kb_order <- c()
for (i in c(1:30)) {
  upstream_3kb_order <- c(upstream_3kb_order, paste0("five-prime-", 3000 - 100*(i-1), "bp"))}
downstream_3kb_order <- c()
for (i in c(1:30)) {
  downstream_3kb_order <- c(downstream_3kb_order, paste0("three-prime-", 100*(i), "bp"))}

# function for the shared theme of plots
shared_theme <- function(my_col, y_breaks) {
  list(theme_pubr(),
    scale_color_manual(values = my_col, labels = c("with CGIs, GC", "w/o CGIs, GC", 
                                                   "with CGIs, missing ratio",  "w/o CGIs, missing ratio")),
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
    guides(linetype=guide_legend(nrow = 2, keyheight = unit(3, "mm")),
           color = guide_legend(nrow = 2, keyheight = unit(3, "mm")))
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
cpg_summary_f = read.csv(paste0(input_path, "CpG.forR.gc_and_indel.tsv"), sep="\t", header=TRUE)
no_cpg_summary_f = read.csv(paste0(input_path,"no_CpG.forR.gc_and_indel.tsv"), sep="\t", header=TRUE)
all_summary_f <- rbind(cpg_summary_f, no_cpg_summary_f)
all_summary_f <- all_summary_f %>% 
  filter(species %in% species_order)
all_summary_f$species <- factor(all_summary_f$species, levels=species_order)
all_summary_f$general_type <- ifelse(all_summary_f$position %in% genic_order, "genic",
                                     ifelse(all_summary_f$position %in% upstream_3kb_order , "5'",
                                            ifelse(all_summary_f$position %in% downstream_3kb_order,"3'", "ELSE")))
all_summary_f <- all_summary_f %>%
  filter(general_type %in% c("genic","5'","3'"))
col_levels=c("True,GC_mean","False,GC_mean","True,missing","False,missing")

fig.4.a.maker <- function(species_s, summary_df) {
  all_summary_f <- summary_df %>%
    filter(species==species_s)
  
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
  
  # (5) melting the dataframe
  melt.utr_cds_introns <- melt(utr_cds_introns,id.vars=c("position","position_rename","type","species","CpG_island"), measure.vars=c("GC_mean","missing"))
  melt.upstream_3kb<- melt(upstream_3kb,id.vars=c("position","position_rename","species","CpG_island"), measure.vars=c("GC_mean","missing"))
  melt.downstream_3kb<- melt(downstream_3kb,id.vars=c("position","position_rename","species","CpG_island"), measure.vars=c("GC_mean","missing"))
  
  # (6) convert CpG_n_variable as factor (for coloring)
  melt.utr_cds_introns$CpG_n_variable <- paste0(melt.utr_cds_introns$CpG_island, ",", melt.utr_cds_introns$variable)
  melt.upstream_3kb$CpG_n_variable <- paste0(melt.upstream_3kb$CpG_island, ",", melt.upstream_3kb$variable)
  melt.downstream_3kb$CpG_n_variable <- paste0(melt.downstream_3kb$CpG_island, ",", melt.downstream_3kb$variable)
  melt.utr_cds_introns$CpG_n_variable <- factor(melt.utr_cds_introns$CpG_n_variable, levels=col_levels)
  melt.upstream_3kb$CpG_n_variable <- factor(melt.upstream_3kb$CpG_n_variable, levels=col_levels)
  melt.downstream_3kb$CpG_n_variable <- factor(melt.downstream_3kb$CpG_n_variable, levels=col_levels)
  my_col2 = c("#d42700","#d4a900","#1a73c7","#8098ad")
  y_breaks = c(0,20,40,60,80)
  
  # (7) generate plots
  gene_body_plot <- ggplot(melt.utr_cds_introns, aes(x =position_rename, y = value, color=CpG_n_variable)) +
    geom_line(aes(group=interaction(CpG_n_variable,type),linetype = type, color=CpG_n_variable), lwd=0.3) +
    geom_point(size=0.2, alpha=0.7) +
    draw_gene_body_lines +
    scale_linetype_manual(name="Type", values = c("solid", "dotted"), labels = c("Exon or upstream/downstream sequence", "Intron"))+
    coord_cartesian(ylim=c(0,80.5), xlim=c(0.5,5.5), expand=FALSE)+
    labs(color="Color", title = "Gene body", linetype="Type") + 
    shared_theme(my_col2, y_breaks)
  
  up3kb_plot <- ggplot(melt.upstream_3kb, aes(x =position_rename, y = value, color=CpG_n_variable)) + 
    geom_line(aes(group=CpG_n_variable, color=CpG_n_variable), lwd=0.3) +
    geom_point(aes(group=CpG_n_variable, color=CpG_n_variable), size=0.2, alpha=0.7) +
    scale_x_continuous(breaks=c(1,15,30), labels=upstream_3kb_labels, name=" ",
                       sec.axis=sec_axis(~. + 0, name=" ", breaks = c(1), labels=c(" "))) + 
    coord_cartesian(ylim=c(0,80.5), xlim=c(0.8,30.5), expand=FALSE)+
    labs(color="Color", title="Upstream", y="Percentage (%)") + 
    shared_theme(my_col2, y_breaks)
  
  down3kb_plot <- ggplot(melt.downstream_3kb, aes(x =position_rename, y = value, color=CpG_n_variable)) + 
    geom_line(aes(group=CpG_n_variable, color=CpG_n_variable), lwd=0.3) +
    geom_point(aes(group=CpG_n_variable, color=CpG_n_variable), size=0.2, alpha=0.7) +
    scale_x_continuous(breaks=c(1,15,30), labels=downstream_3kb_labels, name=" ",
                       sec.axis=sec_axis(~. + 0, name=" ", breaks = c(1), labels=c(" "))) + 
    coord_cartesian(ylim=c(0,80.5), xlim=c(0.8,30.5), expand=FALSE)+
    labs(color="Color", title="Downstream") + 
    shared_theme(my_col2, y_breaks)

  my_list <- list(up3kb_plot + theme(legend.position = "none"),
                  gene_body_plot + theme(legend.position = "none",axis.title.y=element_blank(), axis.text.y=element_blank()),
                  down3kb_plot + theme(legend.position = "none",axis.title.y=element_blank(), axis.text.y=element_blank()))
  fig.4.legend = gtable_filter(ggplotGrob(gene_body_plot), "guide-box")
  fig.4.legend <- arrangeGrob(fig.4.legend)
  fig.4.legend_path <- paste0(output_path, "fig4.a.CGI.legend.png")
  ggsave(fig.4.legend_path, fig.4.legend, width=6, height=0.6,units = "in")
  return(my_list)
}
my_layout <- rbind(c(1,1,1,1,1,2,2,2,2,3,3,3,3))

# generate plots
zf_gene_plot_list <- fig.4.a.maker(species_s="zebra finch", summary_df=all_summary_f)
ZF.plots <- c(zf_gene_plot_list)
fig4.a.zf_plot <- arrangeGrob(grobs=ZF.plots,layout_matrix = my_layout)
fig4.a.zf_filepath <- paste0(output_path, "fig4.a.ZF.CGIplots.png")
ggsave(fig4.a.zf_filepath, fig4.a.zf_plot, width=2.5, height=1.2) 

an_gene_plot_list <- fig.4.a.maker(species_s="anna's hummingbird", summary_df=all_summary_f)
AN.plots <- c(an_gene_plot_list)
fig4.a.an_plot <- arrangeGrob(grobs=AN.plots,layout_matrix = my_layout)
fig4.a.an_filepath <- paste0(output_path, "fig4.a.AN.CGIplots.png")
ggsave(fig4.a.an_filepath, fig4.a.an_plot, width=2.5, height=1.2) 

pl_gene_plot_list <- fig.4.a.maker(species_s="platypus", summary_df=all_summary_f)
PL.plots <- c(pl_gene_plot_list)
fig4.a.pl_plot <- arrangeGrob(grobs=PL.plots,layout_matrix = my_layout)
fig4.a.pl_filepath <- paste0(output_path, "fig4.a.PL.CGIplots.png")
ggsave(fig4.a.pl_filepath, fig4.a.pl_plot, width=2.5, height=1.2)

cl_gene_plot_list <- fig.4.a.maker(species_s="climbing perch", summary_df=all_summary_f)
CL.plots <- c(cl_gene_plot_list)
fig4.a.cl_plot <- arrangeGrob(grobs=CL.plots,layout_matrix = my_layout)
fig4.a.cl_filepath <- paste0(output_path, "fig4.a.CL.CGIplots.png")
ggsave(fig4.a.cl_filepath, fig4.a.cl_plot, width=2.5, height=1.2) 

#################################
######## Fig.4.a (done) #########
#################################

##########################
######## Fig.4.b #########
##########################

# load input dataframe
fig.4.colnames=c("gene","position","species", "gene_biotype","GC_mean","missing_percent","indel","hompol","SNP","softmasked")
fig.4.coltypes=cols("c","c","c","c","n","n","n","n","n","n")
all_gene_record_filepath <- paste0(input_path,"all_type_gene.all.record.tsv")
all_gene_record_df <- vroom(all_gene_record_filepath, delim="\t", col_names=fig.4.colnames,col_types = fig.4.coltypes)

fig.4.b.order <- c(upstream_3kb_order, exon_order, downstream_3kb_order)
fig.4.b.input_df <- all_gene_record_df %>% 
  filter(position %in% fig.4.b.order)
fig.4.b.input_df <- as.data.frame(fig.4.b.input_df)
# more than 10% of missing sequences ==> regarded as missing block
str(fig.4.b.input_df)
fig.4.b.input_df$missing_type <- ifelse(fig.4.b.input_df$missing_percent > 90, "Missing", "Present")
fig.4.b.input_df$num_position <- as.numeric(factor(fig.4.b.input_df$position, levels=fig.4.b.order))
# for each type of gene of each species, calculate GC and missing ratio of each position 
gc_n_missing_summary <- data.frame()
for (species_s in species_order) {
  for (biotype_s in gene_biotype_order) {
    each_species_n_genetype <- fig.4.b.input_df %>%
      filter(species==species_s) %>%
      filter(gene_biotype==biotype_s)
    each_species_n_genetype <- as.data.frame(each_species_n_genetype)
    for (position_s in fig.4.b.order) {
      each_pos_blocks <- subset(each_species_n_genetype, position==position_s)
      # missinb blocks
      missing_pos_blocks <- subset(each_pos_blocks, missing_type=="Missing")
      present_pos_blocks <- subset(each_pos_blocks, missing_type=="Present")
      #GC content (mean, sd)
      missing_gc_mean <- mean(missing_pos_blocks$GC_mean)
      missing_gc_sd <- sd(missing_pos_blocks$GC_mean)
      present_gc_mean <- mean(present_pos_blocks$GC_mean)
      present_gc_sd <- sd(present_pos_blocks$GC_mean)
      # missing ratio
      missing_ratio <- sum(missing_pos_blocks$missing_percent)*100/(length(each_pos_blocks$missing_percent)*100)
      gc_n_missing_summary <- rbind(gc_n_missing_summary, c(biotype_s, species_s, position_s, missing_gc_mean, present_gc_mean,
                                                            missing_gc_sd, present_gc_sd, missing_ratio))
    }}}

# set column names + convert as numeric
colnames(gc_n_missing_summary) <- c("biotype", "species", "position", "missing_gc_mean", "present_gc_mean", "missing_gc_sd", "present_gc_sd", "missing_ratio")
cols.num <- c("missing_gc_mean", "present_gc_mean","missing_gc_sd", "present_gc_sd", "missing_ratio")
gc_n_missing_summary[cols.num] <- sapply(gc_n_missing_summary[cols.num],as.numeric)

gc_n_missing_summary$general_type <- ifelse(gc_n_missing_summary$position %in% upstream_3kb_order, "5'", "3'") 
gc_n_missing_summary$general_type <- factor(gc_n_missing_summary$general_type , levels=c("5'", "3'"))
gc_n_missing_summary$species <- factor(gc_n_missing_summary$species, levels=species_order)
gc_n_missing_summary$biotype <- ifelse(gc_n_missing_summary$biotype == "protein_coding", "protein\ncoding", gc_n_missing_summary$biotype)
gc_n_missing_summary$biotype <- factor(gc_n_missing_summary$biotype, levels=c("protein\ncoding", "lncRNA","snoRNA","snRNA","rRNA","tRNA"))

# numeric positions; classify the genes according to whether they have multiple exons (pcg & lncRNA genes) or not (snoRNA, snRNA, rRNA, tRNA): most of tRNA have less than 3 exons
gc_n_missing_summary$num_position <- as.numeric(factor(gc_n_missing_summary$position, levels=fig.4.b.order))
gc_n_missing_summary$num_position <- ifelse(gc_n_missing_summary$biotype %in% c("protein\ncoding", "lncRNA"),
                             # genes with multiple exons (PCG and lncRNA)
                             ifelse(gc_n_missing_summary$position %in% upstream_3kb_order, 
                                    gc_n_missing_summary$num_position,
                                    ifelse(gc_n_missing_summary$position %in% downstream_3kb_order,
                                           gc_n_missing_summary$num_position + 28, (gc_n_missing_summary$num_position-31)*10 + 36)),
                             # other genes (snoRNA, snRNA, rRNA, tRNA)
                             ifelse(gc_n_missing_summary$position %in% upstream_3kb_order, 
                                    gc_n_missing_summary$num_position,
                                    ifelse(gc_n_missing_summary$position %in% downstream_3kb_order, 
                                           gc_n_missing_summary$num_position + 8, 36)))
# color palette
color_palette <- c("missing_gc_mean"="#d42700", "present_gc_mean"="#9c9c9c", "z_missing_ratio"="#1a73c7")

head(gc_n_missing_summary)

fig.4.b.maker<- function(gc_n_missing_summary,species_s,gene_biotype){
  if(gene_biotype %in% c("rRNA","snRNA","snoRNA","tRNA")){
    gc_n_missing_summary <- subset(gc_n_missing_summary,!gc_n_missing_summary$position %in% c("InternalExon","LastExon"))} else {
      gc_n_missing_summary <- gc_n_missing_summary}
  gc_n_missing_summary$num_position <- as.numeric(factor(gc_n_missing_summary$position, levels=fig.4.b.order))
  gc_n_missing_summary$num_position <- ifelse(gc_n_missing_summary$biotype %in% c("protein\ncoding", "lncRNA"),
                                # PCG and lncRNA
                                ifelse(gc_n_missing_summary$position %in% upstream_3kb_order, 
                                       gc_n_missing_summary$num_position,
                                       ifelse(gc_n_missing_summary$position %in% exon_order,(gc_n_missing_summary$num_position-31)*10 + 36, gc_n_missing_summary$num_position + 28)),
                                # other genes
                                ifelse(gc_n_missing_summary$position %in% upstream_3kb_order, 
                                       gc_n_missing_summary$num_position,
                                       ifelse(gc_n_missing_summary$position %in% downstream_3kb_order, gc_n_missing_summary$num_position + 8, 36)))
  gc_n_missing_summary <- gc_n_missing_summary %>%
    filter(species==species_s) %>%
    filter(gene_biotype==biotype)

  Fig4.d_g.plot <- ggplot(gc_n_missing_summary) + 
    # Average GC content of missing & present blocks
    geom_line(data=subset(gc_n_missing_summary, gc_n_missing_summary$position %in% upstream_3kb_order), 
              aes(x=num_position, y=missing_gc_mean, color="missing_gc_mean"),lwd=0.1) + 
    geom_line(data=subset(gc_n_missing_summary, gc_n_missing_summary$position %in% upstream_3kb_order), 
              aes(x=num_position, y=present_gc_mean, color="present_gc_mean"),lwd=0.1) + 
    geom_line(data=subset(gc_n_missing_summary, gc_n_missing_summary$position %in% downstream_3kb_order), 
              aes(x=num_position, y=missing_gc_mean, color="missing_gc_mean"),lwd=0.1) + 
    geom_line(data=subset(gc_n_missing_summary, gc_n_missing_summary$position %in% downstream_3kb_order), 
              aes(x=num_position, y=present_gc_mean, color="present_gc_mean"),lwd=0.1) + 
    # S.D. of GC content of missing & present blocks
    geom_ribbon(data=subset(gc_n_missing_summary, gc_n_missing_summary$position %in% upstream_3kb_order), 
                aes(x=num_position, ymin=missing_gc_mean-missing_gc_sd, ymax= missing_gc_mean+missing_gc_sd),
                fill="#d42700",alpha=0.15, lwd=0) + 
    geom_ribbon(data=subset(gc_n_missing_summary, gc_n_missing_summary$position %in% upstream_3kb_order), 
                aes(x=num_position, ymin=present_gc_mean-present_gc_sd, ymax= present_gc_mean+present_gc_sd),
                fill="#9c9c9c",alpha=0.15, lwd=0) +
    geom_ribbon(data=subset(gc_n_missing_summary, gc_n_missing_summary$position %in% downstream_3kb_order), 
                aes(x=num_position, ymin=missing_gc_mean-missing_gc_sd, ymax= missing_gc_mean+missing_gc_sd),
                fill="#d42700",alpha=0.15, lwd=0) + 
    geom_ribbon(data=subset(gc_n_missing_summary, gc_n_missing_summary$position %in% downstream_3kb_order), 
                aes(x=num_position, ymin=present_gc_mean-present_gc_sd, ymax= present_gc_mean+present_gc_sd),
                fill="#9c9c9c",alpha=0.15, lwd=0) +
    # Missing ratio of upstream/downstream 3kb regions
    geom_line(aes(x=num_position, y=missing_ratio, color="z_missing_ratio"), 
              stat="identity", lwd=0.2, linetype="solid") + 
    geom_point(aes(x=num_position, y=missing_ratio, color="z_missing_ratio"),
               size=0.1, lwd=0.1) +
    # Missing ratio of exons
    geom_bar(data=subset(gc_n_missing_summary, gc_n_missing_summary$position %in% c(exon_order)), 
             aes(x=num_position-1.5, y=missing_gc_mean), fill="#d42700",stat="identity",lwd=0,
             alpha=0.5,width=3)+
    geom_bar(data=subset(gc_n_missing_summary, gc_n_missing_summary$position %in% c(exon_order)), 
             aes(x=num_position+1.5, y=present_gc_mean), fill="#9c9c9c",stat="identity",lwd=0,
             alpha=0.5,width=3)+
    scale_colour_manual(name=" ",values=color_palette, labels=c("GC-content (Missing)", "GC-content (Present)", "Missing")) + 
    scale_y_continuous(limits =c(0,100), breaks=c(0,50,100))+
    labs(y = "Percentage (%)", fill="Log10(p-value)  ")+
    theme_pubr() + 
    theme(text=element_text(size=6),
          strip.text = element_text(size=7),
          strip.text.x = element_blank(),
          legend.position="bottom",
          strip.background.x = element_blank(),
          legend.key.width = unit(0.6,"cm"),
          legend.key.size = unit(0.2, "cm"),
          legend.key.height = unit(0.2,"cm"),
          axis.line.x.top=element_blank(),
          axis.ticks.length.x.top = unit(0,"cm"),
          axis.line = element_line(size=0.1),
          axis.ticks = element_line(size=0.1),
          legend.box="vertical",
          legend.box.just = "left",
          axis.title.x=element_blank()) + 
    guides(linetype=guide_legend(ncol = 3, keyheight = unit(3, "mm")),
           color = guide_legend(ncol = 3, keyheight = unit(3, "mm")))
  # x axis breaks depending on the gene biotype
  if(gene_biotype %in% c("protein\ncoding", "lncRNA")){
    Fig4.d_g.plot <- Fig4.d_g.plot + scale_x_continuous(breaks=c(1,30, # upstream 3kb 
                                                                 36,46,56, # first, internal, last exons
                                                                 max(gc_n_missing_summary$num_position)-30,max(gc_n_missing_summary$num_position)), #downstream 3kb
                                                        labels=c("-3","TSS      ","F","I","L","      TTS","+3"), expand=c(0,0))}
  else {
    Fig4.d_g.plot <- Fig4.d_g.plot + scale_x_continuous(breaks=c(1,30, # upstream 3kb 
                                                                 36, # internal exon
                                                                 max(gc_n_missing_summary$num_position)-30,max(gc_n_missing_summary$num_position)), #downstream 3kb
                                                        labels=c("-3","TSS      ","E","      TTS","+3"), expand=c(0,0))}
  fig4.d_g.legend = gtable_filter(ggplotGrob(Fig4.d_g.plot), "guide-box")
  fig4.d_g.legend_path = "/disk1/juwan_d1/cat_rerun/1_Main_figures/Fig4/fig4.b.legend.png"
  ggsave(fig4.d_g.legend_path, fig4.d_g.legend, width=6, height=2,units = "in")
  result_plot_list <- list(Fig4.d_g.plot + theme(legend.position = "none", axis.title.y=element_blank()))
  return(result_plot_list)
}
# finch
finch_pcg <- fig.4.b.maker(gc_n_missing_summary=gc_n_missing_summary,species_s="zebra finch",gene_biotype="protein\ncoding")
finch_lncRNA <- fig.4.b.maker(gc_n_missing_summary=gc_n_missing_summary,species_s="zebra finch",gene_biotype="lncRNA")
finch_rRNA <- fig.4.b.maker(gc_n_missing_summary=gc_n_missing_summary,species_s="zebra finch",gene_biotype="rRNA")
finch_snoRNA <- fig.4.b.maker(gc_n_missing_summary=gc_n_missing_summary,species_s="zebra finch",gene_biotype="snoRNA")
finch_snRNA <- fig.4.b.maker(gc_n_missing_summary=gc_n_missing_summary,species_s="zebra finch",gene_biotype="snRNA")
finch_tRNA <- fig.4.b.maker(gc_n_missing_summary=gc_n_missing_summary,species_s="zebra finch",gene_biotype="tRNA")
# anna
anna_pcg <- fig.4.b.maker(gc_n_missing_summary=gc_n_missing_summary,species_s="anna's hummingbird",gene_biotype="protein\ncoding")
anna_lncRNA <- fig.4.b.maker(gc_n_missing_summary=gc_n_missing_summary,species_s="anna's hummingbird",gene_biotype="lncRNA")
anna_rRNA <- fig.4.b.maker(gc_n_missing_summary=gc_n_missing_summary,species_s="anna's hummingbird",gene_biotype="rRNA")
anna_snoRNA <- fig.4.b.maker(gc_n_missing_summary=gc_n_missing_summary,species_s="anna's hummingbird",gene_biotype="snoRNA")
anna_snRNA <- fig.4.b.maker(gc_n_missing_summary=gc_n_missing_summary,species_s="anna's hummingbird",gene_biotype="snRNA")
anna_tRNA <- fig.4.b.maker(gc_n_missing_summary=gc_n_missing_summary,species_s="anna's hummingbird",gene_biotype="tRNA")
# platypus
platypus_pcg <- fig.4.b.maker(gc_n_missing_summary=gc_n_missing_summary,species_s="platypus",gene_biotype="protein\ncoding")
platypus_lncRNA <- fig.4.b.maker(gc_n_missing_summary=gc_n_missing_summary,species_s="platypus",gene_biotype="lncRNA")
platypus_rRNA <- fig.4.b.maker(gc_n_missing_summary=gc_n_missing_summary,species_s="platypus",gene_biotype="rRNA")
platypus_snoRNA <- fig.4.b.maker(gc_n_missing_summary=gc_n_missing_summary,species_s="platypus",gene_biotype="snoRNA")
platypus_snRNA <- fig.4.b.maker(gc_n_missing_summary=gc_n_missing_summary,species_s="platypus",gene_biotype="snRNA")
platypus_tRNA <- fig.4.b.maker(gc_n_missing_summary=gc_n_missing_summary,species_s="platypus",gene_biotype="tRNA")
# perch
perch_pcg <- fig.4.b.maker(gc_n_missing_summary=gc_n_missing_summary,species_s="climbing perch",gene_biotype="protein\ncoding")
perch_lncRNA <- fig.4.b.maker(gc_n_missing_summary=gc_n_missing_summary,species_s="climbing perch",gene_biotype="lncRNA")
perch_rRNA <- fig.4.b.maker(gc_n_missing_summary=gc_n_missing_summary,species_s="climbing perch",gene_biotype="rRNA")
perch_snoRNA <- fig.4.b.maker(gc_n_missing_summary=gc_n_missing_summary,species_s="climbing perch",gene_biotype="snoRNA")
perch_snRNA <- fig.4.b.maker(gc_n_missing_summary=gc_n_missing_summary,species_s="climbing perch",gene_biotype="snRNA")
perch_tRNA <- fig.4.b.maker(gc_n_missing_summary=gc_n_missing_summary,species_s="climbing perch",gene_biotype="tRNA")
# all_genes
all_plots <- c(finch_pcg, finch_lncRNA, finch_rRNA, finch_snoRNA, finch_snRNA, finch_tRNA,
               anna_pcg, anna_lncRNA, anna_rRNA, anna_snoRNA, anna_snRNA, anna_tRNA,
               platypus_pcg, platypus_lncRNA, platypus_rRNA, platypus_snoRNA, platypus_snRNA, platypus_tRNA,
               perch_pcg, perch_lncRNA, perch_rRNA, perch_snoRNA, perch_snRNA, perch_tRNA)
my_layout<-rbind(c(1,2,3,4,5,6),
                 c(7,8,9,10,11,12),
                 c(13,14,15,16,17,18),
                 c(19,20,21,22,23,24))
my_layout<-rbind(c(1,7,13,19),
                 c(2,8,14,20),
                 c(3,9,15,21),
                 c(4,10,16,22),
                 c(5,11,17,23),
                 c(6,12,18,24))

fig.4.b.plots <- arrangeGrob(grobs = all_plots,
                          layout_matrix = my_layout)

fig.4.b.filepath <-"/disk1/juwan_d1/cat_rerun/1_Main_figures/Fig4/fig4.b.png"
ggsave(fig.4.b.filepath, fig.4.b.plots, width=4.5, height=4.0)

#################################
######## Fig.4.b (done) #########
#################################
