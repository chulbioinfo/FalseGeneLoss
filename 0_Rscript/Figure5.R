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
utr_n_cds <- c("five-prime-UTR", "FirstCDS",  "InternalCDS", "LastCDS", "three-prime-UTR")
species_order <- c("zebra finch", "anna's hummingbird", "platypus", "climbing perch")
gene_biotype_order <- c("protein_coding", "lncRNA","rRNA","snoRNA","snRNA","tRNA")

upstream_3kb_order <- c()
for (i in c(1:30)) {
  upstream_3kb_order <- c(upstream_3kb_order, paste0("five-prime-", 3000 - 100*(i-1), "bp"))}
downstream_3kb_order <- c()
for (i in c(1:30)) {
  downstream_3kb_order <- c(downstream_3kb_order, paste0("three-prime-", 100*(i), "bp"))}
upstream_3kb_inside_order <- c()
for (i in c(1:30)) {
  upstream_3kb_inside_order <- c(upstream_3kb_inside_order, paste0("five-inside-", 100*(i), "bp"))}
downstream_3kb_inside_order <- c()
for (i in c(1:30)) {
  downstream_3kb_inside_order <- c(downstream_3kb_inside_order, paste0("three-inside-", 3100-100*(i), "bp"))}

fig5.c_f.order <- c(upstream_3kb_order, upstream_3kb_inside_order, downstream_3kb_inside_order, downstream_3kb_order)

# draw plots
for (species_s in c("zebra finch", "anna's hummingbird")) {
  # (0) select each species
  exon_df <- average_gc_miss_df %>%
    filter(species == species_s) %>%
    filter(position %in% utr_n_cds) %>%
    filter(CpG_island == "protein_coding") 
  exon_df$num_position <- as.numeric(factor(exon_df$position, levels=utr_n_cds))
  fig.5.a_b.path <- paste0("/disk1/juwan_d1/cat_rerun/1_Main_figures/Fig5/fig5.a_b.",species_s,".png")
  png(filename=fig.5.a_b.path, width=3, height=1.8, units="in", res=1000, pointsize = 6)
  par(mar=c(3,3,3,7) + 0.1,cex.axis = 1,cex.lab  = 1)
  # (1) reorder exons --> 5'UTR, first, internal, last, 3'UTR exons
  ordered_exon_df <- exon_df[order(exon_df$num_position), ]
  # (2) GC content 
  plot(ordered_exon_df$num_position, ordered_exon_df$GC_mean, axes=F, 
       ylim=c(min(ordered_exon_df$GC_mean),max(ordered_exon_df$GC_mean)), 
       xlab="", ylab="", type="l",col="black", main="",xlim=c(0.5,5.5))
  axis(2, ylim=c(0,max(ordered_exon_df$GC_mean)),col="black",lwd=2)
  points(ordered_exon_df$num_position,ordered_exon_df$GC_mean,pch=20,col="black")
  mtext(2,text="GC content",line=2)
  # (3) false SNP
  par(new=T)
  plot(ordered_exon_df$num_position, ordered_exon_df$SNP, axes=F, ylim=c(0,max(ordered_exon_df$SNP)), 
       xlab="", ylab="", type="l",lty=3, main="",xlim=c(0.5,5.5),lwd=2,col="blue")
  axis(4, ylim=c(0,max(ordered_exon_df$SNP)),lwd=2,col="blue")
  points(ordered_exon_df$num_position, ordered_exon_df$SNP,pch=20,col="blue")
  mtext(4,text="Freq. of false SNP",col="blue",line=2)
  # (4) false indel
  par(new=T)
  plot(ordered_exon_df$num_position, ordered_exon_df$indel, axes=F, ylim=c(0,max(ordered_exon_df$indel)), xlab="", ylab="", 
       type="l",lty=3, main="",xlim=c(0.5,5.5),lwd=2,col="green")
  axis(4, ylim=c(0,max(ordered_exon_df$indel)),lwd=2,line=3.5,col="green")
  points(ordered_exon_df$num_position, ordered_exon_df$indel,pch=20,col="green")
  mtext(4,text="Freq. of false indel",line=5.5,col="green")
  # (5) draw axis
  axis(1,xlim=c(0.5,5.5),at=c(1,2,3,4,5),
       labels=c("5'","F","I","L","3'")) # pretty(range(ordered_exon_df$num_position),10),
  mtext("Position",side=1,col="black",line=2)
  dev.off()
}

#################################
####### Fig.5.a-b (done) ########
#################################

############################
######## Fig.5.c-f #########
############################

# load input dataframe
all_gene_record_filepath <- paste0(input_path,"all_type_gene.all.record.tsv")
fig.5.c_f.colnames=c("gene","position","species", "gene_biotype","GC_mean","missing_percent","indel","hompol","SNP","softmasked")
fig.5.c_f.coltypes=cols("c","c","c","c","n","n","n","n","n","n")
all_gene_record_df <- vroom(all_gene_record_filepath, delim="\t", col_names=fig.5.c_f.colnames, col_types=fig.5.c_f.coltypes)

fig.5.c_f.order <- c(upstream_3kb_order, upstream_3kb_inside_order, downstream_3kb_inside_order, downstream_3kb_order)
fig.5.c_f.input_df <- all_gene_record_df %>% 
  filter(position %in% fig.5.c_f.order) %>% 
  filter(gene_biotype %in% "protein_coding")
fig.5.c_f.input_df <- as.data.frame(fig.5.c_f.input_df)
fig.5.c_f.input_df$num_position <- as.numeric(factor(fig.5.c_f.input_df$position, levels=fig.5.c_f.order))

#########################
######## (1) SNP ########
#########################
# with 1bp or more false SNP ==> regarded as false SNP block
fig.5.c_f.input_df$snp_type <- ifelse(fig.5.c_f.input_df$SNP >= 1, "SNP", "Normal")
# for each type of gene of each species, calculate GC and missing ratio of each position 
gc_n_SNP_summary <- data.frame()
for (species_s in species_order) {
  each_species_n_genetype <- fig.5.c_f.input_df %>%
    filter(species==species_s)
  each_species_n_genetype <- as.data.frame(each_species_n_genetype)
  for (position_s in fig.5.c_f.order) {
    each_pos_blocks <- subset(each_species_n_genetype, position==position_s)
    # missinb blocks
    snp_pos_blocks <- subset(each_pos_blocks, snp_type=="SNP")
    present_pos_blocks <- subset(each_pos_blocks, snp_type=="Normal")
    #GC content (mean, sd)
    snp_gc_mean <- mean(snp_pos_blocks$GC_mean)
    snp_gc_sd <- sd(snp_pos_blocks$GC_mean)
    normal_gc_mean <- mean(present_pos_blocks$GC_mean)
    normal_gc_sd <- sd(present_pos_blocks$GC_mean)
    # missing ratio
    snp_ratio <- sum(snp_pos_blocks$SNP)*100/(length(each_pos_blocks$SNP)*100)
    gc_n_SNP_summary <- rbind(gc_n_SNP_summary, c(species_s, position_s, snp_gc_mean, normal_gc_mean,snp_gc_sd, normal_gc_sd, snp_ratio))
  }}

# set column names + convert as numeric
colnames(gc_n_SNP_summary) <- c("species", "position", "snp_gc_mean", "normal_gc_mean", "snp_gc_sd", "normal_gc_sd", "snp_ratio")
cols.num <- c("snp_gc_mean", "normal_gc_mean","snp_gc_sd", "normal_gc_sd", "snp_ratio")
head(gc_n_SNP_summary)
gc_n_SNP_summary[cols.num] <- sapply(gc_n_SNP_summary[cols.num],as.numeric)
gc_n_SNP_summary$general_type <- ifelse(gc_n_SNP_summary$position %in% c(upstream_3kb_order, upstream_3kb_inside_order), "5'", "3'") 
gc_n_SNP_summary$general_type <- factor(gc_n_SNP_summary$general_type , levels=c("5'", "3'"))
gc_n_SNP_summary$species <- factor(gc_n_SNP_summary$species, levels=species_order)


###############################
######## (1) SNP (end) ########
###############################

#########################
####### (2) Indel #######
#########################

fig.5.c_f.input_df$indel_type <- ifelse(fig.5.c_f.input_df$indel >= 1, "indel", "Normal")
# for each type of gene of each species, calculate GC and missing ratio of each position 
gc_n_indel_summary <- data.frame()
for (species_s in species_order) {
    each_species_n_genetype <- fig.5.c_f.input_df %>%
      filter(species==species_s) 
    each_species_n_genetype <- as.data.frame(each_species_n_genetype)
    for (position_s in fig.5.c_f.order) {
      each_pos_blocks <- subset(each_species_n_genetype, position==position_s)
      # missinb blocks
      indel_pos_blocks <- subset(each_pos_blocks, indel_type=="indel")
      present_pos_blocks <- subset(each_pos_blocks, indel_type=="Normal")
      #GC content (mean, sd)
      indel_gc_mean <- mean(indel_pos_blocks$GC_mean)
      indel_gc_sd <- sd(indel_pos_blocks$GC_mean)
      normal_gc_mean <- mean(present_pos_blocks$GC_mean)
      normal_gc_sd <- sd(present_pos_blocks$GC_mean)
      # missing ratio
      indel_ratio <- sum(indel_pos_blocks$indel)*100/(length(each_pos_blocks$indel)*100)
      gc_n_indel_summary <- rbind(gc_n_indel_summary, c(species_s, position_s, indel_gc_mean, normal_gc_mean,
                                                    indel_gc_sd, normal_gc_sd, indel_ratio))
    }}

# set column names + convert as numeric
colnames(gc_n_indel_summary) <- c("species", "position", "indel_gc_mean", "normal_gc_mean", "indel_gc_sd", "normal_gc_sd", "indel_ratio")
cols.num <- c("indel_gc_mean", "normal_gc_mean","indel_gc_sd", "normal_gc_sd", "indel_ratio")
gc_n_indel_summary[cols.num] <- sapply(gc_n_indel_summary[cols.num],as.numeric)
gc_n_indel_summary$general_type <- ifelse(gc_n_indel_summary$position %in% c(upstream_3kb_order,upstream_3kb_inside_order), "5'", "3'") 
gc_n_indel_summary$general_type <- factor(gc_n_indel_summary$general_type , levels=c("5'", "3'"))
gc_n_indel_summary$species <- factor(gc_n_indel_summary$species, levels=species_order)

fig.5.c_d.maker<- function(summary_df,species_s){
  summary_df <- summary_df %>%
    filter(species==species_s) 
  summary_df$num_position <- as.numeric(factor(summary_df$position, levels=fig5.c_f.order))
  # color palette
  color_palette <- c("0_snp_gc_mean"="#d42700", "1_normal_gc_mean"="#9c9c9c", "2_snp_ratio"="#1a73c7")
  if (species_s == "climbing perch") {coeff=8} else if (species_s=="anna's hummingbird") {coeff=200} else {coeff=80}
  fig5.c_d.plot <- ggplot(summary_df) + 
    facet_wrap(.~general_type, ncol=2, scales="free_x")+
    # GC content of snp / normal blocks (mean + S.D.))
    geom_line(aes(x=num_position, y=snp_gc_mean, color="0_snp_gc_mean"),lwd=0.3) + 
    geom_line(aes(x=num_position, y=normal_gc_mean, color="1_normal_gc_mean"),lwd=0.3) +
    # snp frequency
    geom_line(aes(x=num_position, y=snp_ratio*coeff, color="2_snp_ratio"), lwd=0.3) + 
    geom_ribbon(aes(x=num_position, ymin=snp_gc_mean-snp_gc_sd, ymax= snp_gc_mean+snp_gc_sd),
                fill="#d42700",linetype="dotted",lwd=0, alpha=0.08) + 
    geom_ribbon(aes(x=num_position, ymin=normal_gc_mean-normal_gc_sd, ymax= normal_gc_mean+normal_gc_sd),
                fill="#9c9c9c",linetype="dotted",lwd=0, alpha=0.08) + 
    geom_vline(data=filter(summary_df, general_type=="5'"), aes(xintercept=30.5), color="grey",linetype="dotted",lwd=0.3, show.legend = FALSE) + 
    geom_vline(data=filter(summary_df, general_type=="3'"),aes(xintercept=90.5), color="grey",linetype="dotted",lwd=0.3, show.legend = FALSE) + 
    scale_colour_manual(name="",values=color_palette, labels=c("%GC(With false SNP)", "%GC(Without false SNP)", "%False SNP frequency")) + 
    scale_y_continuous(limits =c(0,100),name = "GC content (%)", sec.axis = sec_axis(~./coeff, name="False SNP frequency (%)"))+
    scale_x_continuous(breaks=c(1,15,30.5,45,60,61,75,90.5,105,120),
                       labels=c("-3","-1.5","TSS","+1.5","+3","-3","-1.5","TTS","+1.5","+3"),
                       expand=c(0,0)) + 
    labs(y = "Percentage (%)")+
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
          axis.title.x=element_blank())
  fig5.c_d.legend = gtable_filter(ggplotGrob(fig5.c_d.plot), "guide-box")
  fig5.c_d.legend_path = "/disk1/juwan_d1/cat_rerun/1_Main_figures/Fig5/fig5.c_d.legend.png"
  ggsave(fig5.c_d.legend_path, fig5.c_d.legend, width=6, height=2,units = "in")
  fig5.c_d.png <-paste0("/disk1/juwan_d1/cat_rerun/1_Main_figures/Fig5/fig5.c_d.",species_s,".SNP.png")
  ggsave(fig5.c_d.png, fig5.c_d.plot + theme(legend.position="none"), width=3, height=2)
}
fig.5.e_f.maker<- function(summary_df,species_s){
  summary_df <- summary_df %>%
    filter(species==species_s)
  summary_df$num_position <- as.numeric(factor(summary_df$position, levels=fig5.c_f.order))
  print(summary_df$num_position)
  # color palette
  color_palette <- c("0_indel_gc_mean"="#d42700", "1_normal_gc_mean"="#9c9c9c", "2_indel_ratio"="#1a73c7")
  if (species_s == "climbing perch") {coeff=50} else if (species_s=="anna's hummingbird") {coeff=1600} else {coeff=500}
  fig5.e_f.plot <- ggplot(summary_df) + 
    facet_wrap(.~general_type, ncol=2, scales="free_x")+
    # GC content of indel / normal blocks (mean + S.D.))
    geom_line(aes(x=num_position, y=indel_gc_mean, color="0_indel_gc_mean"),lwd=0.3) + 
    geom_line(aes(x=num_position, y=normal_gc_mean, color="1_normal_gc_mean"),lwd=0.3) +
    # indel frequency
    geom_line(aes(x=num_position, y=indel_ratio*coeff, color="2_indel_ratio"), lwd=0.3) + 
    geom_ribbon(aes(x=num_position, ymin=indel_gc_mean-indel_gc_sd, ymax= indel_gc_mean+indel_gc_sd),
                fill="#d42700",linetype="dotted",lwd=0, alpha=0.08) + 
    geom_ribbon(aes(x=num_position, ymin=normal_gc_mean-normal_gc_sd, ymax= normal_gc_mean+normal_gc_sd),
                fill="#9c9c9c",linetype="dotted",lwd=0, alpha=0.08) + 
    geom_vline(data=filter(summary_df, general_type=="5'"), aes(xintercept=30.5), color="grey",linetype="dotted",lwd=0.3, show.legend = FALSE) + 
    geom_vline(data=filter(summary_df, general_type=="3'"),aes(xintercept=90.5), color="grey",linetype="dotted",lwd=0.3, show.legend = FALSE) + 
    scale_colour_manual(name="",values=color_palette, labels=c("%GC(With false indel)", "%GC(Without false indel)", "%False indel frequency")) + 
    scale_y_continuous(limits =c(0,100),name = "GC content (%)", sec.axis = sec_axis(~./coeff, name="False indel frequency (%)"))+
    scale_x_continuous(breaks=c(1,15,30.5,45,60,61,75,90.5,105,120),
                       labels=c("-3","-1.5","TSS","+1.5","+3","-3","-1.5","TTS","+1.5","+3"),
                       expand=c(0,0)) + 
    labs(y = "Percentage (%)")+
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
          axis.title.x=element_blank())
  # x axis breaks depending on the gene biotype
  fig5.e_f.legend = gtable_filter(ggplotGrob(fig5.e_f.plot), "guide-box")
  fig5.e_f.legend_path = "/disk1/juwan_d1/cat_rerun/1_Main_figures/Fig5/fig5.e_f.legend.png"
  ggsave(fig5.e_f.legend_path, fig5.e_f.legend, width=6, height=2,units = "in")
  fig5.e_f.png <-paste0("/disk1/juwan_d1/cat_rerun/1_Main_figures/Fig5/fig5.e_f.",species_s,".indel.png")
  ggsave(fig5.e_f.png, fig5.e_f.plot + theme(legend.position="none"), width=3, height=2)
}

finch_snp_pcg <- fig.5.c_d.maker(summary_df=gc_n_SNP_summary,species_s="zebra finch")
anna_snp_pcg <- fig.5.c_d.maker(summary_df=gc_n_SNP_summary,species_s="anna's hummingbird")
finch_indel_pcg <- fig.5.e_f.maker(summary_df=gc_n_indel_summary,species_s="zebra finch")
finch_indel_pcg <- fig.5.e_f.maker(summary_df=gc_n_indel_summary,species_s="anna's hummingbird")

################################
####### (2) Indel (done) #######
################################