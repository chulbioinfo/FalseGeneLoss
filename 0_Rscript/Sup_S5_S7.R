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
####### Sup.fig S5 #######
##########################

# load input dataframe
input_path <- "/disk1/juwan_d1/cat_rerun/"
all_gene_record_filepath <- paste0(input_path,"all_type_gene.all.record.tsv")
species_order <- c("Zebra finch", "Anna's hummingbird", "Platypus", "Climbing perch")
sup.fig.s5.colnames=c("gene","position","species", "gene_biotype","GC_mean","missing_percent","indel","hompol","SNP","softmasked")
sup.fig.s5.coltypes=cols("c","c","c","c","n","n","n","n","n","n")
all_gene_record_df <- vroom(all_gene_record_filepath, delim="\t", col_names=sup.fig.s5.colnames, col_types=sup.fig.s5.coltypes)

rename_error <- function(input_df) {
  input_df$species <- ifelse(input_df$species=="zebra finch", "Zebra finch", input_df$species)
  input_df$species <- ifelse(input_df$species=="anna's hummingbird", "Anna's hummingbird", input_df$species)
  input_df$species <- ifelse(input_df$species=="platypus", "Platypus", input_df$species)
  input_df$species <- ifelse(input_df$species=="climbing perch", "Climbing perch", input_df$species)
  return(input_df)
}


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


sup.fig.s5_s7.order <- c(upstream_3kb_order, upstream_3kb_inside_order, downstream_3kb_inside_order, downstream_3kb_order)
sup.fig.s5_s7.input_df <- all_gene_record_df %>% 
  filter(position %in% sup.fig.s5_s7.order) %>% 
  filter(gene_biotype %in% "protein_coding")
sup.fig.s5_s7.input_df <- as.data.frame(sup.fig.s5_s7.input_df)
sup.fig.s5_s7.input_df$num_position <- as.numeric(factor(sup.fig.s5_s7.input_df$position, levels=sup.fig.s5_s7.order))
sup.fig.s5_s7.input_df <- rename_error(sup.fig.s5_s7.input_df)

#########################
####### (2) Indel #######
#########################

sup.fig.s5_s7.input_df$indel_type <- ifelse(sup.fig.s5_s7.input_df$indel >= 1, "indel", "Normal")
# for each type of gene of each species, calculate hompol and missing ratio of each position 
hompol_n_indel_summary <- data.frame()
for (species_s in species_order) {
  each_species_n_genetype <- sup.fig.s5_s7.input_df %>%
    filter(species==species_s) 
  each_species_n_genetype <- as.data.frame(each_species_n_genetype)
  for (position_s in sup.fig.s5_s7.order) {
    each_pos_blocks <- subset(each_species_n_genetype, position==position_s)
    # missinb blocks
    indel_pos_blocks <- subset(each_pos_blocks, indel_type=="indel")
    present_pos_blocks <- subset(each_pos_blocks, indel_type=="Normal")
    #hompol content (mean, sd)
    indel_hompol_mean <- mean(indel_pos_blocks$hompol)
    indel_hompol_sd <- sd(indel_pos_blocks$hompol)
    normal_hompol_mean <- mean(present_pos_blocks$hompol)
    normal_hompol_sd <- sd(present_pos_blocks$hompol)
    # missing ratio
    indel_ratio <- sum(indel_pos_blocks$indel)*100/(length(each_pos_blocks$indel)*100)
    hompol_n_indel_summary <- rbind(hompol_n_indel_summary, c(species_s, position_s, indel_hompol_mean, normal_hompol_mean,
                                                      indel_hompol_sd, normal_hompol_sd, indel_ratio))
  }}

# set column names + convert as numeric
colnames(hompol_n_indel_summary) <- c("species", "position", "indel_hompol_mean", "normal_hompol_mean", "indel_hompol_sd", "normal_hompol_sd", "indel_ratio")
cols.num <- c("indel_hompol_mean", "normal_hompol_mean","indel_hompol_sd", "normal_hompol_sd", "indel_ratio")
hompol_n_indel_summary[cols.num] <- sapply(hompol_n_indel_summary[cols.num],as.numeric)
hompol_n_indel_summary$general_type <- ifelse(hompol_n_indel_summary$position %in% c(upstream_3kb_order,upstream_3kb_inside_order), "5'", "3'") 
hompol_n_indel_summary$general_type <- factor(hompol_n_indel_summary$general_type , levels=c("5'", "3'"))
hompol_n_indel_summary$species <- factor(hompol_n_indel_summary$species, levels=species_order)


sup.fig.s5.maker<- function(summary_df,species_s){
  summary_df <- summary_df %>%
    filter(species==species_s)
  summary_df$num_position <- as.numeric(factor(summary_df$position, levels=sup.fig.s5_s7.order))
  print(summary_df$num_position)
  # color palette
  color_palette <- c("0_indel_hompol_mean"="#d42700", "1_normal_hompol_mean"="#9c9c9c", "2_indel_ratio"="#1a73c7")
  if (species_s == "Climbing perch") {coeff=50} else if (species_s=="Anna's hummingbird") {coeff=1600} else {coeff=500}
  sup.fig.s5.plot <- ggplot(summary_df) + 
    facet_wrap(.~general_type, ncol=2, scales="free_x")+
    # hompol content of indel / normal blocks (mean + S.D.))
    geom_line(aes(x=num_position, y=indel_hompol_mean, color="0_indel_hompol_mean"),lwd=0.3) + 
    geom_line(aes(x=num_position, y=normal_hompol_mean, color="1_normal_hompol_mean"),lwd=0.3) +
    # indel frequency
    geom_line(aes(x=num_position, y=indel_ratio*coeff, color="2_indel_ratio"), lwd=0.3) + 
    geom_vline(data=filter(summary_df, general_type=="5'"), aes(xintercept=30.5), color="grey",linetype="dotted",lwd=0.3, show.legend = FALSE) + 
    geom_vline(data=filter(summary_df, general_type=="3'"),aes(xintercept=90.5), color="grey",linetype="dotted",lwd=0.3, show.legend = FALSE) + 
    scale_colour_manual(name="",values=color_palette, labels=c("%homopolymer(With false indel)", "%homopolymer(Without false indel)", "%False indel frequency")) + 
    scale_y_continuous(limits =c(0,100),name = "homopolymer content (%)", sec.axis = sec_axis(~./coeff, name="False indel frequency (%)"))+
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
  sup.fig.s5.legend = gtable_filter(ggplotGrob(sup.fig.s5.plot), "guide-box")
  sup.fig.s5.legend_path = "/disk1/juwan_d1/cat_rerun/1_Main_figures/Sup/sup.fig.S5.legend.png"
  ggsave(sup.fig.s5.legend_path, sup.fig.s5.legend, width=6, height=2,units = "in")
  sup.fig.s5.png <-paste0("/disk1/juwan_d1/cat_rerun/1_Main_figures/Sup/sup.fig.S5.",species_s,".png")
  ggsave(sup.fig.s5.png, sup.fig.s5.plot + theme(legend.position="none"), width=3, height=2)
}

# generate plots
finch_indel_pcg <- sup.fig.s5.maker(summary_df=hompol_n_indel_summary,species_s="Zebra finch")
finch_indel_pcg <- sup.fig.s5.maker(summary_df=hompol_n_indel_summary,species_s="Anna's hummingbird")

#################################
####### Sup.fig S5 (done) #######
#################################


##########################
####### Sup.fig S7 #######
##########################

sup.fig.s5_s7.input_df$missing_type <- ifelse(sup.fig.s5_s7.input_df$missing_percent > 90, "Missing", "Present")
repeat_n_missing_summary <- data.frame()
for (species_s in species_order) {
  each_species_n_genetype <- sup.fig.s5_s7.input_df %>%
    filter(species==species_s) 
  each_species_n_genetype <- as.data.frame(each_species_n_genetype)
  for (position_s in sup.fig.s5_s7.order) {
    each_pos_blocks <- subset(each_species_n_genetype, position==position_s)
    # missinb blocks
    missing_pos_blocks <- subset(each_pos_blocks, missing_type=="Missing")
    present_pos_blocks <- subset(each_pos_blocks, missing_type=="Present")
    #hompol content (mean, sd)
    missing_repeat_mean <- mean(missing_pos_blocks$softmasked)
    missing_repeat_sd <- sd(missing_pos_blocks$softmasked)
    present_repeat_mean <- mean(present_pos_blocks$softmasked)
    present_repeat_sd <- sd(present_pos_blocks$softmasked)
    # missing ratio
    missing_ratio <- sum(missing_pos_blocks$missing_percent)*100/(length(each_pos_blocks$missing_percent)*100)
    repeat_n_missing_summary <- rbind(repeat_n_missing_summary, c(species_s, position_s, missing_repeat_mean, present_repeat_mean,
                                                                  missing_repeat_sd, present_repeat_sd, missing_ratio))
  }}

colnames(repeat_n_missing_summary) <- c("species", "position", "missing_repeat_mean", "present_repeat_mean", "missing_repeat_sd", "present_repeat_sd", "missing_ratio")
cols.num <- c("missing_repeat_mean", "present_repeat_mean","missing_repeat_sd", "present_repeat_sd", "missing_ratio")
repeat_n_missing_summary[cols.num] <- sapply(repeat_n_missing_summary[cols.num],as.numeric)
repeat_n_missing_summary$general_type <- ifelse(repeat_n_missing_summary$position %in% c(upstream_3kb_order,upstream_3kb_inside_order), "5'", "3'") 
repeat_n_missing_summary$general_type <- factor(repeat_n_missing_summary$general_type , levels=c("5'", "3'"))
repeat_n_missing_summary$species <- factor(repeat_n_missing_summary$species, levels=species_order)
repeat_n_missing_summary$header <- paste0(repeat_n_missing_summary$species," (",repeat_n_missing_summary$general_type,")")
repeat_n_missing_summary$header <- factor(repeat_n_missing_summary$header, levels=c("Zebra finch (5')", "Zebra finch (3')", "Anna's hummingbird (5')", "Anna's hummingbird (3')", 
                                                                                    "Platypus (5')", "Platypus (3')", "Climbing perch (5')", "Climbing perch (3')"))

sup.fig.S7.maker<- function(repeat_n_missing_summary){
  repeat_n_missing_summary$num_position <- as.numeric(factor(repeat_n_missing_summary$position, levels=sup.fig.s5_s7.order))
  color_palette <- c("0_missing_repeat_mean"="#d42700", "1_present_repeat_mean"="#9c9c9c", "2_missing_ratio"="#1a73c7")
  sup.fig.S7.plot <- ggplot(repeat_n_missing_summary) + 
    facet_wrap(.~header, ncol=4, scales="free_x") +
    # Repeat content of indel / normal blocks (mean + S.D.))
    geom_line(aes(x=num_position, y=missing_repeat_mean, color="0_missing_repeat_mean", group=species), lwd=0.3) + 
    geom_line(aes(x=num_position, y=present_repeat_mean, color="1_present_repeat_mean", group=species), lwd=0.3) +
    #geom_point(aes(x=num_position, y=missing_repeat_mean, color="0_missing_repeat_mean", group=species), size=0.7) + #shape=species, 
    #geom_point(aes(x=num_position, y=present_repeat_mean, color="1_present_repeat_mean", group=species), size=0.7) + 
    # Missing ratio
    geom_line(aes(x=num_position, y=missing_ratio, color="2_missing_ratio", group=species), lwd=0.3) + 
    #geom_point(aes(x=num_position, y=missing_ratio, color="2_missing_ratio", group=species), size=0.7) +  #shape=species, 
    geom_vline(data=filter(repeat_n_missing_summary, general_type=="5'"), aes(xintercept=30.5), color="grey",linetype="dotted",lwd=0.3, show.legend = FALSE) + 
    geom_vline(data=filter(repeat_n_missing_summary, general_type=="3'"),aes(xintercept=90.5), color="grey",linetype="dotted",lwd=0.3, show.legend = FALSE) + 
    scale_colour_manual(name=" ",values=color_palette, labels=c("repeat content (Missing)", "repeat content (Present)", "Missing ratio")) + 
    scale_y_continuous(limits =c(0,100),name = "Repeat content (%)", sec.axis = sec_axis(~., name="Missing ratio (%)"))+
    scale_x_continuous(breaks=c(1,15,30.5,45,60,61,75,90.5,105,120),
                       labels=c("-3","-1.5","TSS","+1.5","+3","-3","-1.5","TTS","+1.5","+3"),
                       expand=c(0,0)) + 
    labs(y = "Percentage (%)")+
    theme_pubr() + 
    theme(text=element_text(size=6),
          strip.text = element_text(size=7),
          #strip.text.x = element_blank(),
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
  sup.fig.S7.legend = gtable_filter(ggplotGrob(sup.fig.S7.plot), "guide-box")
  sup.fig.S7.legend_path = "/disk1/juwan_d1/cat_rerun/1_Main_figures/Sup/sup.fig.S7.all.legend.png"
  ggsave(sup.fig.S7.legend_path, sup.fig.S7.legend, width=6, height=2,units = "in")
  sup.fig.S7.png <-paste0("/disk1/juwan_d1/cat_rerun/1_Main_figures/Sup/sup.fig.S7.png")
  ggsave(sup.fig.S7.png, sup.fig.S7.plot + theme(legend.position="none"), width=6, height=3)
}

sup.fig.S7.maker(repeat_n_missing_summary)
