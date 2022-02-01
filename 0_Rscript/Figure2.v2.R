library("OmicCircos")
options(stringsAsFactors=FALSE)
set.seed(1234)
trace(circos,edit=T)

# mapping data
chr_df <- read.csv("/disk1/juwan_d1/cat_rerun/zebra_finch/tm/CAT2MISSING_pre_taegut/chromosomeMap/out_circos/chromosomes.size.csv")
gc_df <- read.csv("/disk1/juwan_d1/cat_rerun/zebra_finch/tm/CAT2MISSING_pre_taegut/chromosomeMap/out_circos/chromosomes.gc.csv")
repeat_df <- read.csv("/disk1/juwan_d1/cat_rerun/zebra_finch/tm/CAT2MISSING_pre_taegut/chromosomeMap/out_circos/chromosomes.repeat.csv")
missing_df <- read.csv("/disk1/juwan_d1/cat_rerun/zebra_finch/tm/CAT2MISSING_pre_taegut/chromosomeMap/out_circos/chromosomes.missing_loci.csv")
label_df <- read.csv("/disk1/juwan_d1/cat_rerun/zebra_finch/tm/CAT2MISSING_pre_taegut/chromosomeMap/out_circos/chromosomes.missing_label.csv")
chr_names_df <- read.csv("/disk1/juwan_d1/cat_rerun/zebra_finch/tm/CAT2MISSING_pre_taegut/chromosomeMap/out_circos/chromosomes.name.csv")
gaps_df <- read.csv("/disk1/juwan_d1/cat_rerun/zebra_finch/tm/CAT2MISSING_pre_taegut/chromosomeMap/out_circos/vgp_gap_plot.csv")
alignment_df <- read.csv("/disk1/juwan_d1/cat_rerun/zebra_finch/tm/CAT2MISSING_pre_taegut/chromosomeMap/out_circos/chromosomes.aln.csv")
genes_df <- read.csv("/disk1/juwan_d1/cat_rerun/zebra_finch/tm/CAT2MISSING_pre_taegut/chromosomeMap/out_circos/vgp_genes_plot.csv")

######## plot1 #########
plot1_chr_order <- c("NC_044213","NW_022045306", "NW_022045337",
                     "NC_044211","NC_044214","NW_022045379", "NC_044241","NC_044215","NC_044212","NW_022045339", "NC_044217", "NW_022045287","NC_044219", 
                     "NC_044218","NC_044220","NC_044221","NC_044222","NC_044223","NW_022045299","NC_044224","NC_044216","NC_044225", 
                     "NC_044226","NC_044232","NC_044227","NC_044229","NC_044230","NC_044233","NC_044231","NC_044236",
                     "NC_044238", "NW_022045315","NC_044235")
plot1_chr_order <- c("NC_044213", "NW_022045306", "NW_022045337", "NC_044211", "NC_044214", "NW_022045379",
                     "NC_044241", "NC_044215", "NC_044212", "NW_022045339", "NC_044217", "NW_022045287", "NC_044219",
                     "NC_044218", "NC_044220", "NC_044221", "NC_044222", "NC_044223", "NW_022045299", "NC_044224", "NC_044216",
                     "NC_044225", "NC_044226", "NC_044232", "NC_044227", "NC_044229", "NC_044230")

plot1_chr_df <- chr_df %>% 
  filter(chr %in% plot1_chr_order)
plot1_gc_df <- gc_df %>% 
  filter(chr %in% plot1_chr_order)
plot1_repeat_df <- repeat_df %>% 
  filter(chr %in% plot1_chr_order)
plot1_missing_df <- missing_df %>% 
  filter(chr %in% plot1_chr_order)
plot1_label_df <- label_df %>% 
  filter(chr %in% plot1_chr_order)
plot1_chr_names_df <- chr_names_df %>% 
  filter(chr %in% plot1_chr_order)
plot1_gaps_df <- gaps_df %>% 
  filter(chr %in% plot1_chr_order)
plot1_alignment_df <- alignment_df %>% 
  filter(chr %in% plot1_chr_order)
plot1_genes_df <- genes_df %>% 
  filter(chr %in% plot1_chr_order)

plot1_seg.name <- plot1_chr_order
plot1_chr_db <- segAnglePo(plot1_chr_df,seg=plot1_seg.name)
plot1_chr_colors <- ifelse(grepl("^NC_",plot1_chr_order), "#00c4a7", "#00c4a7")
plot1_line_colors <- c("black")
plot1_gaps_df$num_chrom <- as.numeric(factor(plot1_gaps_df$chr,levels=plot1_chr_order))
ordered_plot1_gaps_df <- plot1_gaps_df[order(plot1_gaps_df$num_chrom),]
length(ordered_plot1_gaps_df$chr)
plot1_gap_colors <- ifelse(grepl("^NC_",ordered_plot1_gaps_df$chr), "#9bc7c0", "#a493cc")

plot1_chr_names_df$value <- paste0(" ", plot1_chr_names_df$value)

Fig2.a.pdf <- "/disk1/juwan_d1/cat_rerun/1_Main_figures/Fig2/fig2.a.pdf"
pdf(Fig2.a.pdf,8,8)
par(mar=c(0.5,0.5,0.5,0.5))
plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="", main="");
circos(cir=plot1_chr_db,R=350,W=200,type="chr", col=plot1_chr_colors, print.chr.lab=FALSE,scale=TRUE,order=plot1_chr_order,lwd=1,cex=5)
circos(cir=plot1_chr_db,R=360,W=0,mapping=plot1_chr_names_df,type="label", side="out",col=c("black"),order=plot1_chr_order, cex=1.2, lwd = 0)
circos(cir=plot1_chr_db,R=315,W=0,mapping=plot1_label_df,type="label", side="in", col=c("black"),order=plot1_chr_order, cex=1.2, lwd = 0)
circos(cir=plot1_chr_db,R=100,W=30,mapping=plot1_alignment_df,type="heatmap_for_aln",col=c("black"),order=plot1_chr_order, B=FALSE, lwd= 1)
circos(cir=plot1_chr_db,R=140,W=15,mapping=plot1_missing_df,type="b3",col=c("black"),order=plot1_chr_order, cex=0.4, lwd= 0)
circos(cir=plot1_chr_db,mapping=plot1_gc_df,R=210,W=60,type="b2", cutoff=0.42, col=c("#fc3f6b","#bfbfbf"),lwd=0,order=plot1_chr_order, B=FALSE,scale=FALSE)
circos(cir=plot1_chr_db,mapping=plot1_repeat_df,R=170, W=60, type="b2",cutoff=0.20, col=c("#5c8dff","#bbc0c4"), lwd=0,order=plot1_chr_order, B=FALSE,scale=FALSE)
circos(cir=plot1_chr_db,mapping=plot1_genes_df, R=150, W=30, type="l",col=c("#ff752b"),order=plot1_chr_order, cex=0.4, lwd= 0.5, B=FALSE,scale = FALSE)
circos(cir=plot1_chr_db,mapping=ordered_plot1_gaps_df,R=335.3,W=32.7,type="b3",col=plot1_gap_colors,lwd=1,order=plot1_chr_order, B=FALSE,scale=FALSE,cutoff=0)
dev.off()

######## plot2 #########
plot2_chr_order <- c("NC_044228", "NW_022045349", "NC_044234", "NC_044237", "NC_044239", "NC_044240", "NC_044242",
                     "NW_022045319", "NC_044243", "NW_022045292", "NW_022045378", "NW_022045341", "NW_022045351", "NW_022045350",
                     "NW_022045362", "NW_022045374", "NW_022045293", "NW_022045327", "NW_022045357")
plot2_chr_order <- c("NC_044233", "NC_044231", "NC_044236", "NC_044238", "NW_022045315", "NC_044235", "NC_044240",
                     "NC_044243", "NC_044234", "NC_044239", "NC_044242", "NW_022045319", "NW_022045292", "NW_022045378",
                     "NC_044237", "NW_022045341", "NW_022045351", "NW_022045357", "NW_022045350", "NW_022045362", 
                     "NC_044228", "NW_022045349", "NW_022045293", "NW_022045374", "NW_022045327")

plot2_chr_df <- chr_df %>% 
  filter(chr %in% plot2_chr_order)
plot2_gc_df <- gc_df %>% 
  filter(chr %in% plot2_chr_order)
plot2_repeat_df <- repeat_df %>% 
  filter(chr %in% plot2_chr_order)
plot2_missing_df <- missing_df %>% 
  filter(chr %in% plot2_chr_order)
plot2_label_df <- label_df %>% 
  filter(chr %in% plot2_chr_order)
plot2_chr_names_df <- chr_names_df %>% 
  filter(chr %in% plot2_chr_order)
plot2_gaps_df <- gaps_df %>% 
  filter(chr %in% plot2_chr_order)
plot2_alignment_df <- alignment_df %>% 
  filter(chr %in% plot2_chr_order)
plot2_genes_df <- genes_df %>% 
  filter(chr %in% plot2_chr_order)

plot2_seg.name <- plot2_chr_order
plot2_chr_db <- segAnglePo(plot2_chr_df,seg=plot2_seg.name)
plot2_chr_colors <- c("#00c4a7", "#00c4a7", "#00c4a7", "#00c4a7", "#a47dff", "#00c4a7",
                      "#00c4a7", "#00c4a7", "#00c4a7", "#00c4a7", "#00c4a7", "#a47dff", 
                      "#00c4a7", "#a47dff", "#00c4a7", "#00c4a7", "#00c4a7", "#00c4a7", 
                      "#00c4a7", "#a47dff", "#00c4a7", "#a47dff", "#a47dff", "#00c4a7", 
                      "#a47dff")
plot2_chr_colors <- c("#00c4a7", "#00c4a7", "#00c4a7", "#00c4a7", "#00c4a7", "#00c4a7",
                      "#00c4a7", "#a47dff", "#00c4a7", "#00c4a7", "#a47dff", "#a47dff", 
                      "#a47dff", "#a47dff", "#00c4a7", "#a47dff", "#a47dff", "#a47dff", 
                      "#a47dff", "#a47dff", "#00c4a7", "#00c4a7", "#a47dff", "#a47dff", 
                      "#a47dff")
plot2_line_colors <- c("black")
plot2_gaps_df$num_chrom <- as.numeric(factor(plot2_gaps_df$chr,levels=plot2_chr_order))
ordered_plot2_gaps_df <- plot2_gaps_df[order(plot2_gaps_df$num_chrom),]
length(ordered_plot2_gaps_df$chr)
plot2_chromosome_level_scaffolds <- c("NC_044233", "NC_044231", "NC_044236", "NC_044238", "NC_044235", "NC_044240",
                                      "NC_044243", "NC_044234", "NC_044239", "NC_044242", "NW_022045292", 
                                      "NC_044237", "NW_022045341", "NW_022045351", "NW_022045357", "NW_022045350", 
                                      "NC_044228", "NW_022045374")
plot2_chromosome_level_scaffolds <- c("NC_044233", "NC_044231", "NC_044236", "NC_044238", "NW_022045315", "NC_044235", 
                                      "NC_044240",  "NC_044234", "NC_044239",
                                       "NC_044237", 
                                     "NC_044228", "NW_022045349")
plot2_gap_colors <- ifelse(ordered_plot2_gaps_df$chr %in% plot2_chromosome_level_scaffolds, "#9bc7c0", "#a493cc")

head(ordered_plot2_gaps_df)

Fig2.b.pdf <- "/disk1/juwan_d1/cat_rerun/1_Main_figures/Fig2/fig2.b.pdf"
pdf(Fig2.b.pdf,8,8)
par(mar=c(0.5,0.5,0.5,0.5))
plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="", main="");
circos(cir=plot2_chr_db,R=350,W=200,type="chr", col=plot2_chr_colors, print.chr.lab=FALSE,scale=TRUE,order=plot2_chr_order,lwd=1,cex=1)
circos(cir=plot2_chr_db,R=360,W=0,mapping=plot2_chr_names_df,type="label", side="out",col=c("black"),order=plot2_chr_order, cex=1.2, lwd = 0)
circos(cir=plot2_chr_db,R=315,W=0,mapping=plot2_label_df,type="label", side="in", col=c("black"),order=plot2_chr_order, cex=1.2, lwd = 0)
circos(cir=plot2_chr_db,R=110,W=30,mapping=plot2_alignment_df,type="heatmap_for_aln",col=c("black"),order=plot2_chr_order, B=FALSE, cex=0.4, lwd= 1)
circos(cir=plot2_chr_db,R=140,W=15,mapping=plot2_missing_df,type="b3",col=c("black"),order=plot2_chr_order, cex=0.4, lwd= 0)
circos(cir=plot2_chr_db,mapping=plot2_gc_df,R=210,W=60,type="b2", cutoff=0.42, col=c("#fc3f6b","#bfbfbf"),lwd=0,order=plot2_chr_order, B=FALSE,scale=FALSE)
circos(cir=plot2_chr_db,mapping=plot2_repeat_df,R=170, W=60, type="b2",cutoff=0.20, col=c("#5c8dff","#bbc0c4"), lwd=0,order=plot2_chr_order, B=FALSE,scale=FALSE)
circos(cir=plot2_chr_db,mapping=plot2_genes_df, R=150, W=30, type="l",col=c("#ff752b"),order=plot2_chr_order, cex=0.4, lwd= 0.5, B=FALSE,scale = FALSE)
circos(cir=plot2_chr_db,mapping=ordered_plot2_gaps_df,R=335.3,W=32.7,type="b3",col=plot2_gap_colors,lwd=1,order=plot2_chr_order, B=FALSE,scale=FALSE,cutoff=0)
dev.off()