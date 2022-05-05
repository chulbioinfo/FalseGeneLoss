library(dplyr)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(grid)
library(gtable)
library(gridExtra)

#################
#### Fig 8.b ####
#################

repeat_test_maker <- function (repeat_type) {
  # load dataframe
  whole_genome_path <- paste0("./input/bTaeGut2.10000.", repeat_type, ".bed") # whole genome of Trio assembly
  chr19_path <- paste0("./input/bTaeGut2.chr19.10000.", repeat_type, ".bed") # chr19 of Trio assembly
  missing_chr19_path <- paste0("./input/bTaeGut2.missing_chr19.10000.", repeat_type, ".bed") # part of chr19 of Trio assembly missing in VGP zebra finch assembly
  
  whole_genome_df <- read.csv(whole_genome_path, sep="\t", header=FALSE)
  chr19_genome_df <- read.csv(chr19_path, sep="\t", header=FALSE)
  missing_chr19_genome_df <- read.csv(missing_chr19_path, sep="\t", header=FALSE)
  
  colnames(whole_genome_df) <- c("chrom","start", "end","overlap")
  colnames(chr19_genome_df) <- c("chrom","start", "end","overlap")
  colnames(missing_chr19_genome_df) <- c("chrom","start", "end","overlap")
  missing_chr19_genome_df$genome_type <- "Missing 2.7 Mbp"
  whole_genome_df$genome_type <- "Average"
  chr19_genome_df$genome_type <- "Chr.19"
  
  whole_vs_chr19 <- rbind(whole_genome_df, chr19_genome_df, missing_chr19_genome_df)
  
  # overlap: repeat content (%)
  whole_vs_chr19$overlap <- whole_vs_chr19$overlap/10000*100
  
  # set y axis limit manually
  y_manual_lim = ifelse(repeat_type %in% c("LTR"), 20, 
                        ifelse(repeat_type == "LINE", 14,
                               ifelse(repeat_type %in% c("Simple_repeat"), 5, 
                                      ifelse(repeat_type=="Low_complexity", 1.5, 1))))
  whole_vs_chr19$genome_type <- factor(whole_vs_chr19$genome_type,levels=c("Average", "Chr.19", "Missing 2.7 Mbp"))
  print(levels(factor(whole_vs_chr19$genome_type)))
  # draw plot
  fig9.b.plot <- ggplot(whole_vs_chr19, aes(x =genome_type, y = overlap, fill = genome_type, group =genome_type, color=genome_type))+
    stat_summary(position = "dodge", fun.data = mean_sd, geom = "errorbar",color="black",lwd=0.1, width=0.5) + 
    stat_summary(position = "dodge", geom = "bar", fun = "mean", lwd=0.1, color="black",width=0.5) +
    scale_fill_manual(values=c("grey", "#5e5e5e", "brown")) + 
    coord_cartesian(ylim=c(0,y_manual_lim),xlim=c(0.5,3.5),expand=FALSE) + 
    labs(y="Repeat content (%)", title=ifelse(grepl("_", repeat_type), gsub("_","\n", repeat_type), paste0(repeat_type,"\n")), fill="Sequences    ") +
    stat_compare_means(method = "anova", label = "p.signif", size=1.5,
                       label.x=2, label.y= y_manual_lim*0.85, 
                       hide.ns=TRUE) + 
    theme_pubr() + 
    theme(text=element_text(size=6),
          title=element_text(size=6),
          axis.text=element_text(size=6),
          axis.text.x=element_blank(),
          axis.title=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position = "bottom",
          legend.key.size = unit(0.6,"line"))
  fig9.b.legend = gtable_filter(ggplotGrob(fig9.b.plot), "guide-box")
  fig9.b.legend_plot <- arrangeGrob(fig9.b.legend)
  fig9.b.filepath <- "./fig9.legend.png"
  ggsave(fig9.b.filepath, fig9.b.legend_plot, width=3, height=1,units = "in")
  fig9.b.plot <- fig9.b.plot + theme(legend.position="none")
  my_list <- list(fig9.b.plot)
  return(my_list)
}

# for each repeat type, draw plots
LINE_plot_list <- repeat_test_maker("LINE")
LTR_plot_list <- repeat_test_maker("LTR")
SINE_plot_list <- repeat_test_maker("SINE")
Simple_repeat_plot_list <- repeat_test_maker("Simple_repeat")
Low_complexity_plot_list <- repeat_test_maker("Low_complexity")

# merge all plots
fig9.plot.list <- c(LINE_plot_list, LTR_plot_list,
                    SINE_plot_list, Simple_repeat_plot_list, Low_complexity_plot_list)
my_layout<-rbind(c(1,2,3,4,5), c(1,2,3,4,5))
fig9.f.plots <- arrangeGrob(grobs = fig9.plot.list, layout_matrix = my_layout)

# save output file
fig9.f.plots_out <-"./fig9.b.png"
ggsave(fig9.f.plots_out, fig9.f.plots, width=4.0, height=1.6)

########################
#### Fig 8.b (done) ####
########################