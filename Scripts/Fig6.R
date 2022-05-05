library("RColorBrewer")
library("ggplot2")
library("ggrepel")
library(plyr)
library(reshape2)
library(tidyr)

############################
######## Fig.6.a-h #########
############################

fgl_names <- c("Completely missing","Exon deletion","Fragmented","Intra-scaffold split",
              "Frameshift","Premature stop codon","Splicing junction disruption", "Ns in coding region")
species_names <- c("Zebra finch", "Platypus", "Anna's hummingbird",  "Climbing perch")
default_colors <- c("#F8766D","#7CAE00","#00BFC4","#C77CFF")
default_colors <- c("#592800", "#b35600",  "#00695c", "#4db8ab")

# convert species names
fgl_name_conversion <- function(input_df) {
  input_df$type <- ifelse(input_df$type == "totallymissing", "Completely missing", input_df$type)
  input_df$type <- ifelse(input_df$type == "exonDeletion", "Exon deletion", input_df$type)
  input_df$type <- ifelse(input_df$type == "fragmented", "Fragmented", input_df$type)
  input_df$type <- ifelse(input_df$type == "intrascaffoldsplit", "Intra-scaffold split", input_df$type)
  input_df$type <- ifelse(input_df$type == "frameshift", "Frameshift", input_df$type)
  input_df$type <- ifelse(input_df$type == "prematurestopcodon", "Premature stop codon", input_df$type)
  input_df$type <- ifelse(input_df$type == "intronexonjunctiondisruption", "Splicing junction disruption", input_df$type)
  input_df$type <- ifelse(input_df$type == "nincodingregion", "Ns in coding region", input_df$type)
  return(input_df)
}

# load input dataframe
zf_path <- "./input/zebra_finch.total.count.tsv"
an_path <- "./input/anna.total.count.tabs"
pl_path <- "./input/platypus.total.count.tabs"
cl_path <- "./input/climbing_perch.total.count.tabs"
zf_df <- read.csv(file=zf_path,header=FALSE,sep="\t",col.names = c("type","count"))
an_df <- read.csv(file=an_path,header=FALSE,sep="\t",col.names = c("type","count"))
pl_df <- read.csv(file=pl_path,header=FALSE,sep="\t",col.names = c("type","count"))
cl_df <- read.csv(file=cl_path,header=FALSE,sep="\t",col.names = c("type","count"))
zf_df$species <- "Zebra finch"
an_df$species <- "Anna's hummingbird"
pl_df$species <- "Platypus"
cl_df$species <- "Climbing perch"
fig.6.a_h.input_df <-fgl_name_conversion(rbind(zf_df, an_df, pl_df, cl_df))
fig.6.a_h.input_df$type <- factor(fig.6.a_h.input_df$type,levels=fgl_names)
fig.6.a_h.input_df$species <- factor(fig.6.a_h.input_df$species,levels =species_names)

# draw plots
fig6.a_h <- ggplot(fig.6.a_h.input_df, aes(x=species, y=count, fill=species)) +
  facet_wrap(~type, scales = "free", ncol=2, dir="v") +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label=paste0(count,"\n")), position = "identity",size=2,vjust = 0, color = "darkslategrey") +
  scale_fill_manual(values=default_colors, labels=species_names) +
  scale_y_continuous(expand = c(0, 0, 0.3 ,0))+
  labs(x="Number of genes affected by false gene loss", y="", fill="Species")+
  theme_classic() +
  theme(text = element_text(size=8),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        strip.text=element_blank(),
        strip.background = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5)) +
  guides(fill = guide_legend(nrow = 1))

# save plots
fig6.a_h.filepath <- "./fig6.a_h.png"
ggsave(fig6.a_h.filepath, fig6.a_h + theme(legend.position = "none"), width=3, height=4.5)
# save legend
fig6.a_h.legend_filepath <- "./fig6.a_h.legend.png"
fig6.a_h.legend = gtable_filter(ggplotGrob(fig6.a_h), "guide-box") 
ggsave(fig6.a_h.legend_filepath, fig6.a_h.legend, width=3, height=1)

###################################
######## Fig.6.a-h (done) #########
###################################


############################
######### Fig.6.i ##########
############################

# input dataframes
zf_ratio_path <- "./input/zebra_finch.total.summary.tsv"
an_ratio_path <- "./input/anna.total.summary.tsv"
pl_ratio_path <- "./input/platypus.total.summary.tsv"
cl_ratio_path <- "./input/climbing_perch.total.summary.tsv"
zf_ratio_df <- read.csv(file=zf_ratio_path,header=TRUE,sep="\t")
an_ratio_df <- read.csv(file=an_ratio_path,header=TRUE,sep="\t")
pl_ratio_df <- read.csv(file=pl_ratio_path,header=TRUE,sep="\t")
cl_ratio_df <- read.csv(file=cl_ratio_path,header=TRUE,sep="\t")
zf_ratio_df$species <- "Zebra finch"
an_ratio_df$species <- "Anna's hummingbird"
pl_ratio_df$species <- "Platypus"
cl_ratio_df$species <- "Climbing perch"
rev_species_names <- rev(c("Zebra finch", "Platypus", "Anna's hummingbird",  "Climbing perch"))
fig6.i_df <- rbind(zf_ratio_df, an_ratio_df, pl_ratio_df, cl_ratio_df)
fig6.i_df$species <- factor(fig6.i_df$species,levels = rev_species_names)
levels(factor(fig6.i_df$species))
fig6.i_df$type <-factor(df$type,levels = c("normal","sequence","structure", "structure_and_sequence"))
# draw plot
fig6.i <- ggplot(fig6.i_df, aes(x=species, y=count, fill=interaction(type, species))) +
  geom_bar(position = "stack",stat="identity") +
  scale_x_discrete(labels=rev_species_names) + 
  theme_classic()+
  theme(text = element_text(size=7),
        legend.position = "none",
        legend.text=element_text(size=7),
        legend.background = element_blank(),
        legend.title = element_blank(),
        axis.title.y = element_blank()) +
  scale_y_continuous(expand = c(0, 0))+
  scale_fill_manual(values=c("#bababa","#4db8ab","#4db8ab","#4db8ab",
                             "#bababa","#00695c","#00695c","#00695c",
                             "#bababa","#b35600","#b35600","#b35600",
                             "#bababa","#592800","#592800","#592800"))+
  labs(x="Species",y="Number of genes", title="Ratio of genes including false gene loss")+
  coord_flip()

fig6.i.filepath <- "./fig6.i.png"
ggsave(fig6.i.filepath, fig6.i, width=5, height=1.5)

###################################
######### Fig.6.i (done) ##########
###################################
