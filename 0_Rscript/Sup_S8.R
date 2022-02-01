library(ggplot2)
#install.packages("devtools")
#library(devtools)
#devtools::install_github("zeehio/facetscales")
#library(facetscales)

##########################
####### Sup.fig.S8 #######
##########################

# basic input variable
species_order <- c("Zebra finch","Anna's hummingbird", "Platypus", "Climbing perch")
error_order <- c("Error in\nprior assembly (FGL)", "Error in\nVGP assembly", "No errors detected", "Filtered out")
# rename directory
rename_error <- function(input_df) {
  input_df$error_type <- ifelse(input_df$error_type=="previous_error", "Error in\nprior assembly (FGL)", input_df$error_type)
  input_df$error_type <- ifelse(input_df$error_type=="vgp_error", "Error in\nVGP assembly", input_df$error_type)
  input_df$error_type <- ifelse(input_df$error_type=="no_errors", "No errors detected", input_df$error_type)
  input_df$error_type <- ifelse(input_df$error_type=="filter_out", "Filtered out", input_df$error_type)
  input_df$FGL_type <- ifelse(input_df$FGL_type=="frameshift", "Frameshift", input_df$FGL_type)
  input_df$FGL_type <- ifelse(input_df$FGL_type=="prematurestopcodon", "Premature stop codon", input_df$FGL_type)
  input_df$FGL_type <- ifelse(input_df$FGL_type=="intronexonjunctiondisruption", "Splicing junction disruption", input_df$FGL_type)
  return(input_df)
}

# load input dataframe
input_path <- "/disk1/juwan_d1/cat_rerun/zebra_finch/tm/CAT2MISSING_pre_taegut/supplementary.count_variants.tabs"
read_mapping_count_df <- read.csv(input_path, header=TRUE, sep="\t")
read_mapping_count_df$species <- factor(read_mapping_count_df$species,levels=species_order)
read_mapping_count_df <- rename_error(input_df=read_mapping_count_df)
read_mapping_count_df$error_type <- factor(read_mapping_count_df$error_type,levels=error_order)
# generate plot
sup.fig.S8.plot <- ggplot(read_mapping_count_df, aes(x=error_type,y=count,fill=error_type, group=species)) + 
  facet_wrap(species~FGL_type, scales="free",ncol=3) + 
  labs(y="Number of variants")+
  scale_x_discrete(labels=error_order)+
  scale_y_continuous(expand=c(0,0)) + 
  scale_fill_manual(values=c("#990000","#009499","#b0b0b0","#474747")) + 
  geom_bar(stat="identity", width = 0.55) +
  labs(fill="") +
  theme_pubr() + 
  theme(text=element_text(size=6),
        axis.text.x=element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position="none",
        legend.key.width = unit(0.3,"cm"),
        legend.key.height = unit(0.3,"cm"),
        axis.line.x.top=element_blank(),
        axis.ticks.length.x.top = unit(0,"cm"),
        axis.line = element_line(size=0.1),
        axis.ticks = element_line(size=0.1),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) 

sup.fig.S8.png <- "/disk1/juwan_d1/cat_rerun/1_Main_figures/Sup/sup.fig.S8.png"
ggsave(sup.fig.S8.png, sup.fig.S8.plot , width=3, height=3)

#################################
####### Sup.fig.S8 (done) #######
#################################