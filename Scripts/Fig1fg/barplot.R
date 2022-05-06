#install.packages("ggplot2")
install.packages("gridExtra")

# Load library
library(ggplot2)
library("gridExtra")

# Set working directory
setwd("C:\\Users\\swear0712\\Desktop\\FGL\\barplot_R\\")

# CpG islands
data_CpG <- read.table("CpGisland.txt",header = T, sep="\t")
attach(data_CpG)
fig_cpg <- ggplot(data_CpG, aes(Species, Missing, fill = factor(Type,levels=c("CpG_islands","Non-CpG_islands"))))+ 
  geom_bar(stat="identity",position = "dodge") +
  scale_x_discrete(limits = c("Zebrafinch","Hummingbird","Platypus","Climbing_perch"))+
  scale_fill_manual(values = c("darkred","darkgrey"),name="Type") +
  ylab("Missing rate (%)") +
  geom_text(aes(label=Missing), vjust=-0.2, position=position_dodge(.9))
detach(data_CpG)


# Repeat
data_repeat <- read.table("repeat.txt",header = T, sep="\t")
attach(data_repeat)
fig_repeat <- ggplot(data_repeat, aes(Species, Missing, fill = factor(Type,levels=c("Repeats","Non-repeats"))))+ 
  geom_bar(stat="identity",position = "dodge") +
  scale_x_discrete(limits = c("Zebrafinch","Hummingbird","Platypus","Climbing_perch"))+
  scale_fill_manual(values = c("darkblue","darkgrey"),name="Type") +
  ylab("Missing rate (%)")+
  geom_text(aes(label=Missing), vjust=-0.2, position=position_dodge(.9))
detach(data_repeat)


grid.arrange(fig_cpg , fig_repeat, ncol = 2, nrow = 1)#, layout_matrix = rbind(c(1,1), c(2,3)))


# Read and Assembly gap
data_read_gap <- read.table("read_assemblygap.txt",header = T, sep="\t")
attach(data_read_gap)
fig_read_gap <- ggplot(data_read_gap, aes(Species, SupportedRate, fill = factor(Type,levels=c("Missing_PreReadDepth","Existing_PreReadDepth", "Missing_PreAsmGap","Existing_PreAsmGap","Missing_PreUnion","Existing_PreUnion"))))+ 
  geom_bar(stat="identity",position = "dodge") +
  scale_x_discrete(limits = c("Platypus","Climbing_perch"))+
  scale_fill_manual(values = c("darkblue","darkgrey","darkblue","darkgrey","darkblue","darkgrey"),name="Type") +
  ylab("Supported rate (%)")
detach(data_read_gap)
fig_read_gap