##########################################################################
## Setup and Run cellex and Cellect tool in linux to get the prioritization csv file.
## https://github.com/perslab/CELLEX 
## https://github.com/perslab/CELLECT

######## For cellex - bubble plot ######## 
library(ggpubr)
dat <- read.csv("prioritization-08112021.csv")
colnames(dat)
names(dat)[ncol(dat)] <- "-log10(Pvalue)"
names(dat)
ggballoonplot(dat, x = "GWAS", y = "Cells",
              fill = "-log10(Pvalue)", 
              color = "lightgray",
              size = "-log10(Pvalue)")+
  gradient_fill(c("white", "grey", "blue", "darkblue")) +  ggtitle("") + xlab("GWAS traits") + ylab("Cell Types") +
  theme(axis.title = element_text(face = "bold", color = "black", size = 16)) +
  theme(axis.text.x = element_text(face = "bold", color = "black", size = 10)) +
  theme(axis.text.y = element_text(face = "bold", color = "black", size = 10))
ggsave("cellex-plot.pdf", width = 5, height = 5, dpi = 600)
##########################################################################