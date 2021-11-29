## Load required libraries
library("Seurat")
library("stringi")
library("data.table")
library("ggpubr")
library("cowplot")

## Read data
local <- readRDS("fasolino_et_al.rds")

## T1D vs Control
cond = "AAB"
disease_order <- c("Control", "AAB", "T1D")
## For T1DvsControl
local_cond <- subset(local, subset = disease_state != cond)
local_cond$disease_state <- droplevels(local_cond$disease_state)
filename = "T1DvsControl_top3_ductal.pdf"
which_cell = "duct_major"
t1 <- '^TPPP3$'; t2 <- '^CXCL6$'; t3 <- '^AC013275.1$';

## Replace t1, t2, and t3 with the below top 3 genes in each cell type while performing analysis for corresponding cell type
# ## acinar: "OLFM4", "INS", "CELA2A"
# ## alpha: "REG1B", "REG3A", "INS"
# ## beta: "IFI6", "INS", "SST"
# ## delta: "REG1B", "SPINK1", "REG3A"
# ## ductal: "TPPP3", "CXCL6", "AC013275.1"

######################### Uncomment this block for T1D vs AAB ###########################
# ## T1D vs AAB
# cond = "Control"
# disease_order <- c("Control", "AAB", "T1D")
# ## For T1DvsControl
# local_cond <- subset(local, subset = disease_state != cond)
# local_cond$disease_state <- droplevels(local_cond$disease_state)
# filename = "T1DvsAAB_top3_ductal.pdf"
# which_cell = "duct_major"
# t1 <- '^JUND$'; t2 <- '^HIST1H1C$'; t3 <- '^MTRNR2L12$';
# ## Replace t1, t2, and t3 with the below top 3 genes in each cell type while performing analysis for corresponding cell type
# ## acinar: "CELA2A", "OLFM4", "SYCN"
# ## alpha: "REG1A", "REG1B", "MTRNR2L12"
# ## beta: "MTRNR2L12", "GNAS", "IFI6"
# ## delta: "MT-ATP6", "CD81", "UBE2M"
# ## ductal: "JUND", "HIST1H1C", "MTRNR2L12"

######################### Uncomment this block for AAB vs Control ###########################
## AAB vs Control
# cond = "T1D"
# disease_order <- c("Control", "AAB", "T1D")
# ## For T1DvsControl
# local_cond <- subset(local, subset = disease_state != cond)
# local_cond$disease_state <- droplevels(local_cond$disease_state)
# filename = "AABvsControl_top3_ductal.pdf"
# which_cell = "duct_major"
# t1 <- '^SNHG25$'; t2 <- '^FXYD2$'; t3 <- '^RPS17$';

## Replace t1, t2, and t3 with the below top 3 genes in each cell type while performing analysis for corresponding cell type
## acinar: "AMY2B", "AMY2A", "MTRNR2L8"
## alpha: "AC099509.1", "PTP4A3", "MTRNR2L12"
## beta: "CTRB2", "CTRC", "PRSS2"
## delta: "GNAS", "MTRNR2L12", "RPS4Y1"
## ductal: "SNHG25", "FXYD2", "RPS17"

##########################################################################################

# ylimit1 <- c(0,10); ylimit2 <- c(0,10); ylimit3 <- c(0,10)

ylimit1 <- c(0,5); ylimit2 <- c(0,5); ylimit3 <- c(0,7)

p <- expr_plot(t1, t2, t3, filename, which_cell, ylimit1, ylimit2, ylimit3)
pdf(filename, width = 8, height = 5)
p
dev.off()

expr_plot <- function(t1, t2, t3, filename, which_cell, ylimit1, ylimit2, ylimit3)
{
  local_cell <- subset(local_cond, subset = cell_label == which_cell)
  gene_var <- c(grep(t1, rownames(local_cell), value = TRUE),
                grep(t2, rownames(local_cell), value = TRUE),
                grep(t3, rownames(local_cell), value = TRUE))
  
  p <- list()
  for (i in gene_var) 
  {
    gene1<- FetchData(local_cell, vars = i)
    colnames(gene1) <- "gene"
    df <- data.frame(local_cell$sample_id, local_cell$disease_state, gene1$gene)
    colnames(df) <- c("Samples", "Type", "Gene")
    # table(df$local_beta.sample_id)
    
    df<-setDT(df)[ , .(Gene = mean(Gene)), by = Samples]
    
    df <- add_column(df, Type = stri_replace_all_regex(str = df$Samples,
                                                       pattern = c('^C.*','^A.*', '^T.*'), 
                                                       replacement = disease_order,
                                                       vectorize_all = FALSE), .after = "Samples") 
    
    p[[i]] <- ggviolin(df, x = "Type", y = "Gene", fill = "Type",
                       palette = c("#00AFBB", "#E7B800"), # "#FC4E07",
                       add=c("boxplot", "jitter"),add.params = list(fill="white"),
                       order = disease_order,
                       shape = "Type", #size = 0.1,
                       ylab = i, xlab = FALSE) #, xlab = "Groups")
  }
  p1 <- ggpar(p[[1]], ylim = ylimit1)
  p2 <- ggpar(p[[2]], ylim = ylimit2)
  p3 <- ggpar(p[[3]], ylim = ylimit3)

  # arrange the three plots in a single row
  prow <- plot_grid( p1 + theme(legend.position="none"),
                     p2 + theme(legend.position="none"),
                     p3 + theme(legend.position="none"),
                     align = 'vh',
                     # labels = c("A", "B", "C"),
                     hjust = -1,
                     nrow = 1
  )
  legend_b <- get_legend(p1 + theme(legend.position="top"))
  
  # add the legend underneath the row we made earlier. Give it 10% of the height
  # of one plot (via rel_heights).
  p <- plot_grid( prow, legend_b, ncol = 1, rel_heights = c(1, .2))
  return(p)
}
