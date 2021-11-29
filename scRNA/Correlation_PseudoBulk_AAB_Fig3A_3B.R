# Load libraries
library(scater)
library(Seurat)
library(tidyverse)
library(cowplot)
library(Matrix.utils)
library(edgeR)
library(dplyr)
library(magrittr)
library(Matrix)
library(purrr)
library(reshape2)
library(S4Vectors)
library(tibble)
library(SingleCellExperiment)
library(pheatmap)
library(apeglm)
library(png)
library(DESeq2)
library(RColorBrewer)

# Read in the dataset
local <- readRDS("fasolino_et_al.rds")

## Remove samples with no GAD levels
## Without HPAP019 and HPAP043
local <- subset(local, subset = sample_id != "AABE")
local <- subset(local, subset = sample_id != "AABA")
## Extract raw counts and metadata to create SingleCellExperiment object
counts <- local@assays$RNA@counts 
metadata <- local@meta.data

# Set up metadata as desired for aggregation and DE analysis 
# Create single cell experiment object
sce <- SingleCellExperiment(assays = list(counts = counts), 
                            colData = metadata)
sce$sample_id <- as.factor(sce$sample_id)
# Identify groups for aggregation of counts
# Subset metadata to only include the cluster and sample IDs to aggregate across
groups <- colData(sce)[, c("cell_label", "sample_id")]
groups$cell_label <- gsub('\\_', '', groups$cell_label)

# Explore the raw counts for the dataset
## Check the assays present
local@assays
## Explore the raw counts for the dataset
local@assays$RNA@counts[1:6,1:6]
## Explore the cellular metadata for the dataset
dim(colData(sce))
head(colData(sce))

# Named vector of cluster names
k_ids <- purrr::set_names(levels(sce$cell_label))

# Total number of clusters
nk <- length(k_ids)
nk

# Named vector of sample names
s_ids <- purrr::set_names(levels(as.factor(sce$sample_id)))

# Total number of samples 
ns <- length(s_ids)
ns
## Determine the number of cells per sample
table(sce$sample_id)

## Turn named vector into a numeric vector of number of cells per sample
n_cells <- as.numeric(table(sce$sample_id))

## Determine how to reorder the samples (rows) of the metadata to match the order of sample names in s_ids vector
m <- match(s_ids, sce$sample_id)

## Create the sample level metadata by combining the reordered metadata with the number of cells corresponding to each sample.
ei <- data.frame(colData(sce)[m, ], 
                 n_cells, row.names = NULL) %>% 
  select(-"cell_label")
ei

# Perform QC if not already performed
dim(sce)

#######################################################
# Calculate qusceality control (QC) metrics using library(scater) with function calculateQCMetrics
sce1 <- perCellQCMetrics(sce)
## The median is a more robust measure relative to the mean which can be influenced by outliers
sce_detected <- !isOutlier(
  metric = sce1$detected,
  nmads = 2, type = "both", log = TRUE)
table(sce_detected)
# Remove outlier cells
sce <- sce[, sce_detected]
dim(sce)
## Remove lowly expressed genes which have less than 10 cells with any counts
sce <- sce[rowSums(counts(sce) > 1) >= 10, ]
dim(sce)
#######################################################
## Count aggregation to sample level

# Aggregate the counts per sample_id and cluster_id

# Subset metadata to only include the cluster and sample IDs to aggregate across
groups <- colData(sce)[, c("cell_label", "sample_id")]
groups$cell_label <- gsub('\\_', '', groups$cell_label)

# Aggregate across cluster-sample groups
pb <- aggregate.Matrix(t(counts(sce)), 
                       groupings = groups, fun = "sum") 
class(pb)

dim(pb)
pb[1:6, 1:6]
# Not every cluster is present in all samples; create a vector that represents how to split samples
splitf <- sapply(stringr::str_split(rownames(pb), 
                                    pattern = "_",  
                                    n = 2), `[`, 1)

# Turn into a list and split the list into components for each cluster and transform, so rows are genes and columns are samples and make rownames as the sample IDs
pb <- split.data.frame(pb, 
                       factor(splitf)) %>%
  lapply(function(u) 
    set_colnames(t(u), 
                 stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")))

# lapply(1:length(pb), function(i) write.csv(pb[[i]],
#                                                file = paste0(names(pb[i]), "_counts.csv"),
#                                                row.names = TRUE))

######### Extract AAB
new_pb <- lapply(1:length(pb), function(i) pb[[i]][, grepl( "AAB", colnames(pb[[i]]))])
names(new_pb) <- names(pb)

new_pb <- lapply(1:length(new_pb), function(i) new_pb[[i]][which(rowSums(new_pb[[i]]) > 10),])
names(new_pb) <- names(pb)
###### Save raw counts ##############
# lapply(1:length(new_pb), function(i) write.csv(new_pb[[i]],
#                                          file = paste0(names(new_pb[i]), "_AAB_counts.csv"),
#                                          row.names = TRUE))
###################
## Correlation analysis on each individual gene
## Without HPAP019 and HPAP043
GAD_Levels <- c(203, 84, 890, 321, 412, 203)

## AAB7 not found in below two samples therefore removed
new_pb$acinarminormhcclassII <- NULL 
new_pb$ductmajor <- NULL
g_corr_val <- NULL
g_pval <- NULL
for (i in 1:length(new_pb)) 
{
  for(j in 1:nrow(new_pb[[i]]))
  {
    g_tmp <- cor.test(new_pb[[i]][j,], GAD_Levels)
    g_corr_val[j] <- g_tmp$estimate[[1]] ## Fixed  
    g_pval[j] <- g_tmp$p.value
  }  
  new_pb[[i]] <- cbind(new_pb[[i]], g_corr_val, g_pval)
}
for (i in 1:length(new_pb))
{
  new_pb[[i]] <- new_pb[[i]][order(new_pb[[i]][,ncol(new_pb[[i]])],decreasing=FALSE),]
}
for (i in 1:length(new_pb))
{
  # new_pb[[i]] <- filter(as.data.frame(new_pb[[i]]), (g_pval < 0.05))
  new_pb[[i]] <- filter(as.data.frame(new_pb[[i]]), (g_pval < 0.05 & g_corr_val > 0.99))
}
for (i in 1:length(new_pb))
{
  new_pb[[i]] <- cbind(rep(names(new_pb)[i],nrow(new_pb[[i]])), rownames(new_pb[[i]]), new_pb[[i]])
}
############ Write out the list to dataframes ##################
# lapply(1:length(new_pb), function(i) write.csv(new_pb[[i]],
#                                                file = paste0(names(new_pb[i]), "_AAB_correlation.csv"),
#                                                row.names = TRUE))
##########################################################################
######## For bar plot
library(ggpubr)
## create a dataframe by binding all the rows of dataframe from the list of dataframes
dat <- bind_rows(new_pb)
colnames(dat)[1] <- "Cell Type"
colnames(dat)[2] <- "Genes"
colnames(dat)
names(dat)[10] <- "-log10(Pvalue)"

## Bar plots
table(dat$`Cell Type`)
df <- data.frame(CellType=c("acinar", "acinar mhcclassII", "alpha", "beta_major", "beta_minor",
                            "delta", "ductal_acinar", "ductal_major", "endothelial",
                            "epsilon", "immune_stellates", "pp", "stellates"),
                 SignificantGenes=c(0, 0, 3, 1473, 6, 50, 1, 0, 0, 2, 8, 0, 4))

df <- df[order(df$SignificantGenes, decreasing = T),]
# Outside/Inside bars
df <- within(df, 
             CellType <- factor(CellType, 
                                  levels=names(sort(table(CellType), 
                                                        decreasing=TRUE))))
df <- data.frame(df,c(rep("pos", nrow(df))))
colnames(df)[3] <- "variable"
ggplot(df, aes(x = reorder(CellType, -SignificantGenes), y = SignificantGenes, fill = variable)) + 
  geom_bar(stat = "identity", fill="steelblue")+
  geom_text(aes(label=SignificantGenes), vjust=-0.3, size=3.5)+
  ggtitle("") + xlab("Cell Types") + ylab("Significant Genes") +
  theme(axis.title = element_text(face = "bold", color = "black", size = 12)) +
  theme(axis.text.x = element_text(face = "bold", color = "black", size = 8)) +
  theme(axis.text.y = element_text(face = "bold", color = "black", size = 10))

ggsave("Bar_plotV4.pdf", width = 10, height = 5, dpi = 600)

beta_cells <- new_pb$betamajor
dim(beta_cells)[1]
variety <- rep( c("HPAP024", "HPAP029", "HPAP038", "HPAP045", "HPAP049", "HPAP050"), each=dim(beta_cells)[1])

Count <- c(log(beta_cells$AABB), log(beta_cells$AABC), log(beta_cells$AABD), 
           log(beta_cells$AABF), log(beta_cells$AABG), log(beta_cells$AABH))
datasome <- data.frame(variety, Count)

# Create a vector named "new_order" containing the desired order
new_order <- with(datasome, reorder(variety , Count, median , na.rm=T))

# Draw the boxplot using this new order
pdf("Box_plotV6.pdf", width = 10, height = 5)
boxplot(datasome$Count ~ new_order , xlab="HPAP AAB+ Donors", ylab="Avg Count" , col="#69b3a2", boxwex=0.4 , main="")
dev.off()
