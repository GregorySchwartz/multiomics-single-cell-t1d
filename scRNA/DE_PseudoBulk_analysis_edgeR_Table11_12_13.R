## load packages
library(ExperimentHub)
library(dplyr)
library(ggplot2)
library(limma)
library(muscat)
library(purrr)
library(SingleCellExperiment)
library(scater)
library(sctransform)

# Read in the dataset
seurat_obj <- readRDS("data/All_Islet_Doublet_Group.rds")
local <- readRDS("/mnt/data0/abhijeet/pilot/scDataLocal/local.rds")
local@meta.data$sample_id <- seurat_obj@meta.data$group
local@meta.data$sample_id <-gsub("\\_","",local@meta.data$sample_id)
rm(seurat_obj)
###############################################################################################
# Extract raw counts and metadata to create SingleCellExperiment object
cts <- local@assays$RNA@counts 
metadat <- local@meta.data

sce_new <- SingleCellExperiment(assays = list(counts = cts), 
                                colData = metadat)
sce_new$sample_id <- as.factor(sce_new$sample_id)

# remove undetected genes
sce_new <- sce_new[rowSums(counts(sce_new) > 0) > 0, ]
dim(sce_new)

# calculate per-cell quality control (QC) metrics
qc <- perCellQCMetrics(sce_new)

# remove cells with few or many detected genes
ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
sce_new <- sce_new[, !ol]
dim(sce_new)

# remove lowly expressed genes ( < 10 counts)
sce_new <- sce_new[rowSums(counts(sce_new) > 1) >= 10, ]
dim(sce_new)

## Log normcounts to calculate log2-transformed normalized expression
# compute sum-factors & normalize
sce_new <- computeLibraryFactors(sce_new)
sce_new <- logNormCounts(sce_new)
colData(sce_new)

assays(sce_new)$vstresiduals <- vst(counts(sce_new), verbosity = FALSE)$y

# sce_new$id <- paste0(sce_new$stim, sce_new$ind)
(sce_new <- prepSCE(sce_new, 
                    kid = "cell_label", # sub population assignments
                    gid = "disease_state",  # group IDs (ctrl/case)
                    sid = "sample_id",   # sample IDs (ctrlid/caseid)
                    drop = TRUE))  # drop all other colData columns

nk <- length(kids <- levels(sce_new$cluster_id))
ns <- length(sids <- levels(sce_new$sample_id))
names(kids) <- kids; names(sids) <- sids

# nb. of cells per cluster-sample
t(table(sce_new$cluster_id, sce_new$sample_id))

# compute UMAP using 1st 20 PCs
sce_new <- runUMAP(sce_new, pca = 20)
# wrapper to prettify reduced dimension plots
.plot_dr <- function(sce_new, dr, col)
  plotReducedDim(sce_new, dimred = dr, colour_by = col) +
  guides(fill = guide_legend(override.aes = list(alpha = 1, size = 3))) +
  theme_minimal() + theme(aspect.ratio = 1)

# downsample to max. 100 cells per cluster
cs_by_k <- split(colnames(sce_new), sce_new$cluster_id)
cs100 <- unlist(sapply(cs_by_k, function(u) 
  sample(u, min(length(u), 100))))

# plot t-SNE & UMAP colored by cluster & group ID
# "TSNE" removed
for (dr in c("UMAP"))
  for (col in c("cluster_id", "group_id"))
    .plot_dr(sce_new[, cs100], dr, col)
pdf("figures/edgeR/UMAP_all.pdf", width = 5, height = 5)
.plot_dr(sce_new[, cs100], dr, col)
dev.off()

pb1 <- aggregateData(sce_new,
                     assay = "counts", fun = "sum",
                     by = c("cluster_id", "sample_id"))
# one sheet per subpopulation
assayNames(pb1)
# pseudobulks for 1st subpopulation
t(head(assay(pb1)))
# Pseudobulk-level multidimensional scaling (MDS) plot
pdf("figures/edgeR/pb_MDS_all.pdf", width = 5, height = 5.5)
(pb_mds <- pbMDS(pb1))
dev.off()

# use very distinctive shaping of groups & change cluster colors
pb_mds1 <- pb_mds + 
scale_shape_manual(values = c(17, 4)) +
 scale_color_manual(values = RColorBrewer::brewer.pal(12, "Set3"))
# change point size & alpha
pb_mds$layers[[1]]$aes_params$size <- 5
pb_mds$layers[[1]]$aes_params$alpha <- 0.6
pdf("figures/edgeR/pb_MDS_all_v1.pdf", width = 5, height = 5.5)
pb_mds
dev.off()

################################################################################
## Run together (All groups and all cell clusters)

# construct design & contrast matrix
ei <- metadata(sce_new)$experiment_info
mm <- model.matrix(~ 0 + ei$group_id)
dimnames(mm) <- list(ei$sample_id, levels(ei$group_id))

################################ T1D-Control ###########################################
contrast <- makeContrasts("T1D-Control", levels = mm)
# run DS analysis
res_new <- pbDS(pb1, method = c("edgeR"), design = mm, contrast = contrast)

# access results table for 1st comparison
tbl <- res_new$table[[1]]
# one data.frame per cluster
names(tbl)

lapply(1:length(tbl), function(i) write.csv(tbl[[i]], 
                                            file = paste0(names(tbl[i]), "_T1D-Control.csv"),
                                            row.names = FALSE))
tbl_fil <- lapply(tbl, function(u) {
  u <- dplyr::filter(u, p_val < 0.05, abs(logFC) > 1)
  dplyr::arrange(u, p_val)
})
## Filtering
lapply(1:length(tbl_fil), function(i) write.csv(tbl_fil[[i]], 
                                            file = paste0(names(tbl_fil[i]), "_T1D-Control_filtered.csv"),
                                            row.names = FALSE))
################################ T1D-AAB ###########################################
contrast <- makeContrasts("T1D-AAB", levels = mm)
# run DS analysis
res_new <- pbDS(pb1, method = c("edgeR"), design = mm, contrast = contrast)

# access results table for 1st comparison
tbl <- res_new$table[[1]]
# one data.frame per cluster
names(tbl)

lapply(1:length(tbl), function(i) write.csv(tbl[[i]], 
                                            file = paste0(names(tbl[i]), "_T1D-AAB.csv"),
                                            row.names = FALSE))
## Filtering
tbl_fil <- lapply(tbl, function(u) {
  u <- dplyr::filter(u, p_val < 0.05, abs(logFC) > 1)
  dplyr::arrange(u, p_val)
})
lapply(1:length(tbl_fil), function(i) write.csv(tbl_fil[[i]], 
                                                file = paste0(names(tbl_fil[i]), "_T1D-AAB_filtered.csv"),
                                                row.names = FALSE))
################################ AAB-Control ###########################################
contrast <- makeContrasts("AAB-Control", levels = mm)
# run DS analysis
res_new <- pbDS(pb1, method = c("edgeR"), design = mm, contrast = contrast)

# access results table for 1st comparison
tbl <- res_new$table[[1]]
# one data.frame per cluster
names(tbl)

lapply(1:length(tbl), function(i) write.csv(tbl[[i]], 
                                            file = paste0(names(tbl[i]), "_AAB-Control.csv"),
                                            row.names = FALSE))
## Filtering
tbl_fil <- lapply(tbl, function(u) {
  u <- dplyr::filter(u, p_val < 0.05, abs(logFC) > 1)
  dplyr::arrange(u, p_val)
})
lapply(1:length(tbl_fil), function(i) write.csv(tbl_fil[[i]], 
                                                file = paste0(names(tbl_fil[i]), "_AAB-Control_filtered.csv"),
                                                row.names = FALSE))
#######################################################################################
