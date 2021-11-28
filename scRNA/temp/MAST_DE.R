setwd("/mnt/alvand/abhijeet/aab/november/maria")

library(MAST)
# Read in the dataset
local <- readRDS("fasolino_et_al.rds")

################## Uncomment any one of the below three ########################
## Create datasets based on different disease states for comparison of cases/control studies 

## 1. T1DvsControl
arg_obj<- subset(local, subset = disease_state != "AAB")

## 2. T1DvsAAB
# arg_obj<- subset(local, subset = disease_state != "Control")

## 3. AABvsControl
# arg_obj <- subset(local, subset = disease_state != "T1D")

## Proceed on any of the above pairwise comparison
## Dropping unused levels from factors
arg_obj$RNA_snn_res.1.2 <- droplevels(arg_obj$RNA_snn_res.1.2)
arg_obj$seurat_clusters <- droplevels(arg_obj$seurat_clusters)
arg_obj$disease <- droplevels(arg_obj$disease)
arg_obj$disease_ontology_term_id <- droplevels(arg_obj$disease_ontology_term_id)
arg_obj$disease_state <- droplevels(arg_obj$disease_state)
arg_obj$cell_label <- droplevels(arg_obj$cell_label)
arg_obj$assay<- droplevels(arg_obj$assay)
arg_obj$ethnicity_ontology_term_id <- droplevels(arg_obj$ethnicity_ontology_term_id)
arg_obj$ethnicity <- droplevels(arg_obj$ethnicity)
arg_obj$sex <- droplevels(arg_obj$sex)
arg_obj$tissue <- droplevels(arg_obj$tissue)
arg_obj$tissue_ontology_term_id <- droplevels(arg_obj$tissue_ontology_term_id)
arg_obj$assay_ontology_term_id <- droplevels(arg_obj$assay_ontology_term_id)
arg_obj$cell_type_ontology_term_id <- droplevels(arg_obj$cell_type_ontology_term_id)
arg_obj$cell_type <- droplevels(arg_obj$cell_type)
arg_obj$development_stage_ontology_term_id <- droplevels(arg_obj$development_stage_ontology_term_id)
arg_obj$development_stage <- droplevels(arg_obj$development_stage)
rm(local)

## Relevel Disease states in order to do comparison
## In T1D vs Control, the argument in relevel needs to be passed as Control
## In T1D vs AAB, pass AAB as argument
## In AAB vs Control, pass Control as argument

arg_obj$disease_state <-relevel(arg_obj$disease_state,"Control") # In order to have more interpretable coefficients, I set the reference level of my genotype condition to be my wild type cells.


## Create a object of T1D vs Control for cell type beta_major
## In case want to perform analysis on only one cell type, lets say acinar, uncomment below two lines
# cts_obj <- subset(arg_obj, subset = cell_label == "acinar")
# arg_obj$disease_state <- droplevels(arg_obj$disease_state)

## Counts- build a matrix
cts <- arg_obj@assays$RNA@counts 
arg_obj@meta.data$well_key <- colnames(arg_obj)

metadat <- data.frame(arg_obj@meta.data$well_key,
                      arg_obj@meta.data$sample_id,
                      arg_obj@meta.data$disease_state) 
names(metadat) <- c("wellKey", "DonorID", "Status")
rownames(metadat) <- metadat[,1]

## Convert counts to log2(counts+1)
log2cts <- log2(cts+1)
log2cts <- Matrix::as.matrix(log2cts)
coldata <- metadat
fData <- data.frame(primerid=rownames(cts))

sca <- MAST::FromMatrix(exprsArray=log2cts, cData=coldata, fData=fData)

cdr2 <- colSums(SummarizedExperiment::assay(sca)>0)
SummarizedExperiment::colData(sca)$ngeneson <- scale(cdr2)
SummarizedExperiment::colData(sca)$Status <-
  factor(SummarizedExperiment::colData(sca)$Status)
SummarizedExperiment::colData(sca)$DonorID <-
  factor(SummarizedExperiment::colData(sca)$DonorID)

## Fit the model
zlmCond <- MAST::zlm(~ ngeneson + Status + (1|DonorID), # Status:DonorID
                                      sca, method='glmer',ebayes = F,
                                      strictConvergence = FALSE)

# summaryCond <- MAST::summary(zlmCond, doLRT='StatusControl')
# 
# summaryDt <- summaryCond$datatable
# fcHurdle <- merge(summaryDt[summaryDt$contrast=='StatusControl'
#                             & summaryDt$component=='logFC', c(1,7,5,6,8)],
#                   summaryDt[summaryDt$contrast=='StatusControl'
#                             & summaryDt$component=='H', c(1,4)],
#                   by = 'primerid')
# fcHurdle <- stats::na.omit(as.data.frame(fcHurdle))



