library(Libra)
library(SingleCellExperiment)
local <- readRDS("fasolino_et_al.rds")

################## Uncomment any one of the below three ########################
## Create datasets based on different disease states for comparison of cases/control studies 

## 1. T1DvsControl
arg_obj<- subset(local, subset = disease_state != "AAB")

## 2. T1DvsAAB
# arg_obj<- subset(local, subset = disease_state != "Control")

## 3. AABvsControl
# arg_obj <- subset(local, subset = disease_state != "T1D")

## Dropping unused levels from factors
arg_obj$RNA_snn_res.1.2 <- droplevels(arg_obj$RNA_snn_res.1.2)
arg_obj$seurat_clusters <- droplevels(arg_obj$seurat_clusters)
arg_obj$disease <- droplevels(arg_obj$disease)
arg_obj$disease_ontology_term_id <- droplevels(arg_obj$disease_ontology_term_id)
arg_obj$disease_state <- droplevels(arg_obj$disease_state)

## In order to have more interpretable coefficients, I set the reference level of my genotype condition to be my wild type cells.
## Relevel Disease states in order to do comparison
## In T1D vs Control, the argument in relevel needs to be passed as Control
## In T1D vs AAB, pass AAB as argument
## In AAB vs Control, pass Control as argument
arg_obj$disease_state <-relevel(arg_obj$disease_state,"Control")

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

## Optionally can remove the below information from metadata as it is not needed
## Uncomment the below lines to remove irrelevant metadata info

# arg_obj$development_stage <- arg_obj$development_stage_ontology_term_id <- NULL
# arg_obj$cell_type_ontology_term_id <- arg_obj$assay_ontology_term_id <- arg_obj$tissue_ontology_term_id <- NULL
# arg_obj$tissue <- arg_obj$sex <- NULL
# arg_obj$RNA_snn_res.1.2 <- NULL
# arg_obj$seurat_clusters <- NULL
# arg_obj$assay <- arg_obj$cell_type <- NULL
# arg_obj$ethnicity <- arg_obj$disease <- NULL
# arg_obj$ethnicity_ontology_term_id <- arg_obj$disease_ontology_term_id <- NULL

arg_obj$disease_state <- factor(arg_obj$disease_state)
arg_obj$cell_label <- as.character(arg_obj$cell_label)
arg_obj$sample_id <- as.character(arg_obj$sample_id)

# Extract raw counts and metadata to create SingleCellExperiment object
cts <- arg_obj@assays$RNA@counts 
metadat <- arg_obj@meta.data

DE = run_de(arg_obj, meta = metadat,
            replicate_col = "sample_id",
            cell_type_col = "cell_label",
            label_col = "disease_state",
            min_cells = 3,
            min_reps = 2,
            min_features = 0,
            de_family = "mixedmodel",
            de_method = "negbinom", #"negbinom",
            de_type = "LRT",
            n_threads = 35
)


## Other approaches tried
## Tried creating a single cell experiment and passed as an argument to run_de function but didn't work
## Uncomment the block
# sce_new <- SingleCellExperiment(assays = list(counts = cts), 
#                                 colData = metadat)
# DE = run_de(sce_new, meta = metadat,
#             replicate_col = "sample_id",
#             cell_type_col = "cell_label",
#             label_col = "disease_state",
#             min_cells = 3,
#             min_reps = 2,
#             min_features = 0,
#             de_family = "mixedmodel",
#             de_method = "negbinom", # "negbinom",
#             # de_type = "LRT",
#             n_threads = 35
#             )

## Tried passing the cts matrix directly but didn't work
# cts_mat <- as.matrix(cts)
# DE = run_de(cts_mat, meta = metadat,
#             replicate_col = "sample_id",
#             cell_type_col = "cell_label",
#             label_col = "disease_state",
#             min_cells = 3,
#             min_reps = 2,
#             min_features = 0,
#             de_family = "mixedmodel",
#             de_method = "negbinom", #"negbinom",
#             de_type = "LRT",
#             n_threads = 35
# )
