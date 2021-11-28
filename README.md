# multiomics-single-cell-t1d
Supplementary code for paper

## Running code
Each script is named after an associated figure in the paper. While each script may be run with the appropriate interpreter (`bash` for `.sh` files, `Rscript` for `.R` files, `python` for `.py` files, `stack` for `.hs` files), input files and folders must be changed to represent the location of the raw data on the user's local system. Many of these scripts use `too-many-cells`, which can be download at https://github.com/GregorySchwartz/too-many-cells with instructions and documentation for additional features.

### scRNA-seq clustering, doublet removal, cell type classification, and doublet identification using DoubletFinder
Seurat v3.1.5 was used for filtering, UMAP generation, and initial clustering. Please refer to the FigS1-S3_Seurat_and_DoubletFinder.Rmd for relevant code to be run in R. 

### Doublet identification using Scrublet 
In addition to DoubletFinder (described above), Scrublet v0.2.1 was also used to identify doublets. Please refer to FigS2A_Scrublet.py script for relevant code to be run using Python.

### Cell Classification 
Garnett was used for cell classification. Please refer to FigS3_Fig1_Garnett.Rmd for relevant code to be run in R. The cell type marker file used in conjunction with the code is found in supplemental materials associated with the paper.

### Label transfer
The label transfer script is a small command line program to use Seurat's label transfer on single-cell data. Usage is:

`Rscript Fig_S4BCEF_label_transfer.R LABELPATH OUTPUTPATH LOGNORMFLAG FILTERANCHORS INPUT REFINPUTS`

Here, `LABELPATH` are the labels for each barcode in a file (see `too-many-cells` documentation), `OUTPUT` is the output label, `LOGNORMFLAG` is whether to use Seurat's `LogNormalize` instead of `SCT`, `FILTERANCHORS` is the number of anchors to use, the `INPUT` is the input scRNA-seq data in matrix market or CSV format, and `REFINPUTS` are a list of reference scRNA-seq data matrices. `get_mat.R` in this repository collects helper functions for this script.

### Creating IMC trees
`stack Fig_5AF_all_IMC_CyTOF_trees_analyze_norm.hs`, (installing `stack`: https://docs.haskellstack.org/en/stable/README/) will generate trees used in the paper. On line `60`, replacing any part of the tuple will cahnge the arguments to `too-many-cells` called. In order, the tuple specifies (see `too-many-cells` documentation): (tree cutting number, a name for the analysis, a path for the input scRNA-seq matrices, a whitelist of included cells, the normalization used to process the data, how many dimensions to drop with LSA, whether to specify a specific node to set as a new root, the labels file containing labels per cell barcode).

### Differential Gene Expression
Differential genes were found using edgeR through TooManyCells with the `-–normalization "NoneNorm"` to invoke edgeR single-cell preprocessing, including normalization and filtering. 

#### Code: 
`too-many-cells differential --prior "" --matrix-path "TotalMed Normalized Matrix" --labels-file "" --normalization "NoneNorm" --nodes "(,)" --labels "([\"\"], [\"\"])" -t 1000000 +RTS -N35`

### Pseudobulk Differential Gene Expression Analysis
Differential genes were found using edgeR through muscat. Please refer to DE_PseudoBulk_analysis_edgeR_Table11_12_13.R for relevant code to be run in R.   

### Correlation analysis for measuring GAD levels
Pseudobulk of samples were aggregated using muscat. Correlation analysis of each individual gene with GAD levels were calculated. A p-value of less than 0.05 and correlation greater than 0.99 was used to filter. Please refer to Correlation_PseudoBulk_AAB_Fig3A_3B.R for the relevant code to be run in R.

### Assessment of common genetic variants associated with T1D using CELLEX and CELLECT
The different expression specificity metrics were calculated using CELLEX tool. The output of CELLEX and GWAS traits were given as input to CELLECT which ranked association of cell types and various disease studies based on p-value. Please refer to cellex_prior_plot_Fig_S8D.R to generate the plot in R.


