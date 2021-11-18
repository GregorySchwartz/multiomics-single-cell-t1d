library(Matrix)
library(Seurat)

readCompressed = function(isCompressed, file) {
  if(isCompressed) {
    return(gzfile(file))
  } else {
    return(file)
  }
}

# Read a 10x matrix in matrix market format from a folder.
read10X = function(path) {
  if (file.exists(paste0(path, "/matrix.mtx.gz"))) {
    matFile = paste0(path, "/matrix.mtx.gz")
    featuresFile = paste0(path, "/features.tsv.gz")
    barcodesFile = paste0(path, "/barcodes.tsv.gz")
    compressed = TRUE
  } else {
    matFile = paste0(path, "/matrix.mtx")
    featuresFile = paste0(path, "/genes.tsv")
    barcodesFile = paste0(path, "/barcodes.tsv")
    compressed = FALSE
  }


  print("Loading files:")
  print(matFile)
  print(featuresFile)
  print(barcodesFile)

  # read in matrix data using the Matrix package
  write("Reading matrix", stderr())
  mat = readMM(readCompressed(compressed, matFile))

  # format cell info
  write("Getting cell info", stderr())
  barcodes = read.csv(readCompressed(compressed, barcodesFile), header = FALSE, stringsAsFactors = FALSE)[,1]

  # format peak info
  write("Formatting features info", stderr())
  features = read.csv(readCompressed(compressed, featuresFile), sep = "\t", header = FALSE, stringsAsFactors = FALSE)[,1]

  write("Assigning index names", stderr())
  row.names(mat) = features
  colnames(mat) = make.names(barcodes, unique = TRUE)

  write("Returning matrix", stderr())
  return(mat)
}

# Read multiple 10x folders.
read10Xs = function(paths) {
  ls = sapply(paths, read10X)
  mat = Reduce(function(acc, x) RowMergeSparseMatrices(acc, x), ls)
  return(mat)
}

# New seurat loading from https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
seuratPipeline = function(mat, cells) {
  # object setup
  mat = CreateSeuratObject(counts = mat)

  labelDf = cells
  rownames(labelDf) = labelDf$item
  mat = AddMetaData(object = mat, metadata = labelDf)

  # MT
  mat = PercentageFeatureSet(mat, pattern = "^MT-", col.name = "percent.mt")

  # preprocessing
  mat = NormalizeData(mat)
  mat = FindVariableFeatures(mat)
  mat = ScaleData(mat, vars.to.regress = c("percent.mt"))

  # dimensionality reduction
  mat = RunPCA(mat)

  # Cluster
  mat = FindNeighbors(mat)
  mat = FindClusters(mat)

  # dimensionality reduction again
  mat = RunUMAP(mat, dims = 1:10)

  return(mat)

}
