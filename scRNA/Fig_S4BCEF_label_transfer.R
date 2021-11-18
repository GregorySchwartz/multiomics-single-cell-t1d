library(Seurat)
library(patchwork)
source("get_mat.R", chdir = TRUE)

# Read a file as matrix market or CSV
read10XOrCsv = function(x) {
  if(dir.exists(x)) {
    return(read10X(x))
  } else {
    return(read.csv(x, row.names = 1, na.strings = "MISSINGVALUE"))
  }
}

# Pipeline or not
toSeurat = function(x, flag, cells) {
  if (flag == 1) {
    sMat = seuratPipeline(x, cells)
    # Remove cells without labels
    sMat = sMat[, colnames(sMat) %in% cells$item]
    return(sMat)
  } else {
    labelDf = cells
    rownames(labelDf) = labelDf$item

    sMat = CreateSeuratObject(counts = x)
    sMat = AddMetaData(object = sMat, metadata = labelDf)
    # Remove cells without labels
    sMat = sMat[,colnames(sMat) %in% cells$item]

    sMat = PercentageFeatureSet(sMat, pattern = "^MT-", col.name = "percent.mt")

    # Transform
    sMat = SCTransform(sMat, vars.to.regress = c("percent.mt"))

    # Reduction
    sMat = RunPCA(sMat)

    # Cluster
    sMat = FindNeighbors(sMat)
    sMat = FindClusters(sMat)

    # Reduction
    sMat = RunUMAP(sMat, reduction = "pca", dims = 1:10)
    return(sMat)
  }
}

# Plot individual objects
plotToFile = function(sMat, output) {
  p = DimPlot(sMat, reduction = "umap", group.by = "label")
  state = sMat$state[1]
  pdf(file.path(output, paste0(as.character(state), "_plot.pdf")), useDingbats = FALSE)
  plot(p)
  dev.off()

  return ()
}

# Format labels
formatLabels = function(xs) {

  return(gsub("-.*$", "", toupper(xs)))

}

# Integrate data
integrateData = function(mats, plotFlag, logNormFlag, filterAnchors, normalizationMethod, output) {
  # Integration
  message("Finding integration anchors")
  integrationFeatures = SelectIntegrationFeatures(object.list = mats)
  
  # For SCT
  if (!(logNormFlag)) {
    mats = PrepSCTIntegration(object.list = mats)
  }
  
  integrationAnchors = FindIntegrationAnchors(object.list = mats, normalization.method = normalizationMethod, anchor.features = integrationFeatures, k.filter = filterAnchors)
  message("Integrating data")
  integrated = IntegrateData(anchorset = integrationAnchors, normalization.method = normalizationMethod)
  # Run the standard workflow for visualization and clustering
  if (logNormFlag) {
    integrated = ScaleData(integrated)
  }
  integrated = RunPCA(integrated)
  integrated = RunUMAP(integrated, reduction = "pca", dims = 1:10)
  if(plotFlag) {
    p1 = DimPlot(integrated, reduction = "umap", group.by = "state")
    p2 = DimPlot(integrated,
      reduction = "umap", group.by = "label", label = TRUE,
      repel = TRUE
    )
    pdf(file.path(output, "integrated_state_plot.pdf"), useDingbats = FALSE)
    plot(p1)
    dev.off()
    pdf(file.path(output, "integrated_label_plot.pdf"), useDingbats = FALSE)
    plot(p2)
    dev.off()
  }
  return(integrated)
}

main = function () {
  args = commandArgs(TRUE)
  labelPath = args[1]
  output = args[2]
  logNormFlag = as.logical(as.numeric(args[3]))
  filterAnchors = as.integer(args[5])
  input = args[5]
  refInputs = args[6:length(args)]

  set.seed(0)
  
  normalizationMethod = if (logNormFlag) {
    "LogNormalize"
  } else {
    "SCT"
  }

  # Make output
  dir.create(file.path(output), recursive = TRUE)

  # Load data
  message("Loading data")
  cells = read.table(file = labelPath, header = TRUE, sep = ",")
  cells$originalItem = cells$item
  cells$item = make.names(cells$item)
  mat = toSeurat(read10XOrCsv(input), logNormFlag, cells)
  refMats = lapply(refInputs, function(x) toSeurat(read10XOrCsv(x), logNormFlag, cells))

  # Plot initial
  sapply(c(mat, unlist(refMats)), function(x) plotToFile(x, output))

  # Integration
  integratedTmp = integrateData(c(list(mat), refMats), TRUE, logNormFlag, filterAnchors, normalizationMethod, output)
  if(length(refMats) > 1) {
    integrated = integrateData(refMats, FALSE, logNormFlag, filterAnchors, normalizationMethod, output)
  } else {
    integrated = refMats[[1]]
  }

  # Query
  message("Finding transfer anchors")
  transferAnchors = FindTransferAnchors(
    reference = integrated, query = mat, k.filter = filterAnchors, normalization.method = normalizationMethod
  )
  message("Predicting labels")
  predictions = TransferData(
    anchorset = transferAnchors, refdata = integrated$label
  )
  message("Adding metadata")
  mat = AddMetaData(mat, metadata = predictions)

  mat@meta.data$predicted.match = mat@meta.data$predicted.id == mat@meta.data$label
  score = table(mat@meta.data$predicted.match)
  write.csv(score, file = file.path(output, "prediction_score.csv"), row.names = FALSE)
  detailedScore = table(mat@meta.data$label, mat@meta.data$predicted.id)
  write.csv(detailedScore, file = file.path(output, "detailed_prediction_score.csv"), row.names = TRUE)

  write.csv(mat@meta.data, file = file.path(output, "predictions.csv"), row.names = TRUE)

}

main()
