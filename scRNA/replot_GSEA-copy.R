## Script by Thomas Kuilman
## path argument: path to output folder of analysis (e.g. PATH/my_analysis.GseaPreranked.1470948568349)
## gene.set argument: name of the gene set (e.g. V$AP1_Q2).
## It is used in a grep command, so multiple matching is possible.
## Also, R regular expressions can be handled, e.g. "IL2[0-9]$"
## Leading "V$" from gene set names are stripped to allow using the grep command.
## In case of multiple grep matches a warning is given and the first option is plotted.
## class.name: the name of the class / variable to which genes have been correlated (e.g. drug-treatment)
## min.score/max.score: when plotting several plots together, set the lower/upper limit

#### modified by Yeqiao 11/21/16
#### This code is incompatible with any command of openning new device... plot in Rstudio and save. 

#### plot a gradient of colors 
showpanel = function(col)
{
  image(z=matrix(1:100, ncol=1), col=col, xaxt="n", yaxt="n" )
}

#### main function
replotGSEA <- function(path, gene.set, class.name, min.score, max.score) {
  if(missing(path)) {
    stop("Path argument is required")
  }
  if (!file.exists(path)) {
    stop("The path folder could not be found. Please change the path")
  }
  if(missing(gene.set)) {
    stop("Gene set argument is required")
  }
  
  ## Load .rnk data
  path.rnk <- list.files(path = file.path(path, "edb"),
                         pattern = ".rnk$", full.names = TRUE)
  gsea.rnk <- read.delim(file = path.rnk, header = FALSE)
  colnames(gsea.rnk) <- c("hgnc.symbol", "metric")
  
  ## Load .edb data
  path.edb <- list.files(path = file.path(path, "edb"),
                         pattern = ".edb$", full.names = TRUE)
  gsea.edb <- read.delim(file = path.edb,
                         header = FALSE, stringsAsFactors = FALSE)
  gsea.edb <- unlist(gsea.edb)
  gsea.metric <- gsea.edb[grep("METRIC=", gsea.edb)]
  gsea.metric <- unlist(strsplit(gsea.metric, " "))
  gsea.metric <- gsea.metric[grep("METRIC=", gsea.metric)]
  gsea.metric <- gsub("METRIC=", "", gsea.metric)
  gsea.edb <- gsea.edb[grep("<DTG", gsea.edb)]
  
  # Select the right gene set
  if (length(gsea.edb) == 0) {
    stop(paste("The gene set name was not found, please provide",
               "a correct name"))
  }
  if (length(grep(paste0(gsub(".\\$(.*$)", "\\1", gene.set), " "), gsea.edb)) > 1) {
    warning(paste("More than 1 gene set matched the gene.set",
                  "argument; the first match is plotted"))
  }
  gsea.edb <- gsea.edb[grep(paste0(gsub(".\\$(.*$)", "\\1", gene.set), " "), gsea.edb)[1]]
  
  # Get template name
  gsea.edb <- gsub(".*TEMPLATE=(.*)", "\\1", gsea.edb)
  gsea.edb <- unlist(strsplit(gsea.edb, " "))
  gsea.template <- gsea.edb[1]
  
  # Get gene set name
  gsea.gene.set <- gsea.edb[2]
  gsea.gene.set <- gsub("GENESET=gene_sets.gmt#", "", gsea.gene.set)
  
  # Get enrichment score
  gsea.enrichment.score <- gsea.edb[3]
  gsea.enrichment.score <- gsub("ES=", "", gsea.enrichment.score)
  
  # Get gene set name
  gsea.normalized.enrichment.score <- gsea.edb[4]
  gsea.normalized.enrichment.score <- gsub("NES=", "",
                                           gsea.normalized.enrichment.score)
  
  # Get nominal p-value
  gsea.p.value <- gsea.edb[5]
  gsea.p.value <- gsub("NP=", "", gsea.p.value)
  gsea.p.value <- as.numeric(gsea.p.value)
  
  # Get FDR
  gsea.fdr <- gsea.edb[6]
  gsea.fdr <- gsub("FDR=", "", gsea.fdr)
  gsea.fdr <- as.numeric(gsea.fdr)
  
  # Get hit indices
  gsea.edb <- gsea.edb[grep("HIT_INDICES=", gsea.edb):length(gsea.edb)]
  gsea.hit.indices <- gsea.edb[seq_len(grep("ES_PROFILE=", gsea.edb) - 1)]
  gsea.hit.indices <- gsub("HIT_INDICES=", "", gsea.hit.indices)
  gsea.hit.indices <- as.integer(gsea.hit.indices)
  
  # Get ES profile
  gsea.edb <- gsea.edb[grep("ES_PROFILE=", gsea.edb):length(gsea.edb)]
  gsea.es.profile <- gsea.edb[seq_len(grep("RANK_AT_ES=", gsea.edb) - 1)]
  gsea.es.profile <- gsub("ES_PROFILE=", "", gsea.es.profile)
  gsea.es.profile <- as.numeric(gsea.es.profile)
  
  
  ## Create GSEA plot
  # Save default for resetting

  def.par <- par(no.readonly = TRUE)
  
  # Create a new device of appropriate size
  dev.new(width = 3, height = 3)
  
  # Create a division of the device
  gsea.layout <- layout(matrix(c(1, 2, 3, 4)), heights = c(1.7, 0.5, 0.2, 2))
  layout.show(gsea.layout)
  
  # Create plots
  par(mar = c(0, 5, 2, 2))
  plot(c(1, gsea.hit.indices, length(gsea.rnk$metric)),
       c(0, gsea.es.profile, 0), type = "l", col = "darkgreen", lwd = 2.5, xaxt = "n",
       xaxs = "i", xlab = "", ylab = "Enrichment score (ES)",
       ylim = c(min.score, max.score),
       main = list(gsea.gene.set, font = 1, cex = 1),
       panel.first = {
         abline(h = seq(round(min(gsea.es.profile), digits = 1),
                        max(gsea.es.profile), 0.1),
                col = "gray95", lty = 2)
         abline(h = 0, col = "gray50", lty = 2)
       })
  plot.coordinates <- par("usr")
  if(gsea.enrichment.score < 0) {
    text(length(gsea.rnk$metric) * 0.01, plot.coordinates[3] * 0.98,
         paste("Nominal p-value:", gsea.p.value, "\nFDR:", gsea.fdr, "\nES:",
               gsea.enrichment.score, "\nNormalized ES:",
               gsea.normalized.enrichment.score), adj = c(0, 0))
  } else {
    text(length(gsea.rnk$metric) * 0.99, plot.coordinates[4] - ((plot.coordinates[4] - plot.coordinates[3]) * 0.03),
         paste("Nominal p-value:", gsea.p.value, "\nFDR:", gsea.fdr, "\nES:",
               gsea.enrichment.score, "\nNormalized ES:",
               gsea.normalized.enrichment.score, "\n"), adj = c(1, 1))
  }
  
  par(mar = c(0, 5, 0, 2))
  plot(0, type = "n", xaxt = "n", xaxs = "i", xlab = "", yaxt = "n",
       ylab = "", xlim = c(1, length(gsea.rnk$metric)))
  abline(v = gsea.hit.indices, lwd = 0.75)
  
  par(mar = c(0, 5, 0, 2))
  
  #### separate the colors for positive and negative fold change
  ncolor.1 = 1000 #sum(gsea.rnk$metric > 0)
  ncolor.2 = 1000 #sum(gsea.rnk$metric <= 0)
  
  #### most positive -> 0
  r.color.1 = colorRampPalette(c("red", "white"))(ncolor.1)
  #### 0 -> most negative 
  r.color.2 = colorRampPalette(c("white", "blue"))(ncolor.2)
  #### put together. the 
  rank.colors = c(r.color.1, r.color.2) 
  
  #### plot color bar
  showpanel(rank.colors)

  box()
  text(length(gsea.rnk$metric) / 2, 0.7,
       labels = ifelse(!missing(class.name), class.name, gsea.template))
  text(length(gsea.rnk$metric) * 0.01, 0.7, "Positive", adj = c(0, NA))
  text(length(gsea.rnk$metric) * 0.99, 0.7, "Negative", adj = c(1, NA))
  
  par(mar = c(5, 5, 0, 2))
  rank.metric <- rle(round(gsea.rnk$metric, digits = 2))
 # plot(gsea.rnk$metric, type = "n", xaxs = "i",
#       xlab = "Rank in ordered gene list", xlim = c(0, length(gsea.rnk$metric)),
#       ylim = c(-1, 1), yaxs = "i",
#       ylab = if(gsea.metric == "None") {"Ranking metric"} else {gsea.metric},
#       panel.first = abline(h = seq(-0.5, 0.5, 0.5), col = "gray95", lty = 2))
  
#  barplot(rank.metric$values, col = "lightgrey", lwd = 0.1, xaxs = "i",
#          xlab = "Rank in ordered gene list", xlim = c(0, length(gsea.rnk$metric)),
#          ylim = c(-1, 1), yaxs = "i", width = rank.metric$lengths, border = NA,
#          ylab = ifelse(gsea.metric == "None", "Ranking metric", gsea.metric), space = 0, add = TRUE)
#  box()
  
  # Reset to default
 # dev.off()
# par(def.par)
}

# Golnaz Sep 6 2021, Fasolino et al revision Nature Metabolism

replotGSEA('Ductal_Acinar_T1D_vs_Control_ImmuneSig.GseaPreranked.1628721199190','DC1.grp', 'na_pos', 0, 0.55)
quartz.save('Ductal_Acinar_T1D_vs_Control_ImmuneSig_DC1.pdf', type = "pdf",dpi=10)#, device = dev.cur())

# Golnaz Sep 6 2021, Fasolino et al revision Nature Metabolism

replotGSEA('Beta_T1D_vs_Control_ImmuneSig.GseaPreranked.1628721185505','DC1.grp', 'na_pos', -0.5, 0.2)
quartz.save('Beta_T1D_vs_Control_ImmuneSig_DC1.pdf', type = "pdf",dpi=10)#, device = dev.cur())
