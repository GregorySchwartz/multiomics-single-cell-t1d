Matrix_path="/mnt/data2/aanchal/data/cyTOF_T1D/formatted_only_proteins_12"
 
too-many-cells make-tree \
-m /mnt/data2/aanchal/data/cyTOF_T1D/merged_file.csv \
--matrix-transpose --normalization "QuantileNorm" \
--prior output_gw/output_hpap_cytof_proteins_12_quantile_norm_7_cut \
--draw-leaf "DrawItem (DrawThresholdContinuous [(\"HLA-DR\", Exact 5), (\"pan-Cytokeratin\", Exact 200)])" \
--dendrogram-output "tree_HLADR-Keratin_5-200.pdf" \
--draw-scale-saturation 10 \
--draw-colors "[\"#e41a1c\", \"#377eb8\", \"#4daf4a\", \"#eaeaea\"]" \
--output output_gw/output_hpap_cytof_proteins_12_quantile_norm_7_cut/dualMarkerOverlay_HLADR-Keratin \
> tmc_HLADR-Keratin.log +RTS -N50
