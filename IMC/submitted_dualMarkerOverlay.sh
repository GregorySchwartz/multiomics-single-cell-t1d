too-many-cells make-tree \
 --matrix-path /mnt/data2/aanchal/data/IMC_T1D/raw_data/imc_fromGW.csv \
--matrix-transpose --normalization "QuantileNorm" \
--prior output_gw/output_hpap_imc_5_cut \
--draw-leaf "DrawItem (DrawThresholdContinuous [(\"HLA.DR\", Exact 2), (\"Keratin\", Exact 10)])" \
--dendrogram-output "tree_HLADR-Keratin_2-10.pdf" \
--draw-scale-saturation 10 \
--draw-colors "[\"#e41a1c\", \"#377eb8\", \"#4daf4a\", \"#eaeaea\"]" \
--output output_gw/output_hpap_imc_5_cut/dualMarkerOverlay/HLADR-Keratin \
> logs/tmc_HLADR-Keratin.log +RTS -N50