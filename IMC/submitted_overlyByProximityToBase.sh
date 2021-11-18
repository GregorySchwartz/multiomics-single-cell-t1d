NODE=257
DISTANCE=20

too-many-cells make-tree \
-j /mnt/data2/aanchal/data/IMC_T1D/raw_data/spatial_information.csv \
--prior output_gw/output_hpap_imc_5_cut \
-o output_gw/output_hpap_imc_5_cut/output_node_${NODE}_proximity_${DISTANCE} \
--dendrogram-output "node_${NODE}_proximity_${DISTANCE}_tree.pdf" \
--labels-output \
--matrix-transpose \
--draw-collection "PieChart" \
--filter-thresholds "(1e-16, 1e-16)" \
--draw-colors "[\"#e41a1c\", \"#cccccc\", \"#377eb8\"]" \
--draw-leaf "DrawItem (DrawProximity ([${NODE}], ${DISTANCE}))" \
--draw-scale-saturation 6 \
> node_${NODE}_proximity_${DISTANCE}.log