NODE=199
BASE="hpap_imc_5_cut"
DISTANCE=20

too-many-cells make-tree \
               -j spatial_information.csv \
               --prior output_${BASE} \
               -o output_${BASE}_node_${NODE}_proximity_${DISTANCE} \
               --dendrogram-output "node_${NODE}_proximity_${DISTANCE}_tree.pdf" \
               --labels-output \
               --matrix-transpose \
               --draw-collection "PieChart" \
               --filter-thresholds "(1e-16, 1e-16)" \
               --draw-colors "[\"#e41a1c\", \"#cccccc\", \"#377eb8\"]" \
               --draw-leaf "DrawItem (DrawProximity ([${NODE}], ${DISTANCE}))" \
               --draw-scale-saturation 6 \
               > node_${NODE}_${BASE}_proximity_${DISTANCE}.csv
