too-many-cells differential \
-m '/mnt/data2/aanchal/data/IMC_T1D/raw_data/_hpap_imc_QuantileNorm_mat.csv' \
--features C.peptide --features CA2 --features CD45 --features Glucagon --features Keratin --features pS6 --features CD4 \
 --features CD8 --features CD3 --features CD14 --features CD11b --features CD20 \
--prior output_gw/output_hpap_imc_5_cut \
--normalization "QuantileNorm" \
--plot-output "output_gw/output_hpap_imc_5_cut/diffProteins_plot_7aQNORM.pdf" \
--nodes "([0], [0])" \
--plot-separate-nodes \
--plot-no-outlier \
--filter-thresholds "(1e-16, 1e-16)" \
--labels-file "output_gw/output_hpap_imc_5_cut/output_node_257_proximity_20/labels.csv" \
--plot-separate-labels \
--labels "([\"Neighbor\"], [])" \
--subsample-groups "0" \
> diff_node_${NODE}_proximity_${DISTANCE}_7aQNORM.log +RTS -N30
