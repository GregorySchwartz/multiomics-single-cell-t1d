#Color by AnnoSpat labels (Label everything other than islets, Acinar, Ductal and Immune  as 'Unknown')
too-many-cells make-tree \
--prior output_gw/output_hpap_cytof_proteins_12_quantile_norm_7_cut \
--labels-file "/mnt/data2/aanchal/data/cyTOF_T1D/labels_ssc/trte_labels_ssc_ELMTUNE(no-1)_cofactor=200_thlist=1_unknownth=75_only_isletsNDuctalAcinarImmune.csv" \
--draw-colors "[ \"#FF0000\", \"#000000\", \"#FF00FF\", \"#00FFFF\", \"#0000FF\", \"#00FF00\", \"#FFFF00\", \"#dd944a\", \"#eaeaea\"  ]" \
--dendrogram-output "tree_AnnoSpat_labels_suffix=TUNE(no-1)_cofactor=200_thlist=1_unknownth=75.pdf" \
--output output_gw/output_hpap_cytof_proteins_12_quantile_norm_7_cut/AnnoSpatOverlay \
> tmc_annospatoverlay.log +RTS -N30
--draw-node-number \
previous --labels-file /mnt/data2/aanchal/data/cyTOF_T1D/labels_ssc/trte_labels_ssc_ELM_4_fixThs3_cofactor=150_only_isletsNDuctalAcinarImmune.csv \




# #Color by Disease state labels (T1D/AAB+/Control)
too-many-cells make-tree \
--prior output_gw/output_hpap_cytof_proteins_12_quantile_norm_7_cut \
--labels-file /mnt/data2/aanchal/data/cyTOF_T1D/labels_state_cytof.csv \
--draw-colors "[ \"#FF0000\", \"#44d43b\", \"#FFFF00\"]" \
--dendrogram-output "tree_diseaseState_labels.pdf" \
--output output_gw/output_hpap_cytof_proteins_12_quantile_norm_7_cut/diseaseState \
> tmc_stateoverlay.log +RTS -N30
#--draw-node-number \

