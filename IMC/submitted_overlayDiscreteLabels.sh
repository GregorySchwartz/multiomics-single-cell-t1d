
# # #Color by AnnoSpat labels (immune subtypes merged and Label everything other than islets, Acinar, Ductal and Immune as 'Unknown')
too-many-cells make-tree \
--prior output_gw/output_hpap_imc_5_cut \
--labels-file /mnt/data2/aanchal/data/IMC_T1D/labels_ssc/trte_labels_ssc_ELM_withImmuneCelltypes_withNegMarkers_d-adaptive2_HPAP026Control_9_only_isletsNDuctalAcinar.csv \
--draw-colors "[ \"#FF0000\", \"#000000\", \"#FF00FF\", \"#00FFFF\", \"#0000FF\", \"#00FF00\", \"#dd944a\", \"#eaeaea\"  ]" \
--dendrogram-output "tree_AnnoSpat_labels_HPAP026Control_9.pdf" \
--output output_gw/output_hpap_imc_5_cut/AnnoSpat \
> logs/tmc_annospatoverlay_HPAP026Control.log +RTS -N30
#--draw-node-number \
#--labels-file /mnt/data2/aanchal/data/IMC_T1D/labels_ssc/trte_labels_ssc_ELM_withImmuneCelltypes_withNegMarkers_d-adaptive2_HPAP026Control_9_only_isletsNDuctalAcinarImmune.csv \

# color by Disease State
too-many-cells make-tree \
--prior output_gw/output_hpap_imc_5_cut \
--labels-file /mnt/data2/aanchal/data/IMC_T1D/raw_data/item_Status_HPAP026shifted2control.csv \
--draw-colors "[ \"#FF0000\", \"#44d43b\", \"#FFFF00\"]" \
--dendrogram-output "tree_diseaseState_labels_HPAP026shifted2control.pdf" \
--output output_gw/output_hpap_imc_5_cut/diseaseState \
> logs/tmc_stateoverlay.log +RTS -N10
#--draw-node-number \
#/mnt/data2/aanchal/data/IMC_T1D/raw_data/labels_state_imc.csv \ --prev labels without shifting hpap026 to control


