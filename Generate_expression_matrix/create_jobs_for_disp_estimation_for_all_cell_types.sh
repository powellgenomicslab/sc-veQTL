new_names=(B_IN B_MEM CD4_NC CD4_ET CD4_SOX4 CD8_NC CD8_ET CD8_S100B DC Mono_C Mono_NC NK_R NK Plasma)

# mle_per_gene_within_ind.R

for i in {0..13}
do

echo ${new_names[i]}

rm -r -f ./${new_names[i]}/

mkdir -p ${new_names[i]}
# mv  ./${new_names[i]} ./old_sct_counts_untransformed

sed 's/ct_name/'${new_names[i]}'/g' mle_per_gene_within_ind.qsub.sh > mle_per_gene_within_ind_${new_names[i]}.qsub.sh

qsub mle_per_gene_within_ind_${new_names[i]}.qsub.sh

done
