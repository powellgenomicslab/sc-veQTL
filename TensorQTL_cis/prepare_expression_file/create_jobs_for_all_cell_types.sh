## Create qsub files for all 14 cell types

new_names=(B_IN B_MEM CD4_NC CD4_ET CD4_SOX4 CD8_NC CD8_ET CD8_S100B DC Mono_C Mono_NC NK_R NK Plasma)
trait=(cell_5_mean_mx cell_5_variance_mx cell_5_residual_mx cell_5_dispersion_mx)

for k in {0..3}
do

echo ${trait[k]}
rm -r -f ./${trait[k]}

mkdir -p ./${trait[k]}

sed 's/generate_expression_file/generate_expression_file_'${trait[k]}'/g' generate_expression_file.qsub.sh > ./${trait[k]}/generate_expression_file.qsub.sh

cd ./${trait[k]}/

for i in {0..13}
do

echo ${new_names[i]} 

# rm -r -f ./${new_names[i]}
# mv  ./${new_names[i]} ./old_sct_counts_untransformed

sed 's/ct_name/'${new_names[i]}'/g' generate_expression_file.qsub.sh > generate_expression_file_${new_names[i]}.qsub.sh

qsub generate_expression_file_${new_names[i]}.qsub.sh

done

cd ../

done
##
