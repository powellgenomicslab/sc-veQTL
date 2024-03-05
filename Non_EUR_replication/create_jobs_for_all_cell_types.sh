## Create.qsub.sh files for all 14 cell types

new_names=(B cDC cM ncM NK PB pDC Progen Prolif T4 T8)

for i in {0..10}
do

echo ${new_names[i]} 

# rm -r -f ./${new_names[i]}
# mv  ./${new_names[i]} ./old_sct_counts_untransformed

sed 's/ct_name/'${new_names[i]}'/g' generate_expression_file_cell_5_dispersion_mx.qsub.sh > generate_expression_file_cell_5_dispersion_mx_${new_names[i]}.qsub.sh

qsub generate_expression_file_${new_names[i]}.qsub.sh

done


##
