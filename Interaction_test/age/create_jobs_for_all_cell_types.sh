####
new_names=(B_IN B_MEM CD4_NC CD4_ET CD4_SOX4 CD8_NC CD8_ET CD8_S100B DC Mono_C Mono_NC NK_R NK Plasma)

# new_names=(B_IN B_MEM)

for i in {1..13}
do

mkdir -p ./${new_names[i]}_cell_5_mean_mx/
mkdir -p ./${new_names[i]}_cell_5_variance_mx/
mkdir -p ./${new_names[i]}_cell_5_dispersion_mx/

rm -r -f ./${new_names[i]}_cell_5_mean_mx/*
rm -r -f ./${new_names[i]}_cell_5_variance_mx/*
rm -r -f ./${new_names[i]}_cell_5_dispersion_mx/*

sed 's/ct_name/'${new_names[i]}'/g' eQTL_OneK1K_cell_5_mean_mx.qsub > eQTL_OneK1K_cell_5_mean_mx_${new_names[i]}.qsub
sed 's/ct_name/'${new_names[i]}'/g' eQTL_OneK1K_cell_5_variance_mx.qsub > eQTL_OneK1K_cell_5_variance_mx_${new_names[i]}.qsub
sed 's/ct_name/'${new_names[i]}'/g' eQTL_OneK1K_cell_5_dispersion_mx.qsub > eQTL_OneK1K_cell_5_dispersion_mx_${new_names[i]}.qsub

qsub eQTL_OneK1K_cell_5_mean_mx_${new_names[i]}.qsub
qsub eQTL_OneK1K_cell_5_variance_mx_${new_names[i]}.qsub
qsub eQTL_OneK1K_cell_5_dispersion_mx_${new_names[i]}.qsub


done

####
