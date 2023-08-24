## Submit qsub scripts for TensorQTL analysis

new_names=(B_IN B_MEM CD4_NC CD4_ET CD4_SOX4 CD8_NC CD8_ET CD8_S100B DC Mono_C Mono_NC NK_R NK Plasma)

estimate=(mean,variance,residual,dispersion)

for k in {0..13}
  do

  for j in {0..3}
  do
  
  # Create directory
  mkdir -p ./${new_names[k]}_cell_5_${estimate[j]}_mx/
  
  # Remove pre-existing results
  # rm -r -f ./${new_names[k]}_cell_5_${estimate[j]}_mx/*

  # Generate job scripts
  sed 's/ct_name/'${new_names[k]}'/g' eQTL_OneK1K_cell_5_${estimate[j]}_mx.qsub > eQTL_OneK1K_cell_5_${estimate[j]}_mx_${new_names[k]}.qsub 

  # Submit
  qsub eQTL_OneK1K_cell_5_${estimate[j]}_mx_${new_names[k]}.qsub

  done

done

##
