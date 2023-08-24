## SGE SETTINGS
#$ -cwd
#$ -S /bin/bash
#$ -q short.q
#$ -pe smp 1
#$ -r yes
#$ -l mem_requested=30G
#$ -N generate_expression_file_cell_5_mean_mx_ct_name
#$ -o stdout_generate_expression_file_cell_5_mean_mx_ct_name
#$ -e stderr_generate_expression_file_cell_5_mean_mx_ct_name
#$ -t 1-22

cd $SGE_O_WORKDIR

i=${SGE_TASK_ID};

Rscript="/share/ScratchGeneral/angxue/software/miniconda3/envs/r_env/bin/Rscript"
dir="/share/ScratchGeneral/angxue/proj/vQTL/TensorQTL/expression/"

mkdir -p ct_name

# touch ./ct_name/chr${i}_generate_expression_file_cell_5_mean_mx_ct_name.log
# > ./ct_name/chr${i}_generate_expression_file_cell_5_mean_mx_ct_name.log

$Rscript --vanilla ${dir}generate_expression_file_cell_5_mean_mx.R ${i} ct_name 





