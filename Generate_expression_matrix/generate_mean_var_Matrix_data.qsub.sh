## SGE SETTINGS
#$ -cwd
#$ -S /bin/bash
#$ -q short.q
#$ -pe smp 1
#$ -r yes
#$ -l mem_requested=80G
#$ -N generate_mean_var
#$ -t 1-14

cd $SGE_O_WORKDIR

i=${SGE_TASK_ID};

# touch ${i}.txt

# module load briglo/R/3.6.0

# Rscript=/share/ScratchGeneral/angxue/software/R-4.0.5/bin/Rscript

Rscript=/share/ScratchGeneral/angxue/software/miniconda3/envs/r_env/lib/R/bin/Rscript

touch generate_mean_var_matrix_cell_type_${i}.log
echo ${JOB_ID} >> generate_mean_var_matrix_cell_type_${i}.log

$Rscript --vanilla generate_mean_var_Matrix_data.R ${i} >> generate_mean_var_matrix_cell_type_${i}.log  

