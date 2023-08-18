## SGE SETTINGS
#$ -cwd
#$ -S /bin/bash
#$ -q short.q
#$ -pe smp 1
#$ -r yes
#$ -l mem_requested=250G
#$ -N extract_RDS
#$ -o stdout_extract_RDS
#$ -e stderr_extract_RDS
#$ -t 1-16

hostname

cd $SGE_O_WORKDIR

i=${SGE_TASK_ID};

Rscript=/share/ScratchGeneral/angxue/software/R-4.0.5/bin/Rscript

touch cell_type_${i}.log

$Rscript --vanilla Extract_RDS_all_cell_types.R ${i} >> cell_type_${i}.log




