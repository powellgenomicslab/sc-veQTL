## SGE SETTINGS
#$ -cwd
#$ -S /bin/bash
#$ -q short.q
#$ -pe smp 1
#$ -r yes
#$ -l mem_requested=140G
#$ -N slingshot_B_cells
#$ -o stdout_slingshot_B_cells
#$ -e stderr_slingshot_B_cells

cd $SGE_O_WORKDIR

#i=${SGE_TASK_ID};

conda activate py37

Rscript413=/share/ScratchGeneral/angxue/software/miniconda3/envs/py37/bin/Rscript

touch slingshot.log

$Rscript413 phate_and_slingshot_B_cells.R >> slingshot.log


