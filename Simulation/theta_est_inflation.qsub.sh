## SGE SETTINGS
#$ -cwd
#$ -S /bin/bash
#$ -q short.q
#$ -pe smp 1
#$ -r yes
#$ -l mem_requested=30G
#$ -N simulation_theta_inflation
#$ -t 1-100

hostname

cd $SGE_O_WORKDIR

i=${SGE_TASK_ID};

Rscript=/share/ScratchGeneral/angxue/software/miniconda3/envs/py37/bin/Rscript

$Rscript --vanilla theta_est_inflation.R ${i}   

