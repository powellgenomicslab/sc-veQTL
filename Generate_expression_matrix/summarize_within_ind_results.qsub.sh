## SGE SETTINGS
#$ -cwd
#$ -S /bin/bash
#$ -q short.q
#$ -pe smp 1
#$ -e stderr_summarize_within_ind_dispersion
#$ -o stdout_summarize_within_ind_dispersion
#$ -r yes
#$ -l mem_requested=20G
#$ -N summarize_within_ind_dispersion
#$ -t 1-22

hostname

cd $SGE_O_WORKDIR

i=${SGE_TASK_ID};

# i=${SGE_TASK_ID};

Rscript=/share/ScratchGeneral/angxue/software/R-4.0.5/bin/Rscript

$Rscript --vanilla summarize_within_ind_results.R ${i}
