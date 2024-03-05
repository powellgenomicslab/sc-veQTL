## SGE SETTINGS
#$ -cwd
#$ -S /bin/bash
#$ -q short.q
#$ -pe smp 1
#$ -e ./ct_name/
#$ -o ./ct_name/
#$ -r yes
#$ -l mem_requested=45G
#$ -N get_MLE_ct_name_parameters
#$ -t 1-980


hostname

cd $SGE_O_WORKDIR

i=${SGE_TASK_ID};

touch ./ct_name/estimate_dispersion_${i}.log

# Rscript=/share/ScratchGeneral/angxue/software/R-4.0.5/bin/Rscript
Rscript=/share/ScratchGeneral/angxue/software/miniconda3/envs/r_env/lib/R/bin/Rscript

$Rscript --vanilla mle_per_gene_within_ind_disp.R ${i} ct_name >> ./ct_name/estimate_dispersion_${i}.log

