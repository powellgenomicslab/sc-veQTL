## SGE SETTINGS
#$ -cwd
#$ -S /bin/bash
#$ -q short.q
#$ -pe smp 1
#$ -r yes
#$ -l mem_requested=30G
#$ -N plain_text_expressoin
#$ -o stdout_plain_text_expressoin
#$ -e stderr_plain_text_expressoin
#$ -t 1-11

cd $SGE_O_WORKDIR

i=${SGE_TASK_ID};

Rscript="/share/ScratchGeneral/angxue/software/miniconda3/envs/r_env/bin/Rscript"
dir="/share/ScratchGeneral/angxue/proj/vQTL/TensorQTL/EAS_replication/expression/"

# mkdir -p ct_name

# touch ./ct_name/chr${i}_plain_text_expressoin.log
# > ./ct_name/chr${i}_plain_text_expressoin.log

$Rscript --vanilla ${dir}plain_text_expressoin_dispersion_mx.R ${i} 





