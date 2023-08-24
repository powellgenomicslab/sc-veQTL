## SGE SETTINGS
#$ -cwd
#$ -S /bin/bash
#$ -q short.q
#$ -pe smp 1
#$ -r yes
#$ -l mem_requested=50G
#$ -N TensorQTL_ct_name
#$ -o ./ct_name_cell_5_dispersion_mx/stdout_eQTL_TensorQTL_ct_name
#$ -e ./ct_name_cell_5_dispersion_mx/stderr_eQTL_TensorQTL_ct_name
#$ -t 1-22

cd $SGE_O_WORKDIR

i=${SGE_TASK_ID};

conda activate py37

mkdir -p ./ct_name_cell_5_dispersion_mx/

touch ./ct_name_cell_5_dispersion_mx/chr${i}_run_TensorQTL_ct_name_cell_5_dispersion_mx.log

python eQTL_OneK1K_cell_5_dispersion_mx.py ${i} ct_name >> ./ct_name_cell_5_dispersion_mx/chr${i}_run_TensorQTL_ct_name_cell_5_dispersion_mx.log 2>&1




