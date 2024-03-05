##############################################################################
# Script information
# Title: Perform TensorQTL for OneK1K
# Author: Angli XUE
# Date: 2022-06-23
# Description: This python script was written to perform TensorQTL analysis on OneK1K
##############################################################################
import pandas as pd
import numpy as np
import torch
import tensorqtl
from tensorqtl import genotypeio, cis, trans
print(f"PyTorch {torch.__version__}")
print(f"Pandas {pd.__version__}")
import os
import sys
from sys import argv
inputArgs = sys.argv[1:]
# Convert the argument to integer
chr_num = int(inputArgs[0])
ct_name = inputArgs[1]

os.makedirs(f"./{ct_name}_cell_5_variance_mx", exist_ok = True)
os.chdir(f"./{ct_name}_cell_5_variance_mx")

# define paths to data
plink_prefix_path = f"/directflow/SCCGGroupShare/projects/angxue/data/onek1k/genotype/plink_chr{chr_num}"
expression_bed = f"/share/ScratchGeneral/angxue/proj/vQTL/TensorQTL/expression/cell_5_variance_mx/{ct_name}/OneK1K_980_samples_{ct_name}_chr{chr_num}.bed.gz"
covariates_file = f"/share/ScratchGeneral/angxue/proj/vQTL/MatrixQTL/round2_more_covar/covariates/cell_5_mean_mx/{ct_name}_covar_peer_factors_PF10.txt"
# interaction_file = f"/share/ScratchGeneral/angxue/proj/vQTL/TensorQTL/interaction/pseudotime/pseudotime_{ct_name}.txt"
prefix = f"OneK1K_{ct_name}"

# load phenotypes and covariates
print("Load phenotypes and covariates")
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)
covariates_df = pd.read_csv(covariates_file, sep="\t", index_col=0).T
# interaction_s = pd.read_csv(interaction_file, sep="\t", index_col=0)
interaction_s = covariates_df.loc[:, ['age']]
covariates_df = covariates_df.drop('age', axis=1)

# PLINK reader for genotypes
pr = genotypeio.PlinkReader(plink_prefix_path)
genotype_df = pr.load_genotypes()
variant_df = pr.bim.set_index("snp")[["chrom", "pos"]]
# Chromosome name changed format
variant_df.chrom = "chr" + variant_df.chrom

## cis-QTL: nominal p-values for all variant-phenotype pairs
# map all cis-associations (results for each chromosome are written to file)

# all genes
# cis.map_nominal(genotype_df, variant_df, phenotype_df, phenotype_pos_df, covariates_df, prefix)

## The lengths must match to compare 
print("Match all the input files")
idx1 = covariates_df.index
idx2 = interaction_s.index
idx3 = genotype_df.columns
idx4 = phenotype_df.columns
idx = idx1.intersection(idx2).intersection(idx3).intersection(idx4)

covariates_df = covariates_df.loc[idx, :]
interaction_s = interaction_s.loc[idx, :]
genotype_df = genotype_df.loc[: ,idx]
phenotype_df = phenotype_df.loc[:, idx]

# interation test
prefix = f"OneK1K_{ct_name}_chr{chr_num}"
cis.map_nominal(genotype_df, variant_df, phenotype_df, phenotype_pos_df, prefix, covariates_df=covariates_df, interaction_df=interaction_s, maf_threshold_interaction=0.05, run_eigenmt=True, output_dir='.', write_top=True, write_stats=True)

# print("Saving the all association pairs")
# inter_df.to_csv(f"{prefix}.gxe_pseudotime_qtl_pairs.chr{chr_num}.csv", sep = "\t")


print("Analaysis finished!")


####
