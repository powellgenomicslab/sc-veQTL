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

os.makedirs(f"./{ct_name}_cell_5_dispersion_mx", exist_ok = True)
os.chdir(f"./{ct_name}_cell_5_dispersion_mx")

# define paths to data
plink_prefix_path = f"/directflow/SCCGGroupShare/projects/angxue/data/onek1k/genotype/plink_chr{chr_num}"
expression_bed = f"/share/ScratchGeneral/angxue/proj/vQTL/TensorQTL/expression/cell_5_dispersion_mx/{ct_name}/OneK1K_980_samples_{ct_name}_chr{chr_num}.bed.gz"
covariates_file = f"/share/ScratchGeneral/angxue/proj/vQTL/MatrixQTL/round2_more_covar/covariates/cell_5_mean_mx/{ct_name}_covar_peer_factors_PF10.txt"
freq_file = f"/directflow/SCCGGroupShare/projects/angxue/data/onek1k/genotype/plink_chr{chr_num}.frqx"
prefix = f"OneK1K_{ct_name}"

# load phenotypes and covariates
print("Load phenotypes and covariates")
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)
covariates_df = pd.read_csv(covariates_file, sep="\t", index_col=0).T

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
covariates_df = covariates_df.loc[phenotype_df.columns, :]
genotype_df = genotype_df.loc[: ,phenotype_df.columns]

# genes on chr22
cis.map_nominal(genotype_df, variant_df,
                phenotype_df.loc[phenotype_pos_df["chr"] == "chr" + str(chr_num)],
                phenotype_pos_df.loc[phenotype_pos_df["chr"] == "chr" + str(chr_num)],
                prefix, covariates_df = covariates_df,
                write_top = True, write_stats = True,
		maf_threshold = 0.05)

# load results
pairs_df = pd.read_parquet(f"{prefix}.cis_qtl_pairs.chr{chr_num}.parquet")
pairs_df.head()

print("Saving the all association pairs")
pairs_df.to_csv(f"{prefix}.cis_qtl_pairs.chr{chr_num}.csv", sep = "\t")

## cis-QTL: empirical p-values for phenotypes
# all genes
# cis_df = cis.map_cis(genotype_df, variant_df, phenotype_df, phenotype_pos_df, covariates_df)
print("Start permutation test")
# genes on a specific chromosome
# This only output the top signal per gene
cis_df = cis.map_cis(genotype_df, variant_df, 
                     phenotype_df.loc[phenotype_pos_df["chr"] == "chr" + str(chr_num)],
                     phenotype_pos_df.loc[phenotype_pos_df["chr"] == "chr" + str(chr_num)],
                     covariates_df = covariates_df, seed = 123456,
                     maf_threshold = 0.05, nperm = 10000)

# calculate chromosome-wide FDR
# need to understand why the author sets lamdba = 0.85
tensorqtl.calculate_qvalues(cis_df, qvalue_lambda = 0.85)

cis_df.head()

# read allele and freq information
# allele_freq_df = pd.read_csv(freq_file, sep="\t", index_col=0)
# allele_freq_df["freq"] = (allele_freq_df["C(HOM A1)"] * 2 + allele_freq_df["C(HET)"]) / ( phenotype_df.shape[1] * 2)
# allele_freq_df = allele_freq_df[["SNP","A1","A2","freq"]]
# allele_freq_df = allele_freq_df.rename(columns = {'SNP':'variant_id'})

# index = cis_df.index
# cis_df = pd.merge(cis_df,allele_freq_df,on="variant_id")
# cis_df.index = index

# save the results
print("Saving the significant association pairs")
cis_df.to_csv(f"{prefix}.sig_cis_qtl_pairs.chr{chr_num}.csv", sep = "\t")

# conditionally independent QTLs
if any(cis_df.qval < 0.05):
	print("Identify conditionally independent QTLs")
	indep_df = cis.map_independent(genotype_df, variant_df, cis_df,
                               phenotype_df, phenotype_pos_df, 
                               covariates_df = covariates_df,
                               maf_threshold = 0.05, nperm = 10000,
                               fdr=0.05, fdr_col="qval")

	# index = indep_df.index
	# indep_df = pd.merge(indep_df,allele_freq_df,on="variant_id")
	# indep_df.index = index

	print("Saving the conditionally independent QTLs")
	indep_df.to_csv(f"{prefix}.independent_cis_qtl_pairs.chr{chr_num}.csv", sep = "\t")

else :
	print("No significant eQTLs with qval < 0.05. Do not identify conditionally independent QTLs.")


print("Analaysis finished!")


####
