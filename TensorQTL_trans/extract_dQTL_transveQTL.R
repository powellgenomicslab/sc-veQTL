####
library(data.table)
library(tidyr)
library(dplyr)
library(stringr)

WD = "/share/ScratchGeneral/angxue/proj/vQTL/TensorQTL/transeQTL/"

sig = read.table("../dispersionQTL_sig_all_14_cell_types_mean03_all_estimates.txt", header = T)

sig$nr_trans_eQTL = NA
sig$nr_trans_veQTL = NA

for(i in 1:nrow(sig)){

gene = sig$phenotype_id[i]
ct_name = sig$Cell_type[i]
snp = sig$variant_id[i]
chr = str_split(snp,":")[[1]][1]
bp = str_split(snp,":")[[1]][2]

## Extract the beta estimates from eQTL, veQTL, residQTL analyses

m = fread(paste0(WD, ct_name, "_cell_5_mean_mx/OneK1K_", ct_name, ".trans_qtl_pairs.chr", chr, ".csv"), header = T)
m = as.data.frame(m)
v = fread(paste0(WD, ct_name, "_cell_5_variance_mx/OneK1K_", ct_name, ".trans_qtl_pairs.chr", chr, ".csv"), header = T)
v = as.data.frame(v)

m2 = m[m$phenotype_id == gene & m$pval < 1e-5, ]
v2 = v[v$phenotype_id == gene & v$pval < 1e-5, ]

# Exclude trans-QTL in MHC region
m2$CHR=str_split(m2$variant_id,":")[[1]][1]
m2$BP=str_split(m2$variant_id,":")[[1]][2]

print(any(m2$CHR == 6 & m2$BP > 28477797 & m2$BP < 33448354))

m2 = m2[!(m2$CHR == 6 & m2$BP > 28477797 & m2$BP < 33448354),]

sig$nr_trans_eQTL[i] = nrow(m2)
sig$nr_trans_veQTL[i] = nrow(v2)

print(i)
}

write.table(sig, "dispersionQTL_with_trans_effects_noMHC.txt", row.names = F, col.names = T, quote = F)

new=read.table("../dispersionQTL_sig_all_14_cell_types_mean03_all_estimates_TSS_rsID.txt",header=T)
new$pair=paste0(new$variant_id,"_",new$phenotype_id)
sig$pair=paste0(sig$variant_id,"_",sig$phenotype_id)
sig2=sig[match(new$pair,sig$pair),]

write.table(sig2, "dispersionQTL_with_trans_effects_50_dQTLs_noMHC.txt", row.names = F, col.names = T, quote = F)

sig2[sig2$nr_trans_eQTL < sig2$nr_trans_veQTL,c("phenotype_id","nr_trans_eQTL","nr_trans_veQTL")]

####
