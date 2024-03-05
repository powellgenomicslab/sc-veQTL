####
library(data.table)

dir = "/share/ScratchGeneral/angxue/proj/vQTL/TensorQTL/scale_gene_wise/"

sig = read.table(paste0(dir, "dispersionQTL_sig_all_14_cell_types_all_estimates_filtered_candidates.txt"), header = T)
sig = sig[sig$candidate, ]

sig2 = sig[sig$Cell_type %in% c("B_IN", "B_MEM"),]


# B_IN
res = c()
data = c()
for(i in 1:22){

tmp = fread(paste0("/share/ScratchGeneral/angxue/proj/vQTL/TensorQTL/interaction/B_cells/B_IN_cell_5_dispersion_mx/OneK1K_B_IN_chr",i,".cis_qtl_top_assoc.txt.gz"),header=T)
tmp = as.data.frame(tmp)

res = rbind(res, tmp)

print(i)
}

res2 = res[res$phenotype_id %in% sig2$phenotype_id,]

# B_MEM
res = c()

for(i in 1:22){

tmp = fread(paste0("/share/ScratchGeneral/angxue/proj/vQTL/TensorQTL/interaction/B_cells/B_MEM_cell_5_dispersion_mx/OneK1K_B_MEM_chr",i,".cis_qtl_top_assoc.txt.gz"),header=T)
tmp = as.data.frame(tmp)

res = rbind(res, tmp)

print(i)
}

res2 = res[res$phenotype_id %in% sig2$phenotype_id,]
res2$Cell_type = "B_MEM"
write.table(res2,"deQTL_interaction_cell_state_results.txt",row.names=F,col.names=T,quote=F)


####
