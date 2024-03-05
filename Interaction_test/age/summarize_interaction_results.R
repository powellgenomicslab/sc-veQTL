####
library(data.table)

new_names <- c("B_IN", "B_MEM", "CD4_NC", "CD4_ET", "CD4_SOX4", "CD8_NC", "CD8_ET",
    "CD8_S100B", "DC", "Mono_C", "Mono_NC", "NK_R", "NK", "Plasma")

dir = "/share/ScratchGeneral/angxue/proj/vQTL/TensorQTL/scale_gene_wise/"

sig = read.table(paste0(dir, "dispersionQTL_sig_all_14_cell_types_all_estimates_filtered_candidates.txt"), header = T)
sig = sig[sig$candidate, ]

# sig2 = sig[sig$Cell_type %in% c("B_IN", "B_MEM"),]
data = c()
for(k in 1:14){

res = c()

for(i in 1:22){
file_name = paste0("/share/ScratchGeneral/angxue/proj/vQTL/TensorQTL/interaction/age/",new_names[k],"_cell_5_dispersion_mx/OneK1K_",new_names[k],"_chr",i,".cis_qtl_top_assoc.txt.gz")
if(file.exists(file_name)){
tmp = fread(file_name, header = T)
tmp = as.data.frame(tmp)

res = rbind(res, tmp)
}
# print(i)
}

sig2 = sig[sig$Cell_type==new_names[k],]
res2 = res[res$phenotype_id %in% sig2$phenotype_id,]
print(new_names[k])
print(sum(res2$pval_adj_bh<0.05))
if(sum(res2$pval_adj_bh<0.05)>0){
res2$Cell_type = new_names[k]
print(res2[res2$pval_adj_bh<0.05,])
data = rbind(data,res2)
}

}

write.table(data, "deQTL_interaction_age_results.txt", row.names=F,col.names=T,quote=F)
####
