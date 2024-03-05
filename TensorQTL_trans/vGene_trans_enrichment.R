####

library(data.table)

new_names <- c("B_IN", "B_MEM", "CD4_NC", "CD4_ET", "CD4_SOX4", "CD8_NC", "CD8_ET",
    "CD8_S100B", "DC", "Mono_C", "Mono_NC", "NK_R", "NK", "Plasma",
    "Erythrocytes","Platelets")

res = data.frame(Cell_type = new_names[1:14], nr_genes = NA, eGene = NA, vGene = NA, prop_V_no_ME = NA, prop_V_with_ME = NA, p_val = NA)

for(i in 1:14){

old = c()

for(k in 1:22){

 tmp = fread(paste0("/share/ScratchGeneral/angxue/proj/vQTL/TensorQTL/",new_names[i],"_cell_5_mean_mx/OneK1K_",new_names[i],".sig_cis_qtl_pairs.chr",k,".csv"), header = T,sep = "\t")
 tmp = as.data.frame(tmp)
 if(ncol(tmp) == 18){tmp = tmp[,-18]}
 old = rbind(old,tmp)
# print(i)

}

old$Cell_type = new_names[i]

new = c()

for(k in 1:22){

 tmp = fread(paste0("/share/ScratchGeneral/angxue/proj/vQTL/TensorQTL/",new_names[i],"_cell_5_variance_mx/OneK1K_",new_names[i],".sig_cis_qtl_pairs.chr",k,".csv"), header = T,sep = "\t")
 tmp = as.data.frame(tmp)
 if(ncol(tmp) == 18){tmp = tmp[,-18]}
 new = rbind(new,tmp)

}
new$Cell_type = new_names[i]

new$qval2 = old$qval

res$nr_genes[i] = nrow(new)
res$eGene[i] = sum(old$qval < 0.05)
res$vGene[i] = sum(new$qval < 0.05)

## Trans-eQTL results
trans=c()

for(k in 1:22){
 tmp = fread(paste0("./", new_names[i], "_cell_5_mean_mx/OneK1K_", new_names[i], ".trans_qtl_pairs.chr",k,".csv"),header=T)
 tmp = as.data.frame(tmp)
 trans = rbind(trans,tmp)
}

# trans = trans[trans$pval<5e-8,]

index1 = which(new$qval<0.05 & new$qval2>0.05)
index2 = which(new$qval<0.05 & new$qval2<0.05)

x1 = sum(new$phenotype_id[index1] %in% trans$phenotype_id)
x2 = sum(new$phenotype_id[index2] %in% trans$phenotype_id)

# x1 = sum(new$variant_id[index1] %in% trans$variant_id)
# x2 = sum(new$variant_id[index2] %in% trans$variant_id)

tres = prop.test(x=c(x1,x2), n=c(length(index1),length(index2)))
tryCatch(prop.test(x=c(x1,x2), n=c(length(index1),length(index2))),warning=function(x)print("Chi-squared approximation may be incorrect"))

res$prop_V_no_ME[i] = as.numeric(tres$estimate[1])
res$prop_V_with_ME[i] = as.numeric(tres$estimate[2])
res$p_val[i] = as.numeric(tres$p.value)

print(new_names[i])

}


####
