####
#! /usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

## Combine all the within-indivdiual variance estimation
library(data.table)

new_names <- c("B_IN", "B_MEM", "CD4_NC", "CD4_ET", "CD4_SOX4", "CD8_NC", "CD8_ET",
    "CD8_S100B", "DC", "Mono_C", "Mono_NC", "NK_R", "NK", "Plasma",
    "Erythrocytes","Platelets")

j <- as.numeric(args[1])

a1=fread(paste0("/share/ScratchGeneral/angxue/proj/vQTL/mean_var_matrix/cell_5/",new_names[j],"_cells_var_mx.txt"),header=T)
a1=as.data.frame(a1)
# Remove those individuals with NAs
ori_name=colnames(a1)
a1=a1[,!is.na(as.numeric(a1[1,]))]

name=colnames(a1)

a=fread(paste0("./",new_names[j],"/",new_names[j],"_within_ind_dispersion_",name[1],".txt"),header=T,check.names=F)
a=as.data.frame(a)

for(i in 2:length(name)){
b=fread(paste0("./",new_names[j],"/",new_names[j],"_within_ind_dispersion_",name[i],".txt"),header=T,check.names=F)
b=as.data.frame(b)

a=cbind(a,b)
}

a=cbind(a[,1],a[,c(name)])
colnames(a)[1]="gene"

# Summarize results matrix
res=as.data.frame(matrix(NA,nrow=nrow(a1),ncol=length(ori_name)+1))
colnames(res)=c("gene",ori_name)
res[,colnames(a)]=a

write.table(res,paste0("./",new_names[j],"_cells_dispersion_mx.txt"),row.names=F,col.names=T,quote=F)

message(new_names[j])
message(paste0(nrow(a)," genes"))


####
