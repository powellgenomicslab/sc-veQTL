#! /usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
k = as.integer(args[1])
ct = args[2]

library(glmGamPoi)
library(fitdistrplus)
library(dplyr)
library(Seurat)

new_names <- c("B_IN", "B_MEM", "CD4_NC", "CD4_ET", "CD4_SOX4", "CD8_NC", "CD8_ET",
    "CD8_S100B", "DC", "Mono_C", "Mono_NC", "NK_R", "NK", "Plasma",
    "Erythrocytes","Platelets")

cell <- c("TCL1A+FCER2+Bcell", "TCL1A-FCER2-Bcell", "CD4+KLRB1-Tcell", "CD4+KLRB1+Tcell",
    "CD4+SOX4+Tcell", "CD8+LTB+Tcell", "CD8+GNLY+NKG7+Tcell", "CD8+S100B+Tcell",
    "Dendriticcell", "MonocyteCD14+", "MonocyteFCGR3A+", "XCL1+NK", "XCL1-NK", "IgJ+Bcell",
    "Erythrocytes","Platelets")

ct_name = as.character(ct)

## Read moment estimators input as priori
# priori=read.table(paste0("/directflow/SCCGGroupShare/projects/angxue/data/onek1k/SCTv2/across_ind_estimates/NB_MOM_estimates_aross_ind_",ct_name,".txt"),header=T)

info=readRDS(paste0("/directflow/SCCGGroupShare/projects/angxue/data/onek1k/",ct_name,"_meta_data.RDS"))
# dim(info)

## Read count matrix
message("Read count matrix")

sct=readRDS(paste0("/directflow/SCCGGroupShare/projects/angxue/data/onek1k/",ct_name,"_sct_counts.RDS"))
# Check data structure
str(sct)

# Convert sct counts to Seurat object
message("Convert sct counts to Seurat object")
ok1k <- CreateSeuratObject(counts = sct, assay = "SCT", project = "onek1k", min.cells = 0, min.features = 0)
# ok1k

info=info[info$UMI %in% rownames(ok1k@meta.data),]

# Add back individual ID from meta.data
ok1k@meta.data$individual=info$individual

# Extract a list of individual ID and their number of cells
a=unlist(table(info$individual))
ind=names(a)

## For each individual, subset the data set
res=as.data.frame(matrix(NA,nrow=nrow(ok1k@assays$SCT@counts),ncol=2))

colnames(res)=c("gene",ind[k])
res$gene = rownames(ok1k@assays$SCT@counts)

message("Estimate the within-individual dispersion")

ss=subset(x = ok1k, subset = individual == ind[kk])
tmp=ss@assays$SCT@counts
tmp2=as.data.frame(tmp)

if(ncol(tmp2)<5){
res[,2]=NA
}else{
fit = glmGamPoi::glm_gp(as.matrix(tmp2), size_factors = FALSE, verbose = T)
res[,2] = fit$overdispersion_shrinkage_list$ql_disp_shrunken
}
write.table(res,paste0("./",ct_name,"/",ct_name,"_within_ind_dispersion_",ind[k],".txt"),row.names=F,col.names=T,quote=F)



#### END
