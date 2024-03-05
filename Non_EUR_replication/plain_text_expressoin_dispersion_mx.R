##############################################################################
# Script information                                                      
# Title: Generate expression files of Perez et al. lupus data for TensorQTL
# Author: Angli XUE
# Date: 2023-12-02
# Description: This R script was written to generate expression files of Perez et al. lupus data for TensorQTL
# We combined the pre-processed pseudo-bulk matrix and gene location files and converted into bed format
##############################################################################

library(data.table, quietly = T)
library(R.utils, quietly = T)

#! /usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# chrNumber <- args[1]
ct <- as.numeric(args[1])

new_names <- c("B", "cDC", "cM", "ncM", "NK", "PB", "pDC", "Progen", "Prolif", "T4", "T8")
ct_name = new_names[ct]

# Read pseudo-bulk matrix
mx = fread(paste0("/directflow/SCCGGroupShare/projects/angxue/data/Perez_et_al_SLE/dispersion_est/",ct_name,"_cells_mean_mx.txt"), header = T)
mx = as.data.frame(mx)
# Remove individuals with all NAs
if(sum(is.na(mx[1, ]))>0){
index = which(is.na(mx[1, ]))
mx = mx[ ,-index]
}
pi0 = rowSums(mx==0)/ncol(mx)

# Read dispersion matrix
dx = fread(paste0("/directflow/SCCGGroupShare/projects/angxue/data/Perez_et_al_SLE/dispersion_est/",ct_name,"_cells_dispersion_mx.txt"), header = T)
dx = as.data.frame(dx)
gene = dx$gene
dx = dx[,-1]
dx = dx[,colnames(mx)]

mx = dx
# Remove genes with pi0 > 0.9
mx = mx[which(pi0 <= 0.9), ]

# log(x+1) and standardization
mx = apply(mx, 2, function(x)log(x+1))
# mx = as.data.frame(scale(mx))
mx = as.data.frame(t(scale(t(mx))))
if(sum(is.na(mx[1, ]))>0){
index = which(is.na(mx[1, ]))
mx = mx[ ,-index]
}

# Modify donor id
name=colnames(mx)
name=as.numeric(name)
index = which(!is.na(name))
new_id = paste0(name[index],"_",name[index])
colnames(mx)[index] = new_id

gene = as.data.frame(gene)
gene = as.data.frame(gene[which(pi0 <= 0.9),])

data = cbind(id=gene[,1],mx)

# Save
write.table(data, paste0(tolower(ct_name), "_cg.expr"), row.names = F, col.names = T, quote = F, sep = "\t")



####
