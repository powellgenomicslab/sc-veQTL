##############################################################################
# Script information                                                      
# Title: Generate expression files of OneK1K for TensorQTL
# Author: Angli XUE
# Date: 2022-06-22
# Description: This R script was written to generate expression files of OneK1K for TensorQTL
# We combined the pre-processed pseudo-bulk matrix and gene location files and converted into bed format
##############################################################################

library(data.table, quietly = T)
library(R.utils, quietly = T)

#! /usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

chrNumber <- args[1]
ct_name <- args[2]

new_names <- c("B_IN", "B_MEM", "CD4_NC", "CD4_ET", "CD4_SOX4", "CD8_NC", "CD8_ET",
    "CD8_S100B", "DC", "Mono_C", "Mono_NC", "NK_R", "NK", "Plasma",
    "Erythrocytes","Platelets")

# Read pseudo-bulk matrix
mx = fread(paste0("/share/ScratchGeneral/angxue/proj/vQTL/mean_var_matrix/cell_5/",ct_name,"_cells_mean_mx.txt"), header = T)
mx = as.data.frame(mx)
# Remove individuals with all NAs
if(sum(is.na(mx[1, ]))>0){
index = which(is.na(mx[1, ]))
mx = mx[ ,-index]
}
pi0 = rowSums(mx==0)/ncol(mx)

mx = fread(paste0("/share/ScratchGeneral/angxue/proj/vQTL/mean_var_matrix/cell_5/",ct_name,"_cells_var_mx.txt"), header = T)
mx = as.data.frame(mx)

# Remove genes with pi0 > 0.9
mx = mx[which(pi0 <= 0.9), ]
mx = mx[,-index]
# log(x+1) and standardization
mx = apply(mx, 2, function(x)log(x+1))
mx = as.data.frame(t(scale(t(mx))))

gene = read.table(paste0("/share/ScratchGeneral/angxue/proj/vQTL/mean_var_matrix/",ct_name,"_cells_gene_list.txt"), header = T)

gene = as.data.frame(gene[which(pi0 <= 0.9),])

loc = fread(paste0("/directflow/SCCGGroupShare/projects/angxue/data/onek1k/Gene_Location_Files/geneloc_chr",chrNumber,".tsv"), header = T)
loc = as.data.frame(loc)

mx = cbind(gene, mx)
loc = loc[ ,c(2:4, 1)]

# Extract the overlapped genes
id = Reduce(intersect, list(loc$geneid, mx$gene))
loc = loc[match(id, loc$geneid), ]
loc$chr = paste0("chr", loc$chr)
mx = mx[match(id, mx$gene), ]

# Combine location and expression
data = cbind(loc, mx[,-1])
data = data[order(data$start, decreasing = F), ]
colnames(data)[1]="#chr"

# Add the TSS position
# Replace the start with TSS location


# The BED file must define the TSS/cis-window center, with start+1 == end
# We used the gene center here
data$start = data$start + as.integer((data$end - data$start)/2)
data$end = data$start + 1


write.table(data, paste0("./",ct_name,"/OneK1K_980_samples_",ct_name,"_chr",chrNumber,".bed"), row.names = F, col.names = T, quote = F, sep = "\t")
gzip(paste0("./",ct_name,"/OneK1K_980_samples_",ct_name,"_chr",chrNumber,".bed"), destname = paste0("./",ct_name,"/OneK1K_980_samples_",ct_name,"_chr",chrNumber,".bed.gz"), overwrite = T)



####
