#### OneK1K data set ####

# Generating the intra-individual mean and variance matrix from SCTransformed count matrix

print(Sys.time())
start <- Sys.time()

#! /usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
print(paste0("Read in argument, ", args))

library(MASS)
setwd("/share/ScratchGeneral/angxue/proj/vQTL/mean_var_matrix/cell_5/")

# Cell type names
new_names <- c("B_IN", "B_MEM", "CD4_NC", "CD4_ET", "CD4_SOX4", "CD8_NC", "CD8_ET", "CD8_S100B", "DC", "Mono_C", "Mono_NC", "NK_R", "NK", "Plasma", "Erythrocytes", "Platelets")
cell <- c("TCL1A+FCER2+Bcell", "TCL1A-FCER2-Bcell", "CD4+KLRB1-Tcell", "CD4+KLRB1+Tcell", "CD4+SOX4+Tcell", "CD8+LTB+Tcell", "CD8+GNLY+NKG7+Tcell", "CD8+S100B+Tcell", "Dendriticcell", "MonocyteCD14+", "MonocyteFCGR3A+", "XCL1+NK", "XCL1-NK", "IgJ+Bcell", "Erythrocytes", "Platelets")

ct_name <- as.character(new_names[as.numeric(args[1])])
print(ct_name)

## Read the data ##
info <- readRDS(paste0("/directflow/SCCGGroupShare/projects/angxue/data/onek1k/", ct_name, "_meta_data.RDS"))
dim(info)

## Number of cells and individuals
print(paste("This data set includes", dim(info)[1], ct_name, "cells from", nrow(table(info$individual)), "individuals"))

print("Number of cells per individual: ")
print(summary(as.numeric(table(as.character(info$individual)))))

## Read sctransformed counts 
print("Start to read in sct counts...")

# The input of “*_sct_counts.RDS” is the exact output from line 78 in https://github.com/powellgenomicslab/PEER_factors/blob/main/1-Extract_datasets/Extract_RDS_all_cell_types.R
sct <- readRDS(paste0("/directflow/SCCGGroupShare/projects/angxue/data/onek1k/", ct_name, "_sct_counts.RDS"))
# Check data structure
str(sct)

## Load packages
library(dplyr)
library(Seurat)

# Convert sct counts to Seurat object
ok1k <- CreateSeuratObject(counts = sct, assay = "SCT", project = "onek1k", min.cells = 0, min.features = 0)
ok1k

# Match Seurat object with previously saved meta.data
# Number might be inconsistent due to filters (e.g., min.cells, min.features)
info <- info[info$UMI %in% rownames(ok1k@meta.data), ]

# Add back individual ID from meta.data 
ok1k@meta.data$individual <- info$individual

# Extract a list of individual ID and their number of cells
a <- unlist(table(info$individual))

# pdf(paste0(ct_name, "_cells_per_indi.pdf"), width = 6.08, height = 4.59)
# hist(a, xlab = "Cells per individual", main = paste0("Number of ", ct_name, " cells per individual"))
# dev.off()

### Identification of highly variable features (feature selection)
ok1k_tmp <- FindVariableFeatures(object = ok1k, selection.method = 'vst', nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(x = VariableFeatures(object = ok1k_tmp), 10)

# plot variable features with and without labels
# plot1 <- VariableFeaturePlot(ok1k_tmp)
# plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

# pdf(paste0(ct_name, "_cells_var_feat.pdf"), width = 12, height = 5)
# plot1 + plot2
# dev.off()

## Check mean-variance relationship within an individual
# To ensure the estimated within-ind variance is accurate, we excluded those individuals with < 5 cells
keep_list <- names(which(a >= 5))
ok1k_tmp <- subset(x = ok1k_tmp, subset = individual %in% keep_list)
print(paste0(length(keep_list), " individuals remain after excluding those with < 5 cells."))
# Update the meta.data (info)
info <- info[info$UMI %in% rownames(ok1k_tmp@meta.data), ]
print(summary(as.numeric(table(as.character(info$individual)))))

# Use the individual with the most cells as matrix format reference
ok1k_max_ind <- subset(x = ok1k_tmp, subset = individual == names(which.max(a)))
ok1k_max_ind <- FindVariableFeatures(object = ok1k_max_ind, selection.method = 'vst', nfeatures = 2000)
hvf.info <- HVFInfo(object = ok1k_max_ind, selection.method = 'vst')

ind <- unique(info$individual)
ind <- ind[order(ind)]
gene <- rownames(hvf.info)
print(paste0(length(ind), " individuals and ", length(gene), " genes included in mean and variance matrices!"))

# write.table(as.data.frame(gene), paste0(ct_name, "_cells_gene_list.txt"), row.names = FALSE, col.names = TRUE, quote = FALSE)

mean_mx <- matrix(NA, nrow = nrow(hvf.info), ncol = length(ind))
var_mx <- matrix(NA, nrow = nrow(hvf.info), ncol = length(ind))
var_std_mx <- matrix(NA, nrow = nrow(hvf.info), ncol = length(ind))
mean_mx <- as.data.frame(mean_mx)
var_mx <- as.data.frame(var_mx)
var_std_mx <- as.data.frame(var_std_mx)

rownames(mean_mx) <- rownames(var_mx) <- rownames(var_std_mx) <- rownames(hvf.info)
colnames(mean_mx) <- colnames(var_mx) <- colnames(var_std_mx) <- ind

for (i in 1:length(ind)) {
  # subset the data within one individual
  ok1k_ind <- subset(x = ok1k_tmp, individual == ind[i])
  # If nr of cell = 1, it will return error because no way to calculate 'variance'
  ok1k_ind <- FindVariableFeatures(object = ok1k_ind, selection.method = 'vst', nfeatures = 2000)
  # Calculate the mean and variance
  hvf.info <- HVFInfo(object = ok1k_ind, selection.method = 'vst')

  mean_mx[, i] <- hvf.info$mean
  var_mx[, i] <- hvf.info$variance

  print(i)
}

# print(paste0(length(gene), " genes included in the main analysis"))
write.table(as.data.frame(gene), paste0(ct_name, "_cells_gene_list.txt"), row.names = FALSE, col.names = TRUE, quote = FALSE)

# Match the column with covariates file
co <- read.table(paste0("/share/ScratchGeneral/seyyaz/onek1k/cell_specific_eQTL_analysis_October19/", cell[which(new_names == ct_name)], "/step1/covariates_chr1.txt"), header = TRUE, check.names = FALSE)
all(colnames(var_mx) %in% colnames(co)[-1])

extra <- colnames(co)[!colnames(co) %in% colnames(var_mx)]
extra <- extra[-1]
tmp <- matrix(NA, nrow = nrow(var_mx), ncol = length(extra))
tmp <- as.data.frame(tmp)
colnames(tmp) <- extra

var_mx <- cbind(var_mx, tmp)
mean_mx <- cbind(mean_mx, tmp)

print(paste0(ncol(co) - 1, " individuals matched with covariate files."))
all(colnames(var_mx) %in% colnames(co)[-1])
all(colnames(co)[-1] %in% colnames(var_mx))

snbuf <- match(colnames(co)[-1], colnames(var_mx))

var_mx <- var_mx[, snbuf]
mean_mx <- mean_mx[, snbuf]

if (all(colnames(var_mx) == colnames(co)[-1])) {
  print("Column match checked!")
}

## Save results
write.table(mean_mx, paste0(ct_name, "_cells_mean_mx.txt"), row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(var_mx, paste0(ct_name, "_cells_var_mx.txt"), row.names = FALSE, col.names = TRUE, quote = FALSE)
# write.table(var_std_mx, "var_std_mx.txt", row.names = FALSE, col.names = TRUE, quote = FALSE)

print("Script ends.")
print(Sys.time())
end <- Sys.time()
print(difftime(end, start, units = "mins"))

#### END ####
