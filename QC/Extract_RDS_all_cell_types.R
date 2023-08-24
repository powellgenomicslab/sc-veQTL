# The code is revised based on /directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/Scripts/Scripts_20201028/Prepare_matrices_for_variance.R from Seyhan, 25-Mar-2021

#! /usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

print(paste0("Read in argument, ",args))

#quit(save="no")

print("Script started!")
library(Seurat)
library(sctransform)

print("Reading in cell_type.RDS...")
df <- readRDS ("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/cell_specific_eQTL_analysis/cell_type.RDS")

df
df$individual[which(df$individual=="870_871" & df$latent=="b1")] <- "966_967"
# Remove two outliers newly identified on May 16th 2022
keep_list=unique(df$individual)
keep_list=keep_list[!keep_list %in% c("88_88","966_967")]
df=subset(x = df, subset = individual %in% keep_list)

df@meta.data$cell <- gsub('\\s+', '', df@meta.data$cell_type)
head(df@meta.data)

# Cell types
cell <- c("TCL1A+FCER2+Bcell", "TCL1A-FCER2-Bcell", "CD4+KLRB1-Tcell", "CD4+KLRB1+Tcell", 
    "CD4+SOX4+Tcell", "CD8+LTB+Tcell", "CD8+GNLY+NKG7+Tcell", "CD8+S100B+Tcell", 
    "Dendriticcell", "MonocyteCD14+", "MonocyteFCGR3A+", "XCL1+NK", "XCL1-NK", "IgJ+Bcell",
    "Erythrocytes","Platelets")
new_names <- c("B_IN", "B_MEM", "CD4_NC", "CD4_ET", "CD4_SOX4", "CD8_NC", "CD8_ET",
    "CD8_S100B", "DC", "Mono_C", "Mono_NC", "NK_R", "NK", "Plasma",
    "Erythrocytes","Platelets")
cells <- data.frame(cell, new_names)

CellsMeta = df@meta.data
head(CellsMeta)

library(dplyr)
CellsMeta <- left_join(CellsMeta, cells, by="cell")
head(CellsMeta)

CellsMetaTrim <- subset(CellsMeta, select = c("new_names"))
head(CellsMetaTrim)

df <- AddMetaData(
     object = df,
     metadata = CellsMetaTrim$new_names,
     col.name = 'new_names'
)
head(x = df@meta.data)

Idents(df) <- "new_names"

## Extract the subset based on cell types
print("Extracting the subset...")

ct_name = as.character(new_names[as.numeric(args[1])])

print(ct_name)

# ct_name = "CD4_NC"
df <- df[, WhichCells(df, ident = ct_name)]
# Generate mean and vairance matrix based on SCT counts
# info=df_meta %>% filter(new_names==ct_name)
# library(dplyr)

# Save the Seurat object 
saveRDS(df, paste0(ct_name,"_QCed_onek1k.RDS"))

# Temp command
# quit(save="no")

# Save SCT counts
sct.data <- (GetAssayData(df, slot = "counts")[, WhichCells(df, ident = ct_name)])
print("Saving sct data...")
saveRDS(sct.data, paste0(ct_name,"_sct_counts.RDS"))
# Why this step?
# The default assay for this RDS file is SCT. So now we change it back to RNA
DefaultAssay(object = df) <- "RNA"
DefaultAssay(object = df)

# Add UMI in the meta data
df@meta.data$UMI=rownames(df@meta.data)

df_meta <- df@meta.data
df_meta <- df_meta %>% filter(new_names==ct_name)

# Save meta data
print("Saving meta data...")
saveRDS(df_meta, paste0(ct_name,"_meta_data.RDS"))

# Normalize the counts
df2 <- NormalizeData(df)

rna.data <- (GetAssayData(df2, slot = "data")[, WhichCells(df, ident = ct_name)])
print("Saving rna data...")
saveRDS(rna.data, paste0(ct_name,"_rna_counts.RDS"))

raw.data <- (GetAssayData(df2, slot = "counts")[, WhichCells(df, ident = ct_name)])
print("Saving raw count data...")
saveRDS(raw.data, paste0(ct_name,"_raw_counts.RDS"))


print("Script ended!")

## END


