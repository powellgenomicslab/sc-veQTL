##############################################################################
# Script information
# Title: Estimate the cell lineage and pseudotime for B cells in OneK1K
# Author: Angli Xue
# Date: 2022-08-25
# Description: This R script was written to estimate the cell lineage and pseudotime for B cells in OneK1K
# This script was modified based on the script from Jose Alquicira Hernandez
##############################################################################

## Import libraries
# Primary
library("tidyverse")
library("dsLib")

# Secondary
library("Seurat")
library("future")
library("phateR")
library("slingshot")
library("dittoSeq")
library("dplyr")

source("plot_functions.R")

# Settings
new_names <- c("B_IN", "B_MEM", "CD4_NC", "CD4_ET", "CD4_SOX4", "CD8_NC", "CD8_ET",
    "CD8_S100B", "DC", "Mono_C", "Mono_NC", "NK_R", "NK", "Plasma",
    "Erythrocytes","Platelets")

cell <- c("TCL1A+FCER2+Bcell", "TCL1A-FCER2-Bcell", "CD4+KLRB1-Tcell", "CD4+KLRB1+Tcell",
    "CD4+SOX4+Tcell", "CD8+LTB+Tcell", "CD8+GNLY+NKG7+Tcell", "CD8+S100B+Tcell",
    "Dendriticcell", "MonocyteCD14+", "MonocyteFCGR3A+", "XCL1+NK", "XCL1-NK", "IgJ+Bcell",
    "Erythrocytes","Platelets")
# Future settings  
options(future.globals.maxSize = 1024^3*100, future.seed = TRUE)
plan(multicore(workers = 2))

# Output settings
output <- set_output("2022-08-31","b_cells_phate")

# count1 = readRDS(paste0("/directflow/SCCGGroupShare/projects/angxue/data/onek1k/B_IN_QCed_onek1k.RDS"))
# count2 = readRDS(paste0("/directflow/SCCGGroupShare/projects/angxue/data/onek1k/B_MEM_QCed_onek1K.RDS"))

# Merge two B cells
# data <- merge(count1, y = count2, project = "onek1k")
# data@meta.data = data@meta.data[,-which(colnames(data@meta.data) %in% c("cell_type"))]

# colnames(data@meta.data)[ncol(data@meta.data)] = "cell_type"

# Read the intermediate b cells
# b_cells = readRDS(here(output, "b_cells_interm.RDS"))
b_cells = readRDS("./results/2022-08-30_b_cells_phate/b_cells_interm.RDS")
b_cells$individual[which(b_cells$individual == "870_871" & b_cells$latent=="b1")] <- "966_967"
keep_list = unique(b_cells$individual)
keep_list = keep_list[!keep_list %in% c("88_88","966_967")]
b_cells = subset(x = b_cells, subset = individual %in% keep_list)

b_cells <- SCTransform(b_cells,
                       vars.to.regress = c("pool", "percent.mt"),
                       conserve.memory = TRUE)

b_cells <- RunPCA(b_cells)
b_cells <- RunUMAP(b_cells, dims = 1:30)

pca <- DimPlot(b_cells, reduction = "pca", group.by = "cell_type")
umap <- DimPlot(b_cells, reduction = "umap", group.by = "cell_type")

ggsave(here(output, "pca.png"), pca, width = 7, height = 5)
ggsave(here(output, "umap.png"), umap, width = 7, height = 5)

pca_pool <- DimPlot(b_cells, reduction = "pca", group.by = "pool") + no_legend()
umap_pool <- DimPlot(b_cells, reduction = "umap", group.by = "pool") + no_legend()

ggsave(here(output, "pca_pool.png"), pca_pool, width = 7, height = 5)
ggsave(here(output, "umap_pool.png"), umap_pool, width = 7, height = 5)

saveRDS(b_cells, here(output, "b_cells.RDS"))

# Step 2
data = b_cells
data$cell_type <- rename_cells2(data$cell_type)
data$cell_type <- factor(data$cell_type, levels = c("B Mem", "B IN"))

# SCTransform
data <- SCTransform(data,
                    vars.to.regress = c("pool", "percent.mt"),
                    do.correct.umi = FALSE,
                    variable.features.n = 500)

data <- RunPCA(data)
data <- RunUMAP(data, dims = 1:30)

data <- FindNeighbors(data, dims = 1:30)
data <- FindClusters(data, resolution = 0.1)

saveRDS(data, here(output, "sct_B_cells_500.RDS"))

# Plot
plot_gene(data, gene = "MS4A1", tag = "all")
plot_gene(data, gene = "CD19", tag = "all")
plot_gene(data, gene = "CD34", tag = "all")
plot_gene(data, gene = "CD3D", tag = "all")
plot_gene(data, gene = "CD3E", tag = "all")
plot_gene(data, gene = "GZMB", tag = "all")
plot_gene(data, gene = "CD8A", tag = "all")

plot_meta(data, var = "seurat_clusters", tag = "all")

Idents(data) <- "seurat_clusters"

markers <- FindAllMarkers(data, only.pos = TRUE)

# This reported error
# Error in h(simpleError(msg, call)) :
#  error in evaluating the argument 'x' in selecting a method for function 'slice': Rle of type 'list' is not supported
# Solved by edit the "slice" as "dplyr::slice"

markers %>% 
  group_by(cluster) %>% 
  arrange(desc(avg_log2FC)) %>% 
  select(-p_val, -p_val_adj) %>% 
  dplyr::slice(1:8) %>% 
  as.data.frame()

# Remove contaminated cells expressing CD34 and T cell markers?
subdata <- data[, !data$seurat_clusters %in% c(4,7)]

subdata <- SCTransform(subdata, 
                       vars.to.regress = c("pool", "percent.mt"), 
                       do.correct.umi = FALSE, 
                       variable.features.n = 500)

subdata <- RunPCA(subdata)
subdata <- RunUMAP(subdata, dims = 1:30)

plot_meta(subdata, var = "cell_type", tag = "sub", 
          dim.1 = 2, dim.2 = 1, format = c("png", "pdf")) 

saveRDS(subdata, here(output, "sct_B_cells_subdata.RDS"))
#   ____________________________________________________________________________
#   Run PHATE                                                               ####

em <- GetAssayData(subdata, slot = "scale.data") %>% t()

p <- ElbowPlot(data, ndims = 50)
ggsave(here(output, "elbow.png"), p, width = 9)

inicio("Run PHATE")
ph <- phate(em, npca = 10)
fin()

saveRDS(ph, file = here(output, "phate_500hvg_10pcs.RDS"))


subdata[["phate"]] <- CreateDimReducObject(ph$embedding, 
                                           key = "PHATE_", 
                                           assay = "SCT")

# plot_meta(subdata, var = "cell_type", reduction.use = "phate", 
#          tag = "sub_500hvg_10pcs")

plot_meta(subdata, var = "cell_type", tag = "sub_500hvg_10pcs", reduction.use = "phate",
          dim.1 = 2, dim.2 = 1, format = c("png", "pdf")) 

#   ____________________________________________________________________________
#   Pseudotime analysis                                                     ####

subdata <- FindNeighbors(subdata, dims = 1:10)
subdata <- FindClusters(subdata, resolution = 0.1)


plot_meta(subdata, 
          var = "seurat_clusters", 
          tag = "sub_500hvg_10pcs_res0.1", 
          reduction.use = "phate",
          dim.1 = 2, dim.2 = 1, format = c("png", "pdf"))


em <- Embeddings(subdata, reduction = "phate")

# Save phate results
p=subdata@reductions$phate@cell.embeddings
p=as.data.frame(p)
p$individual=subdata$individual
saveRDS(p,file=here(output,"phate.RDS"))


inicio("Running pseudotime")
sds <- slingshot(em, clusterLabels = subdata$seurat_clusters, 
                 start.clus = 2, end.clus = 1)
fin()
# This step tooki 4.689 mins 
# ★  Elapsed time: 3.542 hrs

saveRDS(sds, file = here(output, "sds_500hvg_10pcs_0.1_fixed.RDS"))

# How many curves do we get?
print(paste0("Number of curves: ", length(slingCurves(sds))))

i <- 1
curve_1 <- slingCurves(sds)[[i]]
pt <- curve_1$lambda %>% as.data.frame() %>% set_names("pt")
subdata <- AddMetaData(subdata, pt)
curve_1 <- curve_1$s[curve_1$ord, 1:2]
colnames(curve_1) <- c("PHATE_1", "PHATE_2")


p <- plot_meta(subdata, 
               var = "pt", 
               reduction.use = "phate", 
               tag = "curve1", 
               dim.1 = 2, 
               dim.2 = 1, 
               save = FALSE) +
  geom_path(aes(PHATE_2, PHATE_1), data = as.data.frame(curve_1),  size = 1) +
  xlab("PHATE 1") +
  ylab("PHATE 2") +
  scale_color_viridis_c(name = "Pseudotime")

ggsave(here(output, "curve1.png"), p, height = 5.5, width = 8.5)

if(length(slingCurves(sds)) == 2){

i <- 2
curve_2 <- slingCurves(sds)[[i]]
pt <- curve_2$lambda %>% as.data.frame() %>% set_names("pt2")
subdata <- AddMetaData(subdata, pt)

curve_2 <- curve_2$s[curve_2$ord, 1:2]
colnames(curve_2) <- c("PHATE_1", "PHATE_2")

p <- plot_meta(subdata, 
               var = "pt2", 
               reduction.use = "phate", 
               tag = "curve1", 
               dim.1 = 2, 
               dim.2 = 1, 
               save = FALSE) +
  geom_path(aes(PHATE_2, PHATE_1), data = as.data.frame(curve_2),  size = 1) +
  xlab("PHATE 1") +
  ylab("PHATE 2") +
  scale_color_viridis_c(name = "Pseudotime")

ggsave(here(output, "curve2.png"), p, height = 5.5, width = 8.5)
}

subdata@misc[["slingshot"]] <- sds

# Define bins
pt <- slingPseudotime(sds, na = FALSE) %>% 
  as.data.frame() %>% 
  rownames_to_column("barcode") %>% 
  select(barcode, curve = curve_1) %>% 
  as_tibble()


md <- subdata[[]] %>% 
  rownames_to_column("barcode") %>% 
  select(barcode, individual, pool, latent, cell_type) %>% 
  as_tibble()

pt <- full_join(md, pt, by = "barcode")


pt$Q4 <- Hmisc::cut2(pt$curve, 
                     g = 4) %>% 
  `levels<-`("Q" %p% seq_len(4))


pt$Q5 <- Hmisc::cut2(pt$curve, 
                     g = 5) %>% 
  `levels<-`("Q" %p% seq_len(5))


pt$Q6 <- Hmisc::cut2(pt$curve, 
                     g = 6) %>% 
  `levels<-`("Q" %p% seq_len(6))


bins <- pt %>% 
  as.data.frame() %>% 
  column_to_rownames("barcode") %>% 
  select(Q4, Q5, Q6)


subdata <- AddMetaData(subdata, bins)

plot_meta(subdata, 
          var = "Q4", 
          reduction.use = "phate", 
          tag = "pt", 
          dim.1 = 2, 
          dim.2 = 1,
          format = c("png", "pdf"))

plot_meta(subdata, 
          var = "Q5", 
          reduction.use = "phate", 
          tag = "pt", 
          dim.1 = 2, 
          dim.2 = 1,
          format = c("png", "pdf"))

plot_meta(subdata, 
          var = "Q6", 
          reduction.use = "phate", 
          tag = "pt", 
          dim.1 = 2, 
          dim.2 = 1,
          format = c("png", "pdf"))

get_stats <- function(q){
  
  subdata[[]] %>% 
    group_by(individual, !!sym(q)) %>% 
    tally() %>% 
    ungroup() %>% 
    filter(n >= 3) %>% 
    group_by(!!sym(q)) %>% 
    tally()
  
}


map(c("Q4", "Q5", "Q6"), get_stats) %>% 
  map(knitr::kable)

plot_gene(subdata, 
          gene = "IGJ", 
          tag = "phate", 
          reduction.use = "phate",
          dim.1 = 2, 
          dim.2 = 1, 
          order = "increasing")

# END
sessionInfo()






####
