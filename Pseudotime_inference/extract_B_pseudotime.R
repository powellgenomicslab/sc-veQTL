####
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

sds = readRDS("/share/ScratchGeneral/angxue/proj/vQTL/pseudotime/results/2022-08-31_b_cells_phate/sds_500hvg_10pcs_0.1_fixed.RDS")
i <- 1
curve_1 <- slingCurves(sds)[[i]]
pt <- curve_1$lambda %>% as.data.frame() %>% set_names("pt")
data = readRDS("./results/2022-08-31_b_cells_phate/sct_B_cells_subdata.RDS")

meta = data@meta.data
all(rownames(meta) == rownames(pt))
meta$pseudotime = pt$pt

library(dplyr)
res = list()
both = meta %>%group_by(individual) %>% summarize(nr_cell=n(),mean=mean(pseudotime),var=var(pseudotime))
bin = meta[meta$cell_type=="B IN",] %>%group_by(individual) %>% summarize(nr_cell=n(),mean=mean(pseudotime),var=var(pseudotime))
bmem = meta[meta$cell_type=="B Mem",] %>%group_by(individual) %>% summarize(nr_cell=n(),mean=mean(pseudotime),var=var(pseudotime))

res[[1]] = as.data.frame(both)
res[[2]] = as.data.frame(bin)
res[[3]] = as.data.frame(bmem)

saveRDS(res,"pseudotime_per_ind.RDS")


