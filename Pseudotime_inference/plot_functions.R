plot_gene <- function(object, gene, tag = "", reduction.use = "umap", ...){
  p <- dittoDimPlot(object, var = gene, reduction.use = reduction.use, 
                    size = 0.5, ...) + 
    scale_color_viridis_c(name = "Exp") +
    theme(panel.grid = element_blank())
  
  ggsave(here(output, tag %p% "_" %p% gene %p% ".png"), p, 
         height = 5.5, width = 7)
  
}

plot_meta <- function(object, 
                      var, 
                      tag = "", 
                      reduction.use = "umap", 
                      dim.1 = 1, 
                      dim.2 = 2,
                      save = TRUE, 
                      format = "png",
                      height = 5.5, 
                      width = 7, 
                      ...){
  
  red_label <- str_to_upper(reduction.use)
  xlab <- paste(red_label, dim.1)
  ylab <- paste(red_label, dim.2)
  
  
  p <- dittoDimPlot(object, 
                    var = var, 
                    reduction.use = reduction.use, 
                    size = 0.5, 
                    main = "", 
                    xlab = xlab,
                    ylab = ylab, 
                    dim.1 = dim.1,
                    dim.2 = dim.2,
                    ...) + 
    theme_pub(panel.grid = element_blank()) 
  
  if(save){
    for(f in format){
      message("Saving ", f, "...")
      ggsave(here(output, tag %p% "_" %p% var %p% "." %p% f), p, 
             height = height, width = width)      
    }
  }else{
    p
  }
}

rename_cells2 <- function(cell_type){  
  recode(cell_type,
         `CD4+ KLRB1+ T cell`      = "CD4 ET",
         `CD4+ KLRB1- T cell`      = "CD4 NC",
         `CD4+ SOX4+ T cell`       = "CD4 SOX4",
         `CD8+ LTB+ T cell`        = "CD8 NC",
         `CD8+ GNLY+ NKG7+ T cell` = "CD8 ET",
         `CD8+ S100B+ T cell`      = "CD8 S100B",
         `XCL1- NK`                = "NK",
         `XCL1+ NK`                = "NK R",
         `TCL1A- FCER2- B cell`    = "B Mem",
         `TCL1A+ FCER2+ B cell`    = "B IN",
         `IgJ+ B cell`             = "Plasma",
         `Monocyte CD14+`          = "Mono C",
         `Monocyte FCGR3A+`        = "Mono NC",
         `Dendritic cell`          = "DC"
  )
  
}


