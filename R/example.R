if(F)
{## do not run, only for test code ##
    rm(list = ls())
    #library(dplyr)
    #library(rlang)

    source("/home/ningjin/package/CellLabeler_github/R/ude.R")
    source("/home/ningjin/package/CellLabeler_github/R/CellLabeler.R")
    source("/home/ningjin/package/CellLabeler_github/R/celllabeler.object.R")

    ## read a seurat object ##
    if(F){
        library(Seurat)
        seu <- readRDS("/data/public/CellLabeler/OrganData/preprocess_ver2/single_object/cardiac_atrium_tabulasapiens_seurat.RDS")
        seu <- FindVariableFeatures(seu, nfeatures = 2000)
        meta.data <- seu@meta.data
        counts <- LayerData(seu, layer = "counts")
        counts <- counts[VariableFeatures(seu),]
        
        sample.id = as.character(seu@meta.data$sample)
        cluster.id = as.character(seu@meta.data$celltype)
        markers <- readRDS("/data/public/CellLabeler/markerlist_ver2/heart_markers.RDS")

        save(counts,markers, sample.id, cluster.id, file = "/home/ningjin/package/CellLabeler_github/data/exampledata.rda")
    }

    load("/home/ningjin/package/CellLabeler_github/data/exampledata.rda")
    
    ## run celllabeler with counts ##
    out = celllabeler(object = counts, markers = markers, sample.id = sample.id, cluster.id =  cluster.id, num.core = 10)
    out$prediction
    out$ude 
    
    ## run celllabeler with celllabeler object ##
    meta.data = data.frame(sample = sample.id, celltype = cluster.id)
    rownames(meta.data) = colnames(counts)
    obj <- CreateCellLabelerObject(counts = counts, meta.data = meta.data)
    out = celllabeler(object = obj, markers = markers, sample.var = "sample", cluster.var = "celltype")
}