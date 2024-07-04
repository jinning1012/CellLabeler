if(F)
{## do not run, only for test code ##
    rm(list = ls())
    
    data(exampledata)
    ## run celllabeler with counts ##
    out = celllabeler(object = counts, markers = markers, sample.id = sample.id, cluster.id =  cluster.id, num.core = 10)
    out$prediction
    out$ude 
    
    ## create celllabeler object and add ude results
    meta.data = data.frame(sample = sample.id, celltype = cluster.id)
    rownames(meta.data) = colnames(counts)
    obj = CreateCellLabelerObject(counts = counts, meta.data = meta.data, ude.result = out)
    obj = AddMarkers(object = obj, markers = markers, cluster.var = "celltype")

    ## run celllabeler with celllabeler object ##
    meta.data = data.frame(sample = sample.id, celltype = cluster.id)
    rownames(meta.data) = colnames(counts)
    obj <- CreateCellLabelerObject(counts = counts, meta.data = meta.data)
    obj <- celllabeler(object = obj, markers = markers, sample.var = "sample", cluster.var = "celltype")
    obj@ude
    obj@prediction

    ## run celllabeler without inputing markers ##
    obj <- CreateCellLabelerObject(counts = counts, meta.data = meta.data)
    obj <- celllabeler(object = obj, sample.var = "sample", cluster.var = "celltype", num.core = 10)
    obj <- AddMarkers(object = obj, markers = markers, cluster.var = "celltype")
}