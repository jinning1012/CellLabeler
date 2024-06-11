if(F)
{## do not run, only for test code ##
    


    ## ******************************** ##
    ##        load data and run         ##
    ## ******************************** ##
    library(CellLabeler)
    load("data/exampledata.rda")
    ## counts: count matrix; rows are genes and columns are cells
    ## sample.id: character vectors; length = #cells
    ## cluster.id: character vectors; length = #cells
    ## markers: list of markers; names of list specifying cell types
    
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
