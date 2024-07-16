#####################################################################
# Package: CellLabeler
# Version: 0.0.1
# Modified: 2024-6-11 17:22
# Title :  
# Authors: Jin Ning and Shiquan Sun
# Contacts: sqsunsph@xjtu.edu.cn;newlife1012@stu.xjtu.edu.cn
#          Xi'an Jiatong University, Center for Single-Cell Omics and Health
######################################################################

#' Each CellLabeler object has a number of slots which store information. Key slots to access are listed below.
#'
#' @slot counts The raw expression count matrix
#' @slot data Normalized gene matrix
#' @slot meta.data Cell-level metadata
#' @slot ModelFits CellLabeler model fitting results
#' @slot ude CellLabeler detected uniquely differential expression genes
#' @slot prediction CellLabeler prediction dataframe
#' @slot num.core The number of core used in the package
#' 
setClass("CellLabeler", slots=list(
  counts = 'ANY',
  data = 'ANY',
  scale.data = 'ANY',
  meta.data = "data.frame",
  ude = "ANY",
  prediction = "ANY",
  ModelFits = "ANY",
  ModelScores = "ANY",
  num.core = "numeric",
  res.path = "ANY"
))


####################################################
#' Create the CellLabeler object with filtering step
#' @param counts Gene expression count matrix (data.frame), p x n -- p is the number of genes and n is the number of cells
#' @param meta.data Meta data dataframe to store batch and clustering information
#' @param project Character name specifying project name, default is "CellLabeler"
#' @param ude.result A list of results using celllabeler() on counts data
#' @param pct.cells A numeric value indicating at least what proportion of cells the genes expressed on
#' @param min.cells A numeric value indicating at least what number of cells the genes expressed on
#' @param min.features A numeric value indicating at least what number of genes the cells expressed
#' @param min.umi A numeric value indicating at least what total counts of genes the cells expressed
#' @param num.core Number of cores in parallel running
#' @param res.path Path to save the results
#' @return A CellLabeler object with filtered gene expression matrix
#' 
#' @importFrom Matrix rowSums colSums
#' 
#' @examples
#' 
#' data(exampledata)
#' meta.data = data.frame(sample = sample.id, celltype = cluster.id)
#' rownames(meta.data) = colnames(counts)
#' obj = CreateCellLabelerObject(counts = counts, meta.data = meta.data)
#' 
#' @export
#' 
CreateCellLabelerObject = function(counts, meta.data, ude.result = NULL, project = "CellLabeler", pct.cells = 0.005, min.cells = 0, min.features = 0, min.umi = 10, num.core = 1, res.path = NULL){
    ## check data order should consistent
	if(!identical(colnames(counts), rownames(meta.data))){
		stop("The column names of counts and row names of meta.data should be should be matched each other! (counts -- g x n; meta.data -- n x c)")
	}# end fi
    
    ## filtering out lowly expressed genes by percentage = 0.005
    if (pct.cells > 0){
        num.cells = rowSums(counts > 0)
        counts	= counts[which(num.cells >= floor(pct.cells*ncol(counts))),]
    }## end fi

    ## filter genes on the number of cells expressing
    if (min.cells > 0) {
        num.cells = owSums(x = counts > 0)
        counts = counts[which(x = num.cells >= min.cells), ]
    }## end fi
    
    ## Filter based on min.features
    if (min.features > 0) {
        nfeatures = colSums(x = counts > 0)
        counts = counts[, which(x = nfeatures >= min.features)]
        meta.data = meta.data[which(x = nfeatures >= min.features),]
    }## end fi
    
    if (min.umi > 0) {
        total.umi = colSums(x = counts)
        counts = counts[, which(x = total.umi > min.umi)]
        meta.data = meta.data[which(x = total.umi > min.umi),]
    }## end fi
	
    ## inheriting
	object = new(
		Class = "CellLabeler",
		counts = counts,
        meta.data = meta.data,
        data = NULL,
        scale.data = NULL,
        ude = NULL,
        prediction = NULL,
        ModelFits = NULL,
        ModelScores = NULL,
        num.core = num.core,
        res.path = res.path
	)

    ## add ude.results in it
    if(!is.null(ude.result)){
        names_input = names(ude.result)
        if(!all(c("ude","prediction","ModelFits","ModelScores") %in% names_input)){
            stop("## Input UDE results are not correct.")
        }

        object@ude = ude.result$ude
        object@prediction = ude.result$prediction
        object@ModelFits = ude.result$ModelFits
        object@ModelScores = ude.result$ModelScores
    }

	return(object)
}# end function


#' Overview of a Celllabeler object
#' @describeIn CellLabeler-methods Overview of a \code{CellLabeler} object
#' @param object A celllabeler object
#' @return \code{show}: Prints summary to \code{\link[base]{stdout}} and invisibly returns \code{NULL}
#'
#' @importFrom methods show
#'
#' @export
#'
setMethod(f = "show",
  signature = "CellLabeler",
  definition = function(object) {
    cat("## An object of class CellLabeler","\n")
    cat("## Gene number:", nrow(object@counts),"\n")
    cat("## Cell number:", ncol(object@counts),"\n")
    cat("## Meta columns:",paste(colnames(object@meta.data), collapse = ", "), "\n")
    cat('\n')
  })


#' Add the prediction results to meta.data in CellLabelerObject for plot use
#' 
#' @param object A cellLabeler object with predictions
#' @param cluster.var Column name in meta data of the object indicating cluster information
#' @param prediction.var A new defined column name in meta data of the object indicating cluster information
#' 
#' @return A cellLabeler object
#' 
#' @importFrom dplyr left_join
#' 
#' @export
#'  
AddPredictionToMeta = function(object, cluster.var, prediction.var = "prediction"){
    if(!inherits(object,"CellLabeler")){
        stop("Input object should be a CellLabeler object.")
    }

    df1=data.frame(cell = colnames(object@counts), cluster = object@meta.data[,cluster.var], row.names = colnames(object@counts))
    df2=object@prediction
    ## split the prediction if there is any combination
    if(any(grepl(" & ", df2$cluster))){
        df3 = df2[-grep(" & ", df2$cluster),]
        comb_idx = grep(" & ", df2$cluster)
        for(id in comb_idx){
            tmp_df = data.frame(cluster = strsplit(df2[id,"cluster"]," & ")[[1]],
                                prediction = df2[id,"prediction"],
                                score = df2[id,"score"])
            df3 = rbind.data.frame(df3, tmp_df)
        }
    }else{
        df3 = df2
    }
    df = left_join(df1,df3,"cluster")
    object@meta.data[,prediction.var] = df$prediction
    return(object)
}
