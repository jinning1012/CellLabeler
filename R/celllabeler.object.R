#####################################################################
# Package: CellLabeler
# Version: 0.0.1
# Modified: 2024-6-11 17:22
# Title :  
# Authors: Jin Ning and Shiquan Sun
# Contacts: sqsunsph@xjtu.edu.cn;newlife1012@stu.xjtu.edu.cn
#          Xi'an Jiatong University, Department of Biostatistics
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
  num.core = "numeric",
  res.path = "ANY"
))


####################################################
#' Create the CellLabeler object with filtering step
#' @param counts Gene expression count matrix (data.frame), p x n -- p is the number of genes and n is the number of cells
#' @param min_total_counts The minimum counts for each cell for filtering
#' @param percentage The percentage of cells that are expressed for analysis
#' @param meta.data Meta data dataframe to store batch and clustering information
#' @return Returns CellLabeler object with filtered gene expression matrix
#' 
#' @export
CreateCellLabelerObject <- function(counts, meta.data, project = "CellLabeler", pct.cells = 0.005, min.cells = 0, min.features = 0, min.umi = 10, num.core = 1, res.path = NULL){
    ## check data order should consistent
	if(!identical(colnames(counts), rownames(meta.data))){
		stop("The column names of counts and row names of meta.data should be should be matched each other! (counts -- g x n; meta.data -- n x c)")
	}# end fi
    
    ## filtering out lowly expressed genes by percentage = 0.005
    if (pct.cells > 0){
        num.cells <- Matrix::rowSums(counts > 0)
        counts	<- counts[which(num.cells >= floor(pct.cells*ncol(counts))),]
    }## end fi

    ## filter genes on the number of cells expressing
    if (min.cells > 0) {
        num.cells <- Matrix::rowSums(x = counts > 0)
        counts <- counts[which(x = num.cells >= min.cells), ]
    }## end fi
    
    ## Filter based on min.features
    if (min.features > 0) {
        nfeatures <- Matrix::colSums(x = counts > 0)
        counts <- counts[, which(x = nfeatures >= min.features)]
        meta.data <- meta.data[which(x = nfeatures >= min.features),]
    }## end fi
    
    if (min.umi > 0) {
        total.umi <- Matrix::colSums(x = counts)
        counts <- counts[, which(x = total.umi > min.umi)]
        meta.data <- meta.data[which(x = total.umi > min.umi),]
    }## end fi
	
    ## inheriting
	object <- new(
		Class = "CellLabeler",
		counts = counts,
        meta.data = meta.data,
        data = NULL,
        scale.data = NULL,
        ude = NULL,
        prediction = NULL,
        ModelFits = NULL,
        num.core = num.core,
        res.path = res.path
	)
  
	return(object)
}# end function


#' Overview of a Celllabeler object
#' @describeIn CellLabeler-methods Overview of a \code{CellLabeler} object
#'
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
    cat('\n')
  })

