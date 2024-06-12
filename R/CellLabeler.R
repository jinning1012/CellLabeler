#####################################################################
# Package: CellLabeler
# Version: 0.0.1
# Modified: 2024-6-3 08:39
# Title :  
# Authors: Jin Ning and Shiquan Sun
# Contacts: sqsunsph@xjtu.edu.cn;newlife1012@stu.xjtu.edu.cn
#          Xi'an Jiatong University, Department of Biostatistics
######################################################################

#'
#' Fitting the constrained spline model to perform temporal differential expression analysis for time course scRNA-seq data
#' 
#' @param object Gene expression raw count matrix
#' @param sample.id Column names of sample information in metadata. 
#' @param cluster.id Column names of cluster information in metadata. 
#' @param similar.gene Integer value as the number of top expressed genes in considering combining clusters before DE test
#' @param similar.pct Numeric value in [0,1] as the threshold to combine clusters before DE test when either two cluster have enough overlapped highly expressed genes. Default is 0.8.
#' @param lfc Limit testing to genes which show the maximum on average X-fold difference (log-scale) between the interested clusters and the rest others. Default is 0.5.
#' @param min.ccells Minimum number of cells in each cluster required. Default is 10.
#' @param max.genes Maximum number of highly expressed genes in each cluster required. Default is 100.
#' @param mod Character as "combine" or "all", specifying if combining similar clusters or not before DE test. Default is "combine"
#' @param up.thr Numberic value in [0,1] as the proportion of the remaining clusters in which a gene is upregulated in the interested group compared with it. Highter, the more strict. Default 0.9.
#' @param num.core Number of cores in multi-threaded running. Default is 1.
#' 
#' @examples
#' 
#' data(exampledata)
#' res <- celllabeler(object=counts, sample.id = sample.id, cluster.id = cluster.id, markers = markers, num.core = 10)
#' 
#' @export
#' 
celllabeler.default <- function(object, 
                                    markers = NULL,
                                    sample.id = NULL,
                  					cluster.id = NULL,
                  					similar.gene = 20,
                  					similar.pct = 0.8,
                  					lfc = 0.5,
                  					min.ccells = 10,
                                    max.genes = 100,
									mod = "combine",
                                    up.thr = 0.9,
                  					num.core = 1, 
                  					verbose = FALSE) {
    suppressPackageStartupMessages({require("dplyr")})
    ## load in counts ##
    counts = object
    data = Seurat::NormalizeData(counts, verbose = verbose)
    if(mod == "combine")
	{	   
    ##*************************************************##
	##   Combining similar clusters before DE test     ##
	##*************************************************##
    cat("## Combining similar clusters before DE test.\n")
	## compute log2FC for each clusters ##
    clusters = unique(cluster.id)
    TopGenes_bylogFC = parallel::mclapply(clusters, mc.cores = num.core, function(ident.1)
    {
        idx1 = which(cluster.id == ident.1)
        idx2 = which(cluster.id != ident.1)
        avg.ident.1 = log(x = rowMeans(x = exp(matrix(data[,idx1], ncol = length(idx1)))-1) + 1, base = 2) 
        avg.ident.2 = log(x = rowMeans(x = exp(matrix(data[,idx2], ncol = length(idx2)))-1) + 1, base = 2)
        log_foldchange = avg.ident.1 - avg.ident.2
        names(log_foldchange) = rownames(data)
        log_foldchange = log_foldchange[log_foldchange>0]
        output = sort(log_foldchange, decreasing = T)[seq(min(similar.gene, length(log_foldchange)))] %>% names
        return(output)
    })
    names(TopGenes_bylogFC) = clusters

    ## combination idx of either two clusters ##
    comb_mat = combn(clusters,2)
    comb_overlap = lapply(seq(ncol(comb_mat)), function(id) TopGenes_bylogFC[comb_mat[,id]] %>% Reduce(intersect,.) %>% length) %>% unlist
    idx = which(comb_overlap > similar.gene*similar.pct)
    
    adj_matrix <- matrix(0, nrow = length(TopGenes_bylogFC), ncol = length(TopGenes_bylogFC))
    colnames(adj_matrix) = rownames(adj_matrix) = names(TopGenes_bylogFC)
    for(i in idx)
    {
        adj_matrix[comb_mat[1,i], comb_mat[2,i]] <- 1
    }
    adj_matrix <- adj_matrix+t(adj_matrix)-diag(diag(adj_matrix))
    
    suppressPackageStartupMessages({library(igraph)})
    # Convert adjacency matrix to graph #
    g = graph_from_adjacency_matrix(adj_matrix, mode = "undirected")
    components = components(g)
    # membership is a character vector while names being our cluster names and values being the splited group
    membership = components$membership
    ## for each loop, we judge if there still remain more to be combined ##
    cluster.id.ori = cluster.id
    if(max(table(membership)) > 1) 
    {
        message(paste0("## Clusters have been combined before DE test."))
        tobeCombined = names(table(membership))[table(membership) > 1]
        for(id in tobeCombined)
        {
        grouplabel = names(membership)[membership == id]
        cluster.id[cluster.id %in% grouplabel]<- paste(grouplabel,collapse = " & ")
        }
    } ## end combination in graph
    }else{
        cat("## Do not check and combine clusters before DE test.\n")
    } ## end mod fi
    rm(clusters,TopGenes_bylogFC,comb_mat, comb_overlap, adj_matrix, g, components, membership)
    

    ##********************************************************##
	##   Removing clusters with few cells after combination   ##
	##********************************************************##
    clusters = unique(cluster.id)
	clusters.size = table(cluster.id)
    clusters.filter = names(clusters.size)[clusters.size < min.ccells]
    if(length(clusters.filter) > 0)
    {
        for(mm in clusters.filter){
        message(paste0("## Remove ",mm,": #cells =",clusters.size[mm]))
        }
        ## subset ##
        idx.remain = which(cluster.id %in% setdiff(clusters, clusters.filter))
        counts = counts[,idx.remain]
        data = data[,idx.remain]
        sample.id = sample.id[idx.remain]
        cluster.id = cluster.id[idx.remain]
        rm(idx.remain)
    } ## end filter
	
    ##***********************************************************##
	##   Selecting top highly expressed genes for each cluster   ##
	##***********************************************************##
    if(max.genes != Inf)
    {
        clusters = unique(cluster.id) ## update clusters
        TopGenes_bylogFC = parallel::mclapply(clusters, mc.cores = num.core, function(ident.1)
        {
            idx1 = which(cluster.id == ident.1)
            idx2 = which(cluster.id != ident.1)
            avg.ident.1 = log(x = rowMeans(x = exp(matrix(data[,idx1], ncol = length(idx1)))-1) + 1, base = 2) 
            avg.ident.2 = log(x = rowMeans(x = exp(matrix(data[,idx2], ncol = length(idx2)))-1) + 1, base = 2)
            log_foldchange = avg.ident.1 - avg.ident.2
            names(log_foldchange) = rownames(data)
            log_foldchange = log_foldchange[log_foldchange>0]
            output = sort(log_foldchange, decreasing = T)[seq(min(max.genes, length(log_foldchange)))] %>% names
            return(output)
        })
        names(TopGenes_bylogFC) = clusters

        gene.use = unique(unlist(TopGenes_bylogFC))
        counts = counts[gene.use,]
        data = data[gene.use,]
        rm(TopGenes_bylogFC)
        
    } ## end max gene fi

	##########################################################
	cat(paste("## ===== CellLabeler INPUT INFORMATION ====## \n"))
	cat(paste("## number of total cells: ", ncol(data),"\n"))
	cat(paste("## number of total features: ", nrow(data),"\n"))
    cat(paste("## number of cell clusters: ", length(unique(cluster.id)),"\n"))
	cat(paste("## number of cores: ", num.core,"\n"))
	cat(paste("## ========== END INFORMATION ============## \n"))
	cat("\n")
	
	
	##***********************************************************##
	##                      main function                        ##
	##***********************************************************##
    
    ude = FindAllUniqueMarkers(data,sample.id,cluster.id,up.thr,verbose,num.core)
    pct = ComputePCT(counts,cluster.id,allGenes(ude))
    if(is.null(markers)){
        ## no prediction but still computing gene score
        clusters = names(ude)
        deg = lapply(clusters,function(i.celltype) ComputeGeneScore(ude,i.celltype,pct) %>% names)
        names(deg) = clusters
        pred = list(deg = deg)
    }else{
        pred = ComputePrediction(ude,markers,pct)
    }
    

    out = list(ude=pred$deg,prediction = pred$predict, score = pred$scores, ModelFits=ude)
	return(out)
}## end function 



#' CellLabeler: 
#' @param object A CellLabeler object 
#' 
#'
#' @param features Feature names. Default is NULL, then all features will be used.
#' @param sample.var Column names of sample information in metadata. 
#' @param cluster.var Column names of cluster information in metadata. 
#' @param similar.gene Integer value as the number of top expressed genes in considering combining clusters before DE test
#' @param similar.pct Numeric value in [0,1] as the threshold to combine clusters before DE test when either two cluster have enough overlapped highly expressed genes. Default is 0.8.
#' @param lfc Limit testing to genes which show the maximum on average X-fold difference (log-scale) between the interested clusters and the rest others. Default is 0.5.
#' @param min.ccells Minimum number of cells in each cluster required. Default is 10.
#' @param max.genes Maximum number of highly expressed genes in each cluster required. Default is 100.
#' @param mod Character as "combine" or "all", specifying if combining similar clusters or not before DE test. Default is "combine"
#' @param up.thr Numberic value in [0,1] as the proportion of the remaining clusters in which a gene is upregulated in the interested group compared with it. Highter, the more strict. Default 0.9.
#' @param num.core Number of cores in multi-threaded running. Default is 1.
#' 
#' @return object A CellLabeler object
#'
#' @author Jin Ning
#' 
#' @export
#' 
celllabeler.CellLabeler <- function(object, 
                                    features = NULL,
                                    markers = NULL,
                                    sample.var = "sample",
                                    cluster.var = "cluster",
                                    similar.gene = 20,
                                    similar.pct = 0.8,
                                    lfc = 0.5,
                                    max.genes = 100,
                                    min.ccells = 20,
                                    up.thr = 0.9,
                                    mod='combine',
                                    num.core = 1,
                                    verbose = TRUE, ...) {
	
	## parallel parameter setting
	if(num.core == 1){
		if(object@num.core > 1){
            num.core = object@num.core
        }
	}## end fi
	
	## counts
    counts = object@counts
    if(!is.null(features)){
        cat("## Run with specified features ...\n")
        counts = counts[features,]
    }
    ## meta.data
    sample.id = as.character(object@meta.data[,sample.var])
    cluster.id = as.character(object@meta.data[,cluster.var])
	## run main
	results <- celllabeler(object = counts,
                                markers = markers,
                                sample.id = sample.id,
                                cluster.id = cluster.id,
                                similar.gene = similar.gene,
                                similar.pct = similar.pct,
                                lfc = lfc,
                                max.genes = max.genes,
                                min.ccells = min.ccells,
                                up.thr = up.thr,
                                mod=mod,
                                num.core = num.core,
                                verbose = verbose, ... )

	## store back the treated data
	object@ude = results$ude
    object@prediction = results$prediction
    object@ModelFits = results$ModelFits
	return(object)
}## end func


#' Unique marker gene detection and automatic cell type annotation
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname celllabeler
#' @export celllabeler
#'
#' @concept data-access
#'
celllabeler <- function(object, ...) {
	UseMethod(generic = "celllabeler", object = object)
}## end func


##############################################################################################


#########################################
#             CODE END                  #
#########################################