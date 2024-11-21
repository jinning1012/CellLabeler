# suppressPackageStartupMessages({
#     library(dplyr)
#     library(ggplot2)
#     library(rlang)#species
#     library(dplyr)
#     library(magrittr)
#     library(Matrix)
#     library(Rcpp)
# }) 


#' Detect variable genes with observed expression percent lower than a threshold
#' @param counts A raw count gene expression matrix; rows are genes and coumns are cells
#' @param min.pct Numeric value threshold, genes with observed expression percent lower than it would be filtered out
#' 
#' @return A character vectors of variable gene names
#' 
#' @export 
#' 
SelectHVG = function(counts, min.pct = 0.1){
  if(is.null(rownames(counts))){
    stop("Input count matrix is lack of rownames indicating gene names.")
  }
  counts = counts[rowSums(counts)>min.pct,]
  pse_count = rowSums(counts)/sum(counts)
  total_cnt_cell = colSums(counts)
  ## expected lambda_matrix
  lambda_mat = outer(pse_count, total_cnt_cell, FUN = "*")
  expected_pct = rowSums(1-exp(-lambda_mat))/rowSums(counts)
  observed_pct = apply(counts,1,function(x) sum(x>0)/length(x))
  hvgs = names(expected_pct)[(expected_pct-observed_pct)>0.05]
  return(hvgs)
}


#' Update clustering results
#' @param full.snn A snn network constructed from the begining
#' @param normdata A normalized gene expression matrix with rows being genes and columns being samples
#' @param cluster.id A vector of charactors displaying the different clustering labels; names of it should be cell names
#' @param num.de An integer threshold; cluster with DE less than this will not be subclustered
#' @param logfc A numeric threshold for filtering DE
#' @param verbose Bool indicator to print the progress
#' @return A vector of charactors displaying the updated clustering labels
#' @export 
#' 
#' 
UpdateClustering = function(full.snn, normdata, cluster.id, num.de = 20, logfc = 0.1, verbose = T){
    if(is.null(names(cluster.id))){
        stop("The input cluster.id must be vectors with names being cell names!")
    }

    cluster.type = unique(as.character(cluster.id))
    output = c()
    for(i.type in cluster.type){
        if(verbose) cat("Update cluster",i.type,"...\n")
        sub.cell = names(cluster.id)[cluster.id == i.type]
        res.sub.graph = as.WebGraph(x = full.snn[sub.cell, sub.cell])

        # Perform Louvain clustering
        for(resolution in seq(from = 0.1,to = 1,by = 0.1)){
        res.clust = ClustLouvain(object= res.sub.graph, resolution = resolution, verbose = F)
        membership = as.character(res.clust[,1])
        ## judgement ##
        ## 1. cluster proportion
        flag1=length(unique(membership))>1
        cluster_ab = as.numeric(sort(table(membership), decreasing=T))
        flag2 = (as.numeric(cluster_ab[1]/cluster_ab[2]) < 20)

        if(flag1 & flag2){
            ## we test the DE only when cluster proportion fits the request
            ## 2. one versus the rest DE using wilcoxon
            wil_de = presto::wilcoxauc(normdata[,cluster.id == i.type, drop = F], membership)
            wil_de_filtered = lapply(unique(membership), function(i){
            flag_de = list(
                wil_de$group == i,
                wil_de$padj < 0.05,
                wil_de$logFC >= logfc)
            flag_de = Reduce("&",flag_de) 
            #out = sum(flag_de)>round(0.01*nrow(data)) ## if more than 3 DE were found that means this cluster will be remained
            out = sum(flag_de) >= num.de
            return(out) 
            }) %>% unlist

            flag3 = all(wil_de_filtered)
            if(flag1 & flag2 & flag3){
            if(verbose) cat("Finish at Louvain:",resolution,"\n")
            break;
            
            }
        }#end for fi
        }#end for resolution


        if(resolution == 1 & !(flag1 & flag2 & flag3)){
        cat("Reaching the resolution limit before fine clustering.\n")
        fine_cluster = sapply(colnames(normdata)[cluster.id == i.type], function(x){NA})
        }else{
        fine_cluster = paste0(i.type,"--",membership)
        names(fine_cluster) = colnames(normdata)[cluster.id == i.type]
        }
        ## collect for each celltype
        output = c(output, fine_cluster)
    }
    ## order the cells
    output = output[match(colnames(normdata), names(output))]

    return(output)
}



#' Iterative clustering function
#' @param counts Raw count gene expression matrix
#' @param hvg.method A charater displaying the method for HVG identification; 
#' @param pca Bool indicator, if performing PCA before constructed SNN
#' @param npc An integer defining the number of top PCs for downstream analysis
#' @param loop.max An integer defining the maximum of rounds in iterative cluster
#' @param k Number of neighbors
#' @param min.cells The minimum threshold for clustering size
#' @param min.umi For filtering genes
#' @param num.vst Number of features in HighlyVariableGenes
#' @param num.de An integer threshold; cluster with DE less than this will not be subclustered
#' @param logfc A numeric threshold for filtering DE
#' @param verbose Bool indicator to print the progress
#' 
#' @importFrom Seurat FindVariableFeatures NormalizeData ScaleData
#' @return A dataframe of clustering results; rows are cells and columns are rounds of clustering
#' @export 
#' 
#' 
MultipleClustering = function(counts, hvg.method = c("poisson","vst"), pca = T, npc = 20, loop.max = 10, k = 50, min.cells = 100, min.umi = 100,num.hvg = 100, num.vst = 3000, num.de = 20, logfc = 0.1, verbose = T){
    counts = counts[rowSums(counts)>min.umi, , drop = F]
    ## first round of clustering: data preparation and network construction ##
    hvg.method = hvg.method[1]
    if(hvg.method == "poisson"){
    hvgs = SelectHVG(counts)
  
    }else if(hvg.method == "vst"){
    
    ##vst = Seurat::FindVariableFeatures(counts, selection.method = "vst", verbose = F)
    vst = FindVariableFeatures(counts, selection.method = "vst", verbose = F)
    vst = vst %>% arrange(desc(vst.variance.standardized))
    hvgs = rownames(vst)[1:num.vst]
    }

    if(length(hvgs)<num.hvg){
        stop("Too fewer hvgs to perform clustering! Please check the input count matrix.")
    }
    
    ## data preparation ##
    normdata = NormalizeData(counts[hvgs,, drop = F], verbose = F)
    scaledata = ScaleData(normdata, verbose = F)

    ## perform pca ##
    if(pca){
        npc = min(npc,50)
        if(verbose) cat("Performing PCA to obtain the fist",npc,"PCs.\n")
        pca_results = irlba::irlba(t(scaledata), nv = 50, nu = 50, center = FALSE)
        pca_loadings = t(scaledata) %*% pca_results$v[, 1:npc]
        data_pc = pca_loadings
        colnames(data_pc) = paste0("PC",1:ncol(data_pc))
    }else{
        data_pc = t(scaledata)
    }

    ## build neighbors graph ##
    res.graph = FindNeighborCells(object=data_pc, k.param = k, verbose = F)
    full.snn = res.graph$snn
    
    ## clustering steps ##

    ## initial a cluster label ##
    l0 = rep("0",ncol(counts))
    names(l0) = colnames(counts)

    continue_loop = T
    index_loop = 1
    cluster_list = list()
    while(continue_loop){
        cat("## *******************************************\n")
        cat("## Iterative clustering round",index_loop,"...\n")
        
        l1 = UpdateClustering(full.snn, normdata, l0, num.de = num.de, logfc = logfc, verbose = verbose)

        ## 1.if all NA that means no any sub clustering for each of these cluster.id
        if(all(is.na(l1))){
        continue_loop = F
        cat("No valid categories were obtained from this round of clustering.\n")
        
        }else{

        ## transfer label as "1" "2" "3"...
        unique_values = unique(l1[!is.na(l1)])
        l1new = as.character(match(l1, unique_values))
        names(l1new) = names(l1)
        l1new[is.na(l1)] = NA
        cluster_list[[paste0("Round",index_loop)]] = l1new
        rm(unique_values)
        
        ## update loop index and counts for next round
        index_loop = index_loop+1
        if(index_loop > loop.max){
        cat("Hit the maximum iteration limit.\n")
        continue_loop = F ## stop looping for the next round
        }else{
        
        ## selected cells for next round

        ## filtered cell with NA label
        l0 = l1new[!is.na(l1),drop = F]
        
        ## filtered cell with small size
        idx = names(table(l0))[table(l0) < min.cells]
        if(length(idx)>0){
        idx = which(l0 %in% idx)
        l0 = l0[-idx, drop = F]

        }

        if(length(l0) == 0){
            cat("Insufficient genes or cells are available for proceeding to the next clustering round.\n")
            continue_loop = F
        }

        normdata = normdata[, names(l0)]
        
        } ## next clustering round

        } ## 
    }


    ## if the first round of clustering fails to identify valid clustering
    ## we will return the louvain clustering with resolution 1
    if(length(cluster_list)==0){
        print("First round clustering fails to identify valid clustering, we return the louvain clustering with resolution 1.")
        res.clust = ClustLouvain(object=full.snn, resolution = 1, verbose = F)
        res.label = setNames(as.character(res.clust[,1]), rownames(res.clust))
        cluster_list[[1]]= res.label
    }

    ## collect cluster list to data frame ##
    cluster_df = data.frame(Round1 = paste0("Round1--",cluster_list[[1]]), row.names = names(cluster_list[[1]]))
    print(paste0("Total rounds:",length(cluster_list)))
    if(length(cluster_list)>1){
      for(j in 2:length(cluster_list)){
        lab_col = cluster_list[[j]]
        lab_col = lab_col[!is.na(lab_col)]
        new_col = sapply(rownames(cluster_df), function(x){NA})
        new_col[match(names(lab_col), names(new_col))] = paste0("Round",j,"--",lab_col)
        
        ## fill the NA
        new_col[which(is.na(new_col))] = cluster_df[which(is.na(new_col)), j-1]

        cluster_df[,paste0("Round",j)] = new_col
      }
    }
    return(cluster_df)
}


