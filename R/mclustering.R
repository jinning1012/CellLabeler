#####################################################################
# Package: CellLabeler
# Version: 0.0.1
# Modified: 2024-6-3 08:39; 2024-7-16 16:02
# Title :  
# Authors: Jin Ning and Shiquan Sun
# Contacts: sqsunsph@xjtu.edu.cn;newlife1012@stu.xjtu.edu.cn
#          Xi'an Jiatong University, Center for Single-Cell Omics and Health
######################################################################

#' Detect variable genes
#' 
#' @param counts A raw count gene expression matrix; rows are genes and coumns are cells
#' @param min.pct Numeric value threshold, genes with observed expression percent lower than it would be filtered out
#' 
#' @return A character vectors of variable gene names
#' 
#' @importFrom Matrix rowSums colSums
#' 
#' @export 
#' 
select_hvg = function(counts, min.pct = 0.1){
  if(is.null(rownames(counts))){
    stop("Input count matrix if lack of rownames indicating gene names.")
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



#' Detect a more explicit clustering structure with a given clustering labels
#' 
#' @param counts A raw count gene expression matrix; rows are genes and columns are cells
#' @param cluster.id A character vectors defining the original clustering labels
#' @param k Number of the nearest neighbors in knn; default is 50
#' @param num.hvg Minimum number of HVGs detected for the cluster to be further spliting; default is 10
#' @param num.de Minimum number of DEGs detected for each update cluster; default is 10
#' @param logfc LogFC threshold for detecting DEGs
#' @param min.umi Minimum umi of genes in the filtration step
#' @param verbose Bool indicator for printing messages
#' 
#' @return A list of character vectors defining the updated clustering labels and hvgs
#' 
#' @import igraph
#' @importFrom FNN get.knnx
#' @importFrom coop cosine
#' 
#' @export 
#' 
update_clustering = function(counts, cluster.id, k = 50, num.hvg = 10, num.de = 10, logfc = 0.1, min.umi = 100, verbose = T){
  set.seed(1234)
  ## scale the data ##
  counts = counts[rowSums(counts)>min.umi,]
  if(verbose){cat("Dimension of count matrix:",dim(counts),"\n")}
  scaledata = counts %>% Seurat::NormalizeData(verbose = F) %>% 
      Seurat::ScaleData(verbose=F)

  ## remove cells in cluster with size less than 100
  cluster.type = unique(cluster.id)
  if(verbose) cat("Original cluster type:",paste(cluster.type,collapse = ", "),"\n")
  
  output = c()
  hvg_list = list()
  for(i.type in cluster.type){
    if(verbose) cat("Update cluster",i.type,"...\n")
    
    ## compute cosine similarity
    data = counts[,cluster.id == i.type, drop = FALSE]
    hvgs = select_hvg(data)
    hvg_list[[i.type]] = hvgs

    if(verbose) cat("Detect",length(hvgs),"variable genes \n")
    
    if(length(hvgs) < num.hvg){
      ## return NA for cells and record
      fine_cluster = sapply(colnames(counts)[cluster.id == i.type], function(x){NA})
      output = c(output,fine_cluster)
      next
    }

    data = scaledata[hvgs,cluster.id == i.type, drop = FALSE]
    #cosine_similarity = coop::cosine(data)
    cosine_similarity = cosine(data)
    ## knn
    nn.index = get.knnx(cosine_similarity,cosine_similarity,
                        k = k + 1,
                        algorithm = "kd_tree")$nn.index
    nn.index = lapply(1:nrow(nn.index), function(i) nn.index[i, -1])
    adj_matrix = matrix(0, nrow = nrow(cosine_similarity), ncol = ncol(cosine_similarity))
    for (i in seq_along(nn.index)) {
      adj_matrix[i, nn.index[[i]]] = 1
      adj_matrix[nn.index[[i]], i] = 1 
    }
    rownames(adj_matrix) = colnames(adj_matrix) = colnames(data)
    g = igraph::graph_from_adjacency_matrix(adj_matrix, mode = "undirected")

    ## compute weights by jaccard index
    edge_list = as_edgelist(g)
    weights_jac = igraph::similarity(g,method ="jaccard")
    rownames(weights_jac) = colnames(weights_jac) = colnames(data)
    weights = sapply(1:nrow(edge_list), function(i){weights_jac[edge_list[i,1],edge_list[i,2]]})
    E(g)$weight = weights

    # Perform Louvain clustering
    for(resolution in seq(from = 0.1,to = 1,by = 0.1)){
      cl = cluster_louvain(g, resolution = resolution)
      # Inspect the clustering results
      membership = membership(cl)

      ## judgement ##
      ## 1. cluster proportion
      flag1=length(unique(membership))>1
      cluster_ab = as.numeric(sort(table(membership), decreasing=T))
      flag2 = (as.numeric(cluster_ab[1]/cluster_ab[2]) < 20)

      if(flag1 & flag2){
        ## we test the DE only when cluster proportion fits the request
        ## 2. one versus the rest DE using wilcoxon
        wil_de = wilcoxauc(data, membership)
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
      cat("Reaching the most resolution before fine clustering.\n")
      fine_cluster = sapply(colnames(counts)[cluster.id == i.type], function(x){NA})
    }else{
      fine_cluster = paste0(i.type,"--",membership)
      names(fine_cluster) = colnames(counts)[cluster.id == i.type]
    }
    ## collect for each celltype
    output = c(output, fine_cluster)
  }
  ## order the cells
  output = output[match(colnames(counts), names(output))]

  res = list(cluster = output, hvgs = unique(unlist(hvg_list)))
  return(res)
}


#' Perform multiple rounds of clustering by combining HVGs selected, graph construction and Louvain clustering
#' 
#' @param counts A gene by cell raw count gene expression matrix
#' @param loop.max Integer for maximum iteration rounds; default is 10
#' @param k Number of the nearest neighbors in knn; default is 50
#' @param min.cells Minimum size of the cluster to be further spliting; default is 100
#' @param num.hvg Minimum number of HVGs detected for the cluster to be further spliting; default is 10
#' @param num.de Minimum number of DEGs detected for each update cluster; default is 10
#' @param logfc LogFC threshold for detecting DEGs
#' @param min.umi Minimum umi of genes in the filtration step
#' @param verbose Bool indicator for printing messages
#' 
#' @return A dataframe with rownames being cell names and each column defining clustering for each round
#' 
#' @importFrom presto wilcoxauc
#' @export
#' 
#' 
multiple_clustering = function(counts, loop.max = 10, k = 50, min.cells = 100, num.hvg = 10,num.de = 10,logfc = 0.1, min.umi = 100, verbose = T){
  ## initial a cluster label 
  l0 = rep("0",ncol(counts))
  names(l0) = colnames(counts)

  continue_loop = T
  index_loop = 1
  cluster_list = list()
  while(continue_loop){
    cat("## *******************************************\n")
    cat("## Iterative clustering round",index_loop,"...\n")
    
    res = update_clustering(
        counts = counts,
        cluster.id = l0, 
        k = k, 
        num.hvg = num.hvg, 
        num.de = num.de, 
        logfc = logfc,
        min.umi=min.umi,
        verbose)
    l1 = res$cluster
    hvgs = res$hvgs
    
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

      ## update loop index and counts for next round
      index_loop = index_loop+1
      if(index_loop <= loop.max){
        ## filtered cell with NA label
        counts = counts[,which(!is.na(l1)),drop = F]
        l0 = l1new[!is.na(l1),drop = F]
        
        ## filtered cell with small size
        idx = names(table(l0))[table(l0) < min.cells]
        if(length(idx)>0){
          idx = which(l0 %in% idx)
          counts = counts[hvgs,-idx, drop = F]
          l0 = l0[-idx, drop = F]

          if(length(l0) == 0 | length(hvgs)<10){
            cat("Insufficient genes or cells are available for proceeding to the next clustering round.\n")
            continue_loop = F
          }
        }
      }else{
        cat("Hit the maximum iteration limit.\n")
        continue_loop = F ## stop looping for the next round
      }
    }
  }

  ## collect cluster list to data frame ##
  cluster_df = data.frame(Round1 = paste0("Round1--",cluster_list[[1]]), row.names = names(cluster_list[[1]]))
  for(j in 2:length(cluster_list)){
    lab_col = cluster_list[[j]]
    lab_col = lab_col[!is.na(lab_col)]
    new_col = sapply(rownames(cluster_df), function(x){NA})
    new_col[match(names(lab_col), names(new_col))] = paste0("Round",j,"--",lab_col)
    
    ## fill the NA
    new_col[which(is.na(new_col))] = cluster_df[which(is.na(new_col)), j-1]

    cluster_df[,paste0("Round",j)] = new_col
  }
  return(cluster_df)
}