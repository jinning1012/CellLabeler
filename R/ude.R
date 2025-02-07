
#' Permutation of two of all categories
#' @param vec A vector of elements where my permutation performs
#' 
#' @return Returns a list of permutations of two elements
#' 
#' 
perm2 = function(vec){
  out = list()
  flag = 1
  for(i in vec){
    for(j in setdiff(vec,i)){
      out[[flag]] = c(i,j)
      flag = flag+1
    }
  }
  out = do.call(cbind,out)
  return(out)
}

#' Cauchy combination rules
#' @param pvalues A n x m pvalue matrix; n genes, m permutations
#' @param weights A n x m weighs matrix for each genes in each permutation
#' 
#' @export 
#' 
#' @return Returns a n-length combined pvalues 
#' 
CombinePValues = function(pvalues, weights=NULL){
  if(!is.matrix(pvalues)){pvalues = as.matrix(pvalues)}
  ## to avoid extremely values
  pvalues[which(pvalues==0)] = 5.55e-17
  pvalues[which((1-pvalues)<1e-3)] = 0.99
  
  num_pval = ncol(pvalues)
  num_gene = nrow(pvalues)
  if(is.null(weights)){
    weights = matrix(rep(1.0/num_pval, num_pval*num_gene), ncol=num_pval )
  }# end fi
  if( (nrow(weights) != num_gene) || (ncol(weights) != num_pval)){
    stop("the dimensions of weights does not match that of combined pvalues")
  }# end fi
  
  Cstat = tan((0.5 - pvalues)*pi)
  
  wCstat = weights*Cstat
  Cbar = apply(wCstat, 1, sum)
  #combined_pval = 1.0/2.0 - atan(Cbar)/pi
  combined_pval = 1.0 - pcauchy(Cbar)	
  combined_pval[which(combined_pval <= 0)] = 5.55e-17
  return(combined_pval)
}# end func




#' Regress out the sample batch effects
#' @param data.expr An expression matrix to regress the effects of latent.data out
#' of should be the complete expression matrix in genes x cells
#' @param latent.data A matrix or data.frame of latent variables, should be cells
#' x latent variables, the colnames should be the variables to regress
#' @param features.regress An integer vector representing the indices of the
#' genes to run regression on
#' @param model.use Model to use, one of 'linear', 'poisson', or 'negbinom'; pass
#' NULL to simply return data.expr
#' @param use.umi Regress on UMI count data
#' @param verbose Display a progress bar
#' 
#' @importFrom stats as.formula glm lm p.adjust pcauchy predict residuals
#' 
#' @export 
#' 
RegressOutMatrix = function(data.expr,
                             latent.data = NULL,
                             features.regress = NULL,
                             model.use = NULL,
                             use.umi = FALSE,
                             verbose = FALSE) {
  
  ## Do we bypass regression and simply return data.expr?
  bypass = vapply(X = list(latent.data, model.use), FUN = is.null, FUN.VALUE = logical(length = 1L) )
  
  if (any(bypass)) {
    return(data.expr)
  }## end fi
  
  ## Check model.use
  possible.models = c("linear", "poisson", "negbinom")
  if (!model.use %in% possible.models) {
    stop(paste(
      model.use,
      "is not a valid model. Please use one the following:",
      paste0(possible.models, collapse = ", ")
    ))
  }## end fi
  
  ## Check features.regress
  if (is.null(x = features.regress)) {
    features.regress = 1:nrow(x = data.expr)
  }## end fi
  
  if (is.character(x = features.regress)) {
    features.regress = intersect(x = features.regress, y = rownames(x = data.expr))
    if (length(x = features.regress) == 0) {
      stop("Cannot use features that are beyond the scope of data.expr")
    }## end fi
  } else if (max(features.regress) > nrow(x = data.expr)) {
    stop("Cannot use features that are beyond the scope of data.expr")
  }## end fi
  
  ## Check data set dimensions
  if (nrow(x = latent.data) != ncol(x = data.expr)) {
    stop("Uneven number of cells between latent data and expression data")
  }## end fi
  
  use.umi = ifelse(test = model.use != 'linear', yes = TRUE, no = use.umi)
  
  ## Create formula for regression
  vars.to.regress = colnames(x = latent.data)
  fmla = paste('GENE ~', paste(vars.to.regress, collapse = '+'))
  fmla = as.formula(object = fmla)
  if (model.use == "linear") {
    # In this code, we'll repeatedly regress different Y against the same X
    # (latent.data) in order to calculate residuals.  Rather that repeatedly
    # call lm to do this, we'll avoid recalculating the QR decomposition for the
    # latent.data matrix each time by reusing it after calculating it once
    regression.mat = cbind(latent.data, data.expr[1,])
    colnames(regression.mat) = c(colnames(x = latent.data), "GENE")
    #qr = lm(fmla, data = regression.mat, qr = TRUE)$qr
    modelfit = lm(fmla, data = regression.mat, qr = T)
    qr = modelfit$qr
    rm(regression.mat)
  }## end fi
  
  ## Make results matrix
  data.resid = matrix(nrow = nrow(x = data.expr), ncol = ncol(x = data.expr) )
  
  
  for (i in 1:length(x = features.regress)) {
    x = features.regress[i]
    regression.mat = cbind(latent.data, data.expr[x, ])
    colnames(x = regression.mat) = c(vars.to.regress, 'GENE')
    regression.mat = switch(
      EXPR = model.use,
      'linear' = qr.resid(qr = qr, y = data.expr[x,]),
      'poisson' = residuals(object = glm(
        formula = fmla,
        family = 'poisson',
        data = regression.mat),
        type = 'pearson'))
    
    data.resid[i, ] = regression.mat
    

  }## end for
  
  ## start if
  if (use.umi) {
    data.resid = log1p(x = Sweep(x = data.resid, MARGIN = 1, STATS = apply(X = data.resid, MARGIN = 1, FUN = min), FUN = '-' ))
  }##end fi
  
  dimnames(x = data.resid) = dimnames(x = data.expr)
  return(data.resid)
}## end func



#' Fitting the firth logistic regression model to perform unique marker gene analysis for clustered scRNA-seq data
#' 
#' @param data.use Log-normalized gene expression matrix; rows are genes and columns as cells
#' @param sample.id Character vector of sample information in metadata. 
#' @param cluster.id Character vector of cluster information in metadata. 
#' @param upPercent Numberic value in [0,1] as the proportion of the remaining clusters in which a gene is upregulated in the interested group compared with it. Highter, the more strict. Default 0.9.
#' @param num.core Number of cores in multi-threaded running. Default is 1.
#' @param verbose Bool variable to indicate if print the message.
#' 
#' @importFrom utils combn 
#' @import logistf
#' @import foreach
#' @import doParallel
#' @export
#' 
#' 
FindAllUniqueMarkers = function(data.use,
                                 sample.id,
                                 cluster.id,
                                 upPercent = 0.9,
                                 verbose = T,
                                 num.core = 1){
  
  
  ## ******************************************* ##
  ##       function for one comparison           ##
  ## ******************************************* ##
  clusters = unique(cluster.id)
  runOneCompare = function(ident.1, ident.2){
    cells.1 = which(cluster.id == ident.1)
    cells.2 = which(cluster.id == ident.2)
    sub_data.use = data.use[,c(cells.1, cells.2), drop = F]
    sub_group = cluster.id[c(cells.1, cells.2)]
    sub_sample = sample.id[c(cells.1, cells.2)]
    outcome = as.numeric(sub_group == ident.1)
    
    # Step1: decide to remove batch effects
    # 1. escape when there is no batch information
    # 2. escape when each ident contains only one batch
    # modified by ningjin; 2023-9-14 20:06
    if(length(unique(sub_sample))==1){
      resi = sub_data.use
      resi = resi[Matrix::rowSums(resi)!=0, , drop = F] ## remove some genes
    }else{
      latent.data = data.frame(sample = sub_sample)
      rownames(latent.data) = colnames(sub_data.use)
      resi = RegressOutMatrix(data.expr = sub_data.use, 
                                latent.data = latent.data,
                                model.use = 'linear',
                                features.regress = NULL,
                                verbose = F)
      resi = resi[Matrix::rowSums(resi)!=0, ,drop = F]
    }## end batch regression

    ## ********************************** ##
    ##           firth's model            ##
    ## ********************************** ##
    ## modified by ningjin; we change the loop to foreach to catch error if happened

    res_firth = foreach(k = 1:nrow(resi),.errorhandling='pass') %dopar% {
      res_tmp = tryCatch({
        filr = logistf(outcome ~ resi[k,])
        res_tmp = c(beta=filr$coefficients[2], pvalue=filr$prob[2])
      },
      warning = function(war) {return(paste0('logistf warining: escape the gene ',rownames(resi)[k]))},
      error = function(err) {return(paste0('logistf error: escape the gene ',rownames(resi)[k]))},
      finally = {
        # return('other things')
      })
      return(res_tmp)
    }

    ## collect results 
    output_idx = which(unlist(lapply(res_firth, length)) == 2)
    res_firth = res_firth[output_idx] %>% do.call(rbind,.)
    rownames(res_firth) = rownames(resi)[output_idx]
    colnames(res_firth) = c('beta','pvalue')

    return(res_firth)
  }## end function
  
  ## ******************************************************
  ##           run for all group comparison
  ## ******************************************************
  all_res_list = list() 
  ## combination of clusters
  ident.comb = combn(clusters,2)
  ## start loop for comparisons
  all_res_list = parallel::mclapply(1:ncol(ident.comb), mc.cores = num.core, function(i.comb)
  { 
    ident.1 = ident.comb[1,i.comb]
    ident.2 = ident.comb[2,i.comb]
    ## ningjin 20240701
    if(verbose){
      cat(paste0("### Logistf on ",ident.1, " versus ",ident.2,"...\n"))
    }
    output = runOneCompare(ident.1, ident.2) 
    return(output)
  })
  
  res_name = lapply(1:ncol(ident.comb), function(x) paste0(ident.comb[1,x],'_versus_',ident.comb[2,x])) %>% unlist
  names(all_res_list) = res_name
  rm(res_name)
  
  ## for the rest of comparisons
  ident.perm = perm2(clusters)
  for(i in 1:ncol(ident.perm)){
    res_name = paste0(ident.perm[1,i],'_versus_',ident.perm[2,i])
    if(!(res_name %in% names(all_res_list))){
      res_tmp = all_res_list[[paste0(ident.perm[2,i],'_versus_',ident.perm[1,i])]]
      ## inverse the effect
      ## the first column is beta from firth model
      res_tmp[,1] = -as.vector(res_tmp[,1])
      all_res_list[[res_name]] = res_tmp
      rm(res_name)
    }
  }

  ## collect all firth model results
  output = list() # final output 
  for(ident.one in clusters){
    res_name = setdiff(clusters, ident.one)
    res_relative = all_res_list[paste0(ident.one,'_versus_',res_name)]
    ## cauchy combination rule ##
    pvs = lapply(res_relative, function(x) return(x[, 2]))
    shared_genes = lapply(pvs, names) %>% Reduce(intersect,.)
    pvs = lapply(pvs, function(x) x[shared_genes] ) %>% do.call(cbind,.)
    if(ncol(pvs) == 1){
      combined.pv = as.vector(pvs)
    }else{
      combined.pv = CombinePValues(pvs)}
    
    avg_beta = lapply(res_relative, function(x){x[shared_genes,1]}) %>% do.call(cbind,.) %>% rowMeans(.)
    
    de_res = lapply(res_relative, function(x){matrix(x[shared_genes,], nrow = length(shared_genes))}) %>% do.call(cbind,.) %>% as.data.frame()
    colnames(de_res) = lapply(res_name, function(j){paste0(c('beta.','pv.'),j)}) %>% unlist
    rownames(de_res) = shared_genes
    
    de_res[,'gene'] = shared_genes
    de_res[,'avg_beta'] = avg_beta
    de_res[,'combined.pv'] = combined.pv
    de_res[,'p.adjust'] = p.adjust(combined.pv, method = 'BH')
    
    ## check if upregulated ##
    #beta_mat = de_res[,grep('beta.',colnames(de_res))]
    #pv_mat = de_res[,grep('pv.',colnames(de_res))]
    beta_mat = de_res[,paste0("beta.",res_name), drop = F]
    pv_mat = de_res[,paste0("pv.",res_name), drop = F]
    ## modified by ningjin 2023-9-4
    ## we want it upregulated to some extent across groups
    changes = ifelse(apply(beta_mat, 1, function(x) length(which(x>0))/length(x)) >= upPercent , 'upregulated','not upregulated')
    ## ningjin 2024-09-25 10:42
    changes[de_res[,"p.adjust"] >0.05] = "not upregulated"
    de_res$change = changes
    output[[ident.one]] = de_res
  }
  
  ## finally output as list of model fits results
  return(output)
}



#'
#' Compute gene expression percent of single-cell RNA data.
#' 
#' @param counts Gene expression raw count matrix.
#' @param cluster.id Character vector of cluster information in metadata. 
#' @param sub.genes Subset of gene names. Default is NULL, then all genes will be used.
#' 
#' 
#' @export
#' 
ComputePCT = function(counts, cluster.id = NULL, sub.genes = NULL){
  if(!is.null(sub.genes)){
    ## check if sub.genes exceed provided genes
    if(!all(sub.genes %in% rownames(counts))){
      stop(paste0("## Given genes in ComputePCT() exceed the range of all genes from counts data."))
    }
    counts = counts[sub.genes,]
  }
  ngene = nrow(counts)
  ## if clusterid is blank, it means compute gene expression pct across all cells
  if(is.null(cluster.id)){
    out = apply(counts, 1, function(gene){
      length(which(gene > 0))/length(gene)
    })
  }else{
    groups = unique(cluster.id)
    out = matrix(nrow = ngene, ncol = length(groups), dimnames = list(rownames(counts), groups))
    for (ig in seq(length(groups))) {
      group = groups[ig]
      cells.1 = which(cluster.id == group)
      #sub_cnt = matrix(counts[,cells.1], ncol = length(cells.1))
      sub_cnt = counts[,cells.1, drop = F]
      out[,ig] = apply(sub_cnt, 1, function(gene){
        length(which(gene > 0))/length(gene)
      })
    }
  }
  return(out)
}

#'
#' Compute gene score based on combined pvalue and gene expression percents.
#' 
#' @param de_model Results of FindAllUniqueMarkers().
#' @param i.celltype Column names specifying the interested cluster.
#' @param PCT_mat Gene expression percent matrix. Rows are genes and columns are clusters.
#' 
#' @export
#' 
ComputeGeneScore = function(de_model, i.celltype, PCT_mat){
  PCT_mat[which(PCT_mat == 0,arr.ind = T)] = 1e-10
  
  alltypes = names(de_model)
  degs = rownames(de_model[[i.celltype]])[de_model[[i.celltype]]$change == 'upregulated' & de_model[[i.celltype]]$p.adjust < 0.05]
  
  ## Cstats ##
  pvalues = as.matrix(de_model[[i.celltype]][degs,grep('pv.',de_model[[i.celltype]] %>% colnames)])
  pvalues[which(pvalues==0)] = 5.55e-17
  pvalues[which((1-pvalues)<1e-3)] = 0.99
  
  num_pval = ncol(pvalues)
  num_gene = nrow(pvalues)
  weights = matrix(rep(1.0/num_pval, num_pval*num_gene), ncol=num_pval )
  Cstat = tan((0.5 - pvalues)*pi)
  
  ## gene by npval
  wCstat = weights*Cstat
  
  ## pct
  Cbar = matrix(nrow = num_gene, ncol = num_pval, dimnames = list(degs,NULL))
  
  for (j in seq(num_pval)) {
    pct = PCT_mat[degs,i.celltype]/PCT_mat[degs, setdiff(alltypes, i.celltype)[j]]
    Cbar[,j] = wCstat[,j] * pct
  }
  Cbar = Matrix::rowMeans(Cbar)
  output = sort(Cbar, decreasing = T)
  output = output[output > 0]
  return(output)
}


#'
#' Compute gene score and perform predictions.
#' 
#' @param de_model Results of FindAllUniqueMarkers().
#' @param marker_list List of makers, with list names being cell type names.
#' @param PCT_mat Gene expression percent matrix. Rows are genes and columns are clusters.
#' 
#' @importFrom stringr str_to_title
#' 
#' @export
#' 
ComputePrediction = function(de_model, marker_list, PCT_mat){
  alltypes = names(de_model)
  
  # check if the genes from real data is uppercase
  if(rownames(de_model[[1]])[1] == toupper(rownames(de_model[[1]])[1])){
    ## to upper
    marker_list = lapply(marker_list, toupper)
  }else{
    ## to captical first letter
    marker_list = lapply(marker_list, str_to_title)
  }
  
  # init
  score_list = list()
  de_list = list()
  predict_df = data.frame()
  modelmarker_df = data.frame()
  for (i.celltype in alltypes) {
    Cbar = ComputeGeneScore(de_model,i.celltype,PCT_mat)
    Cbar = Cbar[Cbar > 0] # while some statistic will be negative; ningjin 2023-9-5
    #names(Cbar) = toupper(names(Cbar))
    Cbar = sort(Cbar, decreasing = T)
    de = names(Cbar)
    score = lapply(marker_list, function(markers){
      markers = intersect(markers,names(Cbar))
      if(length(markers) == 0){
        out = 0
      }else{
        out = sum(Cbar[markers])
      }
      return(out)
    }) %>% unlist 
    # order the score
    score = sort(score, decreasing = T)


    if(!all(score == 0) ){
      score = (score - min(score))/(max(score) - min(score))
    }
    score_list[[i.celltype]] = score
    ## we record the de genes
    if(is.null(de)){
      de_list[[i.celltype]] = "* No useful DE gene *"
    }else{
      de_list[[i.celltype]] = de
    }
    ## we record the overlapped marker
    if(!all(score == 0)){
      for(im in names(score)[score>0]){
        df = data.frame(
          cluster = i.celltype, 
          celltype = im, 
          ude = paste(intersect(de,marker_list[[im]]), collapse = ","),
          score = score[im])
        rownames(df) = NULL
        modelmarker_df = rbind(modelmarker_df, df)
      }
    }
    
    # predict_df = rbind.data.frame(predict_df,
    #                                data.frame(cluster = i.celltype,
    #                                           prediction = ifelse(all(score == 0),'* Loss of match *',names(score)[1]),
    #                                           score = score[1]))
    
    ## ningjin 2025-1-7
    predict_df = rbind.data.frame(predict_df,
                                   data.frame(cluster = i.celltype,
                                              prediction = ifelse(all(score == 0),'* Loss of match *',names(score)[1])))
    rownames(predict_df) = NULL
  }
  return(list(degs = de_list,
              predict = predict_df,
              scores = score_list,
              markers = modelmarker_df))
}


#'
#' Plot gene score for clusters regarding each possible cell types.
#' 
#' @param res A cellLabeler results with prediction 
#' @return A ggplot displaying prediction scores 
#' @import ggplot2
#' @export
#' 
PlotModelScore = function(res){

  score_list = res$ModelScores
  if(is.null(score_list)){
    stop("## The input results has not been with prediction information.")
  }

  base_df = data.frame()
  selected_types = c()
  for (i in seq(length(score_list))) {
    tmp = data.frame(cluster = names(score_list)[i], prediction = names(score_list[[i]]), score = score_list[[i]])
    base_df = rbind.data.frame(base_df, tmp)
    # we only show important
    selected_types = c(selected_types, names(score_list[[i]])[score_list[[i]] > mean(score_list[[i]])])
  }
  base_df = base_df[base_df$prediction %in% unique(selected_types),]
  rownames(base_df) = NULL

  ## selected score = 1 celltypes
  top1_df = subset(base_df, score==1)
  ## selected score >= 0.9 celltypes
  top2_df = subset(base_df, score>=0.9 & score<1)
  
  p = ggplot(base_df, aes(y = prediction, x = cluster,colour = score, size = score))+
    geom_point(data = top1_df, aes(x = cluster, y = prediction),size=8, shape = 21, colour = 'red',fill=NA)+
    geom_point(data = top2_df, aes(x = cluster, y = prediction),size=8, shape = 21, colour = 'black',fill=NA)+
    scale_size(range = c(1,5),guide = 'none')+
    xlab('Cluster')+ylab('CellLabeler prediction')+
    labs(x = "Cluster", y = "CellLabeler prediction", colour = 'Enrichment score',
      caption = "Dot size also represents the score\nRed circle highlights top 1 score\n Black circle notates score exceed 0.9")+
    geom_point()+
    scale_color_gradient(low = "#70D7DD", high = "#E85C90")+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
          axis.text = element_text(size = 10, colour = 'black'),
          axis.title = element_text(size = 15))
  return(p)
}



#'
#' Collect all genes in the analysis,
#' 
#' @param de_model Results of FindAllUniqueMarkers().
#' @export
#' 
allGenes = function(de_model){
  out = lapply(de_model, rownames) %>% unlist %>% unique()
  return(out)
}





#' Select and rank the ude to build up a list of ude for each clusters
#' mainly used in ude benchmark
#' @param de_model A list of output from FindAllUniqueMarkers()
#' @param counts Raw count gene expression matrix; rows are genes and columns as cells
#' @param cluster.id A vector of cluster labels for cells; same as input in FindAllUniqueMarkers(); if combined, then this cluster.id is the combined cluster label
#' 
#' @return A list of combined statistic for each cluster
#' @importFrom dplyr arrange
#' @export 
#' 
RankAllUniqueMarkers = function(de_model, counts, cluster.id){
  if(!is.list(de_model)){
    stop("TidyUpAllMarkers() input must be a list format.")
  }
  
  ## compute pct mat
  PCT_mat = ComputePCT(counts, cluster.id = cluster.id, sub.genes = allGenes(de_model))
  
  ## compute gene score
  cluster_type = sort(unique(as.character(cluster.id)))
  output_list = list()
  for (i.celltype in cluster_type) {
    Cbar = ComputeGeneScore(de_model,i.celltype,PCT_mat)
    Cbar = Cbar[Cbar > 0] # while some statistic will be negative; ningjin 2023-9-5
    Cbar = sort(Cbar, decreasing = T)
    de = names(Cbar)
    tmp = de_model[[i.celltype]]
    ## record expr pct
    tmp$expr_pct = PCT_mat[tmp$gene,i.celltype]
    ## record gene score
    tmp$genescore = -1
    tmp$genescore[match(de, tmp$gene)] = Cbar
    tmp = tmp %>% arrange(p.adjust, -genescore)
    output_list[[i.celltype]] = tmp
  }
  return(output_list)
}





#' Add markers to celllabeler ude results and perform prediction
#' @param res A celllabeler output list involoving at least ModelFits
#' @param counts A normalized gene expression matrix
#' @param markers Marker list with name being cell type names
#' @param cluster.id Character vectors of cluster information for each cell. 
#' 
#' @return A celllabeler output list involoving predictions
#' 
#' @author Jin ning
#' 
#' @export 
AddMarkers = function(res, counts, markers, cluster.id){
    ## check the input res must involve ModelFits
    if(is.null(res$ModelFits)){
        stop("Please run celllabeler() first to identify marker genes.")
    }

    if(!is.null(res$prediction)){
        cat("There is already prediction results and we then update and the original results.\n")
    }
    
    modelfits = res$ModelFits

    ## combine cluster.id according to ude results
    modelfits_cluster = names(modelfits)
    input_cluster_type = unique(as.character(cluster.id))
    if(!setequal(modelfits_cluster,input_cluster_type)){
        ## we need to do combination
        modelfits_cluster = modelfits_cluster[grep(" & ",modelfits_cluster)]
        for(icomb in modelfits_cluster){
            index_icomb = which(cluster.id %in% strsplit(icomb," & ")[[1]])
            cluster.id[index_icomb] = icomb
        }
    }

    ## here, the input counts are filtered in the steps of CreateCellLabelerObject()
    ## while the input modelfits are run on all raw genes
    ## we filtered out the redundant genes in modelfits

    input_genes = intersect(allGenes(modelfits), rownames(counts))
    modelfits = lapply(modelfits, function(x) x[intersect(input_genes,rownames(x)),])
    
    cat("## Compute expression percent of genes ...\n")
    pct = ComputePCT(counts,cluster.id,input_genes)
    pred = ComputePrediction(modelfits,markers,pct)

    res$prediction = pred$predict
    res$ModelScores = pred$scores
    res$ModelMarkers = pred$markers
    return(res)
}




#' Split the prediction dataframe when there is combination before annotation 
#' @param prediction The dataframe of cell type annotation
#' 
#' @return The dataframe of cell type annotation
#' 
#' @export 

ExpandPrediction = function(prediction){
    comb_idx = which(grepl(" & ",prediction[,1]))
    if(length(comb_idx)>0){
        df1 = prediction[-comb_idx,,drop = F]
        for(idx in comb_idx){
        df = data.frame(
            cluster = strsplit(prediction[idx,1]," & ")[[1]],
            prediction = prediction[idx,2])
        df1 = rbind.data.frame(df1,df)
        }
        prediction = df1
    }
    rownames(prediction)=NULL
    return(prediction)
}




#' Update up-regulated DE genes given a new up.thr percent
#' 
#' @param res A list of results from celllabeler
#' @param up.thr Numeric value in [0,1] specifying the up-regulation threshold
#' @param counts Raw count expression matrix
#' @param cluster.id Vector of charaters specifying clustering label
#' @param markers Marker gene list; if NULL, the orginal prediction will be set empty to avoid misleading
#' @return A list of results in the format of celllabeler output
#' 
#' @export 
AdjustUpregulation = function(res, up.thr, counts, cluster.id, markers = NULL){
    if(is.null(res$ModelFits)){
      stop("The input 'res' should contain 'ModelFits'. \n")
    }

    if(up.thr>1 | up.thr <=0){
      stop("The input 'up.thr' should be in [0,1].")
    }

    modelfits = res$ModelFits
    clusters = names(modelfits)

    new_modelfits = list()
    for(id in clusters){
      de_res = modelfits[[id]]
      res_name = setdiff(clusters, id)
      beta_mat = de_res[,paste0("beta.",res_name), drop = F]
      pv_mat = de_res[,paste0("pv.",res_name), drop = F]
      ## modified by ningjin 2023-9-4
      ## we want it upregulated to some extent across groups
      changes = ifelse(apply(beta_mat, 1, function(x) length(which(x>0))/length(x)) >= up.thr , 'upregulated','not upregulated')
      changes[de_res[,"p.adjust"] >0.05] = "not upregulated"
      de_res$change = changes
      new_modelfits[[id]] = de_res
    }
    
    res$ModelFits = new_modelfits
    new_res = AddMarkers(res=res, counts=counts, markers=markers, cluster.id=cluster.id)
    return(new_res)
}