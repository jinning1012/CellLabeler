
#' Permutation of two of all categories
#' @param vec A vector of elements where my permutation performs
#' 
#' @return Returns a list of permutations of two elements
#' 
#' @export 
#' 
perm2 <- function(vec){
  out <- list()
  flag <- 1
  for(i in vec){
    for(j in setdiff(vec,i)){
      out[[flag]] <- c(i,j)
      flag <- flag+1
    }
  }
  out <- do.call(cbind,out)
  return(out)
}

#' Cauchy combination rules
#' @param pvalues n x m pvalue matrix; n genes, m permutations
#' @param weights
#' 
#' @export 
#' 
#' @return Returns a n-length combined pvalues 
#' 
CombinePValues <- function(pvalues, weights=NULL){
  if(!is.matrix(pvalues)){pvalues <- as.matrix(pvalues)}
  ## to avoid extremely values
  pvalues[which(pvalues==0)] <- 5.55e-17
  pvalues[which((1-pvalues)<1e-3)] <- 0.99
  
  num_pval <- ncol(pvalues)
  num_gene <- nrow(pvalues)
  if(is.null(weights)){
    weights <- matrix(rep(1.0/num_pval, num_pval*num_gene), ncol=num_pval )
  }# end fi
  if( (nrow(weights) != num_gene) || (ncol(weights) != num_pval)){
    stop("the dimensions of weights does not match that of combined pvalues")
  }# end fi
  
  Cstat <- tan((0.5 - pvalues)*pi)
  
  wCstat <- weights*Cstat
  Cbar <- apply(wCstat, 1, sum)
  #combined_pval <- 1.0/2.0 - atan(Cbar)/pi
  combined_pval <- 1.0 - pcauchy(Cbar)	
  combined_pval[which(combined_pval <= 0)] <- 5.55e-17
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
#
#' @export 
#' 
RegressOutMatrix <- function(data.expr,
                             latent.data = NULL,
                             features.regress = NULL,
                             model.use = NULL,
                             use.umi = FALSE,
                             verbose = FALSE) {
  
  ## Do we bypass regression and simply return data.expr?
  bypass <- vapply(X = list(latent.data, model.use), FUN = is.null, FUN.VALUE = logical(length = 1L) )
  
  if (any(bypass)) {
    return(data.expr)
  }## end fi
  
  ## Check model.use
  possible.models <- c("linear", "poisson", "negbinom")
  if (!model.use %in% possible.models) {
    stop(paste(
      model.use,
      "is not a valid model. Please use one the following:",
      paste0(possible.models, collapse = ", ")
    ))
  }## end fi
  
  ## Check features.regress
  if (is.null(x = features.regress)) {
    features.regress <- 1:nrow(x = data.expr)
  }## end fi
  
  if (is.character(x = features.regress)) {
    features.regress <- intersect(x = features.regress, y = rownames(x = data.expr))
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
  
  use.umi <- ifelse(test = model.use != 'linear', yes = TRUE, no = use.umi)
  
  ## Create formula for regression
  vars.to.regress <- colnames(x = latent.data)
  fmla <- paste('GENE ~', paste(vars.to.regress, collapse = '+'))
  fmla <- as.formula(object = fmla)
  if (model.use == "linear") {
    # In this code, we'll repeatedly regress different Y against the same X
    # (latent.data) in order to calculate residuals.  Rather that repeatedly
    # call lm to do this, we'll avoid recalculating the QR decomposition for the
    # latent.data matrix each time by reusing it after calculating it once
    regression.mat <- cbind(latent.data, data.expr[1,])
    colnames(regression.mat) <- c(colnames(x = latent.data), "GENE")
    qr <- lm(fmla, data = regression.mat, qr = TRUE)$qr
    rm(regression.mat)
  }## end fi
  
  ## Make results matrix
  data.resid <- matrix(nrow = nrow(x = data.expr), ncol = ncol(x = data.expr) )
  
  if (verbose) {
    pb <- txtProgressBar(char = '=', style = 3, file = stderr())
  }## end fi
  
  for (i in 1:length(x = features.regress)) {
    x <- features.regress[i]
    regression.mat <- cbind(latent.data, data.expr[x, ])
    colnames(x = regression.mat) <- c(vars.to.regress, 'GENE')
    regression.mat <- switch(
      EXPR = model.use,
      'linear' = qr.resid(qr = qr, y = data.expr[x,]),
      'poisson' = residuals(object = glm(
        formula = fmla,
        family = 'poisson',
        data = regression.mat),
        type = 'pearson'),
      'negbinom' = NBResiduals(
        fmla = fmla,
        regression.mat = regression.mat,
        gene = x))
    
    data.resid[i, ] <- regression.mat
    if (verbose) {
      setTxtProgressBar(pb = pb, value = i / length(x = features.regress))
    }## end fi
  }## end for
  
  if (verbose) {
    close(con = pb)
  }## end fi
  
  ## start if
  if (use.umi) {
    data.resid <- log1p(x = Sweep(x = data.resid, MARGIN = 1, STATS = apply(X = data.resid, MARGIN = 1, FUN = min), FUN = '-' ))
  }##end fi
  
  dimnames(x = data.resid) <- dimnames(x = data.expr)
  return(data.resid)
}## end func


#'
#' Fitting the firth logistic regression model to perform unique marker gene analysis for clustered scRNA-seq data
#' 
#' @param data.use Log-normalized gene expression matrix; rows are genes and columns as cells
#' @param sample.id Character vector of sample information in metadata. 
#' @param cluster.id Character vector of cluster information in metadata. 
#' @param upPercent Numberic value in [0,1] as the proportion of the remaining clusters in which a gene is upregulated in the interested group compared with it. Highter, the more strict. Default 0.9.
#' @param num.core Number of cores in multi-threaded running. Default is 1.
#' @param verbose Bool variable to indicate if print the message.
#' 
#' @export
#' 
#' 
FindAllUniqueMarkers <- function(data.use,
                                 sample.id,
                                 cluster.id,
                                 upPercent = 0.9,
                                 verbose = T,
                                 num.core = 1){
  
  suppressPackageStartupMessages({
    require(logistf)
    require(dplyr)
    require(doParallel)
    #require(utils)
  })
  
  ## ******************************************* ##
  ##       function for one comparison           ##
  ## ******************************************* ##
  clusters = unique(cluster.id)
  runOneCompare <- function(ident.1, ident.2){
    cells.1 = which(cluster.id == ident.1)
    cells.2 = which(cluster.id == ident.2)
    sub_data.use = data.use[,c(cells.1, cells.2)]
    sub_group = cluster.id[c(cells.1, cells.2)]
    sub_sample = sample.id[c(cells.1, cells.2)]
    outcome = as.numeric(sub_group == ident.1)
    
    # Step1: decide to remove batch effects
    # 1. escape when there is no batch information
    # 2. escape when each ident contains only one batch
    # modified by ningjin; 2023-9-14 20:06
    if(length(unique(sub_sample))==1){
      resi <- sub_data.use
      resi <- resi[Matrix::rowSums(resi)!=0, ] ## remove some genes
    }else{
      latent.data = data.frame(sample = sub_sample)
      rownames(latent.data) = colnames(sub_data.use)
      resi <- RegressOutMatrix(data.expr = sub_data.use, 
                                latent.data = latent.data,
                                model.use = 'linear',
                                features.regress = NULL,
                                verbose = F)
      resi <- resi[Matrix::rowSums(resi)!=0, ]
    }## end batch regression
    ## ********************************** ##
    ##           firth's model            ##
    ## ********************************** ##
      ## modified by ningjin; we change the loop to foreach
    ## runSomeGenes() enables to run several genes parallelly while catch errors
    runSomeGenes <- function(index) {
      res_firth_tmp <- foreach(k = index,.errorhandling='pass') %dopar% {
        res_tmp <- tryCatch({
          filr <- logistf(outcome ~ resi[k,])
          res_tmp <- c(beta=filr$coefficients[2], pvalue=filr$prob[2])
        },
        warning = function(war) {return(paste0('logistf warining: escape the gene ',rownames(resi)[k]))},
        error = function(err) {return(paste0('logistf error: escape the gene ',rownames(resi)[k]))},
        finally = {
          # return('other things')
        }) # END tryCatch
        return(res_tmp)
      }# END foreach
      output_idx <- which(unlist(lapply(res_firth_tmp, length)) == 2)
      res_firth_tmp <- res_firth_tmp[output_idx]
      res_firth_tmp <- do.call(rbind,res_firth_tmp)
      rownames(res_firth_tmp) <- rownames(resi)[index[output_idx]]
      colnames(res_firth_tmp) <- c('beta','pvalue')
      return(res_firth_tmp)
    }
    
    len <- nrow(resi)
    res_firth <- Reduce(rbind, 
                        lapply(suppressWarnings(split(seq(len), seq(len / 100))), 
                               function(index) runSomeGenes(index)))
    return(res_firth)
  }## end function
  
  ## ******************************************************
  ##           run for all group comparison
  ## ******************************************************
  all_res_list <- list() 
  ## combination of clusters
  ident.comb <- combn(clusters,2)
  ## start loop for comparisons
  all_res_list <- parallel::mclapply(1:ncol(ident.comb), mc.cores = num.core, function(i.comb)
  { 
    ident.1 <- ident.comb[1,i.comb]
    ident.2 <- ident.comb[2,i.comb]
    ## ningjin 20240701
    if(verbose){
      cat(paste0("### Logistf on ",ident.1, " versus ",ident.1,"...\n"))
    }
    output <- runOneCompare(ident.1, ident.2) 
    return(output)
  })
  
  res_name <- lapply(1:ncol(ident.comb), function(x) paste0(ident.comb[1,x],'_versus_',ident.comb[2,x])) %>% unlist
  names(all_res_list) <- res_name
  rm(res_name)
  
  ## for the rest of comparisons
  ident.perm <- perm2(clusters)
  for(i in 1:ncol(ident.perm)){
    res_name <- paste0(ident.perm[1,i],'_versus_',ident.perm[2,i])
    if(!(res_name %in% names(all_res_list))){
      res_tmp <- all_res_list[[paste0(ident.perm[2,i],'_versus_',ident.perm[1,i])]]
      res_tmp[,1] <- -res_tmp[,1]
      all_res_list[[res_name]] <- res_tmp
      rm(res_name)
    }
  }
  ## collect all firth model results
  output <- list() # final output 
  for(ident.one in clusters){
    res_name <- setdiff(clusters, ident.one)
    res_relative <- all_res_list[paste0(ident.one,'_versus_',res_name)]
    ## cauchy combination rule ##
    pvs <- lapply(res_relative, function(x) return(x[, 2]))
    shared_genes <- lapply(pvs, names) %>% Reduce(intersect,.)
    pvs <- lapply(pvs, function(x) x[shared_genes] ) %>% do.call(cbind,.)
    if(ncol(pvs) == 1){
      combined.pv = as.vector(pvs)
    }else{
      combined.pv = CombinePValues(pvs)}
    
    avg_beta <- lapply(res_relative, function(x){x[shared_genes,1]}) %>% do.call(cbind,.) %>% rowMeans(.)
    
    de_res <- lapply(res_relative, function(x){matrix(x[shared_genes,], nrow = length(shared_genes))}) %>% do.call(cbind,.) %>% as.data.frame()
    colnames(de_res) <- lapply(res_name, function(j){paste0(c('beta.','pv.'),j)}) %>% unlist
    rownames(de_res) <- shared_genes
    
    de_res[,'gene'] <- shared_genes
    de_res[,'avg_beta'] = avg_beta
    de_res[,'combined.pv'] = combined.pv
    de_res[,'p.adjust'] = p.adjust(combined.pv, method = 'BH')
    
    ## check if upregulated ##
    beta_mat <- de_res[,grep('beta.',colnames(de_res))]
    pv_mat <- de_res[,grep('pv.',colnames(de_res))]
    ## modified by ningjin 2023-9-4
    ## we want it upregulated to some extent across groups
    changes <- ifelse(apply(beta_mat, 1, function(x) length(which(x>0))/length(x)) >= upPercent , 'upregulated','not upregulated')
    de_res$change <- changes
    output[[ident.one]] <- de_res
  }
  
  ## computation time
  # time2 <- Sys.time()
  # message(paste0('## Computation time (mins):',round(difftime(time2,time1, units = 'mins'), digits = 3),'.'))
  
  ## finally
  return(output)
}



#'
#' Compute gene expression percent of single-cell RNA data.
#' 
#' @param cnt_data Gene expression raw count matrix.
#' @param cluster.id Character vector of cluster information in metadata. 
#' @param sub.genes Subset of gene names.Default is NULL, then all genes will be used.
#' 
#' 
#' @export
#' 
ComputePCT <- function(cnt_data, cluster.id = NULL, sub.genes = NULL){
  if(!is.null(sub.genes)){
    ## check if sub.genes exceed provided genes
    if(!all(sub.genes %in% rownames(cnt_data))){
      stop(paste0("## Given genes in ComputePCT() exceed the range of all genes from counts data."))
    }
    cnt_data <- cnt_data[sub.genes,]
  }
  ngene = nrow(cnt_data)
  ## if clusterid is blank, it means compute gene expression pct across all cells
  if(is.null(cluster.id)){
    out <- apply(cnt_data, 1, function(gene){
      length(which(gene > 0))/length(gene)
    })
  }else{
    groups <- unique(cluster.id)
    out <- matrix(nrow = ngene, ncol = length(groups), dimnames = list(rownames(cnt_data), groups))
    for (ig in seq(length(groups))) {
      group = groups[ig]
      cells.1 <- which(cluster.id == group)
      sub_cnt <- matrix(cnt_data[,cells.1], ncol = length(cells.1))
      out[,ig] <- apply(sub_cnt, 1, function(gene){
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
ComputeGeneScore <- function(de_model, i.celltype, PCT_mat){
  PCT_mat[which(PCT_mat == 0,arr.ind = T)] <- 1e-10
  
  alltypes <- names(de_model)
  degs <- rownames(de_model[[i.celltype]])[de_model[[i.celltype]]$change == 'upregulated' & de_model[[i.celltype]]$p.adjust < 0.05]
  
  ## Cstats ##
  pvalues <- as.matrix(de_model[[i.celltype]][degs,grep('pv.',de_model[[i.celltype]] %>% colnames)])
  pvalues[which(pvalues==0)] <- 5.55e-17
  pvalues[which((1-pvalues)<1e-3)] <- 0.99
  
  num_pval <- ncol(pvalues)
  num_gene <- nrow(pvalues)
  weights <- matrix(rep(1.0/num_pval, num_pval*num_gene), ncol=num_pval )
  Cstat <- tan((0.5 - pvalues)*pi)
  
  ## gene by npval
  wCstat <- weights*Cstat
  
  ## pct
  Cbar <- matrix(nrow = num_gene, ncol = num_pval, dimnames = list(degs,NULL))
  
  for (j in seq(num_pval)) {
    pct = PCT_mat[degs,i.celltype]/PCT_mat[degs, setdiff(alltypes, i.celltype)[j]]
    Cbar[,j] = wCstat[,j] * pct
  }
  Cbar <- Matrix::rowMeans(Cbar)
  output <- sort(Cbar, decreasing = T)
  output <- output[output > 0]
  return(output)
}


#'
#' Compute gene score and perform predictions.
#' 
#' @param de_model Results of FindAllUniqueMarkers().
#' @param marker_list List of makers, with list names being cell type names.
#' @param PCT_mat Gene expression percent matrix. Rows are genes and columns are clusters.
#' 
#' @export
#' 
ComputePrediction <- function(de_model, marker_list, PCT_mat){
  require(stringr)
  alltypes <- names(de_model)
  
  # check if the genes from real data is uppercase
  if(rownames(de_model[[1]])[1] == toupper(rownames(de_model[[1]])[1])){
    ## to upper
    marker_list <- lapply(marker_list, toupper)
  }else{
    ## to captical first letter
    marker_list <- lapply(marker_list, str_to_title)
  }
  
  # init
  score_list <- list()
  de_list <- list()
  predict_df <- data.frame()
  
  for (i.celltype in alltypes) {
    Cbar <- ComputeGeneScore(de_model,i.celltype,PCT_mat)
    Cbar <- Cbar[Cbar > 0] # while some statistic will be negative; ningjin 2023-9-5
    #names(Cbar) <- toupper(names(Cbar))
    Cbar <- sort(Cbar, decreasing = T)
    de <- names(Cbar)
    score <- lapply(marker_list, function(markers){
      markers <- intersect(markers,names(Cbar))
      if(length(markers) == 0){
        out = 0
      }else{
        out = sum(Cbar[markers])
      }
      return(out)
    }) %>% unlist 
    # order the score
    score <- sort(score, decreasing = T)
    if(!all(score == 0) ){
      score <- (score - min(score))/(max(score) - min(score))
    }
    score_list[[i.celltype]] <- score
    ## we record the de genes
    if(is.null(de)){
      de_list[[i.celltype]] <- "* No useful DE gene *"
    }else{
      de_list[[i.celltype]] <- de
    }

    predict_df <- rbind.data.frame(predict_df,
                                   data.frame(cluster = i.celltype,
                                              prediction = ifelse(all(score == 0),' * Loss of match *',names(score)[1]),
                                              score = score[1]))
    rownames(predict_df) <- NULL
  }
  return(list(predict = predict_df,
              scores = score_list,
              degs = de_list))
}


#'
#' Plot gene score for clusters regarding each possible cell types.
#' 
#' @param model_predict Results of ComputePrediction(). A list contating scores.
#' 
#' @export
#' 
CellLabeler.plot <- function(model_predict){
  require(ggplot2)
  df <- data.frame()
  selected_types <- c()
  score_list <- model_predict$scores
  for (i in seq(length(score_list))) {
    df <- rbind.data.frame(df, data.frame(cluster = names(score_list)[i],
                                          predict = names(score_list[[i]]),
                                          value = score_list[[i]]))
    # we only show important
    selected_types <- c(selected_types,
                        names(score_list[[i]])[score_list[[i]] > mean(score_list[[i]])])
  }
  
  df <- df[df$predict %in% unique(selected_types),]
  
  #rownames(df) <- NULL
  df2 <- df %>% 
    group_by(cluster) %>% 
    summarise(predict = predict[which.max(value)],Value = max(value), Order = 'First') %>% 
    as.data.frame()
  
  ## check score = 0 clusters
  idx <- which(df2$Value == 0)
  if(length(idx) > 0){
    df2$predict[idx] <- '* Loss of match *'
  }
  p <- ggplot(df, aes(y = predict, x = cluster,colour = value, size = value))+
    xlab('Original label')+ylab('CellLabeler prediction')+
    labs(colour = 'Gene Score')+
    geom_point()+
    scale_color_gradient(low = "#161853", high = "#EC255A")+
    scale_size(guide = 'none')+
    geom_point(data = df2, aes(x = cluster, y = predict),size = 10, shape = 21, colour = 'black')+
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
#'
#' 
allGenes <- function(de_model){
  out <- lapply(de_model, rownames) %>% unlist %>% unique()
  return(out)
}

