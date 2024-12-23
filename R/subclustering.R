## subclustering code for jin ning, 
## 2024-7-28 20:56:18
## sqsun

#' we use Louvain method
#' @param object Clustering object
#' @param resolution Resolution parameter
#' @param verbose Print output
#'
#' @rdname ClustLouvain
#' @concept clustering
#' @import dplyr
#'
#' @export
#
ClustLouvain <- function(object,
                       modularity.fxn = 1,
                       resolution = 0.8,
                       algorithm = 1,
                       n.start = 10,
                       n.iter = 10,
                       random.seed = 0,
                       group.singletons = TRUE,
                       temp.file.location = NULL,
                       edge.file.name = NULL,
                       verbose = TRUE,...) {
  
  ### compute graph first
  edge_file <- edge.file.name %||% ''
  if (future::nbrOfWorkers() > 1) {
    clustering.results <- future.apply::future_lapply(
      X = resolution,
      FUN = function(r) {
        ids <- RunModularityClusteringCpp(object, 
                                          modularity.fxn, 
                                          r, 
                                          algorithm, 
                                          n.start, 
                                          n.iter, 
                                          random.seed, 
                                          verbose, 
                                          edge_file)
        names(x = ids) <- colnames(x = object)
        ids <- GroupSingletons(ids = ids, SNN = object, verbose = verbose)
        results <- list(factor(x = ids))
        names(x = results) <- paste0('res.', r)
        return(results)
      })
    clustering.results <- as.data.frame(x = clustering.results)
  }else{
    clustering.results <- data.frame(row.names = colnames(x = object))
    for (r in resolution) {
      ids <- RunModularityClusteringCpp(object, 
                                        modularity.fxn, 
                                        r, 
                                        algorithm, 
                                        n.start, 
                                        n.iter, 
                                        random.seed, 
                                        verbose, 
                                        edge_file)
      names(x = ids) <- colnames(x = object)
      ids <- GroupSingletons(ids = ids, SNN = object, group.singletons = group.singletons, verbose = verbose)
      clustering.results[, paste0("res.", r)] <- factor(x = ids)
    }## end for
  }## fi
  ## return the results
  return(clustering.results)
}## end funcs


####### WebGraph
#' The Neighbor class
#'
#' The Neighbor class is used to store the results of neighbor finding
#' algorithms
#'
#' @slot nn.idx Matrix containing the nearest neighbor indices
#' @slot nn.dist Matrix containing the nearest neighbor distances
#' @slot alg.idx The neighbor finding index (if applicable). E.g. the annoy
#' index
#' @slot alg.info Any information associated with the algorithm that may be
#' needed downstream (e.g. distance metric used with annoy is needed when
#' reading in from stored file).
#' @slot cell.names Names of the cells for which the neighbors have been
#' computed.
#'
#' @name Neighbor-class
#' @rdname Neighbor-class
#' @exportClass Neighbor
#'
Neighbor <- setClass(
  Class = 'Neighbor',
  slots = c(
    nn.idx = 'matrix',
    nn.dist = 'matrix',
    alg.idx = 'ANY',
    alg.info = 'list',
    cell.names = 'character'
  )
)

#' The WebGraph Class
#'
#' The WebGraph class inherits from \code{\link[Matrix:sparseMatrix]{dgCMatrix}}.
#' We do this to enable future expandability of graphs.
#'
#' @slot assay.used Optional name of assay used to generate \code{WebGraph} object
#'
#' @name WebGraph-class
#' @rdname WebGraph-class
#' @exportClass WebGraph
#'
#' @seealso \code{\link[Matrix]{dgCMatrix-class}}
#'
WebGraph <- setClass(
  Class = 'WebGraph',
  contains = "dgCMatrix",
  slots = list(
    assay.used = 'character'  ##'OptionalCharacter'
  )
)


#############################################################################
#' Find neighbor cells
#' @param query Matrix of data to query against object. If missing, defaults to
#' object.
#' @param distance.matrix Boolean value of whether the provided matrix is a
#' distance matrix; note, for objects of class \code{dist}, this parameter will
#' be set automatically
#' @param k.param Defines k for the k-nearest neighbor algorithm
#' @param return.neighbor Return result as \code{\link{Neighbor}} object. Not
#' used with distance matrix input.
#' @param compute.SNN also compute the shared nearest neighbor graph
#' @param prune.SNN Sets the cutoff for acceptable Jaccard index when
#' computing the neighborhood overlap for the SNN construction. Any edges with
#' values less than or equal to this will be set to 0 and removed from the SNN
#' graph. Essentially sets the stringency of pruning (0 --- no pruning, 1 ---
#' prune everything).
#' @param nn.method Method for nearest neighbor finding. Options include: rann,
#' annoy
#' @param annoy.metric Distance metric for annoy. Options include: euclidean,
#' cosine, manhattan, and hamming
#' @param n.trees More trees gives higher precision when using annoy approximate
#' nearest neighbor search
#' @param nn.eps Error bound when performing nearest neighbor seach using RANN;
#' default of 0.0 implies exact nearest neighbor search
#' @param verbose Whether or not to print output to the console
#' @param force.recalc Force recalculation of (S)NN.
#' @param l2.norm Take L2Norm of the data
#' @param cache.index Include cached index in returned Neighbor object
#' (only relevant if return.neighbor = TRUE)
#' @param index Precomputed index. Useful if querying new data against existing
#' index to avoid recomputing.
#'
#' @importFrom RANN nn2
#' @importFrom methods as
#' @importFrom Matrix sparseMatrix
#' @export
#' @concept clustering
#'
FindNeighborCells <- function(object,
                          query = NULL,
                          distance.matrix = FALSE,
                          k.param = 20,
                          return.neighbor = FALSE,
                          compute.SNN = !return.neighbor,
                          prune.SNN = 1/15,
                          nn.method = "annoy",
                          n.trees = 50,
                          annoy.metric = "euclidean",
                          nn.eps = 0,
                          verbose = TRUE,
                          force.recalc = FALSE,
                          l2.norm = FALSE,
                          cache.index = FALSE,
                          index = NULL, ... ) {
  
  ## check

  if (is.null(x = dim(x = object))) {
    warning("Object should have two dimensions, attempting to coerce to matrix", call. = FALSE)
    object <- as.matrix(x = object)
  }## end fi
  
  if (is.null(rownames(x = object))) {
    stop("Please provide rownames (cell names) with the input object")
  }## end fi
  
  n.cells <- nrow(x = object)
  if (n.cells < k.param) {
    warning("k.param set larger than number of cells. Setting k.param to number of cells - 1.", call. = FALSE)
    k.param <- n.cells - 1
  }## end fi
  
  if (l2.norm) {
    object <- L2Norm(mat = object)
    query <- query %iff% L2Norm(mat = query)
  }## end fi
  query <- query %||% object
  # find the k-nearest neighbors for each single cell
  if (!distance.matrix) {
    if (verbose) {
      if (return.neighbor) {
        message("Computing nearest neighbors")
      } else {
        message("Computing nearest neighbor graph")
      }
    }## end fi
    nn.ranked <- NNHelper(data = object,
                          query = query,
                          k = k.param,
                          method = nn.method,
                          n.trees = n.trees,
                          searchtype = "standard",
                          eps = nn.eps,
                          metric = annoy.metric,
                          cache.index = cache.index,
                          index = index)
    
    if (return.neighbor) {
      if (compute.SNN) {
        warning("The SNN graph is not computed if return.neighbor is TRUE.", call. = FALSE)
      }
      return(nn.ranked)
    }## end fi
    
    nn.ranked <- Indices(object = nn.ranked)
  } else {
    if (verbose) {
      message("Building SNN based on a provided distance matrix")
    }## end fi
    
    knn.mat <- matrix(data = 0, ncol = k.param, nrow = n.cells)
    knd.mat <- knn.mat
    for (i in 1:n.cells) {
      knn.mat[i, ] <- order(object[i, ])[1:k.param]
      knd.mat[i, ] <- object[i, knn.mat[i, ]]
    }
    nn.ranked <- knn.mat[, 1:k.param]
  }## end fi
  
  ## convert nn.ranked into a WebGraph
  j <- as.numeric(x = t(x = nn.ranked))
  i <- ((1:length(x = j)) - 1) %/% k.param + 1
  nn.matrix <- as(object = sparseMatrix(i = i, j = j, x = 1, dims = c(nrow(x = object), nrow(x = object))), Class = "WebGraph")
  rownames(x = nn.matrix) <- rownames(x = object)
  colnames(x = nn.matrix) <- rownames(x = object)
  neighbor.graphs <- list(nn = nn.matrix)
  
  if(compute.SNN) {
    if (verbose) {
      message("Computing SNN")
    }## end fi
    
    snn.matrix <- ComputeSNN(nn_ranked = nn.ranked, prune = prune.SNN)
    
    rownames(x = snn.matrix) <- rownames(x = object)
    colnames(x = snn.matrix) <- rownames(x = object)
    snn.matrix <- as.WebGraph(x = snn.matrix)
    neighbor.graphs[["snn"]] <- snn.matrix
  }## end fi
  return(neighbor.graphs)
}## end funcs


# Internal helper function to dispatch to various neighbor finding methods
#
# @param data Input data
# @param query Data to query against data
# @param k Number of nearest neighbors to compute
# @param method Nearest neighbor method to use: "rann", "annoy"
# @param cache.index Store algorithm index with results for reuse
# @param ... additional parameters to specific neighbor finding method
#
#' @importFrom methods new
#

NNHelper <- function(data, query = data, k, method, cache.index = FALSE, ...) {
  args <- as.list(x = sys.frame(which = sys.nframe()))
  args <- c(args, list(...))
  results <- (
    switch(
      EXPR = method,
      "rann" = {
        args <- args[intersect(x = names(x = args), y = names(x = formals(fun = nn2)))]
        do.call(what = 'nn2', args = args)
      },
      "annoy" = {
        args <- args[intersect(x = names(x = args), y = names(x = formals(fun = AnnoyNN)))]
        do.call(what = 'AnnoyNN', args = args)
      },
      stop("Invalid method. Please choose one of 'rann', 'annoy'")
    ) )
  ##
  n.ob <- new(Class = 'Neighbor',
              nn.idx = results$nn.idx,
              nn.dist = results$nn.dists,
              alg.info = results$alg.info %||% list(),
              cell.names = rownames(x = query))
  
  if (isTRUE(x = cache.index) && !is.null(x = results$idx)) {
    slot(object = n.ob, name = "alg.idx") <- results$idx
  }## end fi
  return(n.ob)
}## end func




# Group single cells that make up their own cluster in with the cluster they are
# most connected to.
#
# @param ids Named vector of cluster ids
# @param SNN SNN graph used in clustering
# @param group.singletons Group singletons into nearest cluster. If FALSE, assign all singletons to
# a "singleton" group
#
# @return Returns Seurat object with all singletons merged with most connected cluster
#
GroupSingletons <- function(ids, SNN, group.singletons = TRUE, verbose = TRUE) {
  # identify singletons
  singletons <- c()
  singletons <- names(x = which(x = table(ids) == 1))
  singletons <- intersect(x = unique(x = ids), singletons)
  if (!group.singletons) {
    ids[which(ids %in% singletons)] <- "singleton"
    return(ids)
  }
  # calculate connectivity of singletons to other clusters, add singleton
  # to cluster it is most connected to
  cluster_names <- as.character(x = unique(x = ids))
  cluster_names <- setdiff(x = cluster_names, y = singletons)
  connectivity <- vector(mode = "numeric", length = length(x = cluster_names))
  names(x = connectivity) <- cluster_names
  new.ids <- ids
  for (i in singletons) {
    i.cells <- names(which(ids == i))
    for (j in cluster_names) {
      j.cells <- names(which(ids == j))
      subSNN <- SNN[i.cells, j.cells]
      set.seed(1) # to match previous behavior, random seed being set in WhichCells
      if (is.object(x = subSNN)) {
        connectivity[j] <- sum(subSNN) / (nrow(x = subSNN) * ncol(x = subSNN))
      } else {
        connectivity[j] <- mean(x = subSNN)
      }
    }
    m <- max(connectivity, na.rm = T)
    mi <- which(x = connectivity == m, arr.ind = TRUE)
    closest_cluster <- sample(x = names(x = connectivity[mi]), 1)
    ids[i.cells] <- closest_cluster
  }
  if (length(x = singletons) > 0 && verbose) {
    message(paste(
      length(x = singletons),
      "singletons identified.",
      length(x = unique(x = ids)),
      "final clusters."
    ))
  }
  return(ids)
}## end func

#' Run annoy
#' @param data Data to build the index with
#' @param query A set of data to be queried against data
#' @param metric Distance metric; can be one of "euclidean", "cosine", "manhattan","hamming"
#' @param n.trees More trees gives higher precision when querying
#' @param k Number of neighbors
#' @param search.k During the query it will inspect up to search_k nodes which gives you a run-time tradeoff between better accuracy and speed.
#' @param include.distance Include the corresponding distances
#' @param index optional index object, will be recomputed if not provided
#' 
#' @import dplyr
#' 
#' 
AnnoyNN <- function(data,
                    query = data,
                    metric = "euclidean",
                    n.trees = 50,
                    k,
                    search.k = -1,
                    include.distance = TRUE,
                    index = NULL) {
  idx <- index %||% AnnoyBuildIndex(
    data = data,
    metric = metric,
    n.trees = n.trees)
  nn <- AnnoySearch(
    index = idx,
    query = query,
    k = k,
    search.k = search.k,
    include.distance = include.distance)
  nn$idx <- idx
  nn$alg.info <- list(metric = metric, ndim = ncol(x = data))
  return(nn)
}## end func

# Build the annoy index
#
# @param data Data to build the index with
# @param metric Distance metric; can be one of "euclidean", "cosine", "manhattan",
# "hamming"
# @param n.trees More trees gives higher precision when querying
#
#' @importFrom RcppAnnoy AnnoyEuclidean AnnoyAngular AnnoyManhattan AnnoyHamming
#
AnnoyBuildIndex <- function(data, metric = "euclidean", n.trees = 50) {
  f <- ncol(x = data)
  a <- switch(
    EXPR = metric,
    "euclidean" =  new(Class = RcppAnnoy::AnnoyEuclidean, f),
    "cosine" = new(Class = RcppAnnoy::AnnoyAngular, f),
    "manhattan" = new(Class = RcppAnnoy::AnnoyManhattan, f),
    "hamming" = new(Class = RcppAnnoy::AnnoyHamming, f),
    stop ("Invalid metric")
  )
  for (ii in seq(nrow(x = data))) {
    a$addItem(ii - 1, data[ii, ])
  }
  a$build(n.trees)
  return(a)
}## end func

#' Search an Annoy approximate nearest neighbor index
#' @param Annoy index, built with AnnoyBuildIndex
#' @param query A set of data to be queried against the index
#' @param k Number of neighbors
#' @param search.k During the query it will inspect up to search_k nodes which gives you a run-time tradeoff between better accuracy and speed.
#' @param include.distance Include the corresponding distances in the result
#' 
#' @return A list with 'nn.idx' (for each element in 'query', the index of the nearest k elements in the index) and 'nn.dists' (the distances of the nearestk elements)
#
#' @importFrom future plan
#' @importFrom future.apply future_lapply
#
AnnoySearch <- function(index, query, k, search.k = -1, include.distance = TRUE) {
  n <- nrow(x = query)
  idx <- matrix(nrow = n,  ncol = k)
  dist <- matrix(nrow = n, ncol = k)
  convert <- methods::is(index, "Rcpp_AnnoyAngular")
  if (!inherits(x = future::plan(), what = "multicore")) {
    oplan <- future::plan(strategy = "sequential")
    on.exit(future::plan(oplan), add = TRUE)
  }
  res <- future.apply::future_lapply(X = 1:n, FUN = function(x) {
    res <- index$getNNsByVectorList(query[x, ], k, search.k, include.distance)
    # Convert from Angular to Cosine distance
    if (convert) {
      res$dist <- 0.5 * (res$dist * res$dist)
    }
    list(res$item + 1, res$distance)
  })
  for (i in 1:n) {
    idx[i, ] <- res[[i]][[1]]
    if (include.distance) {
      dist[i, ] <- res[[i]][[2]]
    }
  }
  return(list(nn.idx = idx, nn.dists = dist))
}## end func




#' As web graph functions
#' @export 
#' 
#'
as.WebGraph <- function(x, ...) {
  UseMethod(generic = 'as.WebGraph', object = x)
}## end func

#' As web graph functions when given a neighbor
as.WebGraph.Neighbor <- function(x, weighted = TRUE, ...) {
  # CheckDots(...)
  j <- as.integer(x = Indices(object = x) - 1)
  i <- as.integer(x = rep(x = (1:nrow(x = x)) - 1, times = ncol(x = x)))
  vals <- if (weighted) {
    as.vector(x = Distances(object = x))
  } else {
    1
  }
  graph <- new(
    Class = "dgTMatrix",
    i = i,
    j = j,
    x = vals,
    Dim = as.integer(x = c(nrow(x = x), nrow(x = x)))
  )
  colnames(x = graph) <- rownames(x = graph) <- Cells(x = x)
  graph <- as.WebGraph.Matrix(x = graph)
  return(graph)
}## end func



#' As webGraph matrix
#' @export 
#' 
#' 
as.WebGraph.Matrix <- function(x, ...) {
  # CheckDots(...)
  x <- as.sparse(x = x)
  if (is.null(x = rownames(x = x))) {
    stop("Please provide rownames to the matrix before converting to a WebGraph.")
  }
  if (is.null(x = colnames(x = x))) {
    stop("Please provide colnames to the matrix before converting to a WebGraph.")
  }
  return(as(object = x, Class = "WebGraph"))
}## end func


