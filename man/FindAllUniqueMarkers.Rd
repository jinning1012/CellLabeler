% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ude.R
\name{FindAllUniqueMarkers}
\alias{FindAllUniqueMarkers}
\title{Fitting the firth logistic regression model to perform unique marker gene analysis for clustered scRNA-seq data}
\usage{
FindAllUniqueMarkers(
  data.use,
  sample.id,
  cluster.id,
  upPercent = 0.9,
  verbose = T,
  num.core = 1
)
}
\arguments{
\item{data.use}{Log-normalized gene expression matrix; rows are genes and columns as cells}

\item{sample.id}{Character vector of sample information in metadata.}

\item{cluster.id}{Character vector of cluster information in metadata.}

\item{upPercent}{Numberic value in \link{0,1} as the proportion of the remaining clusters in which a gene is upregulated in the interested group compared with it. Highter, the more strict. Default 0.9.}

\item{verbose}{Bool variable to indicate if print the message.}

\item{num.core}{Number of cores in multi-threaded running. Default is 1.}
}
\description{
Fitting the firth logistic regression model to perform unique marker gene analysis for clustered scRNA-seq data
}
