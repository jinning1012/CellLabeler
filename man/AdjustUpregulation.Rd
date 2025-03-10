% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ude.R
\name{AdjustUpregulation}
\alias{AdjustUpregulation}
\title{Update up-regulated DE genes given a new up.thr percent}
\usage{
AdjustUpregulation(res, up.thr, counts, cluster.id, markers = NULL)
}
\arguments{
\item{res}{A list of results from celllabeler}

\item{up.thr}{Numeric value in \link{0,1} specifying the up-regulation threshold}

\item{counts}{Raw count expression matrix}

\item{cluster.id}{Vector of charaters specifying clustering label}

\item{markers}{Marker gene list; if NULL, the orginal prediction will be set empty to avoid misleading}
}
\value{
A list of results in the format of celllabeler output
}
\description{
Update up-regulated DE genes given a new up.thr percent
}
