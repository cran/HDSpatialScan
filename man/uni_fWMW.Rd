% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wmw_fonc.R
\name{uni_fWMW}
\alias{uni_fWMW}
\title{Index for the NPFSS scan procedure (univariate functional case)}
\usage{
uni_fWMW(signs, matrix_clusters)
}
\arguments{
\item{signs}{numeric matrix. Matrix of signs of the data, the rows correspond to the sites (or the individuals) and each column represents an observation time.}

\item{matrix_clusters}{numeric matrix. Matrix in which each column represents a potential cluster. It is the result of the "clusters" function.}
}
\value{
numeric vector.
}
\description{
This function returns the index we want to maximize on the set of potential clusters, for each potential cluster
}
