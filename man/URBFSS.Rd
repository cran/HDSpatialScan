% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scan_fonc.R
\name{URBFSS}
\alias{URBFSS}
\title{URBFSS scan procedure}
\usage{
URBFSS(
  data,
  MC = 999,
  typeI = 0.05,
  nbCPU = 1,
  times = NULL,
  initialization,
  permutations
)
}
\arguments{
\item{data}{matrix. Matrix of the data, the rows correspond to the sites (or to the individuals if the observations are by individuals and not by sites) and each column represents an observation time. The times must be the same for each site/individual.}

\item{MC}{numeric. Number of Monte-Carlo permutations to evaluate the statistical significance of the clusters. By default: 999.}

\item{typeI}{numeric. The desired type I error. A cluster will be evaluated as significant if its associated p-value is less than typeI. By default 0.05.}

\item{nbCPU}{numeric. Number of CPU. If nbCPU > 1 parallelization is done. By default: 1.}

\item{times}{numeric. Times of observation of the data. By default NULL.}

\item{initialization}{list. Initialization for the scan procedure (see \code{\link{InitScan}} for more details).}

\item{permutations}{matrix. Indices of permutations of the data.}
}
\value{
An object of class ResScanOutputUniFunct.
}
\description{
This function computes the URBFSS (Univariate Rank-Based Functional scan statistic).
}
\seealso{
\code{\link{MRBFSS}} which is the multivariate version of the URBFSS
}
