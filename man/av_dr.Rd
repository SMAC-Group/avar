% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/av_processes.R
\name{av_dr}
\alias{av_dr}
\title{Calculate Theoretical Allan Variance for Drift Process}
\usage{
av_dr(delta, n)
}
\arguments{
\item{delta}{A \code{double} value for the noise parameter \eqn{\delta}{delta}.}

\item{n}{An \code{integer} value for the size of the cluster.}
}
\value{
A \code{double} indicating the theoretical allan variance for the drift
process.
}
\description{
This function allows us to calculate the theoretical allan variance for
drift process.
}
\note{
This function is based on the calculation of the theoretical allan variance
for drift process raised in "Allan Variance of Time Series Models for
Measurement Data" by Nien Fan Zhang, 2008, Metrologia, 45(5): 549.
This calculation is fundamental and necessary for the study in "A Study of the Allan Variance
for Constant-Mean Non-Stationary Processes" by Xu et al., 2017, IEEE Signal Processing Letters, 24(8): 1257–1260.
}
\examples{
av1 = av_dr(delta = 1, n = 5)
av2 = av_dr(delta = 2, n = 8)
}
