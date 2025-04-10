% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LPfitEM.R
\name{LPfitEM}
\alias{LPfitEM}
\title{Estimating a lognormal-Pareto mixture via the ECME algorithm}
\usage{
LPfitEM(y, eps, maxiter, qxmin0 = 0.5, nboot = 0)
}
\arguments{
\item{y}{numerical vector: random sample from the mixture.}

\item{eps}{non-negative scalar: tolerance for the stopping rule.}

\item{maxiter}{non-negative integer: maximum number of iterations of the ECME algorithm.}

\item{qxmin0}{scalar, 0 < qxmin0 < 1: quantile level used for determining the starting value of xmin. Defaults to 0.5.}

\item{nboot}{non-negative integer: number of bootstrap replications used for estimating the standard errors. If omitted, no standard errors are computed.}
}
\value{
A list with the following elements:

pars: estimated parameters (p, alpha, mu, sigma, xmin).

loglik: maximized log-likelihood.

thRank: estimated rank of xmin.

niter: number of iterations.

npareto: estimated number of Pareto observations.

postProb: matrix of posterior probabilities.

bootEst: matrix of estimated parameters at each bootstrap replication.

bootstd: bootstrap standard errors of the estimators.
}
\description{
This function fits a lognormal-Pareto mixture by means of the ECME algorithm.
}
\details{
Estimation of a lognormal-Pareto mixture via the ECME algorithm. Standard errors are computed via non-parametric bootstrap.
}
\examples{
ysim <- rLnormParMix(100,.9,0,1,5,1)
mixFit <- LPfitEM(ysim,eps=1e-10,maxiter=1000,nboot=0)


}
