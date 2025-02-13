#' ECME-based testing for a Pareto tail
#'
#' This function draws a bootstrap sample from the null (lognormal) distribution and computes the test for the null hypothesis of
#' a pure lognormal distribution versus the alternative of a lognormal-Pareto
#' mixture, where the parameters of the latter are estimated by means of the
#' ECME algorithm. To be only called from ParallelTestEM.
#' @param x list: sequence of integers 1,...,K, where K is the mumber of datasets. Set x = 1 in case
#' of a single dataset.
#' @param n sample size.
#' @param muNull log-expectation value under the null hypothesis.
#' @param sigmaNull log-standard deviation under the null hypothesis.
#' @return A list with the following elements:
#'
#' LR: observed value of the llr test.
#' @keywords mixture; profile likelihood; log-likelihood ratio test.
#' @export
#' @import stats
#' @examples
#' n = 100
#' muNull = mean(log(TN2016))
#' sigmaNull = sd(log(TN2016))
#' res = LPtestEM(1,n,muNull,sigmaNull)

LPtestEM <- function(x,n,muNull,sigmaNull)
{
  ysim <- sort(rlnorm(n,muNull,sigmaNull))
  est0 <- c(mean(log(ysim)),sd(log(ysim)))
  ell0 <- sum(log(dlnorm(ysim,est0[1],est0[2])))
  temp1 <- LPfitEM(ysim,eps=1e-12,maxiter=1000)
  ell1 <- temp1$loglik
  LR <- pmax(0,2*(ell1-ell0)) # sometimes slightly negative, as ell1 is computed numerically
  results <- list(res=LR)
  results
}
