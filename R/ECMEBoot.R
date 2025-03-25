#' Bootstrap standard errors for the MLEs of a lognormal-Pareto mixture
#'
#' This function draws a bootstrap sample and uses it to estimate the parameters of a lognormal-Pareto mixture distribution. Since this is typically called by LPfitEM, see the help of LPfitEM for examples.
#' @param x list: sequence of integers 1,...,K, where K is the mumber of datasets. Set x = 1 in case
#' of a single dataset.
#' @param y numerical vector: observed sample.
#' @param eps non-negative scalar: starting value of the log-expectation of the lognormal distribution on the log scale.
#' @param maxiter non-negative integer: starting value of the log-variance of the lognormal distribution on the log scale.
#' @return Estimated parameters obtained from a bootstrap sample.
#' @keywords mixture; ECME algorithm.
#' @details At each bootstrap replication, the mixture is estimated via the ECME algorithm. The function is typically called by LPfitEM.
#' @export

ECMEBoot = function(x,y,eps,maxiter)
{
  samSiz <- length(y)
  indici = sample(samSiz, samSiz, replace = TRUE)
  yboot = sort(y[indici])
  temp <- LPfitEM(yboot,eps,maxiter)
  resBest <- temp$pars
  results <- list(res=resBest)
}
