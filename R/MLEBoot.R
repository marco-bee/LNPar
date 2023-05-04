#' Bootstrap standard errors for the estimators of a lognormal-Pareto mixture
#'
#' This function computes non-parametric bootstrap standard errors for the estimators of a lognormal-Pareto mixture distribution. The implementation is based on parallel computing.
#' @param y numerical vector: random sample from the mixture.
#' @param minRank integer: minimum possible rank of the threshold.
#' @param p0 (0<p0<1): starting value of the mixing weight.
#' @param alpha0 non-negative scalar: starting value of the Pareto shape parameter.
#' @param mu0 scalar: starting value of the expectation of the lognormal distribution on the log scale.
#' @param sigma0 non-negative scalar: starting value of the standard deviation of the lognormal distribution on the log scale.
#' @return Bootstrap standard errors of the estimators.
#' @keywords mixture; profile likelihood; EM algorithm.
#' @details At each bootstrap replication, the mixture is estimated with thresholds equal to ys(minRank), ys(minRank+1),..., ys(n),
#' where n is the sample size and ys is the sample in ascending order. The function is typically called by LPfit (see the examples below).
#' @export
#' @examples
#' resBoot <- MLEBootPar(y,100,20,.5,1.5,0,1)
#'
#' # Typical use via LPfit
#'
#' resFit <- LPfit(y,90,500)
#' parsStd <- resFit$bootstd
#' @references  Bee, M. (2022), “On discriminating between lognormal and Pareto tail: a mixture-based approach”,
#' Advances in Data Analysis and Classification, https://doi.org/10.1007/s11634-022-00497-4

MLEBoot = function(x,y,minRank,p0,alpha0,mu0,Psi0)
{
  samSiz <- length(y)
  indici = sample(samSiz, samSiz, replace = TRUE)
  yboot = sort(y[indici])
  th <- yboot[minRank:(samSiz-1)]
  nthresh <- length(th)
  resMat <- matrix(0,nthresh,5)
  paretoObs <- cbind(th,matrix(0,nthresh,2))
  for (i in 1:nthresh)
  {
    a <- th[i]
    Res <- par_logn_mix_known(yboot, p0, a, alpha0, mu0, sqrt(Psi0))
    resMat[i,] <- c(Res$prior,Res$alpha,Res$mu,Res$sigma,Res$loglik)
  }
  indice <- which.max(resMat[,5])
  xminhat <- th[indice]
  temp <- par_logn_mix_known(y, p0, xminhat, alpha0, mean(y), sd(y))
  resBest <- c(temp$prior,temp$alpha,temp$mu,temp$sigma,temp$loglik)
  results <- list(res=resBest)
}
