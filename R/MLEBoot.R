#' Bootstrap standard errors for the estimators of a lognormal-Pareto mixture
#'
#' This function draws a bootstrap sample and uses it to estimate the parameters of a lognormal-Pareto mixture distribution. Since this is typically called by LPfit, see the help of LPfit for examples.
#' @param x list: sequence of integers 1,...,K, where K is the mumber of datasets. Set x = 1 in case
#' of a single dataset.
#' @param y numerical vector: observed sample.
#' @param minRank positive integer: minimum possible rank of the threshold.
#' @param takeOut integer: minimum number of observations above the threshold (see Details).
#' @param p0 (0<p0<1): starting value of the mixing weight.
#' @param alpha0 non-negative scalar: starting value of the Pareto shape parameter.
#' @param mu0 scalar: starting value of the log-expectation of the lognormal distribution on the log scale.
#' @param Psi0 non-negative scalar: starting value of the log-variance of the lognormal distribution on the log scale.
#' @return Estimated parameters obtained from a bootstrap sample.
#' @keywords mixture; profile likelihood; EM algorithm.
#' @details At each bootstrap replication, the mixture is estimated with thresholds equal to ys(minRank), ys(minRank+1),..., ys(n-takeOut),
#' where n is the sample size and ys is the sample in ascending order.
#' If the estimated prior probability is larger than pimin, the prior probability is set equal to 1, and a pure lognormal is estimated via MLE.
#' The function is typically called by LPfit (see the example below).
#' @export
#' @examples
#' mixFit <- LPfit(TN2016,90,5,10)
#' @references  Bee, M. (2022), “On discriminating between lognormal and Pareto tail: a mixture-based approach”,
#' Advances in Data Analysis and Classification, https://doi.org/10.1007/s11634-022-00497-4

MLEBoot = function(x,y,minRank,takeOut,p0,alpha0,mu0,Psi0)
{
  samSiz <- length(y)
  indici = sample(samSiz, samSiz, replace = TRUE)
  yboot = sort(y[indici])
  th <- yboot[minRank:(samSiz-takeOut)]
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
  if (is.nan(temp$prior[1]) == T || temp$prior[1] == 1 || is.nan(temp$prior[2]) == T)
  {
    if (is.nan(temp$prior[1]) == T)
      prior[1] = 1
    alpha <- NA				# M step: alpha
    mu <- mean(log(y))
    sigma <- ((N-1)/N)*sd(log(y))
    loglik <- sum(log(dlnorm(y,mu,sigma)))                  # evaluate log-likelihood function
    resBest <- c(prior,alpha,mu,sigma,loglik)
    change = 0
    break
  }
  results <- list(res=resBest)
}
