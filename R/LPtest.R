#' Testing for a Pareto tail
#'
#' This function allows you to compute the bootstrap test for the null hypothesis of a pure lognormal distribution versus the alternative of a lognormal-Pareto mixture.
#' @param n sample size
#' @param muNull lognormal expected value under the null hypothesis.
#' @param sigmaNull lognormal standard deviation under the null hypothesis.
#' @param minRank minimum possible rank of the threshold.
#' @param obsTest value of the test statistics computed with the data under analysis.
#' @param nboot number of bootstrap replications.
#' @return A list with the following elements:
#' 
#' lr: simulated values of the llr test under the null hypothesis.
#' 
#' pval: p-value of the test.
#' @keywords mixture; profile likelihood; log-likelihood ratio test.
#' @export
#' @examples
#' resTest <- LPtest(100,0,1,0.67,100,90)

LPtest <- function(n,muNull,sigmaNull,obsTest,nboot,minRank)
{
lr <- rep(0,nboot)
for (i in 1:nboot)
{
  ysim <- sort(rlnorm(n,muNull,sigmaNull))
  th <- ysim[minRank:(n-1)]
  est0 <- c(mean(log(ysim)),sd(log(ysim)))
  ell0 <- sum(log(dlnorm(ysim,est0[1],est0[2])))
  xmin0 <- quantile(ysim,.45)
  # p0 <- length(ysim[ysim<xmin0])/n
  p0 <- .75
  #  alpha0 <- length(ysim[ysim>xmin0]) / (sum(log(ysim[ysim>xmin0]/xmin0)))
  alpha0 <- 3
  mu0 <- mean(log(ysim))+1
  sigma0 <- sd(log(ysim))
  resMat <- matrix(0,n-minRank,5)
  
  for (s in minRank:(n-1))
  {
    a <- ysim[s]
    Res <- par_logn_mix_known(ysim, p0, a, alpha0, mu0, sigma0)
    resMat[s-minRank+1,] <- c(Res$prior,Res$alpha,Res$mu,Res$sigma,Res$loglik)
    temp <- Res$post[,1]
    post1 <- temp[temp<1]
  }
  indice <- which.max(resMat[,5]) # rank (starting form the smallest observations) of xminhat (corresponding to the largest likelihood)
  xminhat <- th[indice]
  temp1 <- par_logn_mix_known(ysim, p0, xminhat, alpha0, mu0, sigma0)
  ell1 <- temp1$loglik
  lr[i] <- pmax(0,2*(ell1-ell0)) # sometimes slightly negative, as ell1 is computed numerically
}
obsp <- length(lr[lr>obsTest])/nboot
results <- list("lr" = lr, "pval" = obsp)
results
}