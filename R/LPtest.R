#' Profile-based testing for a Pareto tail
#'
#' This function draws a bootstrap sample from the null (lognormal) distribution and computes the test for the null hypothesis of
#' a pure lognormal distribution versus the alternative of a lognormal-Pareto
#' mixture, where the parameters of the latter are estimated via maximum profile
#' likelihood. To be only called from ParallelTest. Estimation unde rthe alternative is perfromed
#'
#' @param x list: sequence of integers 1,...,K, where K is the mumber of datasets. Set x = 1 in case
#' of a single dataset.
#' @param n sample size.
#' @param muNull lognormal expected value under the null hypothesis.
#' @param sigmaNull lognormal standard deviation under the null hypothesis.
#' @param minRank minimum possible rank of the threshold.
#' @return A list with the following elements:
#'
#' LR: observed value of the llr test.
#' @export
#' @import stats
#' @examples
#' n = 100
#' muNull = mean(log(TN2016))
#' sigmaNull = sd(log(TN2016))
#' minRank = 90
#' res = LPtest(1,n,muNull,sigmaNull,minRank)
#' @references{
#'   \insertRef{bee24a}{LNPar}
#' }
#'

LPtest <- function(x,n,muNull,sigmaNull,minRank)
{
  ysim <- sort(rlnorm(n,muNull,sigmaNull))
  th <- ysim[minRank:(n-1)]
  est0 <- c(mean(log(ysim)),sd(log(ysim)))
  ell0 <- sum(log(dlnorm(ysim,est0[1],est0[2])))
  xmin0 <- quantile(ysim,.45)
  p0 <- .75
  alpha0 <- 3
  mu0 <- mean(log(ysim))+1
  sigma0 <- sd(log(ysim))
  th <- ysim[minRank:(n-1)]
  nthresh <- length(th)
  resMat <- matrix(0,nthresh,5)

  for (k in 1:nthresh)
  {
    a <- th[k]
    Res <- par_logn_mix_known(ysim, p0, a, alpha0, mu0, sigma0)
    resMat[k,] <- c(Res$prior,Res$alpha,Res$mu,Res$sigma,Res$loglik)
    temp <- Res$post[,1]
    post1 <- temp[temp<1]
  }
  indice <- which.max(resMat[,5]) # rank (starting form the smallest observations) of xminhat (corresponding to the largest likelihood)
  xminhat <- th[indice]
  temp1 <- par_logn_mix_known(ysim, p0, xminhat, alpha0, mu0, sigma0)
  ell1 <- temp1$loglik
  LR <- pmax(0,2*(ell1-ell0)) # sometimes slightly negative, as ell1 is computed numerically
  results <- list(res=LR)
  results
}
