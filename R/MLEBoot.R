#' Bootstrap standard errors for the estimators of a lognormal-Pareto mixture
#'
#' This function computes non-parametric bootstrap standard errors for the estimators of a lognormal-Pareto mixture distribution.
#' @param y numerical vector: random sample from the mixture.
#' @param nboot number of bootstrap replications.
#' @param nthresh minimum possible rank of the threshold.
#' @param p0 starting value of the mixing weight.
#' @param alpha0 starting value of the Pareto shape parameter.
#' @param mu0 starting value of the expectation of the lognormal distribution on the log scale.
#' @param mu0 starting value of the standard deviation of the lognormal distribution on the log scale.
#' @return bootstrap standard errors of the estimators.
#' @keywords mixture; profile likelihood; EM algorithm.
#' @details At each bootstrap replication, the mixture is estimated with thresholds equal to ys(n-nthresh), ys(n-nthresh+1),..., ys(n),
#' where n is the sample size and ys is the sample in ascending order. The function is typically called by LPfit (see the examples below).
#' @export
#' @examples
#' resBoot <- MLEBoot(y,100,20,.5,1.5,0,1)
#' 
#' # Now call MLEBoot from LPfit
#' 
#' resFit <- LPfit(y,90,500)
#' parsStd <- resFit$bootstd
#' @references  Bee, M. (2020), “On discriminating between lognormal and Pareto tail: a mixture-based approach”, DEM working paper 2020/9. https://www.economia.unitn.it/alfresco/download/workspace/SpacesStore/06b9cd83-9754-4b9e-8e7c-7954d43395b7/DEM2020_09.pdf

MLEBoot = function(y,nboot,nthresh,p0,alpha0,mu0,Psi0)
{
  require(tcltk)
  pb <- tkProgressBar("test progress bar", "Some information in %", 0, nboot, 50)
  samSiz <- length(y)
  indice = rep(0,nboot)
  xminhat = rep(0,nboot)
  resBest = matrix(0,nboot,5)
  for (j in 1:nboot)
  {
    indici = sample(samSiz, samSiz, replace = TRUE)
    yboot = sort(y[indici])
    th <- yboot[(n-nthresh):n] 
    nthresh <- length(th)
    resMat <- matrix(0,nthresh,5)
    paretoObs <- cbind(th,matrix(0,nthresh,2))
    for (i in 1:nthresh)
    {
      a <- th[i]
      Res <- par_logn_mix_known(yboot, p0, a, alpha0, mu0, sqrt(Psi0))
      resMat[i,] <- c(Res$prior,Res$alpha,Res$mu,Res$sigma,Res$loglik)
    }
    indice[j] <- which.max(resMat[,5])
    xminhat[j] <- th[indice[j]]
    temp <- par_logn_mix_known(y, p0, xminhat[j], alpha0, mean(y), sd(y))
    resBest[j,] <- c(temp$prior,temp$alpha,temp$mu,temp$sigma,temp$loglik)
    write(c(resBest[j,],indice[j],xminhat[j]),paste('pars_boot.txt',sep=''),ncolumns = 7, append=T)
    info <- sprintf("%d%% done", round(j/(nboot/100)))
    setTkProgressBar(pb, j, sprintf("test (%s)", info), info)
  }
  close(pb)
  varcov = cov(resBest[,1:4])
  stddev = sqrt(diag(varcov))
  results <- list(std=stddev)
}