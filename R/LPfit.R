#' Estimating a lognormal-Pareto mixture by maximizing the profile log-likelihood
#'
#' This function fits a lognormal-Pareto mixture by maximizing the profile log-likelihood.
#' @param y numerical vector: random sample from the mixture.
#' @param minRank integer: minimum possible rank of the threshold.
#' @param nboot number of bootstrap replications used for estimating the standard errors. If omitted, no standard errors are computed.
#' @return A list with the following elements:
#'
#' xmin: estimated threshold.
#'
#' prior: estimated mixing weight.
#'
#' postProb: matrix of posterior probabilities.
#'
#' alpha: estimated Pareto shape parameter.
#'
#' mu: estimated expectation of the lognormal distribution on the lognormal scale.
#'
#' sigma: estimated standard deviation of the lognormal distribution on the lognormal scale.
#'
#' loglik: maximized log-likelihood.
#'
#' nit: number of iterations.
#'
#' npareto: estimated number of Pareto observations.
#'
#' bootstd: bootstrap standard errors of the estimators.
#' @details Estimation is implemented as in Bee (2022). As of standard errors, at each bootstrap replication the mixture is estimated with thresholds equal to ys(minRank), ys(minRank+1),..., ys(n),
#' where n is the sample size and ys is the sample sorted in ascending order. The latter procedure is implemented via parallel computing.
#' If the algorithm does not converge in 1000 iterations, a message is displayed.
#' @keywords mixture; profile likelihood.
#' @export
#' @examples
#' mixFit <- LPfit(TN2016,90,0)
#' @references{
#'   \insertRef{bee22}{LNPar}
#' }
#'
#'
#' @importFrom Rdpack reprompt

LPfit <- function(y,minRank,nboot)
{
  ys <- sort(y)
  n <- length(ys)

  # initial values

  a0 <- ys[length(ys)-minRank]
  p0 <- length(ys[ys<a0])/n
  alpha0 <- length(ys[ys>a0]) / (sum(log(ys[ys>a0]/a0)))
  mu0 <- mean(log(ys))+1
  Psi0 <- var(log(ys))
#  th <- ys
  th <- ys[minRank:(n-1)]
  nthresh <- length(th)
  resMat <- matrix(0,nthresh,5)
  paretoObs <- cbind(th,matrix(0,nthresh,2))
  # for (i in minRank:(n-1))
  for (i in 1:nthresh)
  {
    a <- th[i]
    Res <- par_logn_mix_known(ys, p0, a, alpha0, mean(ys), sd(ys))
    resMat[i,] <- c(Res$prior,Res$alpha,Res$mu,Res$sigma,Res$loglik)
  }
  indice <- which.max(resMat[,5])
  xminhat <- th[indice]
  # xminhat <- ys[indice+minRank]
  resBest <- par_logn_mix_known(ys, p0, xminhat, alpha0, mean(ys), sd(ys))
  npareto <- n * (1-resBest$prior)
  prior <- resBest$prior
  postProb <- resBest$post
  alpha <- resBest$alpha
  mu <- resBest$mu
  sigma <- resBest$sigma
  loglik <- resBest$loglik
  nit <- resBest$nit
  if (nboot==0)
  {
    results <- list(xmin=xminhat,prior=prior,postProb=postProb,alpha=alpha,mu=as.double(mu),sigma=as.vector(sigma),loglik=loglik,nit=nit,npareto=npareto)
    return(results)
  }
  else
  {
    nreps.list <- sapply(1:nboot, list)
    chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
    if (nzchar(chk) && chk == "TRUE") {
      n.cores <- 2L
    } else {
      n.cores <- parallel::detectCores()
    }
    clust <- parallel::makeCluster(n.cores)
    BootMat = matrix(0,nboot,5)
    temp <- parallel::parLapply(clust,nreps.list, MLEBoot,ys,minRank,p0,alpha0,mean(ys),var(ys))
    parallel::stopCluster(cl=clust)
    for (i in 1:nboot)
    {
      BootMat[i,] = as.vector(unlist(temp[[i]]))
    }
    varcov = cov(BootMat)
    stddev = sqrt(diag(varcov))
    results <- list(xmin=xminhat,prior=prior,postProb=postProb,alpha=alpha,mu=as.double(mu),sigma=as.vector(sigma),loglik=loglik,nit=nit,npareto=npareto,bootEst=BootMat[,1:4],bootStd=stddev)
    return(results)
  }
}
