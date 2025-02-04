#' Estimating a lognormal-Pareto mixture via the ECME algorithm
#'
#' This function fits a lognormal-Pareto mixture by means of the ECME algorithm.
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
#' @details Estimation of a lognormal-Pareto mixture via the ECME algorithm.
#' @keywords mixture; ECME algorithm.
#' @export
#' @examples
#' mixFit <- LPfitEM(TN2016,0)
#'
#'
#' @importFrom Rdpack reprompt

LPfitEM <- function(y,eps,maxiter,qxmin0)
{
  pars <- matrix(0,maxiter,5)
  loglik <- rep(0,maxiter)
  ys <- sort(y)
  n <- length(ys)

  # initial values

  xmin0 <- quantile(ys,qxmin0)
  p0 <- length(ys[ys<xmin0])/n
  alpha0 <- length(ys[ys>xmin0]) / (sum(log(ys[ys>xmin0]/xmin0)))
  mu0 <- mean(log(ys))+1
  sigma0 <- sd(log(ys))
  parold = c(p0,alpha0,mu0,sigma0,xmin0)
  post_p <- matrix(0,n,2)		# open matrix for posterior probabilities
  y1 <- ys[ys<=xmin0]
  y2 <- ys[ys>xmin0]
  f2 <- rep(0,n)
  f1 <- rep(0,n)
  change = 100
  xmin = xmin0
  mu <- mu0
  sigma <- sigma0
  alpha <- alpha0
  p <- p0
  nit <- 1

  while (change > eps && nit <= maxiter)
  {
    y1 <- ys[ys<xmin]
    y2 <- ys[ys>=xmin]
    n1 <- length(y1)
    n2 <- length(y2)
    f1 <- rep(0,n)
    f2 <- rep(0,n)
    post_p <- matrix(0,n,2)
    f1 <- dlnorm(ys,mu,sigma)
    f2[(n1+1):n] <- alpha * (xmin/y2)^alpha * (1/y2)
    f <- p * f1 + (1-p) * f2
    post_p[1:n1,1] <- 1
    post_p[1:n1,2] <- 0
    post_p[(n1+1):n,1] <- (p * f1[(n1+1):n]) / f[(n1+1):n]
    post_p[(n1+1):n,2] <- ((1-p) * f2[(n1+1):n]) / f[(n1+1):n]
    if (p>=1-1/n)
    {
      p = 1
      alpha <- NA
      xmin = NA
      mu <- mean(log(ys))
      sigma <- sqrt(((n - 1)/n)) * sd(log(ys))
      parsBestNA <- c(p[1], alpha, mu, sigma, xmin)
      loglik <- sum(log(dlnorm(ys, mu, sigma)))
      post_p[,1] <- 1
      post_p[,2] <- 0
      change = 0
      break
    }

    p <- mean(post_p[,1])              # CM step: prior probabilities
    alpha <- n*(1-p) / (sum(t(post_p[(n1+1):n,2]) %*% log(y2/xmin)))
    mu <- (1/(n*p)) * (log(ys) %*% post_p[,1])
    sigma <- sqrt((1/(n*p)) * log(ys)^2 %*% post_p[,1] - mu^2)

    # CM step 2

    xmin <- optimize(ll_lnormparmix,c(0,ys),p,mu,sigma,alpha,ys,maximum=TRUE)$maximum

    pars[nit,] <- c(p, alpha, mu, sigma, xmin)
    loglik[nit] <- ll_lnormparmix(xmin,p,mu,sigma,alpha,ys)

    change = max(abs(pars[nit,]-parold))
    parold = pars[nit,]
    if (change > eps && nit == maxiter)
    {
      indice = which.max(loglik)
      parsb = pars[indice,]
      parsBestNA <- parsb
      rmin = parsb[5]
      max_loglik = loglik[indice]
      rankECME = length(ys[ys<=rmin])
    }
    else
    {
      parsBestNA <- pars[nit,]
      rmin = pars[5]
      max_loglik = loglik[nit]
      rankECME = length(ys[ys<=rmin])
    }
    nit = nit + 1
  }
  out <- list(pars = parsBestNA , loglik = max_loglik, rankECME = rankECME,
              niter = nit - 1, post_p = post_p)
  return(out)
}

