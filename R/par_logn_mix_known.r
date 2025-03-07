#' Estimate the parameters of a lognormal-Pareto density, assuming a known threshold
#'
#' This function estimates the parameters of a Pareto and a lognormal density, assuming a known threshold.
#' @param y non-negative numerical vector: random sample from the mixture.
#' @param prior1 scalar (0<prior1<1): starting value of the prior probability.
#' @param th positive scalar: threshold.
#' @param alpha non-negative scalar: starting value of the Pareto shape parameter.
#' @param mu scalar: starting value of the lognormal parameter mu.
#' @param sigma positive scalar: starting value of the lognormal parameter sigma.
#' @return A list with the following elements:
#'
#' xmin: estimated threshold.
#'
#' prior: estimated mixing weight.
#'
#' post: matrix of posterior probabilities.
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
#' @export
#' @examples
#' mixFit <- par_logn_mix_known(TN2016, .5, 4700, 3, 7, 1.2)

par_logn_mix_known <- function(y, prior1, th, alpha, mu, sigma)
{
y <- sort(y)
N <- length(y)                                     # number of observations
eps <- 1e-10                                   # convergence criterion
change <- 1                           # initial test value for convergence
maxiter <- 1000                             # maximum number of iterations
nit <- 0                                    # initialize iteration counter
param <- c(prior1, alpha, mu, sigma)            # arrange all parameters in a vector
post_p <- matrix(0, N, 2)		# open matrix for posterior probabilities
y1 <- y[y<th]
y2 <- y[y>=th]
prior <- c(prior1,1-prior1)
n1 <- length(y1)
n2 <- length(y2)
f2 <- rep(0,N)
f1 <- rep(0,N)

if (n1==0) # handle separately the case of no observations below the threshold
{
  while (change > eps && nit <= maxiter)             # start iterations
  {
    parold <- param                               # store old parameter values
    f1 <- dlnorm(y,mu,sigma)    # evaluate
    f2 <- alpha * (th^alpha) / (y2^(alpha + 1))
    f <- prior[1] * f1 + prior[2] * f2       # component densities and mixture density
    post_p[,1] <- (prior[1] * f1) / f
    post_p[,2] <- (prior[2] * f2) / f
    prior <- colMeans(post_p)              # M step: prior probabilities
    if (is.nan(prior[1]) == T || prior[1] == 0 || is.nan(prior[2]) == T)
    {
      alpha <- N / (sum(log(y/th)))				# M step: alpha
      loglik <- sum(log(alpha * (th^alpha) / (y2^(alpha + 1))))                  # evaluate log-likelihood function
      change = 0
      break
    }
    if (prior[1] > 0)
    {
      alpha <- N*prior[2] / (sum(t(post_p[,2]) %*% log(y2/th)))				# M step: alpha
      mu <- (1/(N*prior[1])) * (log(y) %*% post_p[,1])
      sigma <- sqrt((1/(N*prior[1])) * log(y)^2 %*% post_p[,1] - mu^2)
      f1 <- dlnorm(y,mu,sigma)    # evaluate
      f2 <- alpha * (th^alpha) / (y2^(alpha + 1))
      f <- prior[1] * f1 + prior[2] * f2       # component densities and mixture density
      loglik <- sum(log(f))                  # evaluate log-likelihood function
      param <- c(prior[1], alpha, mu, sigma)              # arrange all parameters in a vector
      change <- max(abs(param - parold))           # test value for convergence
      nit <- nit + 1                               # increase iteration counter
    }
  }
  if (nit >= maxiter)
  {
    cat(c("algorithm did not converge in ", maxiter,  " iterations","\n"))   # warning message if convergence is not reached
  }
}

if (n1>0)
{
  while (change > eps && nit <= maxiter)             # start iterations
  {
    parold <- param                               # store old parameter values
    f1 <- dlnorm(y,mu,sigma)    # evaluate
    f2[(n1+1):N] <- alpha * (th^alpha) / (y2^(alpha + 1))  # component
    f <- prior[1] * f1 + prior[2] * f2       # densities and mixture density
    post_p[1:n1,1] <- 1         # E step
    post_p[(n1+1):N,1] <- (prior[1] * f1[(n1+1):N]) / f[(n1+1):N]
    post_p[(n1+1):N,2] <- (prior[2] * f2[(n1+1):N]) / f[(n1+1):N]
    prior <- colMeans(post_p)              # M step: prior probabilities
    if (is.nan(prior[1]) == T || prior[1] == 1 || is.nan(prior[2]) == T)
    {
      if (is.nan(prior[1]) == T)
        prior[1] = NA
      alpha <- NA				# M step: alpha
      mu <- mean(log(y))
      sigma <- sd(log(y))
      loglik <- sum(log(dlnorm(y,mu,sigma)))                  # evaluate log-likelihood function
      change = 0
      break
    }

    if (prior[1] < 1)
    {
      alpha <- N*prior[2] / (sum(t(post_p[(n1+1):N,2]) %*% log(y2/th)))				# M step: alpha
      mu <- (1/(N*prior[1])) * (log(y) %*% post_p[,1])
      sigma <- sqrt((1/(N*prior[1])) * log(y)^2 %*% post_p[,1] - mu^2)
      f1 <- dlnorm(y,mu,sigma)    # evaluate
      f2[(n1+1):N] <- alpha * (th^alpha) / (y2^(alpha + 1))  # component
      f <- prior[1] * f1 + prior[2] * f2       # densities and mixture density
      loglik <- sum(log(f))                  # evaluate log-likelihood function
      param <- c(prior[1], alpha, mu, sigma)              # arrange all parameters in a vector
    }
    change <- max(abs(param - parold))           # test value for convergence
    nit <- nit + 1                               # increase iteration counter
  }
  if (nit >= maxiter)
  {
    cat(c("algorithm did not converge in ", maxiter,  " iterations","\n"))   # warning message if convergence is not reached
  }
}
results <- list("prior" = prior[1], "post" = post_p, "alpha " = alpha, "mu " = mu, "sigma " = sigma, "loglik" = loglik, "nit" = nit)
results
}
