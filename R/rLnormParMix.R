#' Random number simulation for a mixture of a lognormal and a Pareto r.v.
#'
#' This function simulates random numbers for a mixture of a lognormal and a Pareto r.v.
#' @param n positive integer: number of simulated random numbers.
#' @param pi scalar, 0 < pi < 1: mixing weight.
#' @param mu scalar: expected value of the lognormal distribution on the log scale.
#' @param sigma positive scalar: standard deviation of the lognormal distribution on the log scale.
#' @param xmin positive scalar: threshold.
#' @param alpha non-negative scalar: Pareto shape parameter.
#' @return n iid random numbers from the lognormal-Pareto distribution.
#' @keywords mixture.
#' @export
#' @examples
#' ySim <- rLnormParMix(100,.5,0,1,4,1.5)

rLnormParMix = function(n,pi,mu,sigma,xmin,alpha)
{
  p <- rbinom(n,1,pi)
  y1 <- rlnorm(sum(p),mu,sigma)
  y2 <- rpareto(n-sum(p), xmin, alpha)
  if (sum(p)==n)
    y <- y1
  if (sum(p)==0)
    y <- y2
  y <- c(y1,y2)
  return(y)
}
