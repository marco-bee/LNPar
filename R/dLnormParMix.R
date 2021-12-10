#' density of a mixture of a lognormal and a Pareto r.v.   
#'
#' This function computes the density of a mixture of a lognormal and a Pareto r.v.
#' @param x numerical vector: values where the density has to be evaluated. 
#' @param p scalar, 0 < p < 1: mixing weight.
#' @param mu scalar: expected value of the lognormal distribution on the log scale.
#' @param sigma scalar: standard deviation of the lognormal distribution on the log scale.
#' @param xmin scalar: threshold.
#' @param alpha non-negative scalar: Pareto shape parameter.
#' @return Density of the lognormal-Pareto distribution evaluated at x.
#' @keywords mixture.
#' @export
#' @examples
#' mixDens <- dLnormParMix(5,.5,0,1,4,1.5)

dLnormParMix <- function(x,pi,mu,sigma,xmin,alpha)
  
{
  source('dpareto.r')
  y1 <- dlnorm(x,mu,sigma)
  y2 <- dpareto(x,xmin, alpha)
  y <- pi * y1 + (1-pi) * y2
  return(y)
}