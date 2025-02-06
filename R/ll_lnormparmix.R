#' Log-likelihood with respect to xmin
#'
#' This function evaluates the log-likelihood function with respect to xmin for a mixture of a lognormal and a Pareto r.v.,
#' assuming to know the numerical values of all the other parameters.
#' @param x positive scalar: value of xmin where the function is evaluated.
#' @param pi scalar, 0 < pi < 1: mixing weight.
#' @param mu scalar: expected value of the lognormal distribution on the log scale.
#' @param sigma positive scalar: standard deviation of the lognormal distribution on the log scale.
#' @param alpha non-negative scalar: Pareto shape parameter.
#' @param y (nx1) vector: random sample from the mixture.
#' @return ll numerical value of the log-likelihood function.
#' @keywords mixture.
#' @export
#' @examples
#' y <- rLnormParMix(100,.5,0,1,4,1.5)
#' llMix <- ll_lnormparmix(x,pi,mu,sigma,alpha,y(3,.5,0,1,1.5,y)

ll_lnormparmix <- function(x,pi,mu,sigma,alpha,y)
{
  ll = sum(log(dLnormParMix(y,pi,mu,sigma,x,alpha)))
  return(ll)
}
