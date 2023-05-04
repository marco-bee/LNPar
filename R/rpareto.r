#' Random number generation for a Pareto r.v.
#'
#' This function simulates random numbers for a Pareto r.v.
#' @param n positive integer: number of simulated random numbers.
#' @param xmin positive scalar: Pareto scale parameter.
#' @param alpha non-negative scalar: Pareto shape parameter.
#' @return n iid random numbers from the Pareto distribution.
#' @keywords random number simulation.
#' @export
#' @examples
#' ySim <- rpareto(5,4,1.5)

rpareto <- function(n, xmin, alpha)
{
x = xmin * (1 - runif(n))^(-1/alpha)
return(x)
}
