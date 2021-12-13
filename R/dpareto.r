#' density of a Pareto r.v.
#'
#' This function evaluates the density of a Pareto r.v.
#' @param x numerical vector: values where the density has to be evaluated.
#' @param xmin scalar: Pareto scale parameter.
#' @param alpha non-negative scalar: Pareto shape parameter.
#' @return Density of the Pareto distribution evaluated at x.
#' @keywords density.
#' @export
#' @examples
#' parDens <- dPareto(5,4,1.5)

dpareto <- function(x, xmin, alpha)
{
y = matrix(0,length(x),1)
x1 = x[x>=xmin]
y[x<xmin] = 0
y[x>=xmin] = alpha * (xmin^alpha) / (x1^(alpha+1))
return(y)
}
