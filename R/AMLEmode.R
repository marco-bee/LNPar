#' Approximating the mode of a multivariate empirical distribution
#'
#' This function approximates the mode of a multivariate empirical distribution by means of:
#' 1. the sample mean
#' 2. the product of the maxima of the univariate kernel densities estimated using the marginals
#' 3. the maximum of the multivariate kernel density
#' 4. the maximum of the product of the univariate kernel densities
#'
#' @param ABCsam (m x k) matrix: ABC sample, where m is the ABC sample size and k is the
#' number of parameters.
#' @return A list containing the 4 approximate modes.
#' @details The bandwidth is estimated via smoothed cross-validation
#' @keywords dynamic mixture; approximate maximum likelihood.
#' @export
#' @examples
#' res <- AMLEmode(ABCsam)

AMLEmode <- function(ABCsam)
{
  require(ks)
  N = nrow(ABCsam) # ABC sample size
  qq = ncol(ABCsam)

  # 1. sample mean

  AMLE1 = colMeans(ABCsam)

  # 2. univariate kd

  AMLE2 = rep(0,qq)
  for (i in 1:qq)
  {
    vettore = ABCsam[,i]
    bw = hscv(vettore)
    f_eps = kde(vettore,bw,eval.points=vettore)
    AMLE2[i] = vettore[which.max(f_eps$estimate)]
  }

  # 3. multivariate kde

  bw = Hscv(ABCsam)
  f_eps = kde(ABCsam,bw,eval.points=ABCsam)
  AMLE3 = ABCsam[which.max(f_eps$estimate),]

  # 4. maximum of the product of the univariate kernel densities

  f_eps_g = matrix(0,N,qq)
  for (i in 1:qq)
  {
    vettore = ABCsam[,i]
    bw = hscv(vettore)
    f_eps = kde(vettore,bw,eval.points=vettore)
    f_eps_g[,i] = f_eps$estimate
  }
  F_eps = apply(f_eps_g,1,prod)
  AMLE4 = ABCsam[which.max(F_eps),]
  res = list(SampleMean = AMLE1, UnivariateKD = AMLE2, MultivariateKD = AMLE3, ProductUnivariateKD = AMLE4)
  return(res)
  }
