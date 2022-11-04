#' Creating bootstrap samples from observed data and computing MLEs
#'
#' This function creates bootstrap samples of input data and fits a dynamic mixture via
#' standard maximum likelihood.
#' @param y numerical vector: observed random sample from the mixture.
#' @param nboot number of bootstrap replications
#' @return A list with the following elements:
#'
#' MLE: maximum likelihood estimates obtained from each bootstrap sample.
#'
#' errors: number of times the MLE algorithm breaks down.
#' @details MLEs are computed by means of the optim function. When it breaks
#' down, the sample is discarded and a new one is generated. The function keeps
#' track of the number of times this happens.
#' @keywords dynamic mixture; MLE; non-parametric bootstrap.
#' @export
#' @examples
#' bootMLEs <- Bootsamples(TN2016,100)

Bootsamples = function(y,nboot)
{
  n = length(y)
  MLE = matrix(0,nboot,6)
  i = 1
  j = 0
  while (i <= nboot)
  {
    tryCatch({
      yboot = sample(y,n,replace=TRUE)
      mediana = median(yboot)
      y1 = yboot[yboot<mediana]
      y2 = yboot[yboot>=mediana]
      mu0 = MASS::fitdistr(y1,'lognormal')$estimate[1]
      sigma0 = MASS::fitdistr(y1,'lognormal')$estimate[2]
      xi0 = evir::gpd(y2,mediana)$par.ests['xi']
      beta0 = evir::gpd(y2,mediana)$par.ests['beta']
      muc0 = quantile(yboot,.25)
      tau0 = log(sd(yboot)/2)
      x0Lik = as.numeric(c(muc0,tau0,mu0,sigma0,xi0,beta0))
      res <- optim(x0Lik,dynloglik, gr=NULL,yboot,method='L-BFGS-B',lower=c(-Inf,.01,-Inf,.05,10^-10,.1),upper=c(Inf,Inf,Inf,10,Inf,100),control=list(fnscale=-1))
      MLE[i,] <- res$par # muc, tau, mu, sigma, xi, beta
    },
    error = function(e) {'errore'}) -> condizione
    if(class(condizione) != "character") i = i + 1
    if(class(condizione) == "character") j = j + 1
    }
  return(list(MLE=MLE,errors=j))
}
