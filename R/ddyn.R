#' Density of a Lognormal-GPD dynamic mixture
#'
#' This function evaluates the density of a Lognormal-GPD dynamic mixture.
#' @param x non-negative vector: points where the function is evaluated.
#' @param pars (6 by 1) numerical vector: values of CA1, CA2, meanlog, sdlog, xi, beta.
#' @return density of the lognormal-GPD mixture evaluated at x.
#' @keywords dynamic mixture.
#' @export
#' @examples
#' dLNPar <- ddyn(x,pars)

ddyn <- function(x,pars)
{
library(evir)
library(pracma)
muc <- pars[1]
tau <- pars[2]
mu <- pars[3]
sigma <- pars[4]
xi <- pars[5]
beta <- pars[6]
f <- function(x) (dgpd(x, xi, mu=0, beta)-dlnorm(x,mu,sigma)) * atan((x-muc)/tau)
p <- pcauchy(x,muc,tau)
temp <- (1-p) * dlnorm(x,mu,sigma) + p * dgpd(x, xi, mu=0, beta)
I <- NULL
I1 <- 10
i <- 1
while (abs(I1) > .0001)
{
  temp1 <- pracma::quadinf(f,i-1,i)
  I1 <- temp1$Q
  I <- c(I,I1)
  i <- i + 1
  if (i > 2000)
  {
    print(c(muc,tau,mu,sigma,xi,beta))
    print(c(i,I1))
  }
}
Z <- 1+(1/pi) * sum(I)
f <- temp/Z
return(f)
}
