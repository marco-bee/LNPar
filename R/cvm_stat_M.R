#' Computing the Cramér - von Mises distance between two samples
#'
#' This function Computing the Cramér - von Mises distance between two samples.
#' @param vec1 (n x 1) vector: first sample.
#' @param vec2 (n x 1) vector: second sample.
#' @param power power of the integrand. Equal to 2 when using the L2 distance
#' @return out scalar: Cramér - von Mises distance between the two samples
#' @details This function computes the Cramér - von Mises distance
#' between two samples. See Drovandi and Frazier (2022, p. 7).
#' @keywords dynamic mixture; simulation.
#' @export
#' @examples
#' out = cvm_stat_M(runif(100),rnorm(100),2)
#' @references{
#'   \insertRef{dro22}{LNPar}
#' }
#'
#'
#' @importFrom Rdpack reprompt

cvm_stat_M <- function(vec1,vec2,power)
{
# compute CvM distance

n1 = length(vec1)
n2 = length(vec2)
n = n1+n2
joint_sample = c(vec1,vec2)
ee = c((1/n1)*rep(1,n1),rep(0,n2))
ff = c(rep(0,n1),(1/n2)*rep(1,n2))
tempp = sort(joint_sample,index.return=TRUE);
temp = tempp$x
ind = tempp$ix
d = joint_sample[ind]
e = ee[ind]
f = ff[ind]
out = 0
Ecur = 0
Fcur = 0
height = 0
for (i in 1:(n-1))
{
  Ecur = Ecur + e[i]
  Fcur = Fcur + f[i]
  height = abs(Fcur-Ecur)
  if (d[i] != d[i+1])
    out = out + height * power
}
return(out)
}
