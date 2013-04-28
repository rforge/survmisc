\name{quantileSurv}
\alias{quantileSurv}
\title{Quantiles for \code{Surv} object}
\usage{
  quantileSurv(s, q = c(25, 50, 75), alpha = 0.05)
}
\arguments{
  \item{s}{A \code{Surv} object}

  \item{q}{Vector of quantiles (expressed as percentage)}

  \item{alpha}{Significance level \eqn{\alpha}{alpha}}
}
\value{
  A matrix with quantile, and upper and lower confidence
  intervals. \cr Intervals are calculated from
  \eqn{\sigma}{sigma} which is: \deqn{\sigma (t) = \sqrt{
  \frac{Var[\hat{S}(t)]}{\hat{S}^2(t)}} }{ ( Var[S(t)] /
  (S(t)^2) )^0.5} The intervals given are: \item{linear}{
  \deqn{ \hat{S}(t) \pm Z_{1- \alpha} \sigma (t)
  \hat{S}(t)}{ S(t)+- Z(1-alpha) sigma(t) S(t)} Where
  \eqn{\hat{S}(t) }{S(t)} is the Kaplan-Meier survival
  estimate. \cr } \item{log transform}{ \deqn{ [
  \hat{S}(t)^{\frac{1}{\theta}}, \hat{S}(t)^{\theta} ] }{
  [S(t)^(1/theta), S(t)^theta]} Where \eqn{\theta}{theta}
  is: \deqn{ \exp{ \frac{Z_{1- \alpha} \sigma (t)}{
  \log{\hat{S}(t)}}} }{ exp ( Z(1-alpha)sigma(t) /
  log(S(t)) )} } \item{arcsin-sqrt}{ Arcsine-square root
  transform. \cr \cr Upper: \deqn{ \sin^2(\max[0,
  \arcsin{\sqrt{\hat{S}(t)}} - \frac{Z_{1-
  \alpha}\sigma(t)}{2} \sqrt{
  \frac{\hat{S}(t)}{1-\hat{S}(t)}}]) }{ sin^2( max[0,
  arcsin S(t)^0.5 - Z(1-alpha)sigma(t)/2
  (S(t)/1-S(t))^0.5])} Lower: \deqn{
  \sin^2(\min[\frac{\pi}{2}, \arcsin{\sqrt{\hat{S}(t)}} +
  \frac{Z_{1- \alpha}\sigma(t)}{2} \sqrt{
  \frac{\hat{S}(t)}{1-\hat{S}(t)}}]) }{ sin^2( min[pi/2,
  arcsin S(t)^0.5 + Z(1-alpha)sigma(t)/2
  (S(t)/1-S(t))^0.5])} }
}
\description{
  Quantiles for \code{Surv} object
}
\examples{
data(bmt, package="KMsurv")
b1 <- bmt[bmt$group==1, ] # ALL patients
s1 <- Surv(time=b1$t2, event=b1$d3)
quantileSurv(s1)
}
\references{
  Strawderman RL, Parzen MI, Wells MT 1993 Accurate
  Confidence Limits for Quantiles under Random Censoring.
  \emph{Biometrics} 1993 \bold{43}(4):1399--415.
  \href{http://www.jstor.org/stable/2533506}{JSTOR}
}

