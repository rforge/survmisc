\name{compNSurv}
\alias{compNSurv}
\title{Compare more than two survival curves}
\usage{
  compNSurv(t, n, e, n1, e1, FHp, FHq, scores = "",
    round1 = 5)
}
\arguments{
  \item{t}{Time}

  \item{n}{No. at risk (overall)}

  \item{e}{No. events (overall)}

  \item{n1}{No. at risk per group}

  \item{e1}{No. of events per group}

  \item{FHp}{\eqn{p} for Fleming-Harrington test}

  \item{FHq}{\eqn{q} for Fleming-Harrington test}

  \item{scores}{scores; if no value is given this is set as
  \eqn{1,2,...,K} the no. of predictors}

  \item{round1}{No. digits to which to round (for display)}
}
\value{
  A list with 3 elements. The first are 2 data frames. In
  each, there is one row for each of the weights, as given
  above. \describe{ \item{lrTests}{log-rank family of
  tests. \cr Columns give the value of the test statistic,
  \eqn{Q}, degrees of freedom (\eqn{= K-1} for \eqn{K}
  groups) and corresponding p value the from chi-squared
  distribution.} \item{trendTests}{Tests for trend. \cr
  Columns give the value of \eqn{Z} and corresponding p
  value from the normal distribution.} \item{scores}{A
  vector of the scores used in calculating the trendTests
  above.} }
}
\description{
  Compare more than two survival curves. Inputs are vectors
  corresponding to observations at a set of discrete time
  points for right censored data.
}
\details{
  \describe{ \item{lrTests}{log-rank family of tests: \cr
  These are given by the general expression: \cr \cr \deqn{
  Q = \sum{ W_i (e_i - \hat{e}_i)}^T \sum{ W_i \hat{V_i}
  W_i}^{-1} \sum{ W_i (e_i - \hat{e}_i)} }{ Q = [sum
  W(e-E)]^T [sum WVW]^-1 [sum W(e-E)] } Where \eqn{W} is
  the weight, given as in \code{comp2Surv},
  \eqn{\hat{e}}{E} is the no. of expected events for that
  time and \eqn{\hat{V}} is the variance-covariance matrix
  given by \code{covMatSurv}. \cr \cr If there are \eqn{K}
  groups, then \eqn{K-1} are selected (arbitrary). Likewise
  the corresponding variance-covariance matrix is reduced
  to the appropriate \eqn{K-1 \times K-1} dimensions.
  \eqn{Q} is distributed as chi-squared with \eqn{K-1}
  degree of freedom. \cr \cr }

  \item{trendTests}{Tests for trend. \cr \cr These are
  designed to detect ordered differences in survival
  curves, that is, for at least one group: \deqn{S_1(t)
  \geq S_2(t) \geq ... \geq S_K(t) \quad t \leq \tau}{
  S1(t) >= S2(t) >= ... >= SK(t) for t <= tau} where
  \eqn{\tau}{tau} is the largest \eqn{t} where all groups
  have at least one subject at risk. \cr The null
  hypothesis is that \deqn{S_1(t) = S_2(t) = ... = S_K(t)
  \quad t \leq \tau}{ S1(t) = S2(t) = ... = SK(t) for t <=
  tau} \cr Scores used to construct the test are typically
  \eqn{s = 1,2,...,K} but may be given as a vector
  representing a numeric characteristic of the group. \cr
  They are calculated by finding \deqn{ Z_j(t_i) =
  \sum_{t_i \leq \tau} W(t_i)[e_{ji} - n_{ji}
  \frac{e_i}{n_i}], \quad j=1,2,...,K}{ Z[t(i)] = SUM
  W[t(i)] [e[j](i) - n[j](i).e(i)/n(i) ]} \cr The test
  statistic is \deqn{Z = \frac{ \sum_{j=1}^K
  s_jZ_j(\tau)}{\sqrt{\sum_{j=1}^K \sum_{g=1}^K s_js_g
  \sigma_{jg}}} }{ Z = SUM(j=1...K) s[j]Z[j] / SUM(j=1..K)
  SUM(g=1..K) s[j]s[g]sigma[jg]} where \eqn{\sigma}{sigma}
  is the the appropriate element in the variance-covariance
  matrix (as in \code{comp2Surv}). \cr \cr If ordering is
  present, the statistic \eqn{Z} will be greater than the
  upper \eqn{\alpha}{alpha}th percentile of a standard
  normal distribution. \cr \cr } }
}
\examples{
data(tneBMT)
compNSurv(t=tneBMT$t,n=tneBMT$n,e=tneBMT$e,
 n1=as.matrix(tneBMT[,grep("n_",colnames(tneBMT))]),
 e1=as.matrix(tneBMT[,grep("e_",colnames(tneBMT))]),
 FHp=1,FHq=1)
}
\references{
  Tarone RE, Ware J 1977 On Distribution-Free Tests for
  Equality of Survival Distributions.
  \emph{Biometrika};\bold{64}(1):156--60.
  \href{http://www.jstor.org/stable/2335790}{JSTOR}
}
\seealso{
  \code{\link{comp2Surv}} for weights

  Calls \code{\link{covMatSurv}}

  Called by \code{\link{compSurvfit}}
}

