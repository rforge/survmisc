\name{comp2Surv}
\alias{comp2Surv}
\title{Compare two survival curves}
\usage{
  comp2Surv(n, e, n1, e1, FHp = 1, FHq = 1, lim = 10000,
    round1 = 5)
}
\arguments{
  \item{n}{No. at risk (overall)}

  \item{e}{No. events (overall)}

  \item{n1}{No. at risk in group 1 (arbitrary; either
  group)}

  \item{e1}{No. events in group 1 (arbitrary; either
  group)}

  \item{FHp}{\eqn{p} for Fleming-Harrington test}

  \item{FHq}{\eqn{q} for Fleming-Harrington test}

  \item{lim}{limit; used the calculating estimate of
  \eqn{Pr[\sup|B(t)|>x]}{Pr[sup|B(x)]} for supremum tests.
  Higher values will be more accuate, but this is
  asymptotic for values of \eqn{>1e4}.}

  \item{round1}{No. digits to which to round (for display)}
}
\value{
  A list of 2 data frames. In each, there is one row for
  each of the weights, as given above. \describe{
  \item{lrTests}{log-rank family of tests. \cr Columns give
  the value of the test statistic, \eqn{Q}, degrees of
  freedom (=1 for 2 groups) and corresponding p value from
  chi-squared distribution.} \item{supTests}{Supremum
  (Renyi) family of tests. \cr Columns give the value of
  the test statistic \eqn{Q} and corresponding p value.} }
}
\description{
  Compare two survival curves. \cr \cr Inputs are vectors
  corresponding to observations at a set of discrete time
  points for right censored data.
}
\details{
  \describe{ \item{lrTests}{log-rank family of tests: \cr
  These are given by the general expression: \cr \cr \deqn{
  Q = \sum{ W_i (e_i - \hat{e}_i)}^T \sum{ W_i \hat{V_i}
  W_i}^{-1} \sum{ W_i (e_i - \hat{e}_i)} }{ Q = [SUM
  W(e-E)]^T [SUM WVW]^-1 [SUM W(e-E)] } Where \eqn{W} is
  the weight, given below, \eqn{e} is the no. of events,
  \eqn{\hat{e}}{E} is the no. of expected events for that
  time and \eqn{\hat{V}} is the variance-covariance matrix
  given by \code{covMatSurv}. \cr The sum is taken to the
  largest observed survival time (censored observations are
  excluded). \cr \cr For 2 groups, this can also be written
  as: \deqn{ Q = \frac{ \sum W_i (e1_i - n1_i
  \frac{e_i}{n_i})}{ \sqrt{\sum W_i^2
  \frac{n1_i}{n_i}(1-\frac{n1_1i}{n_i})
  (\frac{n_i-e_i}{n_i-1}) e_i} } }{ Z = SUM W [e1 -
  n1.(e/n)] / (SUM W^2.e1/e.(1-(n1/n)).(n-e/n-1).e)^0.5 }
  where \eqn{n1_1i}{n1} refers to the no. at risk in group
  1. \cr \cr The weights are given as follows: \itemize{
  \item Log-rank \eqn{W = 1} \item Gehan-Breslow
  generalized Wilcoxon \eqn{W = n}, the no. at risk \item
  Tarone-Ware \eqn{W = \sqrt{n}}{W = n^0.5} \item
  Peto-Peto, \eqn{W = \bar{S}(t)}{W = S(t)}, a modified
  estimator of survival function given by:
  \deqn{\bar{S}(t)=\prod{1 - \frac{e_i}{n_i+1}}}{ S(t) =
  CUMPROD ( 1- e/(n+1))} \item modified Peto-Peto (by
  Andersen) \eqn{W = \bar{S}(t) \frac{n}{n+1} }{S(t) n/n+1}
  \item Fleming-Harrington at \eqn{t_0}{t0} \eqn{W=1} at
  and thereafter \deqn{ W = \hat{S}(t_{i-1})^p
  [1-\hat{S}(t_{i-1})]^q}{ W = S(t(i-1))^p [1 - S(i-1) ]^q
  } Here \eqn{\hat{S}} is the Kaplan-Meier (product-limit)
  estimator. \cr Note that both \eqn{p} and \eqn{q} need to
  be \eqn{>= 0} \cr }}

  \item{supTests}{Supremum (Renyi) family of tests. \cr
  These are designed to detect differences in survival
  curves which cross. That is, an early difference in
  survival in favor of one group is balanced by a later
  reversal. \cr The same weights as above are used. They
  are calculated by finding \deqn{ Z(t_i) = \sum_{t_k \leq
  t_i} W(t_k)[e1_k - n1_k\frac{e_k}{n_k}], \quad
  i=1,2,...,k}{ Z(t[i]) = SUM W(t[k]) [ e1[k] -
  n1[k]e[k]/n[k] ]} (which is similar to the numerator used
  to find \eqn{Q} in the log-rank test for 2 groups above)
  \cr and it's variance: \deqn{ \sigma^2(\tau) = \sum_{t_k
  \leq \tau} W(t_k)^2 \frac{n1_k n2_k (n_k-e_k) e_k}{n_k^2
  (n_k-1)} }{ simga^2(tau) = SUM(k=1,2,...,tau) W(t[k]) [
  n1[k].n2[k].(n[k]-e[k]).e[k] / n[k]^2.(n[k]-1) ] } where
  \eqn{\tau}{tau} is the largest \eqn{t} where both groups
  have at least one subject at risk. \cr \cr Then
  calculate: \deqn{ Q = \frac{ \sup{|Z(t)|}}{\sigma(\tau)}
  ,t<\tau }{ Q = sup( |Z(t)| ) / sigma(tau), t<tau} \cr
  When the null hypothesis is true, the distribution of
  \eqn{Q} is approximately \deqn{Q \sim \sup{|B(x)|, \quad
  0 \leq x \leq 1}}{ Q ~ sup( |B(x)|, 0 <= x <= 1)} \cr For
  a standard Brownian motion (Wiener) process: \deqn{
  Pr[\sup|B(t)|>x] = 1-\frac{4}{\pi} \sum_{k=0}^{\infty}
  \frac{(-1)^k}{2k+1} \exp{ \frac{-\pi^2(2k+1)^2}{8x^2}}}{
  Pr[sup|B(t)|>x] = 1 - 4/pi SUM (-1)^k/2k+1 exp(-pi^2
  (2k+1)^2/x^2)} \cr } }
}
\note{
  Fleming-Harrington weights: \itemize{ \item \eqn{p = q =
  0} gives the log-rank test, i.e. \eqn{W=1} \item
  \eqn{p=1, q=0} gives a version of the
  Mann-Whitney-Wilcoxon test (tests if populations
  ditributions are identical) \item \eqn{q=0, p>0} gives
  more weight to differences early on \item \eqn{q>0, p=0}
  gives more weight to differences later on }
}
\examples{
data(tneKidney)
comp2Surv(n=tneKidney$n,e=tneKidney$e,n1=tneKidney$n_1,e1=tneKidney$e_1)
}
\references{
  Peto R, Peto J 1972 Asymptotically Efficient Rank
  Invariant Test Procedures. \emph{J Royal Statistical
  Society} \bold{135}(2):186--207.
  \href{http://www.jstor.org/stable/2344317}{JSTOR}

  Fleming TR, Harrington DP, O'Sullivan M 1987 Supremum
  Versions of the Log-Rank and Generalized Wilcoxon
  Statistics. \emph{J American Statistical Association}
  \bold{82}(397):312--20.
  \href{http://www.jstor.org/stable/2289169}{JSTOR}

  Billingsly P 2009 \emph{Convergence of Probability
  Measures.} New York: John Wiley & Sons.
  \href{http://books.google.com/books/about/Convergence_of_Probability_Measures.html?id=GzjbezrsrFcC}{Google
  Books}
}
\seealso{
  Calls \code{\link{covMatSurv}}

  Called by \code{\link{compSurvfit}}
}

