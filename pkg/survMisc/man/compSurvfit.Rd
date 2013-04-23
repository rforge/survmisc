\name{compSurvfit}
\alias{compSurv}
\alias{compSurvfit}
\title{Compare two survival curves}
\usage{
  compSurvfit(s, FHp = 1, FHq = 1, round1 = 5, scores = "",
    lim = 10000)
}
\arguments{
  \item{s}{A \code{survfit} object}

  \item{FHp}{\eqn{p} for Fleming-Harrington test}

  \item{FHq}{\eqn{q} for Fleming-Harrington test}

  \item{round1}{No. digits to which to round (for display)}

  \item{scores}{scores (used if >2 groups) if no value is
  given this is set as \eqn{1,2,...,K} the no. of
  predictors}

  \item{lim}{limit (used if 2 groups) limit to which to
  calculate \eqn{Pr[\sup|B(t)|>x]}{Pr[sup|B(x)]})}
}
\value{
  A list with the following elements: \describe{
  \item{tne}{Time, no. at risk and no. events, overall and
  per predictor.  One row for each time where at least one
  event occurs} \item{tneP}{A list of dataframes, one for
  each predictor.  Each shows time, no. at risk and no.
  events.} \item{tests}{Tests as returned by
  \code{comp2Surv} or \code{compNSurv} depending on whether
  there are 2 or more than 2 predictors.} }
}
\description{
  Compare survival curves.
}
\examples{
data(kidney,package="KMsurv")
s <- survfit(Surv(time=time, event=delta) ~ type, data=kidney )
compSurvfit(s)
data(gastric)
s1 <- survfit(Surv(time=time,event=event) ~ group, data=gastric)
compSurvfit(s1)
data(bmt,package="KMsurv")
b1 <- bmt[bmt$group==1,]
s2 <- survfit(Surv(time=t2, event=d3) ~ group, data=bmt)
compSurvfit(s2)
}
\references{
  Gehan A. A Generalized Wilcoxon Test for Comparing
  Arbitrarily Singly-Censored Samples. Biometrika 1965 Jun.
  52(1/2):203--23.
  \href{http://www.jstor.org/stable/2333825}{JSTOR}

  Tarone RE, Ware J 1977 On Distribution-Free Tests for
  Equality of Survival Distributions.
  \emph{Biometrika};\bold{64}(1):156--60.
  \href{http://www.jstor.org/stable/2335790}{JSTOR}
}
\seealso{
  \code{\link{tneSurvByGroup}}

  \code{\link{comp2Surv}}

  \code{\link{compNSurv}}
}

