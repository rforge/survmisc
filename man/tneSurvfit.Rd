\name{tneSurvfit}
\alias{tneSurvfit}
\title{Time, No. at risk, No. events, No. expected (by predictor)}
\usage{
  tneSurvfit(s = s)
}
\arguments{
  \item{s}{A \code{survfit} object}
}
\value{
  A \code{data frame} with the following columns:
  \item{t}{time} \item{n}{no. at risk (total)} \item{e}{no.
  events (total)} \item{p}{predictor} \item{np}{no. at risk
  (by predictor)} \item{Ep}{no. events expected (by
  predictor)} \item{e_Ep}{no. events minus no. events
  expected}
}
\description{
  Gives time, no. at risk and no. events for a
  \code{survfit} object describing right censored data with
  one predictor. Also per-predictor estimates of no. at
  risk, no. expected and no. events - no. expected. \cr \cr
  No. events expected (per predictor) is given by:
  \deqn{\frac{e_i(n[p]_i)}{n_i}}{ e(i)n[p](i) / n(i)} where
  \eqn{n[p]_i} is the no. at risk for the predictor.
}
\examples{
data(kidney,package="KMsurv")
s <- survfit(Surv(time=time, event=delta) ~ type, data=kidney)
tneSurvfit(s)
s <- survfit(Surv(time=time, event=delta) ~ 1, data=kidney)
tneSurvfit(s)
}

