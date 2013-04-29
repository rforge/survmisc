\docType{data}
\name{tneKidney}
\alias{tneKidney}
\title{Time, no. at risk, no. events for kidney data}
\format{A data frame with 16 rows and 7 columns
Columns are:
\describe{
 \item{t}{time}
 \item{n_1}{no. at risk in group 1}
 \item{e_1}{no. events in group 1}
 \item{n_2}{no. at risk in group 2}
 \item{e_2}{no. events in group 2}
 \item{n}{no. at risk overall}
 \item{e}{no. events overall}
}}
\source{
  Klein J, Moeschberger 2003 \emph{Survival Analysis}, 2nd
  edition. New York: Springer. Example 7.9, pg 224.
}
\description{
  Data on the recurrence times to infection, at the point
  of insertion of the catheter, for kidney patients using
  portable dialysis equipment. Generated from
  \code{data("kidney",package="KMsurv")}. Catheters may be
  removed for reasons other than infection, in which case
  the observation is censored. Each patient has exactly 2
  observations.
}
\references{
  McGilchrist CA, Aisbett CW 1991 Regression with Frailty
  in Survival Analysis. \emph{Biometrics}
  \bold{47}(2):461--6.
  \href{http://www.jstor.org/stable/2532138}{JSTOR}
}
\seealso{
  \code{kidney} in package \code{survival}
}
\keyword{datasets}

