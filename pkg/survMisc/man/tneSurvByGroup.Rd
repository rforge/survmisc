\name{tneSurvByGroup}
\alias{tneSurvByGroup}
\title{Time, No. at risk, No. events by group}
\usage{
  tneSurvByGroup(t, e, p, np, onlyEvents = FALSE)
}
\arguments{
  \item{t}{Times}

  \item{e}{No. events (overall)}

  \item{p}{Predictor}

  \item{np}{No. events for this predictor}

  \item{onlyEvents}{if \code{TRUE} shows only times at
  which at least one event occurred. Otherwise shows
  \emph{all} times recorded (including those censored).}
}
\value{
  A list with one element corresponding to each value of
  the predictor \eqn{p}. \cr \cr Each element is a matrix,
  with one row for each observation. \cr The columns in the
  matrix are: \item{t}{time} \item{n}{no. at risk}
  \item{e}{no. events}
}
\description{
  Gives time, no. at risk and no. events for survival data
  by group
}
\examples{
data(kidney,package="KMsurv")
s <- survfit(Surv(time=time, event=delta) ~ type, data=kidney )
df2 <- tneSurvfit(s)
tneSurvByGroup(t=df2$t,e=df2$e,p=df2$p,np=df2$np)
}
\references{
  Example is from: Klein J, Moeschberger M 2003
  \emph{Survival Analysis}, 2nd edition. New York:
  Springer. Example 7.2, pg 210.
}

