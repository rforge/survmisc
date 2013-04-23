##' @name compSurvfit
##' @export
##' @aliases compSurv
##' @include tneSurv.R
##' @include tneSurvByGroup.R
##' @include comp2Surv.R
##' @include compNSurv.R
##' @title Compare two survival curves
##' @description Compare survival curves.
##' @param s A \code{survfit} object
##' @param FHp \eqn{p} for Fleming-Harrington test
##' @param FHq \eqn{q} for Fleming-Harrington test
##' @param round1 No. digits to which to round (for display)
##' @param scores scores (used if >2 groups) if no value is given this is set as \eqn{1,2,...,K} the no. of predictors
##' @param lim limit (used if 2 groups) limit to which to calculate \eqn{Pr[\sup|B(t)|>x]}{Pr[sup|B(x)]})
##' @return A list with the following elements: \describe{
##' \item{tne}{Time, no. at risk and no. events, overall and per predictor.
##'  One row for each time where at least one event occurs}
##' \item{tneP}{A list of dataframes, one for each predictor.
##'  Each shows time, no. at risk and no. events.}
##' \item{tests}{Tests as returned by \code{comp2Surv} or \code{compNSurv}
##'  depending on whether there are 2 or more than 2 predictors.}
##' }
##' @seealso \code{\link{tneSurvByGroup}}
##' @seealso \code{\link{comp2Surv}}
##' @seealso \code{\link{compNSurv}}
##' @references Gehan A.
##' A Generalized Wilcoxon Test for Comparing Arbitrarily Singly-Censored Samples.
##' Biometrika 1965 Jun. 52(1/2):203--23.
##' \href{http://www.jstor.org/stable/2333825}{JSTOR}
##' @references Tarone RE, Ware J 1977
##' On Distribution-Free Tests for Equality of Survival Distributions.
##' \emph{Biometrika};\bold{64}(1):156--60.
##' \href{http://www.jstor.org/stable/2335790}{JSTOR}
##' @examples
##' data(kidney,package="KMsurv")
##' s <- survfit(Surv(time=time, event=delta) ~ type, data=kidney )
##' compSurvfit(s)
##' data(gastric)
##' s1 <- survfit(Surv(time=time,event=event) ~ group, data=gastric)
##' compSurvfit(s1)
##' data(bmt,package="KMsurv")
##' b1 <- bmt[bmt$group==1,]
##' s2 <- survfit(Surv(time=t2, event=d3) ~ group, data=bmt)
##' compSurvfit(s2)
compSurvfit <- function(s, FHp=1, FHq=1, round1=5, scores="", lim=1e4) {
    if (!class(s)=="survfit") stop("Only applies to object of class 'survfit'")
    if (FHp<0|FHq<0) stop("Values for p and q for Fleming-Harrington tests must be >=0")
### get time + no. events from survival formula
    df2 <- tneSurvfit(s)
    p1 <- sort(unique(df2$p))
###
### time, number (at risk), events by group
    tne1 <- tneSurvByGroup(t=df2$t,e=df2$e,p=df2$p,np=df2$np)
### merge all observations based on time
    m1 <- data.table::data.table(Reduce(function(...) merge(..., by="t", all=TRUE), tne1))
### for 'n' carry last observation back (to fill in missing values in first rows)
    data.table::set(m1,j=grep("n",colnames(m1)), value=zoo::na.locf(m1[,grep("n",colnames(m1)),with=FALSE], fromLast=TRUE))
### for 'n' carry last observation forward
#    data.table::set(m1,j=grep("n",colnames(m1)), value=zoo::na.locf(m1[,grep("n",colnames(m1)),with=FALSE]))
### for remaining 'n' and 'e', replace NA with zero
### for 'n' this will be elements in the tail of the vector
 for (j in seq_len(ncol(m1))){
     data.table::set(m1,which(is.na(m1[[j]])),j,0)
    }
### names
    n1 <- c('n_','e_')
    n1 <- as.vector(outer(n1,p1,paste,sep=""))
    data.table::setnames(m1,old=names(m1),new=c("t",n1))
### remove times with no events
    m1 <- m1[!rowSums(m1[,grep("e",colnames(m1)),with=FALSE])==0, ]
### initialize constants (prevent warning from R CMD check)
    n <- e <- NULL
### make no. at risk (total) per time period *inefficient method*
    m1[, n := rowSums(m1[,grep("n",colnames(m1)),with=FALSE]) ]
### total events per time period
    m1[, e := rowSums(m1[,grep("e",colnames(m1)),with=FALSE]) ]
### 2 groups only:
    if (length(p1)==2){
        res1 <- comp2Surv(n=m1$n,e=m1$e,n1=m1$n_1,e1=m1$e_1,
                          FHp=1,FHq=1,round1=round1)
    } else {
        res1 <- compNSurv(t=m1$t,n=m1$n,e=m1$e,
                          n1=as.matrix(m1[,grep("n_",colnames(m1)),with=FALSE]),
                          e1=as.matrix(m1[,grep("e_",colnames(m1)),with=FALSE]),
                          FHp=FHp,FHq=FHq,round1=round1)
    }
###
    res <- list(
        tne = m1,
        tneP = tne1,
        tests = res1)
    return(res)
}
