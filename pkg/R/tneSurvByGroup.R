##' @name tneSurvByGroup
##' @export
##' @include tneSurvfit.R
##' @title Time, No. at risk, No. events by group
##' @description Gives time, no. at risk and no. events for survival data by group
##' @param t Times
##' @param e No. events (overall)
##' @param p Predictor
##' @param np No. events for this predictor
##' @param onlyEvents if \code{TRUE} shows only times at which at least one event occurred. Otherwise shows \emph{all} times recorded (including those censored).
##' @return A list with one element corresponding to each value of the predictor \eqn{p}.
##' \cr \cr
##' Each element is a matrix, with one row for each observation.
##' \cr
##' The columns in the matrix are:
##' \item{t}{time}
##' \item{n}{no. at risk}
##' \item{e}{no. events}
##' @references Example is from:
##' Klein J, Moeschberger M 2003
##' \emph{Survival Analysis}, 2nd edition.
##' New York: Springer.
##' Example 7.2, pg 210.
##' @examples
##' data(kidney,package="KMsurv")
##' s <- survfit(Surv(time=time, event=delta) ~ type, data=kidney )
##' df2 <- tneSurvfit(s)
##' tneSurvByGroup(t=df2$t,e=df2$e,p=df2$p,np=df2$np)
tneSurvByGroup <- function(t,e,p,np,onlyEvents=FALSE){
    if(!isTRUE( all.equal(length(t),length(e),length(p),length(np)) )) stop ("All vectors must be of equal length")
    if(!isTRUE(all(vapply(c(t,e,p,np),FUN=is.numeric,FUN.VALUE=TRUE)==TRUE))) stop("All vectors must be numeric")
    p1 <- sort(unique(p))
### subset and aggregate
    subAgg <- function(i){
        sub1 <- data.frame(t=t[p==p1[i]],
                           np=np[p==p1[i]],
                           e=e[p==p1[i]])
        res1 <- as.matrix(
            stats::aggregate(
                sub1, by=list(sub1$t),
                FUN=identity)[2:4])
        if (is.list(res1)){
### no names yet, so refer to columns by no.
### take max for t and np
            res1[,1:2] <- lapply(res1[,1:2], FUN=max)
### take sum for e
            res1[,3] <- lapply(res1[,3], FUN=sum)
            res1 <- matrix(unlist(res1),nrow=nrow(res1),ncol=3)
        }
        colnames(res1) <- c("t","n","e")
        if (onlyEvents) res1 <- res1[!res1[,3]==0, ]
        return(res1)
    }
    res <- sapply(X=1:length(p1),FUN=subAgg)
    names(res) <- paste("pred_",p1,sep="")
    return(res)
}
