##' @name tneSurv
##' @export
##' @title Time, No. at risk, No. events
##' @description Gives time, no. at risk and no. events for a \code{Surv} object describing right censored data.
##' @param s A \code{Surv} object
##' @param onlyEvents if \code{TRUE} shows only times at which at least one event occurred.
##' \cr
##' Otherwise shows \emph{all} times recorded (including those censored).
##' @return A data frame with columns:
##' \item{t}{time}
##' \item{n}{no. at risk}
##' \item{e}{no. events}
##' @examples
##' df0 <- data.frame(t=c(1,1,2,3,5,8,13,21),
##'                   e=rep(c(0,1),4))
##' s1 <- Surv(df0$t,df0$e,type="right")
##' tneSurv(s1)
tneSurv <- function(s,onlyEvents=TRUE){
    if(!class(s)=="Surv") stop("Only applies to class 'Surv'")
    if(!attr(s,which="type")=="right") warning("Applies to right censored data")
### no. subjects leaving study in each time period
    r1 <- as.matrix(stats::aggregate(s[,"time"],
                                     by=list(s[,"time"]),
                                     FUN=length))
### no. events in each time period
    r2 <- as.matrix(stats::aggregate(s[,"status"],
                                     by=list(s[,"time"]),
                                     FUN=sum))
### join by time
    m1 <- merge(r1,r2,by=1)
    nr1 <- nrow(s)
### make no. at risk
    m1$n <- if(nrow(m1)==1){
        1} else {
            c(nrow(s),nrow(s)-(cumsum(m1[1:(nrow(m1)-1),2])))
        }
### reorder
    m1 <- m1[,c(1,4,3)]
    colnames(m1) <- c("t","n","e")
### (optional) show only times where events occurred
    if (onlyEvents) m1 <- m1[!m1[,3]==0, ]
    return(m1)
}
