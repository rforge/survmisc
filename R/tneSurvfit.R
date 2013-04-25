##' @name tneSurvfit
##' @export
##' @title Time, No. at risk, No. events, No. expected (by predictor)
##' @description Gives time, no. at risk and no. events for a \code{survfit} object describing right censored data with one predictor. Also per-predictor estimates of no. at risk, no. expected and no. events - no. expected.
##' \cr \cr
##' No. events expected (per predictor) is given by:
##' \deqn{\frac{e_i(n[p]_i)}{n_i}}{
##'  e(i)n[p](i) / n(i)}
##' where \eqn{n[p]_i} is the no. at risk for the predictor.
##' @param s A \code{survfit} object
##' @return A \code{data frame} with the following columns:
##'  \item{t}{time}
##'  \item{n}{no. at risk (total)}
##'  \item{e}{no. events (total)}
##'  \item{p}{predictor}
##'  \item{np}{no. at risk (by predictor)}
##'  \item{Ep}{no. events expected (by predictor)}
##'  \item{e_Ep}{no. events minus no. events expected}
##' @examples
##' data(kidney,package="KMsurv")
##' s <- survfit(Surv(time=time, event=delta) ~ type, data=kidney)
##' tneSurvfit(s)
##' s <- survfit(Surv(time=time, event=delta) ~ 1, data=kidney)
##' tneSurvfit(s)
tneSurvfit <- function(s=s){
    if (!class(s)=="survfit") stop("Only applies to object of class 'survfit'")
    if (length(eval(s$call[[2]])[[3]])>1) stop("Only applies to survfit formulas with one predictor")
### get location to evaluate variables (environment or data frame)
    if (is.null(s$call$data)){
        loc1 <- environment(eval(parse(text=as.character(s$call[2]))))
    } else {
        loc1 <- eval(s$call$data)
    }
### length of data frame (no. rows)
    l1 <- length(get(ls(loc1),loc1))
### hold predictor, time,  event, no. at risk (total)
### no. at risk (per predictor), no. Expected events (per predictor)
### no. events - no. Expected events (per predictor)
    df2 <- data.frame(matrix(0,nrow=l1,ncol=7))
    colnames(df2) <- c("t","n","e","p","np","Ep","e_Ep")
### get names of time, event and predictor...
### for name of time:
### get formula: eval(s$call[[2]]),
### then get LHS (Surv object): eval(s$call[[2]])[2]
### then remove trailing (): eval(s$call[[2]])[2][[1]]
### then get time variable
    t1 <- as.character(eval(s$call[[2]])[2][[1]][2])
    df2$t <- get(t1,loc1)
    e1 <- as.character(eval(s$call[[2]])[2][[1]][3])
    df2$e <- with(loc1, get(e1))
### for predictor:
### get RHS of formula: eval(s$call[[2]])[3]
### then remove trailing (): eval(s$call[[2]])[2][[1]]
    p1 <-  as.character(eval(s$call[[2]])[3][[1]])
### if one predictor (intercept-only) model:
    if (p1==1){
        df2$p <- rep(1,l1)
        } else {
            df2$p <-  with(loc1, get(p1))
            }
    df2 <- df2[order(df2$t, decreasing=FALSE), ]
    rownames(df2) <- NULL
### make no. at risk (total)
    df2[1,"n"] <- nrow(df2)
    df2[2:nrow(df2),"n"] <- nrow(df2)-cumsum(df2[,"e"][-nrow(df2)])
### make no. at risk (per predictor)
    p1 <- sort(unique(df2$p))
    for (i in 1:length(p1)) {
### total no. at risk
        s1 <- sum(df2$p==p1[i])
        df2[which(df2$p==p1[i]),"np"] <- seq(s1,1)
    }
### make no. expected events (per predictor)
    df2[,"Ep"] <- (df2[,"e"]*df2[,"np"]) / df2[,"n"]
### make events - expected
    df2[,"e_Ep"] <- df2[,"e"]-df2[,"Ep"]
    return(df2)
}
