##' @name compNSurv
##' @include covMatSurv.R
##' @include tneSurvfit.R
##' @export
##' @title Compare more than two survival curves
##' @description Compare more than two survival curves.
##' Inputs are vectors corresponding to observations at a set of discrete time points for right censored data.
##' @param t Time
##' @param n No. at risk (overall)
##' @param e No. events (overall)
##' @param n1 No. at risk per group
##' @param e1 No. of events per group
##' @param FHp \eqn{p} for Fleming-Harrington test
##' @param FHq \eqn{q} for Fleming-Harrington test
##' @param scores scores; if no value is given this is set as \eqn{1,2,...,K} the no. of predictors
##' @param round1 No. digits to which to round (for display)
##' @return A list with 3 elements. The first are 2 data frames. In each, there is one row for each of the weights, as given above. \describe{

##' \item{lrTests}{log-rank family of tests. \cr
##' Columns give the value of the test statistic, \eqn{Q}, degrees of freedom (\eqn{= K-1} for \eqn{K} groups) and corresponding p value the from chi-squared distribution.}
##' \item{trendTests}{Tests for trend. \cr
##' Columns give the value of \eqn{Z} and corresponding p value from the normal distribution.}
##' \item{scores}{A vector of the scores used in calculating the trendTests above.}
##' }
##' @details \describe{
##'  \item{lrTests}{log-rank family of tests:
##' \cr
##'  These are given by the general expression:
##' \cr \cr
##' \deqn{ Q = \sum{ W_i (e_i - \hat{e}_i)}^T \sum{ W_i \hat{V_i} W_i}^{-1} \sum{ W_i (e_i - \hat{e}_i)} }{
##'  Q = [sum W(e-E)]^T [sum WVW]^-1 [sum W(e-E)] }
##'  Where \eqn{W} is the weight, given as in \code{comp2Surv}, \eqn{\hat{e}}{E} is the no. of expected events for that time and \eqn{\hat{V}} is the variance-covariance matrix given by \code{covMatSurv}.
##' \cr \cr
##' If there are \eqn{K} groups, then \eqn{K-1} are selected (arbitrary). Likewise the corresponding variance-covariance matrix is reduced to the appropriate \eqn{K-1 \times K-1} dimensions.
##' \eqn{Q} is distributed as chi-squared with \eqn{K-1} degree of freedom.
##' \cr \cr
##' }
##'
##' \item{trendTests}{Tests for trend.
##' \cr \cr
##' These are designed to detect ordered differences in survival curves,
##' that is, for at least one group:
##' \deqn{S_1(t) \geq S_2(t) \geq ... \geq S_K(t) \quad t \leq \tau}{
##'  S1(t) >= S2(t) >= ... >= SK(t) for t <= tau}
##' where \eqn{\tau}{tau} is the largest \eqn{t} where all groups have at least one subject at risk.
##' \cr
##' The null hypothesis is that
##' \deqn{S_1(t) = S_2(t) = ... = S_K(t) \quad t \leq \tau}{
##' S1(t) = S2(t) = ... = SK(t) for t <= tau}
##' \cr
##' Scores used to construct the test are typically \eqn{s = 1,2,...,K} but may be given as a vector representing a numeric characteristic of the group.
##' \cr
##' They are calculated by finding
##' \deqn{ Z_j(t_i) = \sum_{t_i \leq \tau} W(t_i)[e_{ji} - n_{ji} \frac{e_i}{n_i}], \quad j=1,2,...,K}{
##' Z[t(i)] = SUM W[t(i)] [e[j](i) - n[j](i).e(i)/n(i) ]}
##' \cr
##' The test statistic is
##' \deqn{Z = \frac{ \sum_{j=1}^K s_jZ_j(\tau)}{\sqrt{\sum_{j=1}^K \sum_{g=1}^K s_js_g \sigma_{jg}}} }{
##' Z = SUM(j=1...K) s[j]Z[j] / SUM(j=1..K) SUM(g=1..K) s[j]s[g]sigma[jg]}
##' where \eqn{\sigma}{sigma} is the the appropriate element in the variance-covariance matrix (as in \code{comp2Surv}).
##' \cr \cr
##' If ordering is present, the statistic \eqn{Z} will be greater than the upper \eqn{\alpha}{alpha}th percentile of a standard normal distribution.
##' \cr \cr
##'   }
##' }
##' @seealso \code{\link{comp2Surv}} for weights
##' @seealso Calls \code{\link{covMatSurv}}
##' @seealso Called by \code{\link{compSurvfit}}
##' @references Tarone RE, Ware J 1977
##' On Distribution-Free Tests for Equality of Survival Distributions.
##' \emph{Biometrika};\bold{64}(1):156--60.
##' \href{http://www.jstor.org/stable/2335790}{JSTOR}
##' @examples
##' data(tneBMT)
##' compNSurv(t=tneBMT$t,n=tneBMT$n,e=tneBMT$e,
##'  n1=as.matrix(tneBMT[,grep("n_",colnames(tneBMT))]),
##'  e1=as.matrix(tneBMT[,grep("e_",colnames(tneBMT))]),
##'  FHp=1,FHq=1)
compNSurv <- function (t,n,e,n1,e1,
                       FHp,FHq,scores="",round1=5){
        if(!isTRUE( all.equal(length(n),length(e),length(n1),length(e1)) )) stop ("All vectors must be of equal length")
        if(!isTRUE(all(vapply(c(n,e,n1,e1),FUN=is.numeric,FUN.VALUE=TRUE)==TRUE))) stop("All vectors must be numeric")
### events observed minus expected
    eME <- e1-(n1*e/n)
### degrees of freedom
    df1 <- ifelse(is.null(dim(eME)),1,ncol(eME))-1
### make covariance matrix (n groups)
    v1 <- covMatSurv(t,n,e,n1)
    cov1 <- rowSums(v1, dims=2)
### display chisq, degrees of freedom and rounded result
    dis <- function(chi1,df=df1,rounded=round1){
        c(chi1,df1, round(1-stats::pchisq(chi1,df),digits=rounded))
    }
### WEIGHTS
### log-rank, weight = 1
    eME1 <- colSums(eME)
    chi1 <- eME1[1:df1] %*% solve(cov1[1:df1,1:df1]) %*% matrix((eME1[1:df1]),nrow=df1,ncol=1)
    LR1 <- dis(chi1)
### Gehan-Breslow generalized Wilcoxon, weight = n
    v2 <- sweep(v1,MARGIN=3,STATS=n^2,FUN="*")
    covG <- rowSums(v2,dims=2)
    eMEG <- colSums(sweep(eME,1,n,"*"))
    chi1 <- eMEG[1:df1] %*% solve(covG[1:df1,1:(df1)]) %*% matrix((eMEG[1:df1]),df1,1)
    GB1 <- dis(chi1)
### Tarone-Ware, weight = sqrt(n)
    v2 <- sweep(v1,MARGIN=3,STATS=n,FUN="*")
    covTw <- rowSums(v2,dims=2)
    eMETw <- colSums(sweep(eME,1,sqrt(n),"*"))
    chi1 <- eMETw[1:df1] %*% solve(covTw[1:df1,1:df1]) %*% matrix((eMETw[1:df1]),df1,1)
    TW1 <- dis(chi1)
### Peto-Peto, weight = S(t) = modified estimator of survival function
    S1 <- cumprod( (1- (e/(n+1))) )
    v2 <- sweep(v1,MARGIN=3,STATS=S1^2,FUN="*")
    covPP <- rowSums(v2,dims=2)
    eMEPP <- colSums(sweep(eME,1,S1,"*"))
    chi1 <- eMEPP[1:df1] %*% solve(covPP[1:df1,1:df1]) %*% matrix((eMEPP[1:df1]),df1,1)
    PP1 <- dis(chi1)
### modified Peto-Peto (by Andersen), weight = S(t)n/n+1
    S2 <- (S1*n)/(n+1)
    v2 <- sweep(v1,MARGIN=3,STATS=S2^2,FUN="*")
    covmPP <- rowSums(v2,dims=2)
    eMEmPP <- colSums(sweep(eME,1,S2,"*"))
    chi1 <- eMEmPP[1:df1] %*% solve(covmPP[1:df1,1:df1]) %*% matrix((eMEmPP[1:df1]),df1,1)
    mPP1 <- dis(chi1)
### Fleming-Harrington; weight 1st element is 1 as depends on [i-1]
    S1 <- cumprod( (n-e)/n )
    S1 <- c(1,S1[1:(length(S1)-1)])
    FH1 <- S1^FHp*((1-S1)^FHq)
    v2 <- sweep(v1,MARGIN=3,STATS=FH1^2,FUN="*")
    covFH <- rowSums(v2,dims=2)
    eMEFH <- colSums(sweep(eME,1,FH1,"*"))
    chi1 <- eMEFH[1:df1] %*% solve(covFH[1:df1,1:df1]) %*% matrix((eMEFH[1:df1]),df1,1)
    FH1 <- dis(chi1)
### results
    FHn <- paste("Flem~-Harr~ with p=",FHp,", q=",FHq,sep="")
    res1 <- rbind(LR1,GB1,TW1,PP1,mPP1,FH1)
    dimnames(res1) <- list(c(
        "Log-rank",
        "Gehan-Breslow (mod~ Wilcoxon)",
        "Tarone-Ware",
        "Peto-Peto",
        "Mod~ Peto-Peto (Andersen)",
        FHn), c(
            "ChiSq", "df", "p"))
###
### Trend tests
    suppressWarnings(if (scores=="") scores <- 1:(df1+1) )
### no. predictors (from degrees of freedom, above)
    lp1 <- df1+1
### display chisq, degrees of freedom and rounded result
    dis2 <- function(Z1,rounded=round1){
        c(Z1,round(1-stats::pnorm(Z1),digits=rounded))
      }
### log-rank weights
    sfun <- function(i,j) scores[i]*scores[j]*cov1[i,j]
    p2 <- as.vector(sapply(lp1, function (x) rep(x,lp1)))
    s1 <- sum(mapply(sfun, rep(seq(lp1),lp1), p2))
### eME1 = obs - exp, from log-rank test above
    Z1 <- sum(scores*eME1)/sqrt(s1)
    LRT1 <- dis2(Z1)
### Gehan-Breslow generalized Wilcoxon, weight = n
    sfun <- function(i,j) scores[i]*scores[j]*covG[i,j]
    s1 <- sum(mapply(sfun, rep(seq(lp1),lp1), p2))
    Z1 <- sum(scores*eMEG)/sqrt(s1)
    GBT1 <- dis2(Z1)
### Tarone-Ware, weight = sqrt(n)
    sfun <- function(i,j) scores[i]*scores[j]*covTw[i,j]
    s1 <- sum(mapply(sfun, rep(seq(lp1),lp1), p2))
    Z1 <- sum(scores*eMETw)/sqrt(s1)
    TWT1 <- dis2(Z1)
### Peto-Peto, weight = S(t) = modified estimator of survival function
    sfun <- function(i,j) scores[i]*scores[j]*covPP[i,j]
    s1 <- sum(mapply(sfun, rep(seq(lp1),lp1), p2))
    Z1 <- sum(scores*eMEPP)/sqrt(s1)
    PPT1 <- dis2(Z1)
### modified Peto-Peto (by Andersen), weight = S(t)n/n+1
    sfun <- function(i,j) scores[i]*scores[j]*covmPP[i,j]
    s1 <- sum(mapply(sfun, rep(seq(lp1),lp1), p2))
    Z1 <- sum(scores*eMEmPP)/sqrt(s1)
    mPPT1 <- dis2(Z1)
### Fleming-Harrington; weight 1st element is 1 as depends on [i-1]
    sfun <- function(i,j) scores[i]*scores[j]*covFH[i,j]
    s1 <- sum(mapply(sfun, rep(seq(lp1),lp1), p2))
    Z1 <- sum(scores*eMEFH)/sqrt(s1)
    FHT1 <- dis2(Z1)
### result
    FHnT <- paste("Trend F-H with p=",FHp,", q=",FHq,sep="")
    res2 <- rbind(LRT1,GBT1,TWT1,PPT1,mPPT1,FHT1)
    dimnames(res2) <- list(c(
        "Log-rank",
        "Gehan-Breslow (mod~ Wilcoxon)",
        "Tarone-Ware",
        "Peto-Peto",
        "Mod~ Peto-Peto (Andersen)",
        FHnT), c(
            "Z", "p"))
### final
    return( list(lrTests=res1,
                 trendTests=res2,
                 scores=scores))
}
