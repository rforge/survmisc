##' @name comp2Surv
##' @export
##' @title Compare two survival curves
##' @description Compare two survival curves.
##' \cr \cr
##' Inputs are vectors corresponding to observations at a set of discrete time points for right censored data.
##' @param n No. at risk (overall)
##' @param e No. events (overall)
##' @param n1 No. at risk in group 1 (arbitrary; either group)
##' @param e1 No. events in group 1 (arbitrary; either group)
##' @param FHp \eqn{p} for Fleming-Harrington test
##' @param FHq \eqn{q} for Fleming-Harrington test
##' @param lim limit; used the calculating estimate of
##' \eqn{Pr[\sup|B(t)|>x]}{Pr[sup|B(x)]} for supremum tests.
##' Higher values will be more accuate,
##' but this is asymptotic for values of \eqn{>1e4}.
##' @param round1 No. digits to which to round (for display)
##' @return A list of 2 data frames. In each, there is one row for each of the weights, as given above. \describe{
##' \item{lrTests}{log-rank family of tests. \cr
##' Columns give the value of the test statistic, \eqn{Q}, degrees of freedom (=1 for 2 groups) and corresponding p value from chi-squared distribution.}
##' \item{supTests}{Supremum (Renyi) family of tests. \cr
##' Columns give the value of the test statistic \eqn{Q}  and corresponding p value.}
##' }
##' @details \describe{
##'  \item{lrTests}{log-rank family of tests:
##'  \cr
##'  These are given by the general expression:
##' \cr \cr
##' \deqn{ Q = \sum{ W_i (e_i - \hat{e}_i)}^T \sum{ W_i \hat{V_i} W_i}^{-1} \sum{ W_i (e_i - \hat{e}_i)} }{
##'  Q = [SUM W(e-E)]^T [SUM WVW]^-1 [SUM W(e-E)] }
##'  Where \eqn{W} is the weight, given below, \eqn{e} is the no. of events,
##' \eqn{\hat{e}}{E} is the no. of expected events for that time
##' and \eqn{\hat{V}} is the variance-covariance matrix
##' given by \code{covMatSurv}.
##' \cr
##' The sum is taken to the largest observed survival time (censored observations are excluded).
##' \cr \cr
##' For 2 groups, this can also be written as:
##' \deqn{ Q = \frac{ \sum W_i (e1_i - n1_i \frac{e_i}{n_i})}{
##'  \sqrt{\sum W_i^2 \frac{n1_i}{n_i}(1-\frac{n1_1i}{n_i}) (\frac{n_i-e_i}{n_i-1}) e_i} } }{
##'  Z = SUM W [e1 - n1.(e/n)] / (SUM W^2.e1/e.(1-(n1/n)).(n-e/n-1).e)^0.5 }
##' where \eqn{n1_1i}{n1} refers to the no. at risk in group 1.
##' \cr \cr
##' The weights are given as follows: \itemize{
##'  \item  Log-rank \eqn{W = 1}
##'  \item Gehan-Breslow generalized Wilcoxon \eqn{W = n}, the no. at risk
##'  \item Tarone-Ware \eqn{W = \sqrt{n}}{W = n^0.5}
##'  \item Peto-Peto, \eqn{W = \bar{S}(t)}{W = S(t)}, a  modified estimator of survival function given by:
##' \deqn{\bar{S}(t)=\prod{1 - \frac{e_i}{n_i+1}}}{
##'  S(t) = CUMPROD ( 1- e/(n+1))}
##'  \item modified Peto-Peto (by Andersen) \eqn{W = \bar{S}(t) \frac{n}{n+1} }{S(t) n/n+1}
##'  \item Fleming-Harrington at \eqn{t_0}{t0} \eqn{W=1} at and thereafter
##'  \deqn{ W = \hat{S}(t_{i-1})^p [1-\hat{S}(t_{i-1})]^q}{
##'  W = S(t(i-1))^p [1 - S(i-1) ]^q }
##'  Here \eqn{\hat{S}} is the Kaplan-Meier (product-limit) estimator.
##' \cr
##' Note that both \eqn{p} and \eqn{q} need to be \eqn{>= 0}
##' \cr
##' }}
##'
##' \item{supTests}{Supremum (Renyi) family of tests.
##' \cr
##' These are designed to detect differences in survival curves which cross. That is, an early difference in survival in favor of one group is balanced by a later reversal.
##' \cr
##' The same weights as above are used.
##' They are calculated by finding
##' \deqn{ Z(t_i) = \sum_{t_k \leq t_i} W(t_k)[e1_k - n1_k\frac{e_k}{n_k}], \quad i=1,2,...,k}{
##' Z(t[i]) = SUM W(t[k]) [ e1[k] - n1[k]e[k]/n[k] ]}
##' (which is similar to the numerator used to find \eqn{Q} in the log-rank test for 2 groups above) \cr
##' and it's variance:
##' \deqn{ \sigma^2(\tau) = \sum_{t_k \leq \tau} W(t_k)^2 \frac{n1_k n2_k (n_k-e_k) e_k}{n_k^2 (n_k-1)} }{
##' simga^2(tau) = SUM(k=1,2,...,tau) W(t[k]) [ n1[k].n2[k].(n[k]-e[k]).e[k] / n[k]^2.(n[k]-1) ] }
##' where \eqn{\tau}{tau} is the largest \eqn{t} where both groups have at least one subject at risk.
##' \cr \cr
##' Then calculate:
##' \deqn{ Q = \frac{ \sup{|Z(t)|}}{\sigma(\tau)} ,t<\tau }{
##' Q = sup( |Z(t)| ) / sigma(tau), t<tau}
##' \cr
##' When the null hypothesis is true, the distribution of \eqn{Q} is approximately
##' \deqn{Q \sim \sup{|B(x)|, \quad 0 \leq x \leq 1}}{
##' Q ~ sup( |B(x)|, 0 <= x <= 1)}
##' \cr
##' For a standard Brownian motion (Wiener) process:
##' \deqn{ Pr[\sup|B(t)|>x] = 1-\frac{4}{\pi} \sum_{k=0}^{\infty} \frac{(-1)^k}{2k+1} \exp{ \frac{-\pi^2(2k+1)^2}{8x^2}}}{
##' Pr[sup|B(t)|>x] = 1 - 4/pi SUM (-1)^k/2k+1 exp(-pi^2 (2k+1)^2/x^2)}
##' \cr
##'   }
##'  }
##' @note Fleming-Harrington weights: \itemize{
##' \item \eqn{p = q = 0} gives the log-rank test, i.e. \eqn{W=1}
##' \item \eqn{p=1, q=0} gives a version of the Mann-Whitney-Wilcoxon test (tests if populations ditributions are identical)
##' \item \eqn{q=0, p>0} gives more weight to differences early on
##' \item \eqn{q>0, p=0} gives more weight to differences later on
##' }
##'
##' @seealso Calls \code{\link{covMatSurv}}
##' @seealso Called by \code{\link{compSurvfit}}
##' @references Peto R, Peto J 1972
##' Asymptotically Efficient Rank Invariant Test Procedures.
##' \emph{J Royal Statistical Society} \bold{135}(2):186--207.
##' \href{http://www.jstor.org/stable/2344317}{JSTOR}
##' @references Fleming TR, Harrington DP, O'Sullivan M 1987
##' Supremum Versions of the Log-Rank and Generalized Wilcoxon Statistics.
##' \emph{J  American Statistical Association} \bold{82}(397):312--20.
##' \href{http://www.jstor.org/stable/2289169}{JSTOR}
##' @references Billingsly P 2009
##' \emph{Convergence of Probability Measures.}
##' New York: John Wiley & Sons.
##' \href{http://books.google.com/books/about/Convergence_of_Probability_Measures.html?id=GzjbezrsrFcC}{Google Books}
##' @examples
##' data(tneKidney)
##' comp2Surv(n=tneKidney$n,e=tneKidney$e,n1=tneKidney$n_1,e1=tneKidney$e_1)
comp2Surv <- function(n,e,n1,e1,
                      FHp=1,FHq=1,lim=1e4,round1=5){
    if (FHp<0|FHq<0) stop("Values for p and q for Fleming-Harrington tests must be >=0")
    if(!isTRUE( all.equal(length(n),length(e),length(n1),length(e1)) )) stop ("All vectors must be of equal length")
    if(!isTRUE(all(vapply(c(n,e,n1,e1),FUN=is.numeric,FUN.VALUE=TRUE)==TRUE))) stop("All vectors must be numeric")
###
### make observed - expected for group 1
    eME1 <- e1-(n1*e/n)
### variance (2 groups)
    var1 <- (n1/n)*(1-(n1/n))*((n-e)/(n-1))*e
### display chisq, degrees of freedom and rounded result
    dis1 <- function(chi1,df1=1,rounded=round1){
        c(chi1,df1, round(1-stats::pchisq(chi1,df1),digits=rounded))
    }
###
### WEIGHTS
### log-rank, weight = 1
    chi1 <- sum(eME1)^2/sum(var1)
    LR1 <- dis1(chi1)
### Gehan-Breslow generalized Wilcoxon, weight = n
    chi1 <- sum(n*eME1)^2/sum(n^2*var1)
    GB1 <- dis1(chi1)
### Tarone-Ware, weight = sqrt(n)
    chi1 <- sum(sqrt(n)*eME1)^2/ sum(sqrt(n)^2*var1)
    TW1 <- dis1(chi1)
### Peto-Peto, weight = S(t) = modified estimator of survival function
    S1 <- cumprod( (1- (e/(n+1))) )
    chi1 <- sum(S1*eME1)^2/sum(S1^2*var1)
    PP1 <- dis1(chi1)
### modified Peto-Peto (by Andersen), weight = S(t)n/n+1
    S2 <- (S1*n)/(n+1)
    chi1 <- sum(S2*eME1)^2/sum(S2^2*var1)
    mPP1 <- dis1(chi1)
### Fleming-Harrington; weight 1st element is 1 as depends on [i-1]
    S3 <- cumprod( (n-e)/n )
    S3 <- c(1,S3[1:(length(S3)-1)])
    FHw <- S3^FHp*((1-S3)^FHq)
    chi1 <- sum(FHw*eME1)^2/sum(FHw^2*var1)
    FH1 <- dis1(chi1)
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
### Renyi statistics (analagous to 2-sample Kolmogorov-Smirnov test)
### aka supremum tests
###
### Probability of supremum of absolute value of
### standard Brownian motion process B(t)
    probSupBr <- function(y, limit=lim){
        k <- seq(from=0,to=limit,by=1)
        res <- 1-( (4/pi)* sum(  ((-1)^k) /(2*k+1) * exp( (-(pi^2)*(2*k+1)^2) / (8*y^2) )  ) )
        return(res)
    }
    Z1 <- cumsum(e1-(n1*(e/n)))
### display chisq, degrees of freedom and rounded result
    dis2 <- function(Q1,rounded=round1){
        c(Q1,round(probSupBr(Q1),digits=rounded))
      }
    Q1 <- max(abs(Z1))/sqrt(sum(var1))
### log-rank weights
    RenLR1 <- dis2(Q1)
### Gehan-Breslow generalized Wilcoxon, weight = n
    Z1 <- cumsum( n* (e1-(n1*(e/n))) )
    Q1 <- max(abs(Z1))/sqrt(sum(n^2*var1))
    RenGB1 <- dis2(Q1)
### Tarone-Ware, weight = sqrt(n)
    Z1 <- cumsum( sqrt(n)* (e1-(n1*(e/n))) )
    Q1 <- max(abs(Z1))/sqrt(sum(n*var1))
    RenTW1 <- dis2(Q1)
### Peto-Peto, weight = S(t) = modified estimator of survival function
    S1 <- cumprod( (n-e)/n )
    Z1 <- cumsum( S1 *(e1-(n1*(e/n))) )
    Q1 <- max(abs(Z1))/sqrt(sum(S1^2*var1))
    RenPP1 <- dis2(Q1)
### modified Peto-Peto (by Andersen), weight = S(t)n/n+1
    Zfun <- function(i)
    Z1 <- cumsum( S2*(e1-(n1*(e/n))) )
    Q1 <- max(abs(Z1))/sqrt(sum(S2^2*var1))
    RenmPP1 <- dis2(Q1)
### Fleming-Harrington; weight 1st element is 1 as depends on [i-1]
    Zfun <- function(i)
    Z1 <- cumsum( FHw *(e1-(n1*(e/n))) )
    Q1 <- max(abs(Z1))/sqrt(sum(FHw^2*var1))
    RenFH1 <- dis2(Q1)
### result
    FHnR <- paste("Renyi ",FHn,sep="")
    res2 <- rbind(RenLR1,RenGB1,RenTW1,RenPP1,RenmPP1,RenFH1)
    dimnames(res2) <- list(c(
        "Log-rank",
        "Gehan-Breslow (mod~ Wilcoxon)",
        "Tarone-Ware",
        "Peto-Peto",
        "Mod~ Peto-Peto (Andersen)",
        FHnR), c(
            "Q","p"))
### final
    return( list(lrTests=res1,
                 supTests=res2))
}
