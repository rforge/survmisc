\docType{data}
\name{gastric}
\alias{gastric}
\title{gastric cancer trial data}
\format{A data frame with 90 rows and 3 columns
\cr
Columns are:
\describe{
 \item{time}{time in days}
 \item{event}{death}
 \item{group}{\describe{
   \item{0}{chemotherapy}
   \item{1}{chemotherapy + radiotherapy}
   }}
 }}
\source{
  Klein J, Moeschberger 2003 \emph{Survival Analysis}, 2nd
  edition. New York: Springer. Example 7.9, pg 224.
}
\description{
  gastric cancer trial data
}
\details{
  Data from a trial of locally unresectable gastic cancer.
  \cr Patients (45 in each group) were randomized to one of
  two groups: chemotheapy vs chemotherapy + radiotherapy.
}
\references{
  Gastrointestinal Tumor Study Group 1982 A comparison of
  combination chemotherapy and combined modality therapy
  for locally advanced gastric carcinoma.
  \emph{Cancer};\bold{49}(9):1771--7.
  \href{http://www.ncbi.nlm.nih.gov/pubmed/6176313}{Pubmed}.

  Stablein DM, Koutrouvelis IA 1985 A two-sample test
  sensitive to crossing hazards in uncensored and singly
  censored data. \emph{Biometrics};\bold{41}(3):643--52.
  \href{http://www.jstor.org/stable/2531284}{JSTOR}.
}
\keyword{datasets}

