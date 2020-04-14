\name{survData}
\alias{survData}
\docType{data}
\title{
Toy dataset for survival analysis.
}
\description{
Small dataset without any left trunction for survival analysis.
}
\usage{data("survData")}
\format{
  A data frame with 300 observations on the following 17 variables.
  \describe{
    \item{\code{time}}{a numeric vector}
    \item{\code{status}}{a numeric vector}
    \item{\code{ImpCov1}}{a numeric vector}
    \item{\code{ImpCov2}}{a numeric vector}
    \item{\code{ImpCov3}}{a numeric vector}
    \item{\code{ImpCov4}}{a numeric vector}
    \item{\code{ImpCov5}}{a numeric vector}
    \item{\code{Cov6}}{a numeric vector}
    \item{\code{Cov7}}{a numeric vector}
    \item{\code{Cov8}}{a numeric vector}
    \item{\code{Cov9}}{a numeric vector}
    \item{\code{Cov10}}{a numeric vector}
    \item{\code{Cov11}}{a numeric vector}
    \item{\code{Cov12}}{a numeric vector}
    \item{\code{Cov13}}{a numeric vector}
    \item{\code{Cov14}}{a numeric vector}
    \item{\code{Cov15}}{a numeric vector}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(survData)
# test <- OncoCast(data=survData,formula = Surv(time,status)~.,
# method = "LASSO",runs = 30,
# save = FALSE,nonPenCol = NULL,cores =2)
}
\keyword{datasets}
