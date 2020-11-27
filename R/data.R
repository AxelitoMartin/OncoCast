#' Toy dataset for survival analysis without left truncation
#'
#' @format A data frame with 300 observations and 17 variables
#' \describe{
#'   \item{time}{Time to event continuous variable.}
#'   \item{status}{Binary variable specifying if an event occured}
#'   \item{ImpCov's}{Binary features that were used to generate the survival outcome.}
#'   \item{Cov's}{Binary noise features.}
#' }
"survData"

#' Toy dataset for survival analysis with left truncation
#'
#' @format A data frame with 1000 observations and 113 variables
#' \describe{
#'   \item{time1}{Delayed entry time continuous variable.}
#'   \item{time2}{Time to event continuous variable.}
#'   \item{status}{Binary variable specifying if an event occured}
#'   \item{ImpCov's}{Binary features that were used to generate the survival outcome.}
#'   \item{Cov's}{Binary noise features.}
#' }
"survData.LT"

#' Toy dataset for gaussian analysis
#'
#' @format A data frame with 300 observations and 51 variables
#' \describe{
#'   \item{y}{Normally distributed outcome.}
#'   \item{Cov's}{Covariates of interest.}
#' }
"normData"

#' Toy dataset for binary analysis
#'
#' @format A data frame with 300 observations and 51 variables
#' \describe{
#'   \item{y}{Binomial distributed outcome.}
#'   \item{Cov's}{Covariates of interest.}
#' }
"binData"

#' appTest_data.csv: Test data used for the Rshiny application, showing the format to be used for importing data to be analysed via OncoCast
#'
#' @format A data frame with 200 observations and 53 variables
#' \describe{
#'   \item{time}{Time to event continuous variable.}
#'   \item{status}{Binary variable specifying if an event occured}
#'   \item{ContCov's}{Continuous features.}
#'   \item{BinCov's}{Binary features.}
#' }


#' appTest_valdata.csv: Test data used for the Rshiny application, showing the format to be used for importing data to validate an OncoCast run.
#'
#' @format A data frame with 100 observations and 53 variables
#' \describe{
#'   \item{time}{Time to event continuous variable.}
#'   \item{status}{Binary variable specifying if an event occured}
#'   \item{ContCov's}{Continuous features.}
#'   \item{BinCov's}{Binary features.}
#' }
