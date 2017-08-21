#' SLWeightedScreen: A collection of observation weight-respecting screening algorithms for SuperLearner
#'
#' The SLWeightedScreen package provides feature selection or 'screen' functions
#' that plug into the R package for the ensemble learner
#' \code{\link[SuperLearner]{SuperLearner}}. Specifically, the screening
#' algorithms in this package utilize observation weights
#' (\code{obsWeights}) in their calculations.
#'
#' @section Filter methods:
#'
#' \itemize{
#'  \item{Correlation:}{\itemize{
#'  	\item{\code{\link{screen.wgtd.corP}}}
#' 		\item{\code{\link{screen.wgtd.corRank}}}
#' 	}}
#'  \item{Comparison test:}{\itemize{
#'  	\item{\code{\link{screen.wgtd.ttestP}}}
#' 		\item{\code{\link{screen.wgtd.ttestRank}}}
#' 	}}
#' }
#'
#' @section Embedded methods:
#'
#' \itemize{
#'  \item{Regularized/penalized regression:}{\itemize{
#'  	\item{\code{\link{screen.wgtd.lasso}}}
#' 		\item{\code{\link{screen.wgtd.elasticnet}}}
#' 	}}
#' }
#'
#'
#' @docType package
#' @name SLWeightedScreen
NULL