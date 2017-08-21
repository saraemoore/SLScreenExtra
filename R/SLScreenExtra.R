#' SLScreenExtra: A collection of screening algorithms for SuperLearner
#'
#' The SLScreenExtra package provides feature selection or \code{screen.*()}
#' functions that plug into the R package for the ensemble learner
#' \code{\link[SuperLearner]{SuperLearner}}. Some of the screening
#' algorithms (\code{screen.wgtd.*()}) in this package utilize observation
#' weights (\code{obsWeights}) in their calculations.
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
#' @name SLScreenExtra
NULL