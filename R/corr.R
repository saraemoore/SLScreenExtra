#' Weighted correlation screening algorithm
#'
#' Performs feature selection according to the ranking of weighted
#' correlation coefficient estimates. Implemented via
#' \code{\link[wCorr]{weightedCorr}}.
#'
#' @param Y Outcome (numeric vector). See \code{\link[SuperLearner]{SuperLearner}}
#' for specifics.
#' @param X Predictor variable(s) (data.frame or matrix). See
#' \code{\link[SuperLearner]{SuperLearner}} for specifics.
#' @param family Error distribution to be used in the model:
#' \code{\link[stats]{gaussian}} or \code{\link[stats]{binomial}}.
#' Currently unused. See \code{\link[SuperLearner]{SuperLearner}}
#' for specifics.
#' @param obsWeights Optional numeric vector of observation weights. See
#' \code{\link[SuperLearner]{SuperLearner}} for specifics.
#' @param id Cluster identification variable. Currently unused.
#' @param method Which correlation coefficient to compute. Currently accepts
#' \code{"pearson"} or \code{"spearman"}.
#' @param k Minimum number of features to select.
#' @param ... Currently unused.
#' @return A logical vector with length equal to \code{ncol(X)}
#' @export
#' @importFrom wCorr weightedCorr
#' @importFrom stats var
#' @examples
#' # based on examples in SuperLearner package
#' set.seed(1)
#' n <- 100
#' p <- 20
#' X <- matrix(rnorm(n*p), nrow = n, ncol = p)
#' X <- data.frame(X)
#' Y <- X[, 1] + sqrt(abs(X[, 2] * X[, 3])) + X[, 2] - X[, 3] + rnorm(n)
#' obsWeights <- 1/runif(n)
#' screen.wgtd.corRank(Y, X, gaussian(), obsWeights, seq(n), k = 3)
#'
#' screen.wgtd.corRank3 <- function(..., k = 3){
#'     screen.wgtd.corRank(..., k = k)
#' }
#'
#' library(SuperLearner)
#' sl = SuperLearner(Y, X, family = gaussian(), cvControl = list(V = 2),
#'                   obsWeights = obsWeights,
#'                   SL.library = list(c("SL.glm", "All"),
#'                                     c("SL.glm.interaction", "screen.wgtd.corRank3")))
#' sl
#' sl$whichScreen
screen.wgtd.corRank <- function(Y, X, family, obsWeights, id, method = "pearson", k = 2, ...) {
    if(!method%in%c("pearson", "spearman")) {
        stop("Correlation method ", method, " not supported by screen.wgtd.corRank")
    }
    if(is.null(obsWeights)) {
        obsWeights <- rep(1, length(Y))
    }
    listCor <- apply(X, 2, function(x, Y, method, obsWeights) {
    	# if x is homogenous, bump it to the bottom of the list
        ifelse(var(x) <= 0, 0, wCorr::weightedCorr(x, Y, method = method, weights = obsWeights))
    }, Y = Y, method = method, obsWeights = obsWeights)

    whichVariable <- (rank(-abs(listCor)) <= k)
    return(whichVariable)
}

#' Weighted correlation screening algorithm
#'
#' Performs feature selection according to the ranking of P-values
#' returned from weighted correlations. Implemented via
#' \code{\link[weights]{wtd.cor}}.
#'
#' @inheritParams screen.wgtd.corRank
#' @param method Which correlation coefficient to compute. Currently only
#' accepts \code{"pearson"}.
#' @param minPvalue To pass the screen, resulting P-values must not exceed this
#' number.
#' @param k Minimum number of features to select. Only used
#' if less than this number of features are selected using \code{minPvalue}.
#' @param ... Passed to \code{\link[weights]{wtd.cor}}. These arguments control
#' bootstrapping of P-values and standard errors as well as forced scaling of
#' weights.
#' @inherit screen.wgtd.corRank return
#' @export
#' @importFrom weights wtd.cor
#' @importFrom stats var
#' @examples
#' # based on example in SuperLearner package
#' set.seed(1)
#' n <- 100
#' p <- 20
#' X <- matrix(rnorm(n*p), nrow = n, ncol = p)
#' X <- data.frame(X)
#' Y <- X[, 1] + sqrt(abs(X[, 2] * X[, 3])) + X[, 2] - X[, 3] + rnorm(n)
#' obsWeights <- 1/runif(n)
#' screen.wgtd.corP(Y, X, gaussian(), obsWeights, seq(n), minPvalue = 0.000001)
#'
#' screen.wgtd.corP01 <- function(..., minPvalue = 0.01){
#'     screen.wgtd.corP(..., minPvalue = minPvalue)
#' }
#'
#' library(SuperLearner)
#' sl = SuperLearner(Y, X, family = gaussian(), cvControl = list(V = 2),
#'                   obsWeights = obsWeights,
#'                   SL.library = list(c("SL.glm", "All"),
#'                                     c("SL.glm.interaction", "screen.wgtd.corP01")))
#' sl
#' sl$whichScreen
screen.wgtd.corP <- function(Y, X, family, obsWeights, id, method = "pearson", minPvalue = 0.1, k = 2, ...) {
	if(method!="pearson") {
        stop("Correlation method ", method, " not supported by screen.wgtd.corRank")
    }
	listP <- apply(X, 2, function(x, Y, obsWeights) {
    	# if x is homogeneous, bump it to the bottom of the list
		ifelse(var(x) <= 0, 1, weights::wtd.cor(x, y = Y, weight = obsWeights, ...)[, "p.value"])
	}, Y = Y, obsWeights = obsWeights)
	whichVariable <- (listP <= minPvalue)
	if (sum(whichVariable) < k) {
	  warning('number of variables with p value less than minPvalue is less than k')
		whichVariable[rank(listP) <= k] <- TRUE
	}
	return(whichVariable)
}
