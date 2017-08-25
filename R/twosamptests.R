#' Weighted t-test screening algorithm
#'
#' Performs feature selection according to the ranking of t statistics
#' returned from weighted t-tests. Implemented via
#' \code{\link[weights]{wtd.t.test}}.
#'
#' @param Y Outcome (numeric vector). See \code{\link[SuperLearner]{SuperLearner}}
#' for specifics.
#' @param X Predictor variable(s) (data.frame or matrix). See
#' \code{\link[SuperLearner]{SuperLearner}} for specifics.
#' @param family Error distribution to be used in the model:
#' \code{\link[stats]{gaussian}} or \code{\link[stats]{binomial}}.
#' See \code{\link[SuperLearner]{SuperLearner}} for specifics.
#' @param obsWeights Optional numeric vector of observation weights. See
#' \code{\link[SuperLearner]{SuperLearner}} for specifics.
#' @param id Cluster identification variable. Currently unused.
#' @param minscreen Minimum number of features to select (aka rank). Only used
#' if less than this number of features are selected using \code{minPvalue}.
#' @param ... Passed to \code{\link[weights]{wtd.t.test}}. These arguments
#' control bootstrapping of P-values and standard errors as well as forced
#' scaling of weights.
#' @return A logical vector with length equal to \code{ncol(X)}
#' @export
#' @importFrom weights wtd.t.test
#' @importFrom stats var
#' @examples
#' # based on example in SuperLearner package
#' set.seed(1)
#' n <- 100
#' p <- 20
#' X <- matrix(rnorm(n*p), nrow = n, ncol = p)
#' X <- data.frame(X)
#' Y <- rbinom(n, 1, plogis(.2*X[, 1] + .1*X[, 2] - .2*X[, 3] + .1*X[, 3]*X[, 4] - .2*abs(X[, 4])))
#' obsWeights <- 1/runif(n)
#' screen.wgtd.ttestRank(Y, X, binomial(), obsWeights, seq(n), minscreen = 4)
#'
#' screen.wgtd.ttestRank4 <- function(..., minscreen = 4){
#'     screen.wgtd.ttestRank(..., minscreen = minscreen)
#' }
#'
#' library(SuperLearner)
#' sl = SuperLearner(Y, X, family = binomial(), cvControl = list(V = 2),
#'                   obsWeights = obsWeights,
#'                   SL.library = list(c("SL.lm", "All"),
#'                                     c("SL.lm", "screen.wgtd.ttestRank4")))
#' sl
#' sl$whichScreen
screen.wgtd.ttestRank <- function(Y, X, family, obsWeights, id, minscreen = 2, ...) {
	if(family$family == "gaussian") {
		stop('t-test screening undefined for gaussian family, look at screen.wgtd.corP or screen.wgtd.corRank')
	}
	if(family$family == "binomial") {
		listT <- apply(X, 2, function(x, Y, obsWeights) {
			if(var(x) <= 0)
				return(0)
			x <- split(x, as.factor(Y))
			# no variability in Y:
			if(length(x) < 2)
				return(0) # this should probably instead be NA but SL might choke
			obsWeights <- split(obsWeights, as.factor(Y))
	    	weights::wtd.t.test(x[[1]], y = x[[2]],
	    						weight = obsWeights[[1]], weighty = obsWeights[[2]],
	    						samedata = FALSE, alternative="two.tailed", ...)$coefficients[["t.value"]]
		}, Y = Y, obsWeights = obsWeights)
	}
	whichVariable <- (rank(-abs(listT)) <= minscreen)
	return(whichVariable)
}

#' Weighted t-test screening algorithm
#'
#' Performs feature selection according to the ranking of P-values
#' returned from weighted t-tests. Implemented via
#' \code{\link[weights]{wtd.t.test}}.
#'
#' @param Y Outcome (numeric vector). See \code{\link[SuperLearner]{SuperLearner}}
#' for specifics.
#' @param X Predictor variable(s) (data.frame or matrix). See
#' \code{\link[SuperLearner]{SuperLearner}} for specifics.
#' @param family Error distribution to be used in the model:
#' \code{\link[stats]{gaussian}} or \code{\link[stats]{binomial}}.
#' See \code{\link[SuperLearner]{SuperLearner}} for specifics.
#' @param obsWeights Optional numeric vector of observation weights. See
#' \code{\link[SuperLearner]{SuperLearner}} for specifics.
#' @param id Cluster identification variable. Currently unused.
#' @param minPvalue To pass the screen, resulting P-values must not exceed this
#' number.
#' @param minscreen Minimum number of features to select (aka rank). Only used
#' if less than this number of features are selected using \code{minPvalue}.
#' @param ... Passed to \code{\link[weights]{wtd.t.test}}. These arguments
#' control bootstrapping of P-values and standard errors as well as forced
#' scaling of weights.
#' @return A logical vector with length equal to \code{ncol(X)}
#' @export
#' @importFrom weights wtd.t.test
#' @importFrom stats var
#' @examples
#' # based on example in SuperLearner package
#' set.seed(1)
#' n <- 100
#' p <- 20
#' X <- matrix(rnorm(n*p), nrow = n, ncol = p)
#' X <- data.frame(X)
#' Y <- rbinom(n, 1, plogis(.2*X[, 1] + .1*X[, 2] - .2*X[, 3] + .1*X[, 3]*X[, 4] - .2*abs(X[, 4])))
#' obsWeights <- 1/runif(n)
#' screen.wgtd.ttestP(Y, X, binomial(), obsWeights, seq(n), minPvalue = 0.05)
#'
#' screen.wgtd.ttestP05 <- function(..., minPvalue = 0.05){
#'     screen.wgtd.ttestP(..., minPvalue = minPvalue)
#' }
#'
#' library(SuperLearner)
#' sl = SuperLearner(Y, X, family = binomial(), cvControl = list(V = 2),
#'                   obsWeights = obsWeights,
#'                   SL.library = list(c("SL.lm", "All"),
#'                                     c("SL.lm", "screen.wgtd.ttestP05")))
#' sl
#' sl$whichScreen
screen.wgtd.ttestP <- function(Y, X, family, obsWeights, id, minPvalue = 0.1, minscreen = 2, ...) {
	if(family$family == "gaussian") {
		stop('t-test screening undefined for gaussian family, look at screen.wgtd.corP or screen.wgtd.corRank')
	}
	if(family$family == "binomial") {
		listP <- apply(X, 2, function(x, Y, obsWeights) {
			if(var(x) <= 0)
				return(0)
			x <- split(x, as.factor(Y))
			# no variability in Y:
			if(length(x) < 2)
				return(0) # this should probably instead be NA but SL might choke
			obsWeights <- split(obsWeights, as.factor(Y))
	    	weights::wtd.t.test(x[[1]], y = x[[2]],
	    						weight = obsWeights[[1]], weighty = obsWeights[[2]],
	    						samedata = FALSE, alternative="two.tailed", ...)$coefficients[["p.value"]]
		}, Y = Y, obsWeights = obsWeights)
	}
	whichVariable <- (listP <= minPvalue)
	if (sum(whichVariable) < minscreen) {
		warning('number of variables with p value less than minPvalue is less than minscreen')
		whichVariable[rank(listP) <= minscreen] <- TRUE
	}
	return(whichVariable)
}