#' Weighted t-test screening algorithm
#'
#' Performs feature selection according to the ranking of t statistics
#' returned from weighted t-tests.
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
	    						...)$coefficients[["t.value"]]
		}, Y = Y, obsWeights = obsWeights)
	}
	whichVariable <- (rank(-abs(listT)) <= minscreen)
	return(whichVariable)
}

#' Weighted t-test screening algorithm
#'
#' Performs feature selection according to the ranking of P-values
#' returned from weighted t-tests.
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
	    						...)$coefficients[["p.value"]]
		}, Y = Y, obsWeights = obsWeights)
	}
	whichVariable <- (listP <= minPvalue)
	if (sum(whichVariable) < minscreen) {
		warning('number of variables with p value less than minPvalue is less than minscreen')
		whichVariable[rank(listP) <= minscreen] <- TRUE
	}
	return(whichVariable)
}
