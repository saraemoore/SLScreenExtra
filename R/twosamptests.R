#' Weighted t-test screening algorithm
#'
#' Performs feature selection according to the ranking of t statistics
#' or P-values returned from weighted t-tests. Implemented via
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
#' @param selector A string corresponding to a subset selecting function
#' implemented in the FSelector package. One of:
#' \code{\link[FSelector]{cutoff.k}} or \code{\link[FSelector]{cutoff.k.percent}}.
#' Ignored if \code{minP} is non-\code{NULL} and at least \code{k} features
#' have P-values at or below \code{minP}. Default: \code{"cutoff.k"}.
#' @param k Numeric. Minimum number or proportion of features to select.
#' Passed through to the \code{selector}. For \code{\link[FSelector]{cutoff.k}},
#' this is an integer indicating the number of features to keep from \code{X}.
#' For \code{\link[FSelector]{cutoff.k.percent}}, this is instead the proportion
#' of features to keep. Ignored if \code{minP} is non-\code{NULL} and at least
#' \code{k} features have P-values at or below \code{minP}.
#' @param minP Numeric. To pass the screen, resulting P-values must not exceed this
#' number. Ignored if \code{NULL} (default) or if fewer than \code{k} features
#' have P-values at or below this value.
#' @param ... Passed to \code{\link[weights]{wtd.t.test}}. These arguments
#' control bootstrapping of P-values and standard errors as well as forced
#' scaling of weights.
#' @return A logical vector with length equal to \code{ncol(X)}
#' @export
#' @importFrom weights wtd.t.test
#' @importFrom stats var
#' @importFrom methods is
#' @importFrom FSelector cutoff.k cutoff.k.percent
#' @rdname screen.wgtd.ttest
#' @examples
#' # based on example in SuperLearner package
#' set.seed(1)
#' n <- 100
#' p <- 20
#' X <- matrix(rnorm(n*p), nrow = n, ncol = p)
#' X <- data.frame(X)
#' Y <- rbinom(n, 1, plogis(.2*X[, 1] + .1*X[, 2] - .2*X[, 3] + .1*X[, 3]*X[, 4] - .2*abs(X[, 4])))
#' obsWeights <- 1/runif(n)
#' screen.wgtd.ttest(Y, X, binomial(), obsWeights, seq(n), k = 4)
#'
#' screen.wgtd.ttest4 <- function(..., k = 4){
#'     screen.wgtd.ttest(..., k = k)
#' }
#'
#' library(SuperLearner)
#' sl = SuperLearner(Y, X, family = binomial(), cvControl = list(V = 2),
#'                   obsWeights = obsWeights,
#'                   SL.library = list(c("SL.lm", "All"),
#'                                     c("SL.lm", "screen.wgtd.ttest4")))
#' sl
#' sl$whichScreen
screen.wgtd.ttest <- function(Y, X, family, obsWeights, id,
							  selector = c("cutoff.k", "cutoff.k.percent"),
							  k = switch(selector,
							 			 cutoff.k = ceiling(0.5*ncol(X)),
							 			 cutoff.k.percent = 0.5,
							 			 NULL),
							  minP = NULL,
							  ...) {

    stat <- "t.value"
    if(!is.null(minP)) {
    	stat <- "p.value"
    }
    default_val <- ifelse(stat == "t.value", 0, 1)
    selector <- match.arg(selector)
    X <- as.data.frame(X)

	if(!is(family, "family")) {
        stop("screen.wgtd.ttest(): please supply a 'family' of gaussian() or binomial().")
    }
    filter_res <- rep(default_val, ncol(X))
	if(family$family == "binomial") {
		filter_res <- apply(X, 2, function(x, Y, obsWeights) {
			if(var(x) <= 0)
				return(default_val)
			x <- split(x, as.factor(Y))
			# no variability in Y:
			if(length(x) < 2)
				return(default_val) # this should probably instead be NA but SL might choke
			obsWeights <- split(obsWeights, as.factor(Y))
	    	weights::wtd.t.test(x[[1]], y = x[[2]],
	    						weight = obsWeights[[1]], weighty = obsWeights[[2]],
	    						samedata = FALSE, alternative = "two.tailed", ...)$coefficients[[stat]]
		}, Y = Y, obsWeights = obsWeights)
	} else {
        stop("screen.wgtd.ttest(): family '", family$family, " not supported.")
    }

	whichScreen <- rep(FALSE, ncol(X))
    if(stat == "p.value") {
	    whichScreen <- (filter_res <= minP)
	}
	if(sum(whichScreen) < k) {
		if(stat == "p.value") {
			warning("Number of variables with P-value less than ", minP,
					" is ", sum(whichScreen), ", which is less than ", k, ".")
			# flip p vals so larger == more 'important'
			filter_res <- 1 - filter_res
		}
		# selector_f <- match.fun(selector)
	    # match.fun has scoping problems :(
	    selector_f <- get(selector)
		subset <- selector_f(attrs = as.data.frame(abs(filter_res)), k = k)
		whichScreen <- colnames(X) %in% subset
	}
	return(whichScreen)
}

#' Weighted t-test screening algorithm
#'
#' Convenience wrapper to perform feature selection according to the ranking of
#' P-values returned from weighted t-tests. Implemented via
#' \code{\link[weights]{wtd.t.test}}.
#'
#' @inheritParams screen.wgtd.ttest
#' @inherit screen.wgtd.ttest return
#' @param k Numeric. Minimum number or proportion of features to select.
#' Passed through to the \code{selector}. For \code{\link[FSelector]{cutoff.k}},
#' this is an integer indicating the number of features to keep from \code{X}.
#' For \code{\link[FSelector]{cutoff.k.percent}}, this is instead the proportion
#' of features to keep. Default: \code{1}. Ignored if at least \code{k} features
#' have P-values at or below \code{minP}.
#' @param minP Numeric. To pass the screen, resulting P-values must not exceed
#' this number. Default is \code{0.1}. Ignored if fewer than \code{k} features
#' have P-values at or below this value.
#' @export
#' @examples
#' # based on example in SuperLearner package
#' set.seed(1)
#' n <- 100
#' p <- 20
#' X <- matrix(rnorm(n*p), nrow = n, ncol = p)
#' X <- data.frame(X)
#' Y <- rbinom(n, 1, plogis(.2*X[, 1] + .1*X[, 2] - .2*X[, 3] + .1*X[, 3]*X[, 4] - .2*abs(X[, 4])))
#' obsWeights <- 1/runif(n)
#' screen.wgtd.ttestP(Y, X, binomial(), obsWeights, seq(n), minP = 0.05)
#'
#' screen.wgtd.ttestP05 <- function(..., minP = 0.05){
#'     screen.wgtd.ttestP(..., minP = minP)
#' }
#'
#' library(SuperLearner)
#' sl = SuperLearner(Y, X, family = binomial(), cvControl = list(V = 2),
#'                   obsWeights = obsWeights,
#'                   SL.library = list(c("SL.lm", "All"),
#'                                     c("SL.lm", "screen.wgtd.ttestP05")))
#' sl
#' sl$whichScreen
screen.wgtd.ttestP <- function(..., selector = "cutoff.k", k = 1, minP = 0.1) {
	screen.wgtd.ttest(..., selector = selector, k = k, minP = minP)
}