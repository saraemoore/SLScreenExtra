#' Screening algorithms implemented in the FSelector package
#'
#' A SuperLearner-compatible interface to functions in the
#' \code{\link[FSelector]{FSelector}} package.
#'
#' @param Y Outcome (numeric vector). See \code{\link[SuperLearner]{SuperLearner}}
#' for specifics.
#' @param X Predictor variable(s) (data.frame or matrix). See
#' \code{\link[SuperLearner]{SuperLearner}} for specifics.
#' @param family Error distribution to be used in the model:
#' \code{\link[stats]{gaussian}} or \code{\link[stats]{binomial}}.
#' Currently unused. See \code{\link[SuperLearner]{SuperLearner}}
#' for specifics.
#' @param obsWeights Optional numeric vector of observation weights. Currently
#' unused.
#' @param id Cluster identification variable. Currently unused.
#' @param filter A string corresponding to a feature ranking or selecting function
#' implemented in the FSelector package. One of: \code{\link[FSelector]{cfs}},
#' \code{\link[FSelector]{chi.squared}}, \code{\link[FSelector]{consistency}},
#' \code{\link[FSelector]{gain.ratio}}, \code{\link[FSelector]{information.gain}},
#' \code{\link[FSelector]{linear.correlation}}, \code{\link[FSelector]{oneR}},
#' \code{\link[FSelector]{random.forest.importance}},
#' \code{\link[FSelector]{rank.correlation}}, \code{\link[FSelector]{relief}},
#' \code{\link[FSelector]{symmetrical.uncertainty}}. Default: \code{"cfs"}.
#' Note that "filter" is a misnomer in the case of embedded feature selection
#' methods such as \code{\link[FSelector]{random.forest.importance}}.
#' @param filter_params A named list of tuning parameter arguments specific to
#' the chosen \code{filter}. Default of \code{NULL} should be used when either
#' 1) the chosen \code{filter} does not utilize tuning parameter(s) or
#' 2) the default values should be retained.
#' @param selector A string corresponding to a subset selecting function
#' implemented in the FSelector package. One of:
#' \code{\link[FSelector]{cutoff.biggest.diff}},
#' \code{\link[FSelector]{cutoff.k}}, \code{\link[FSelector]{cutoff.k.percent}},
#' or \code{"all"}. Note that \code{"all"} is a not a function but indicates
#' pass-thru should be performed in the case of a \code{filter} which selects
#' rather than ranks features. Default: \code{"cutoff.biggest.diff"}.
#' @param k Passed through to the \code{selector} in the case where \code{selector} is
#' \code{\link[FSelector]{cutoff.k}} or \code{\link[FSelector]{cutoff.k.percent}}.
#' Otherwise, should remain NULL (the default). For \code{\link[FSelector]{cutoff.k}},
#' this is an integer indicating the number of features to keep from \code{X}.
#' For \code{\link[FSelector]{cutoff.k.percent}}, this is instead the proportion
#' of features to keep.
#' @param verbose Should debugging messages be printed? Default: \code{FALSE}.
#' @param ... Currently unused.
#' @return A logical vector with length equal to \code{ncol(X)}.
#' @import FSelector
#' @importFrom methods is
screen.FSelector <- function(Y, X, family, obsWeights, id,
							 filter = c("cfs", "chi.squared", "consistency", "gain.ratio", "information.gain",
							 			"linear.correlation", "oneR", "random.forest.importance",
							 			"rank.correlation", "relief", "symmetrical.uncertainty"),
							 filter_params = NULL,
							 selector = c("cutoff.biggest.diff", "cutoff.k", "cutoff.k.percent", "all"),
							 k = switch(selector,
							 			cutoff.k = ceiling(0.5*ncol(X)),
							 			cutoff.k.percent = 0.5,
							 			NULL),
							 verbose = FALSE,
							 ...) {

    filter <- match.arg(filter)
    selector <- match.arg(selector)
    if(!is(family, "family")) {
    	stop("screen.FSelector(): please supply a 'family' of gaussian() or binomial().")
    }

    df <- NULL
    if (family$family == "gaussian" & !(filter %in% c("chi.squared", "consistency", "oneR",
    												 "gain.ratio", "information.gain",
    												 "symmetrical.uncertainty"))) {
        df <- data.frame(X, Y = Y)
    } else if (family$family == "binomial" & !grepl("correlation$", filter)) {
        df <- data.frame(X, Y = as.factor(Y))
    } else {
    	stop("screen.FSelector(): family '", family$family,
    		 "' not supported in combination with filter '", filter, "'.")
    }
    y_name <- colnames(df)[ncol(df)]

	filter_f <- match.fun(filter)
	filter_res <- do.call(filter_f,
						  c(list(formula = as.simple.formula(".", y_name), data = df),
						  	filter_params))
	if(verbose) {
		print(filter_res)
	}

	subset <- colnames(df)[-ncol(df)]
	if(selector != "all") {
		selector_f <- match.fun(selector)
		# c(list()) trick so k disappears if NULL
		subset <- do.call(selector_f, c(list(attrs = filter_res), k = k))
	}
	if(verbose) {
		f <- as.simple.formula(subset, y_name)
		print(f)
	}

	whichScreen <- colnames(df) %in% subset
	# remove Y
	return(whichScreen[-length(whichScreen)])
}

#' Correlation Feature Selection (CFS) screening algorithm
#'
#' CFS (Hall, 1999) utilizes \code{\link[FSelector]{best.first.search}} to find
#' columns of \code{X} correlated with \code{Y} but not with one another (i.e.,
#' not redundant). CFS, combined with a search algorithm, does not rank features
#' and therefore does not allow for specification of either the number of
#' features to be chosen (\code{k}) or the method by which they should be chosen
#' (\code{selector}).
#'
#' @inheritParams screen.FSelector
#' @inherit screen.FSelector return
#' @references \url{http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.37.4643}
#' @export
#' @examples
#' data(iris)
#' Y <- as.numeric(iris$Species=="setosa")
#' X <- iris[,-which(colnames(iris)=="Species")]
#' screen.FSelector.cfs(Y, X, binomial())
#'
#' data(mtcars)
#' Y <- mtcars$mpg
#' X <- mtcars[,-which(colnames(mtcars)=="mpg")]
#' screen.FSelector.cfs(Y, X, gaussian())
#'
#' # based on examples in SuperLearner package
#' set.seed(1)
#' n <- 100
#' p <- 20
#' X <- matrix(rnorm(n*p), nrow = n, ncol = p)
#' X <- data.frame(X)
#' Y <- X[, 1] + sqrt(abs(X[, 2] * X[, 3])) + X[, 2] - X[, 3] + rnorm(n)
#'
#' library(SuperLearner)
#' sl = SuperLearner(Y, X, family = gaussian(), cvControl = list(V = 2),
#'                   SL.library = list(c("SL.glm", "All"),
#'                                     c("SL.glm.interaction", "screen.FSelector.cfs")))
#' sl
#' sl$whichScreen
screen.FSelector.cfs <- function(Y, X, family,
							 	 verbose = FALSE,
							 	 ...) {
	screen.FSelector(Y, X, family, filter = "cfs", filter_params = NULL,
					 selector = "all", k = NULL, verbose = verbose, ...)
}

#' Chi-squared test screening algorithm
#'
#' Cramer's V, derived from Pearson's chi-squared statistic, is used to
#' find columns of \code{X} which are associated with \code{Y}. Implemented
#' for \code{binomial()} family only and designed to be used with binary or
#' categorical \code{X}. Continuous \code{X} will be discretized by
#' \code{\link{FSelector}} and \code{\link[RWeka]{Discretize}} using the MDL
#' method (Fayyad & Irani, 1993).
#'
#' @inheritParams screen.FSelector
#' @inherit screen.FSelector return
#' @references \url{http://hdl.handle.net/2014/35171}
#' @export
#' @examples
#' data(iris)
#' Y <- as.numeric(iris$Species=="setosa")
#' X <- iris[,-which(colnames(iris)=="Species")]
#' screen.FSelector.chi.squared(Y, X, binomial(), selector = "cutoff.k.percent", k = 0.5)
#'
#' # based on example in SuperLearner package
#' set.seed(1)
#' n <- 100
#' p <- 20
#' X <- matrix(rnorm(n*p), nrow = n, ncol = p)
#' X <- data.frame(X)
#' Y <- rbinom(n, 1, plogis(.2*X[, 1] + .1*X[, 2] - .2*X[, 3] + .1*X[, 3]*X[, 4] - .2*abs(X[, 4])))
#'
#' library(SuperLearner)
#' sl = SuperLearner(Y, X, family = binomial(), cvControl = list(V = 2),
#'                   SL.library = list(c("SL.lm", "All"),
#'                                     c("SL.lm", "screen.FSelector.chi.squared")))
#' sl
#' sl$whichScreen
screen.FSelector.chi.squared <- function(Y, X, family,
										 selector = c("cutoff.biggest.diff", "cutoff.k", "cutoff.k.percent"),
										 k = switch(selector,
							 						cutoff.k = ceiling(0.5*ncol(X)),
							 						cutoff.k.percent = 0.5,
							 						NULL),
							 			 verbose = FALSE,
							 			 ...) {
    selector <- match.arg(selector)
	screen.FSelector(Y, X, family, filter = "chi.squared", filter_params = NULL,
					 selector = selector, k = k, verbose = verbose, ...)
}

#' Consistency screening algorithm
#'
#' The \code{\link[FSelector]{consistency}} algorithm utilizes
#' \code{\link[FSelector]{best.first.search}} to find the optimal combination of
#' columns of \code{X} to minimize 'inconsistency' with the outcome \code{Y}.
#' Implemented for \code{binomial()} family only and designed to be used with
#' binary or categorical \code{X}. Continuous \code{X} will be discretized by
#' \code{\link{FSelector}} and \code{\link[RWeka]{Discretize}} using the MDL
#' method (Fayyad & Irani, 1993). Search algorithms do not rank features and
#' therefore this algorithm does not allow for specification of either the
#' number of features to be chosen (\code{k}) or the method by which they
#' should be chosen (\code{selector}).
#'
#' @inheritParams screen.FSelector
#' @inherit screen.FSelector return
#' @references \url{http://hdl.handle.net/2014/35171}
#' @export
#' @examples
#' data(iris)
#' Y <- as.numeric(iris$Species=="setosa")
#' X <- iris[,-which(colnames(iris)=="Species")]
#' screen.FSelector.consistency(Y, X, binomial())
#'
#' # based on example in SuperLearner package
#' set.seed(1)
#' n <- 100
#' p <- 20
#' X <- matrix(rnorm(n*p), nrow = n, ncol = p)
#' X <- data.frame(X)
#' Y <- rbinom(n, 1, plogis(.2*X[, 1] + .1*X[, 2] - .2*X[, 3] + .1*X[, 3]*X[, 4] - .2*abs(X[, 4])))
#'
#' library(SuperLearner)
#' sl = SuperLearner(Y, X, family = binomial(), cvControl = list(V = 2),
#'                   SL.library = list(c("SL.lm", "All"),
#'                                     c("SL.lm", "screen.FSelector.consistency")))
#' sl
#' sl$whichScreen
screen.FSelector.consistency <- function(Y, X, family,
							 			 verbose = FALSE,
							 			 ...) {
	screen.FSelector(Y, X, family, filter = "consistency", filter_params = NULL,
					 selector = "all", k = NULL, verbose = verbose, ...)
}

#' Entropy-based screening algorithms
#'
#' Information gain, gain ratio, and symmetrical uncertainty scores are
#' calculated from the Shannon entropy of \code{X} and \code{Y}. Information
#' gain (\code{\link[FSelector]{information.gain}}) (or, equivalently, mutual
#' information) is a measure of entropy reduction achieved by the feature with
#' regard to the outcome, \code{Y}. The information gain ratio
#' (\code{\link[FSelector]{gain.ratio}}) is a normalized version of information
#' gain, normalized by the entropy of the feature. Symmetrical uncertainty
#' (\code{\link[FSelector]{symmetrical.uncertainty}}) is a normalized and
#' bias-corrected version of information gain. Implemented for \code{binomial()}
#' family only and designed to be used with binary or categorical \code{X}.
#' Continuous \code{X} will be discretized by \code{\link{FSelector}} and
#' \code{\link[RWeka]{Discretize}} using the MDL method (Fayyad & Irani, 1993).
#'
#' @param filter Character string. One of: \code{"symmetrical.uncertainty"} (default),
#' \code{"gain.ratio"}, or \code{"information.gain"}
#' @param unit Unit in which entropy is measured by
#' \code{\link[entropy]{entropy}}. Character string. One of: \code{"log"}
#' (default), \code{"log2"}, or \code{"log10"}.
#' @inheritParams screen.FSelector
#' @inherit screen.FSelector return
#' @references \url{http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.37.4643}
#' \url{http://hdl.handle.net/2014/35171}
#' @export
#' @examples
#' data(iris)
#' Y <- as.numeric(iris$Species=="setosa")
#' X <- iris[,-which(colnames(iris)=="Species")]
#' screen.FSelector.entropy(Y, X, binomial(), selector = "cutoff.k.percent", k = 0.75)
#'
#' # based on example in SuperLearner package
#' set.seed(1)
#' n <- 100
#' p <- 20
#' X <- matrix(rnorm(n*p), nrow = n, ncol = p)
#' X <- data.frame(X)
#' Y <- rbinom(n, 1, plogis(.2*X[, 1] + .1*X[, 2] - .2*X[, 3] + .1*X[, 3]*X[, 4] - .2*abs(X[, 4])))
#'
#' library(SuperLearner)
#' sl = SuperLearner(Y, X, family = binomial(), cvControl = list(V = 2),
#'                   SL.library = list(c("SL.lm", "All"),
#'                                     c("SL.lm", "screen.FSelector.entropy")))
#' sl
#' sl$whichScreen
screen.FSelector.entropy <- function(Y, X, family,
									 filter = c("symmetrical.uncertainty", "gain.ratio", "information.gain"),
									 unit = formals(information.gain)$unit,
									 selector = c("cutoff.biggest.diff", "cutoff.k", "cutoff.k.percent"),
									 k = switch(selector,
							 					cutoff.k = ceiling(0.5*ncol(X)),
							 					cutoff.k.percent = 0.5,
							 					NULL),
							 		 verbose = FALSE,
							 		 ...) {
    filter <- match.arg(filter)
    selector <- match.arg(selector)
	screen.FSelector(Y, X, family,
					 filter = filter, filter_params = list(unit = unit),
					 selector = selector, k = k,
					 verbose = verbose, ...)
}

#' Correlation screening algorithm
#'
#' Features are ranked and selected based on the strength of their correlation
#' with the outcome \code{Y}. Implemented for \code{gaussian()} family only.
#'
#' @param filter Character string. One of \code{"linear.correlation"} (Pearson)
#' or \code{"rank.correlation"} (Spearman).
#' @inheritParams screen.FSelector
#' @inherit screen.FSelector return
#' @export
#' @examples
#' data(mtcars)
#' Y <- mtcars$mpg
#' X <- mtcars[,-which(colnames(mtcars)=="mpg")]
#' screen.FSelector.correlation(Y, X, gaussian(), filter = "rank.correlation")
#'
#' # based on examples in SuperLearner package
#' set.seed(1)
#' n <- 100
#' p <- 20
#' X <- matrix(rnorm(n*p), nrow = n, ncol = p)
#' X <- data.frame(X)
#' Y <- X[, 1] + sqrt(abs(X[, 2] * X[, 3])) + X[, 2] - X[, 3] + rnorm(n)
#'
#' library(SuperLearner)
#' sl = SuperLearner(Y, X, family = gaussian(), cvControl = list(V = 2),
#'                   SL.library = list(c("SL.glm", "All"),
#'                                     c("SL.glm", "screen.FSelector.correlation")))
#' sl
#' sl$whichScreen
screen.FSelector.correlation <- function(Y, X, family,
								 filter = c("linear.correlation", "rank.correlation"),
								 selector = c("cutoff.biggest.diff", "cutoff.k", "cutoff.k.percent"),
								 k = switch(selector,
							 				cutoff.k = ceiling(0.5*ncol(X)),
							 				cutoff.k.percent = 0.5,
							 				NULL),
							 	 verbose = FALSE,
							 	 ...) {
    filter <- match.arg(filter)
    selector <- match.arg(selector)
	screen.FSelector(Y, X, family,
					 filter = filter, filter_params = NULL,
					 selector = selector, k = k,
					 verbose = verbose, ...)
}

#' OneR screening algorithm
#'
#' AKA "1R" or "One Rule" (Holte, 1993). Features ranked and selected by the
#' accuracy of their predictions of the outcome \code{Y} based on a single
#' rule (or single-level decision tree). Implemented for \code{binomial()}
#' family only and designed to be used with binary or categorical \code{X}.
#' Continuous \code{X} will be discretized by \code{\link{FSelector}} and
#' \code{\link[RWeka]{Discretize}} using the MDL method (Fayyad & Irani,
#' 1993).
#'
#' @inheritParams screen.FSelector
#' @inherit screen.FSelector return
#' @references \url{https://doi.org/10.1023/A:1022631118932}
#' \url{http://hdl.handle.net/2014/35171}
#' @export
#' @examples
#' data(iris)
#' Y <- as.numeric(iris$Species=="setosa")
#' X <- iris[,-which(colnames(iris)=="Species")]
#' screen.FSelector.oneR(Y, X, binomial(), selector = "cutoff.k.percent", k = 0.5)
#'
#' # based on example in SuperLearner package
#' set.seed(1)
#' n <- 100
#' p <- 20
#' X <- matrix(rnorm(n*p), nrow = n, ncol = p)
#' X <- data.frame(X)
#' Y <- rbinom(n, 1, plogis(.2*X[, 1] + .1*X[, 2] - .2*X[, 3] + .1*X[, 3]*X[, 4] - .2*abs(X[, 4])))
#'
#' library(SuperLearner)
#' sl = SuperLearner(Y, X, family = binomial(), cvControl = list(V = 2),
#'                   SL.library = list(c("SL.lm", "All"),
#'                                     c("SL.lm", "screen.FSelector.oneR")))
#' sl
#' sl$whichScreen
screen.FSelector.oneR <- function(Y, X, family,
								  selector = c("cutoff.biggest.diff", "cutoff.k", "cutoff.k.percent"),
								  k = switch(selector,
								  			 cutoff.k = ceiling(0.5*ncol(X)),
								  			 cutoff.k.percent = 0.5,
								  			 NULL),
								  verbose = FALSE,
								  ...) {
    selector <- match.arg(selector)
	screen.FSelector(Y, X, family,
					 filter = "oneR", filter_params = NULL,
					 selector = selector, k = k,
					 verbose = verbose, ...)
}

#' Random Forest screening algorithm
#'
#' The \code{\link[FSelector]{random.forest.importance}} algorithm uses
#' \code{\link[randomForest]{randomForest}} (with \code{ntree = 1000}) to
#' estimate the specified type of importance for each column of \code{X}.
#'
#' @param type Importance type. Integer: \code{1}, indicating mean decrease in
#' accuracy (for \code{binomial()} family) or percent increase in mean squared
#' error (for \code{gaussian()} family) when comparing predictions using
#' the original variable versus a permuted version of the variable (column of
#' \code{X}), or \code{2}, indicating the increase in
#' node purity achieved by splitting on that column of \code{X} (for
#' \code{binomial()} family, measured by Gini index; for \code{gaussian()},
#' measured by residual sum of squares). For default value, see
#' \code{\link[FSelector]{random.forest.importance}}.
#' @inheritParams screen.FSelector
#' @inherit screen.FSelector return
#' @export
#' @examples
#' data(iris)
#' Y <- as.numeric(iris$Species=="setosa")
#' X <- iris[,-which(colnames(iris)=="Species")]
#' screen.FSelector.random.forest.importance(Y, X, binomial(), selector = "cutoff.k.percent", k = 0.75)
#'
#' data(mtcars)
#' Y <- mtcars$mpg
#' X <- mtcars[,-which(colnames(mtcars)=="mpg")]
#' screen.FSelector.random.forest.importance(Y, X, gaussian(), type = 2)
#'
#' # based on examples in SuperLearner package
#' set.seed(1)
#' n <- 100
#' p <- 20
#' X <- matrix(rnorm(n*p), nrow = n, ncol = p)
#' X <- data.frame(X)
#' Y <- X[, 1] + sqrt(abs(X[, 2] * X[, 3])) + X[, 2] - X[, 3] + rnorm(n)
#'
#' library(SuperLearner)
#' sl = SuperLearner(Y, X, family = gaussian(), cvControl = list(V = 2),
#'                   SL.library = list(c("SL.glm", "All"),
#'                                     c("SL.glm", "screen.FSelector.random.forest.importance")))
#' sl
#' sl$whichScreen
screen.FSelector.random.forest.importance <- function(Y, X, family,
													  type = formals(random.forest.importance)$importance.type,
													  selector = c("cutoff.biggest.diff", "cutoff.k", "cutoff.k.percent"),
													  k = switch(selector,
													  			 cutoff.k = ceiling(0.5*ncol(X)),
													  			 cutoff.k.percent = 0.5,
													  			 NULL),
													  verbose = FALSE,
													  ...) {
    selector <- match.arg(selector)
	screen.FSelector(Y, X, family,
					 filter = "random.forest.importance",
					 filter_params = list(importance.type = type),
					 selector = selector, k = k,
					 verbose = verbose, ...)
}

#' RReliefF screening algorithm
#'
#' The \code{\link[FSelector]{relief}} algorithm implements the RReliefF
#' (Robnik-Sikonja & Kononenko, 1997) feature quality estimation algorithm,
#' an extension to ReliefF (Kononenko, 1994) and Relief (Kira & Rendell, 1992)
#' algorithms. RReliefF is compatible with both classification and regression
#' problems and is well-suited to \code{X} with strong associations between
#' features.
#'
#' @param neighbours.count Number of neighboring observations to find for each
#' observation sampled from \code{X}
#' @param sample.size Number of observations to sample from \code{X}
#' @inheritParams screen.FSelector
#' @inherit screen.FSelector return
#' @references \url{https://www.aaai.org/Library/AAAI/1992/aaai92-020.php},
#' \url{https://doi.org/10.1007/3-540-57868-4_57},
#' \url{http://dl.acm.org/citation.cfm?id=645526.657141}
#' @export
#' @examples
#' data(iris)
#' Y <- as.numeric(iris$Species=="setosa")
#' X <- iris[,-which(colnames(iris)=="Species")]
#' screen.FSelector.relief(Y, X, binomial(), selector = "cutoff.k.percent", k = 0.75)
#'
#' data(mtcars)
#' Y <- mtcars$mpg
#' X <- mtcars[,-which(colnames(mtcars)=="mpg")]
#' screen.FSelector.relief(Y, X, gaussian(), neighbours.count = 3, sample.size = 15)
#'
#' # based on examples in SuperLearner package
#' set.seed(1)
#' n <- 100
#' p <- 20
#' X <- matrix(rnorm(n*p), nrow = n, ncol = p)
#' X <- data.frame(X)
#' Y <- X[, 1] + sqrt(abs(X[, 2] * X[, 3])) + X[, 2] - X[, 3] + rnorm(n)
#'
#' library(SuperLearner)
#' sl = SuperLearner(Y, X, family = gaussian(), cvControl = list(V = 2),
#'                   SL.library = list(c("SL.glm", "All"),
#'                                     c("SL.glm", "screen.FSelector.relief")))
#' sl
#' sl$whichScreen
screen.FSelector.relief <- function(Y, X, family,
									neighbours.count = formals(relief)$neighbours.count,
									sample.size = formals(relief)$sample.size,
									selector = c("cutoff.biggest.diff", "cutoff.k", "cutoff.k.percent"),
									k = switch(selector,
											   cutoff.k = ceiling(0.5*ncol(X)),
											   cutoff.k.percent = 0.5,
											   NULL),
									verbose = FALSE,
									...) {
    selector <- match.arg(selector)
	screen.FSelector(Y, X, family,
					 filter = "relief",
					 filter_params = list(neighbours.count = neighbours.count,
					 					  sample.size = sample.size),
					 selector = selector, k = k,
					 verbose = verbose, ...)
}
