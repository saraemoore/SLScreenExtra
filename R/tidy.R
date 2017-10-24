checkColName = function(var, df) {
	if(var %in% colnames(df)) {
		warning("Column name exists. Prefixing with a dot.")
		checkColName(paste0(".", var), df)
	} else {
		return(var)
	}
}

#' Move rownames to a column for a matrix, then convert to data.frame
#'
#' Inspired by \code{\link[tibble]{rownames_to_column}} but not a proper
#' extension of it (yet).
#'
#' @param mat A matrix.
#' @param var Name of new column to use for rownames. If the column exists, a
#' period will be prepended until a unique column name is found.
#' @param stringsAsFactors Logical. When casting to \code{\link{data.frame}},
#' should character vectors be converted to factors?
rownames_to_column.matrix = function(mat, var = "rowname", stringsAsFactors = FALSE) {
	rn = rownames(mat)
	rownames(mat) = NULL
	df = as.data.frame(mat, stringsAsFactors = stringsAsFactors)

	# check if var is already a column and warn if so
	var = checkColName(var, df)
	df[,var] = rn

	if(stringsAsFactors&is.character(df[,var])) {
		df[,var] = factor(df[,var])
	}

	return(df)
}

#' Tidying method(s) for a SuperLearner fit object
#'
#' This method extends \code{\link[broom]{tidy}} to tidy the results from a
#' \code{\link[SuperLearner]{SuperLearner}} fit (screening, prediction, or
#' both) into a summary.
#'
#' @name SuperLearner_tidiers
#' @rdname SuperLearner_tidiers
#'
#' @param x object of class \code{\link[SuperLearner]{SuperLearner}}
#' @param algorithm one of "prediction" (the default), "screening", or "both"
#' (where "both" indicates that information on both "prediction" and
#' "screening" should be reported).
#' @param stringsAsFactors Set to \code{FALSE} by default. Experimental feature
#'  -- note that setting this to \code{TRUE} may produce warnings/errors
#' (and/or may have no effect).
#' @param ... passed through to internal functions handling summaries of
#' feature selection and/or prediction models (depending on the supplied value
#' for \code{algorithm}). See details (below) for specific arguments accepted.
#' @return A \code{data.frame} without rownames. Column names included depend
#' on the supplied value for \code{algorithm}.
#' @details This method can be used to summarize information related to the
#' screening algorithm(s), prediction algorithm(s), or both.
#'
#' @section Optional argument(s):
#'
#' If \code{algorithm} is set to "screening" or "both", the optional argument
#' \code{includeAll} can be set to \code{TRUE} to include results from the
#' pass-thru screening algorithm, \code{"All"}. By default, \code{includeAll}
#' is set to \code{FALSE} and these results are excluded.
#'
#' @section Resulting \code{data.frame}:
#'
#' Columns in resulting \code{data.frame} depend upon selected \code{algorithm}:
#' \describe{
#'     \item{"prediction"}{
#'         One row per element in \code{SL.library}. Five columns:
#'         \describe{
#'             \item{"estimate"}{Coefficient estimate}
#'             \item{"cvRisk"}{Estimate of cross-validated risk.}
#'             \item{"screener"}{Screening algorithm}
#'             \item{"predictor"}{Prediction algorithm}
#'             \item{"discrete"}{Logical. Is this the discrete SuperLearner?}
#'         }
#'     }
#'     \item{"screening"}{Number of rows equal to the product of
#'         [number columns in \code{X}] and
#'         [number of unique screening algorithms in \code{SL.library}, not
#'         including \code{"All"} (by default)].
#'         Three columns:
#'         \describe{
#'             \item{"screener"}{Screening algorithm}
#'             \item{"term"}{Column of \code{X}}
#'             \item{"selected"}{Logical. Did the screening algorithm select
#'                   this column of \code{X}?}
#'         }
#'     }
#'     \item{"both"}{Number of rows equal to the product of
#'         [number columns in \code{X}] and
#'         [number of elements in \code{SL.library}, not including any elements
#'         where the screening algorithm was set to \code{"All"} (by default)].
#'         Seven columns constituting the union of the columns returned
#'         for "prediction" and "screening" (described above).
#'     }
#' }
#' @seealso \code{\link{tidy.CV.SuperLearner}}
#' @importFrom SuperLearner SuperLearner
#' @importFrom broom tidy
#' @export
#' @examples
#' # based on an example in the SuperLearner package
#' set.seed(1)
#' n <- 100
#' p <- 10
#' X <- matrix(rnorm(n*p), nrow = n, ncol = p)
#' X <- data.frame(X)
#' Y <- rbinom(n, 1, plogis(.2*X[, 1] + .1*X[, 2] - .2*X[, 3] + .1*X[, 3]*X[, 4] - .2*abs(X[, 4])))
#'
#' library(SuperLearner)
#' sl = SuperLearner(Y, X, family = binomial(), cvControl = list(V = 2),
#'                   SL.library = list(c("SL.mean", "screen.wgtd.corP"),
#'                                     c("SL.mean", "screen.wgtd.ttest"),
#'                                     c("SL.glm", "screen.wgtd.corP"),
#'                                     c("SL.glm", "screen.wgtd.ttest")))
#'
#' library(broom)
#' tidy(sl)
#' tidy(sl, algorithm = "screening")
#' tidy(sl, algorithm = "both")
#'
tidy.SuperLearner = function(x, algorithm = c("prediction", "screening", "both"),
							 stringsAsFactors = FALSE, ...) {
    algorithm <- match.arg(algorithm)

    switch(algorithm,
    	screening = tidyFeatures.SuperLearner(x, includePred = FALSE,
    										  stringsAsFactors = stringsAsFactors, ...),
    	prediction = tidyModels.SuperLearner(x,
    										 stringsAsFactors = stringsAsFactors, ...),
    	both = tidyFeatures.SuperLearner(x, includePred = TRUE,
    									 stringsAsFactors = stringsAsFactors, ...))
}

#' Tidying method(s) for a CV.SuperLearner fit object
#'
#' This method extends \code{\link[broom]{tidy}} to tidy the results from a
#' \code{\link[SuperLearner]{CV.SuperLearner}} fit (screening, prediction, or
#' both) into a summary.
#'
#' @name CV.SuperLearner_tidiers
#' @rdname CV.SuperLearner_tidiers
#'
#' @param x object of class \code{\link[SuperLearner]{CV.SuperLearner}}
#' @param ... Passed through to \code{\link{tidy.SuperLearner}}. Optional
#' arguments include \code{algorithm}, \code{includeAll}, and
#' \code{stringsAsFactors}. See \code{\link{tidy.SuperLearner}} for
#' documentation.
#' @return A \code{data.frame} without rownames. Column names included depend
#' on the supplied value for the optional \code{\link{tidy.SuperLearner}}
#' argument \code{algorithm}. See \code{\link{tidy.SuperLearner}} for more
#' details on returned \code{data.frame}. Note that the resulting
#' \code{data.frame} from \code{tidy.CV.SuperLearner} will contain one
#' additional column, however: "fold," indicating the (outer) SuperLearner
#' cross-validation fold number.
#' @seealso \code{\link{tidy.SuperLearner}}
#' @importFrom SuperLearner CV.SuperLearner
#' @importFrom broom tidy
#' @importFrom dplyr bind_rows mutate %>%
#' @export
#' @examples
#' # based on an example in the SuperLearner package
#' set.seed(1)
#' n <- 100
#' p <- 10
#' X <- matrix(rnorm(n*p), nrow = n, ncol = p)
#' X <- data.frame(X)
#' Y <- rbinom(n, 1, plogis(.2*X[, 1] + .1*X[, 2] - .2*X[, 3] + .1*X[, 3]*X[, 4] - .2*abs(X[, 4])))
#'
#' library(SuperLearner)
#' cvsl = CV.SuperLearner(Y, X, family = binomial(),
#'                        SL.library = list(c("SL.mean", "screen.FSelector.oneR"),
#'                                          c("SL.mean", "screen.wgtd.ttest"),
#'                                          c("SL.glm", "screen.FSelector.oneR"),
#'                                          c("SL.glm", "screen.wgtd.ttest")),
#'                        cvControl = list(V = 2),
#'                        innerCvControl = list(list(V = 2)))
#'
#' library(broom)
#' tidy(cvsl)
#' tidy(cvsl, algorithm = "screening")
#' tidy(cvsl, algorithm = "both")
#'
tidy.CV.SuperLearner = function(x, ...) {
	resByFold = sapply(x$AllSL, tidy, ..., simplify = FALSE)
	bind_rows(resByFold, .id = "fold") %>%
		mutate(fold = as.numeric(fold))
}

tidyModels <- function (x, ...) {
	UseMethod("tidyModels", x)
}

#' @importFrom dplyr mutate
#' @importFrom SuperLearner SuperLearner
tidyModels.SuperLearner <- function (x, stringsAsFactors = FALSE, ...) {

	res = data.frame(estimate = x$coef,
					 cvRisk = x$cvRisk,
					 screener = x$SL.library$screenAlgorithm[x$SL.library$library$rowScreen],
					 predictor = x$SL.library$library$predAlgorithm,
					 stringsAsFactors = stringsAsFactors)

	# indicate discrete SuperLearner:
	# (there could be ties)
	if(length(which(x$cvRisk <= min(x$cvRisk))) > 1) {
		warning("Multiple prediction algorithms yielded identical cross-validated risk estimates, ",
				"therefore multiple discrete SuperLearners will be indicated.")
	}
	res = mutate(res, discrete = cvRisk <= min(cvRisk))

	return(res)
}

#' @importFrom dplyr bind_rows mutate %>%
#' @importFrom SuperLearner CV.SuperLearner
tidyModels.CV.SuperLearner <- function (x, ...) {
	resByFold = sapply(x$AllSL, tidyModels, ..., simplify = FALSE)
	bind_rows(resByFold, .id = "fold") %>%
		mutate(fold = as.numeric(fold))
}

tidyFeatures <- function (x, ...) {
	UseMethod("tidyFeatures", x)
}

#' @importFrom tidyr gather
#' @importFrom dplyr full_join filter
#' @importFrom SuperLearner SuperLearner
tidyFeatures.SuperLearner = function(x, includePred = FALSE, includeAll = FALSE, stringsAsFactors = FALSE, ...) {

	res = x$whichScreen
	colnames(res) = x$varNames
	# rownames are probably duplicated --
	#     first move rownames to a column, then make into data.frame
	res = rownames_to_column.matrix(res, "screener", stringsAsFactors = stringsAsFactors)
	res = gather(res, term, selected, -screener)

	if(includePred) {
		coefByLibAlg = tidyModels(x, stringsAsFactors = stringsAsFactors)
		res = full_join(res, coefByLibAlg, by = "screener")
	}

	if(!includeAll) {
		res = filter(res, screener!="All")
		if(nrow(res)==0) {
			warning("No features will be returned. If no screening algorithm was used, try includeAll = TRUE.")
		}
	}

	# no rownames on tidy() output
	rownames(res) = NULL
	return(res)
}

#' @importFrom dplyr bind_rows mutate %>%
#' @importFrom SuperLearner CV.SuperLearner
tidyFeatures.CV.SuperLearner = function(x, ...) {

	resByFold = sapply(x$AllSL, tidyFeatures, ..., simplify = FALSE)
	bind_rows(resByFold, .id = "fold") %>%
		mutate(fold = as.numeric(fold))
}
