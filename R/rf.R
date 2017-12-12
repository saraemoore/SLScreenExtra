#' "Best of both worlds" Random Forest screening algorithm
#'
#' Customizability of \code{\link[SuperLearner]{screen.randomForest}} combined
#' with the \code{\link[FSelector]{cutoff}} selectors of
#' \code{\link[FSelector]{FSelector}}.
#' @inheritParams screen.FSelector
#' @inherit screen.FSelector return
#' @param selector A string corresponding to a subset selecting function
#' implemented in the FSelector package. One of:
#' \code{\link[FSelector]{cutoff.biggest.diff}} (default),
#' \code{\link[FSelector]{cutoff.k}}, or
#' \code{\link[FSelector]{cutoff.k.percent}}.
#' @param nTree Integer. Number of trees. Default: 1000.
#' @param mTry Integer. Number of columns of \code{X} sampled at each split.
#' Default: square root (\code{gaussian()} family) or one third
#' (\code{binomial()} family) of total number of features, rounded down.
#' @param nodeSize Integer. Minimum number of observations in terminal nodes.
#' Default: 5 (\code{gaussian()} family) or 1 (\code{binomial()} family).
#' @param importanceType Importance type. \code{"permutation"} (default) indicates
#' mean decrease in accuracy (for \code{binomial()} family) or percent increase
#' in mean squared error (for \code{gaussian()} family) when comparing
#' predictions using the original variable versus a permuted version of the
#' variable (column of \code{X}). \code{"impurity"} indicates increase in
#' node purity achieved by splitting on that column of \code{X} (for
#' \code{binomial()} family, measured by Gini index; for \code{gaussian()},
#' measured by residual sum of squares). See
#' \code{\link[randomForest]{randomForest}} for more details, where
#' \code{"permutation"} corresponds to \code{type = 1} and \code{"impurity"}
#' corresponds to \code{type = 2}.
#' @param maxNodes Maximum number of terminal nodes allowed in a tree. Default
#' (\code{NULL}) indicates that trees should be grown to maximum possible size.
#' See \code{\link[randomForest]{randomForest}} for more details.
#' @importFrom randomForest randomForest importance
#' @importFrom methods is
#' @importFrom FSelector cutoff.biggest.diff cutoff.k cutoff.k.percent
#' @export
#' @examples
#' data(iris)
#' Y <- as.numeric(iris$Species=="setosa")
#' X <- iris[,-which(colnames(iris)=="Species")]
#' screen.randomForest.imp(Y, X, binomial(), selector = "cutoff.k.percent", k = 0.75)
#'
#' data(mtcars)
#' Y <- mtcars$mpg
#' X <- mtcars[,-which(colnames(mtcars)=="mpg")]
#' screen.randomForest.imp(Y, X, gaussian(), importanceType = "impurity")
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
#'                                     c("SL.glm", "screen.randomForest.imp")))
#' sl
#' sl$whichScreen
screen.randomForest.imp = function (Y, X, family, obsWeights, id,
                                    selector = c("cutoff.biggest.diff", "cutoff.k", "cutoff.k.percent"),
                                    k = switch(selector,
                                               cutoff.k = ceiling(0.5*ncol(X)),
                                               cutoff.k.percent = 0.5,
                                               NULL),
                                    nTree = 1000,
                                    mTry = ifelse(family$family == "gaussian", floor(sqrt(ncol(X))), max(floor(ncol(X)/3), 1)),
                                    nodeSize = ifelse(family$family == "gaussian", 5, 1),
                                    importanceType = c("permutation", "impurity"),
                                    maxNodes = NULL,
                                    verbose = FALSE,
                                    ...)
{
    importanceType <- match.arg(importanceType)
    importanceType <- ifelse(importanceType == "permutation", 1, 2)

    selector <- match.arg(selector)
    X <- as.data.frame(X)
    if(!is(family, "family")) {
        stop("screen.randomForest.imp(): please supply a 'family' of gaussian() or binomial().")
    } else if(!(family$family %in% c("binomial", "gaussian"))) {
        stop("screen.randomForest.imp(): family '", family$family, " not supported.")
    }

    if (family$family == "binomial") {
        Y <- as.factor(Y)
        # importance type 1 is MeanDecreaseAccuracy for classification
        # importance type 1 is %IncMSE for regression
    }
    rf_fit <- randomForest(X, y = Y,
                           ntree = nTree, mtry = mTry,
                           nodesize = nodeSize, maxnodes = maxNodes,
                           importance = TRUE, keep.forest = FALSE)

    ## previously:
    ## returned all variables passing a certain accuracy criteria
    ## rather than top k variables (or whatever)
    ## however, danger of returning zero features here.
    # imp.cutoff = ifelse(family$family == "gaussian", 1, 0.005)
    # return(as.vector(importance(rf_fit, type = importanceType) > imp.cutoff))

    filter_res = as.data.frame(randomForest::importance(rf_fit, type = importanceType))
    # selector_f <- match.fun(selector)
    # match.fun has scoping problems
    selector_f <- get(selector)
    # c(list()) trick so k disappears if NULL
    subset <- do.call(selector_f, c(list(attrs = filter_res), k = k))
    if(verbose) {
        f <- as.simple.formula(subset, "Y")
        print(f)
    }

    return(colnames(X) %in% subset)
}


#' Screen features via a fast implementation of Random Forest
#'
#' Speed up \code{\link[SuperLearner]{screen.randomForest}} or
#' \code{\link{screen.randomForest.imp}}. Uses the
#' \code{\link[FSelector]{cutoff}} selectors.
#' @inheritParams screen.randomForest.imp
#' @inherit screen.randomForest.imp return
#' @param importanceType Importance type. \code{"permutation"} (default) indicates
#' mean decrease in accuracy (for \code{binomial()} family) or percent increase
#' in mean squared error (for \code{gaussian()} family) when comparing
#' predictions using the original variable versus a permuted version of the
#' variable (column of \code{X}). \code{"impurity"} indicates increase in
#' node purity achieved by splitting on that column of \code{X} (for
#' \code{binomial()} family, measured by Gini index; for \code{gaussian()},
#' measured by variance of the responses). See
#' \code{\link[ranger]{ranger}} for more details.
#' @param scalePermutationImportance Scale permutation importance by standard
#' error. Ignored if \code{importanceType = "impurity"}. See
#' \code{\link[ranger]{ranger}} for more details.
#' @param probabilityTrees Logical. If family is \code{binomial()} and
#' \code{probabilityTrees} is FALSE (the default), classification trees are
#' grown. If family is \code{binomial()} and
#' \code{probabilityTrees} is TRUE (the default), probability trees are
#' grown (Malley et al., 2012). Ignored if family is \code{gaussian()}, for
#' which regression trees are always grown. See \code{\link[ranger]{ranger}}
#' for more details.
#' @param numThreads Number of threads. Default: 1.
#' @importFrom ranger ranger
#' @importFrom methods is
#' @importFrom FSelector cutoff.biggest.diff cutoff.k cutoff.k.percent
#' @references \url{http://dx.doi.org/10.18637/jss.v077.i01}
#' \url{http://dx.doi.org/10.1023/A:1010933404324}
#' \url{http://dx.doi.org/10.3414/ME00-01-0052}
#' @export
#' @examples
#' data(iris)
#' Y <- as.numeric(iris$Species=="setosa")
#' X <- iris[,-which(colnames(iris)=="Species")]
#' screen.ranger(Y, X, binomial(), selector = "cutoff.k.percent", k = 0.75)
#'
#' data(mtcars)
#' Y <- mtcars$mpg
#' X <- mtcars[,-which(colnames(mtcars)=="mpg")]
#' screen.ranger(Y, X, gaussian(), importanceType = "impurity")
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
#'                                     c("SL.glm", "screen.ranger")))
#' sl
#' sl$whichScreen
screen.ranger <- function(Y, X, family,
                          selector = c("cutoff.biggest.diff", "cutoff.k", "cutoff.k.percent"),
                          k = switch(selector,
                                     cutoff.k = ceiling(0.5*ncol(X)),
                                     cutoff.k.percent = 0.5,
                                     NULL),
                          nTree = 1000,
                          mTry = ifelse(family$family == "gaussian", floor(sqrt(ncol(X))), max(floor(ncol(X)/3), 1)),
                          nodeSize = ifelse(family$family == "gaussian", 5, 1),
                          importanceType = c("permutation", "impurity"),
                          scalePermutationImportance = TRUE,
                          probabilityTrees = FALSE,
                          numThreads = 1,
                          verbose = FALSE,
                          ...) {

    importanceType <- match.arg(importanceType)
    selector <- match.arg(selector)
    if(!is(family, "family")) {
        stop("screen.ranger(): please supply a 'family' of gaussian() or binomial().")
    }

    df <- NULL
    if (family$family == "gaussian") {
        df <- data.frame(X, Y = Y)
        probability <- FALSE
    } else if (family$family == "binomial") {
        df <- data.frame(X, Y = as.factor(Y))
    } else {
        stop("screen.ranger(): family '", family$family, "' not supported.")
    }
    y_name <- colnames(df)[ncol(df)]

    rf_fit <- ranger(data = df, dependent.variable.name = y_name,
                     num.trees = nTree, mtry = mTry,
                     importance = importanceType, probability = probabilityTrees,
                     write.forest = FALSE, min.node.size = nodeSize,
                     scale.permutation.importance = scalePermutationImportance,
                     num.threads = numThreads, verbose = verbose)

    filter_res = as.data.frame(rf_fit$variable.importance)
    # selector_f <- match.fun(selector)
    # match.fun has scoping problems
    selector_f <- get(selector)
    # c(list()) trick so k disappears if NULL
    subset <- do.call(selector_f, c(list(attrs = filter_res), k = k))

    whichScreen <- colnames(df) %in% subset
    # remove Y
    return(whichScreen[-length(whichScreen)])
}
