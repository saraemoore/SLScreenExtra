#' Non-linear regression screening algorithm
#'
#' Performs feature selection via "Multivariate Adaptive Regression Splines"/
#' "Fast MARS" using \code{\link[earth]{earth}}'s implementation.
#' @inheritParams screen.randomForest.imp
#' @inherit screen.randomForest.imp return
#' @param obsWeights Passed on via \code{\link[earth]{earth}}'s \code{weights}
#' argument. Expect slower computation time if \code{obsWeights} are provided.
#' @param importanceType Variable importance criterion. One of:
#' \code{"nsubsets"} ("number of subsets"), \code{"rss"}, or \code{"gcv"}.
#' @param degree Maximum degree of interaction. Default: 2. 1 would indicate
#' no interaction terms should be used.
#' @param penalty Generalized Cross Validation (GCV) penalty per knot.
#' Default: 3.
#' @param kForward Maximum number of terms created by the forward pass
#' (including the intercept). Default: twice the number of features (in
#' \code{X}) plus one OR 21 -- whichever is greater.
#' @param pMethod Pruning method. Default: \code{"cv"}: select the number of
#' terms yielding the maximum mean out-of-fold R-Squared over the
#' cross-validated model fits (CVRSq). See \code{\link[earth]{earth}} for other
#' possible values.
#' @param nFold Number of cross-validation folds. Must be \code{>0} if
#' \code{pmethod = "cv"}. Default: 5.
#' @param ... Currently unused.
#' @export
#' @importFrom earth earth evimp
#' @importFrom methods is
#' @importFrom stats gaussian binomial
#' @importFrom FSelector cutoff.biggest.diff cutoff.k cutoff.k.percent
#' @examples
#' data(iris)
#' Y <- as.numeric(iris$Species=="setosa")
#' X <- iris[,-which(colnames(iris)=="Species")]
#' screen.earth(Y, X, binomial(), selector = "cutoff.k.percent", k = 0.75)
#'
#' data(mtcars)
#' Y <- mtcars$mpg
#' X <- mtcars[,-which(colnames(mtcars)=="mpg")]
#' screen.earth(Y, X, gaussian(), importanceType = "rss")
#'
#' # based on examples in SuperLearner package
#' set.seed(1)
#' n <- 250
#' p <- 20
#' X <- matrix(rnorm(n*p), nrow = n, ncol = p)
#' X <- data.frame(X)
#' Y <- X[, 1] + sqrt(abs(X[, 2] * X[, 3])) + X[, 2] - X[, 3] + rnorm(n)
#'
#' library(SuperLearner)
#' sl = SuperLearner(Y, X, family = gaussian(), cvControl = list(V = 2),
#'                   SL.library = list(c("SL.glm", "All"),
#'                                     c("SL.glm", "screen.earth")))
#' sl
#' sl$whichScreen
screen.earth <- function (Y, X, family, obsWeights, id,
                          selector = c("cutoff.biggest.diff", "cutoff.k", "cutoff.k.percent"),
                          k = switch(selector,
                                     cutoff.k = ceiling(0.5*ncol(X)),
                                     cutoff.k.percent = 0.5,
                                     NULL),
                          importanceType = c("nsubsets", "rss", "gcv"),
                          degree = 2, penalty = 3, kForward = max(21, 2 * ncol(X) + 1),
                          pMethod = "cv", nFold = 5, ...)
{
    # guarantee column names
    X <- as.data.frame(X)

    selector <- match.arg(selector)
    importanceType <- match.arg(importanceType)
    if(missing(obsWeights)) {
        obsWeights <- NULL
    }

    if(!is(family, "family")) {
        stop("screen.earth(): please supply a 'family' of gaussian() or binomial().")
    } else if(!(family$family %in% c("binomial", "gaussian"))) {
        stop("screen.earth(): family '", family$family, " not supported.")
    }

    earth_fit <- earth::earth(x = X, y = Y, weights = obsWeights,
                              glm = list(family = match.fun(family$family)),
                              degree = degree, nk = kForward,
                              penalty = penalty, pmethod = pMethod,
                              nfold = nFold)

    earth_imp <- earth::evimp(earth_fit, trim = FALSE)
    # note that all measures are potentially subject to ties.
    # usually they will agree, but not always.
    # bigger is better for all three criteria.
    filter_res <- earth_imp[,importanceType, drop = FALSE]
    rownames(filter_res) <- colnames(X)[earth_imp[,"col"]]

    # selector_f <- match.fun(selector)
    # match.fun has scoping problems
    selector_f <- get(selector)
    # c(list()) trick so k disappears if NULL
    subset <- do.call(selector_f, c(list(attrs = as.data.frame(filter_res)), k = k))
    # filter_res <- earth_imp[rev(order(filter_res)), "col"]

    # whichVariable <- rep(FALSE, ncol(X))
    # kthVar <- min(k, nrow(filter_res))
    # whichVariable[filter_res[seq(kthVar), "col"]] <- TRUE
    # return(whichVariable)

    return(colnames(X) %in% subset)

}

#' Non-linear regression screening algorithm
#'
#' Performs feature selection via "Multivariate Adaptive Regression Splines"/
#' "Fast MARS" using \code{\link[earth]{earth}}'s implementation.
#' @param ... Arguments passed on to \code{\link{screen.earth}}.
#' @param pMethod Pruning method. Default: \code{"backward"}.
#' @param nFold Number of cross-validation folds. Default: 0 (cross-validation
#' disabled).
#' @export
#' @examples
#' data(iris)
#' Y <- as.numeric(iris$Species=="setosa")
#' X <- iris[,-which(colnames(iris)=="Species")]
#' screen.earth.backwardprune(Y, X, binomial())
#'
#' data(mtcars)
#' Y <- mtcars$mpg
#' X <- mtcars[,-which(colnames(mtcars)=="mpg")]
#' screen.earth.backwardprune(Y, X, gaussian())
#'
#' # based on examples in SuperLearner package
#' set.seed(1)
#' n <- 250
#' p <- 20
#' X <- matrix(rnorm(n*p), nrow = n, ncol = p)
#' X <- data.frame(X)
#' Y <- X[, 1] + sqrt(abs(X[, 2] * X[, 3])) + X[, 2] - X[, 3] + rnorm(n)
#'
#' library(SuperLearner)
#' sl = SuperLearner(Y, X, family = gaussian(), cvControl = list(V = 2),
#'                   SL.library = list(c("SL.glm", "All"),
#'                                     c("SL.glm", "screen.earth.backwardprune")))
#' sl
#' sl$whichScreen
screen.earth.backwardprune = function(..., pMethod = "backward", nFold = 0){
    screen.earth(..., pMethod = pMethod, nFold = nFold)
}
