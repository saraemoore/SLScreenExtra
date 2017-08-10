#' Weighted regularized regression screening algorithm
#'
#' Performs feature selection via \code{\link[glmnet]{cv.glmnet}}.
#'
#' @param Y Outcome (numeric vector). See \code{\link[SuperLearner]{SuperLearner}}
#' for specifics.
#' @param X Predictor variable(s) (data.frame or matrix). See
#' \code{\link[SuperLearner]{SuperLearner}} for specifics.
#' @param family Error distribution to be used in the model: \code{"gaussian"}
#' or \code{"binomial"}.  See \code{\link[SuperLearner]{SuperLearner}} for
#' specifics.
#' @param obsWeights Optional numeric vector of observation weights. See
#' \code{\link[SuperLearner]{SuperLearner}} for specifics.
#' @param id Cluster identification variable. Currently unused.
#' @param alpha The elasticnet mixing parameter. Default for
#' \code{screen.wgtd.elasticnet} is (arbitrarily) \code{0.5}. Forced to
#' \code{1} for \code{screen.wgtd.lasso} and \code{0} for
#' \code{screen.wgtd.ridge}. See \code{\link[glmnet]{glmnet}} for specifics.
#' @param minscreen Minimum number of features to select (aka rank). Only used
#' if less than this number of features are selected using \code{minPvalue}.
#' @param nfolds Number of cross-validation folds to use when choosing optimal
#' \code{lambda}. Default is \code{10}. See \code{\link[glmnet]{cv.glmnet}} for
#' specifics.
#' @param nlambda Number of \code{lambda} values to try. Default is \code{100}.
#' See \code{\link[glmnet]{glmnet}} for specifics.
#' @param ... Currently unused.
#' @return A logical vector with length equal to \code{ncol(X)}
#' @export
#' @importFrom glmnet cv.glmnet coef.glmnet
#' @importFrom stats coef model.matrix
screen.wgtd.elasticnet <- function (Y, X, family, obsWeights, id, alpha = 0.5, minscreen = 2, nfolds = 10, nlambda = 100, ...) {
    if (!is.matrix(X)) {
        X <- model.matrix(~-1 + ., X)
    }
    fitCV <- glmnet::cv.glmnet(x = X, y = Y, weights = obsWeights,
                               lambda = NULL, type.measure = "deviance",
                               nfolds = nfolds, family = family$family,
                               alpha = alpha, nlambda = nlambda)
    # remove intercept:
    whichVariable <- (as.numeric(coef(fitCV$glmnet.fit, s = fitCV$lambda.min))[-1] != 0)
    if (sum(whichVariable) < minscreen) {
        warning("fewer than minscreen variables passed the glmnet screen, increased lambda to allow minscreen variables")
        sumCoef <- apply(as.matrix(fitCV$glmnet.fit$beta), 2, function(x) sum((x != 0)))
        newCut <- which.max(sumCoef >= minscreen)
        whichVariable <- (as.matrix(fitCV$glmnet.fit$beta)[, newCut] != 0)
    }
    return(whichVariable)
}

#' @rdname screen.wgtd.elasticnet
#' @export
screen.wgtd.lasso <- function (Y, X, family, obsWeights, id, ...) {
    screen.wgtd.elasticnet(Y, X, family, obsWeights, id, alpha = 1, ...)
}

#' @rdname screen.wgtd.elasticnet
#' @export
screen.wgtd.ridge <- function (Y, X, family, obsWeights, id, ...) {
    screen.wgtd.elasticnet(Y, X, family, obsWeights, id, alpha = 0, ...)
}
