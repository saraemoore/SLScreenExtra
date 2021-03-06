% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fselector.R
\name{screen.FSelector.relief}
\alias{screen.FSelector.relief}
\title{RReliefF screening algorithm}
\usage{
screen.FSelector.relief(
  Y,
  X,
  family,
  neighbours.count = formals(relief)$neighbours.count,
  sample.size = formals(relief)$sample.size,
  selector = c("cutoff.biggest.diff", "cutoff.k", "cutoff.k.percent"),
  k = switch(selector, cutoff.k = ceiling(0.5 * ncol(X)), cutoff.k.percent = 0.5, NULL),
  verbose = FALSE,
  ...
)
}
\arguments{
\item{Y}{Outcome (numeric vector). See \code{\link[SuperLearner]{SuperLearner}}
for specifics.}

\item{X}{Predictor variable(s) (data.frame or matrix). See
\code{\link[SuperLearner]{SuperLearner}} for specifics.}

\item{family}{Error distribution to be used in the model:
\code{\link[stats]{gaussian}} or \code{\link[stats]{binomial}}.
Currently unused. See \code{\link[SuperLearner]{SuperLearner}}
for specifics.}

\item{neighbours.count}{Number of neighboring observations to find for each
observation sampled from \code{X}}

\item{sample.size}{Number of observations to sample from \code{X}}

\item{selector}{A string corresponding to a subset selecting function
implemented in the FSelector package. One of:
\code{\link[FSelector]{cutoff.biggest.diff}},
\code{\link[FSelector]{cutoff.k}}, \code{\link[FSelector]{cutoff.k.percent}},
or \code{"all"}. Note that \code{"all"} is a not a function but indicates
pass-thru should be performed in the case of a \code{filter} which selects
rather than ranks features. Default: \code{"cutoff.biggest.diff"}.}

\item{k}{Passed through to the \code{selector} in the case where \code{selector} is
\code{\link[FSelector]{cutoff.k}} or \code{\link[FSelector]{cutoff.k.percent}}.
Otherwise, should remain NULL (the default). For \code{\link[FSelector]{cutoff.k}},
this is an integer indicating the number of features to keep from \code{X}.
For \code{\link[FSelector]{cutoff.k.percent}}, this is instead the proportion
of features to keep.}

\item{verbose}{Should debugging messages be printed? Default: \code{FALSE}.}

\item{...}{Currently unused.}
}
\value{
A logical vector with length equal to \code{ncol(X)}.
}
\description{
The \code{\link[FSelector]{relief}} algorithm implements the RReliefF
(Robnik-Sikonja & Kononenko, 1997) feature quality estimation algorithm,
an extension to ReliefF (Kononenko, 1994) and Relief (Kira & Rendell, 1992)
algorithms. RReliefF is compatible with both classification and regression
problems and is well-suited to \code{X} with strong associations between
features.
}
\examples{
data(iris)
Y <- as.numeric(iris$Species=="setosa")
X <- iris[,-which(colnames(iris)=="Species")]
screen.FSelector.relief(Y, X, binomial(), selector = "cutoff.k.percent", k = 0.75)

data(mtcars)
Y <- mtcars$mpg
X <- mtcars[,-which(colnames(mtcars)=="mpg")]
screen.FSelector.relief(Y, X, gaussian(), neighbours.count = 3, sample.size = 15)

# based on examples in SuperLearner package
set.seed(1)
n <- 100
p <- 20
X <- matrix(rnorm(n*p), nrow = n, ncol = p)
X <- data.frame(X)
Y <- X[, 1] + sqrt(abs(X[, 2] * X[, 3])) + X[, 2] - X[, 3] + rnorm(n)

library(SuperLearner)
sl = SuperLearner(Y, X, family = gaussian(), cvControl = list(V = 2),
                  SL.library = list(c("SL.glm", "All"),
                                    c("SL.glm", "screen.FSelector.relief")))
sl
sl$whichScreen
}
\references{
\url{https://www.aaai.org/Library/AAAI/1992/aaai92-020.php},
\url{https://doi.org/10.1007/3-540-57868-4_57},
\url{http://dl.acm.org/citation.cfm?id=645526.657141}
}
