\name{lagsarlmtree}
\alias{lagsarlmtree}
\alias{print.lagsarlmtree}
\alias{plot.lagsarlmtree}
\alias{coef.lagsarlmtree}
\alias{logLik.lagsarlmtree}
\alias{predict.lagsarlmtree}
\alias{sctest.lagsarlmtree}
\alias{impacts.lagsarlmtree}

\title{Spatial Lag Model Trees}

\description{
  Model-based recursive partitioning based on linear regression
  adjusting for a (global) spatial simultaneous autoregressive lag.
}

\usage{
lagsarlmtree(formula, data, listw = NULL, method = "eigen",
  zero.policy = NULL, interval = NULL, control = list(),
  rhowystart = NULL, abstol = 0.001, maxit = 100, 
  dfsplit = TRUE, verbose = FALSE, plot = FALSE, ...)
}

\arguments{
  \item{formula}{formula specifying the response variable and regressors
    and partitioning variables, respectively. For details see below.}
  \item{data}{data.frame to be used for estimating the model tree.}
  \item{listw}{a weights object for the spatial lag part of the model.}
  \item{method}{"eigen" (default) - the Jacobian is computed as the product 
of (1 - rho*eigenvalue) using \code{eigenw}, and "spam" or "Matrix_J" for strictly symmetric weights lists of styles "B" and "C", or made symmetric by similarity (Ord, 1975, Appendix C) if possible for styles "W" and "S", using code from the spam or Matrix packages to calculate the determinant; \dQuote{Matrix} and \dQuote{spam_update} provide updating Cholesky decomposition methods; "LU" provides an alternative sparse matrix decomposition approach. In addition, there are "Chebyshev" and Monte Carlo "MC" approximate log-determinant methods; the Smirnov/Anselin (2009) trace approximation is available as "moments". Three methods: "SE_classic", "SE_whichMin", and "SE_interp" are provided experimentally, the first to attempt to emulate the behaviour of Spatial Econometrics toolbox ML fitting functions. All use grids of log determinant values, and the latter two attempt to ameliorate some features of "SE_classic".
}
   \item{zero.policy}{default NULL, use global option value; if TRUE assign zero to the lagged value of zones without 
neighbours, if FALSE (default) assign NA - causing \code{lagsarlm()} to terminate with an error}
   \item{interval}{default is NULL, search interval for autoregressive parameter}
   \item{control}{list of extra control arguments - see \code{\link[spdep]{lagsarlm}}}
   \item{rhowystart}{numeric. A vector of length \code{nrow(data)}, to be used as
    an offset in estimation of the first tree. \code{NULL} by default,
    which results in an initialization with the root model (without partitioning).}
  \item{abstol}{numeric. The convergence criterion used for estimation of the model.
    When the difference in log-likelihoods of the model from two consecutive
    iterations is smaller than \code{abstol}, estimation of the model tree has
    converged.} 
  \item{maxit}{numeric. The maximum number of iterations to be performed in
    estimation of the model tree.}
  \item{dfsplit}{logical or numeric. \code{as.integer(dfsplit)} is the
    degrees of freedom per selected split employed when extracting
    the log-likelihood.}
  \item{verbose}{Should the log-likelihood value of the estimated
    model be printed for every iteration of the estimation?}
  \item{plot}{Should the tree be plotted at every iteration of the estimation?
    Note that selecting this option slows down execution of the function.}
  \item{\dots}{Additional arguments to be passed to \code{lmtree()}.
    See \code{\link[partykit]{mob_control}} documentation for details.}
}
  
\details{
  Spatial lag trees learn a tree where each terminal node is associated with
  different regression coefficients while adjusting for a (global) spatial
  simultaneous autoregressive lag. This allows for detection of subgroup-specific
  coefficients with respect to selected covariates, while adjusting for spatial
  correlations in the data. The estimation algorithm iterates between (1) estimation
  of the tree given an offset of the spatial lag effect, and (2) estimation of the
  spatial lag model given the tree structure.

  The code is still under development and might change in future versions.
}

\value{
The function returns a list with the following objects:
  \item{formula}{The formula as specified with the \code{formula} argument.}
  \item{call}{the matched call.}
  \item{tree}{The final \code{lmtree}.}
  \item{lagsarlm}{The final \code{lagsarlm} model.}
  \item{data}{The dataset specified with the \code{data} argument
    including added auxiliary variables \code{.rhowy} and \code{.tree}
    from the last iteration.}
  \item{nobs}{Number of observations.}
  \item{loglik}{The log-likelihood value of the last iteration.}
  \item{df}{Degrees of freedom.}
  \item{dfsplit}{degrees of freedom per selected split as specified with the \code{dfsplit}
    argument.}
  \item{iterations}{The number of iterations used to estimate the \code{lagsarlmtree}.} 
  \item{maxit}{The maximum number of iterations specified with the \code{maxit} argument.}
  \item{rhowystart}{Offset in estimation of the first tree as specified in the \code{rhowystart} 
    argument.}
  \item{abstol}{The prespecified value for the change in log-likelihood to evaluate
    convergence, as specified with the \code{abstol} argument.} 
  \item{listw}{The \code{listw} object used.}
  \item{mob.control}{A list containing control parameters passed to
    \code{lmtree()}, as specified with \dots.}
}

\seealso{\code{\link[stats]{lm}}, \code{\link[spdep]{lagsarlm}}, \code{\link[partykit]{lmtree}}}

\references{ 
  Wagner M, Zeileis A (2019).
    Heterogeneity and Spatial Dependence of Regional Growth in the EU: A Recursive Partitioning Approach.
    German Economic Review, \bold{20}(1), 67--82.
    \doi{10.1111/geer.12146}
    \url{https://eeecon.uibk.ac.at/~zeileis/papers/Wagner+Zeileis-2019.pdf}
}

\examples{
## data and spatial weights
data("GrowthNUTS2", package = "lagsarlmtree")
data("WeightsNUTS2", package = "lagsarlmtree")

## spatial lag model tree
system.time(tr <- lagsarlmtree(ggdpcap ~ gdpcap0 + shgfcf + shsh + shsm |
    gdpcap0 + accessrail + accessroad + capital + regboarder + regcoast + regobj1 + cee + piigs,
  data = GrowthNUTS2, listw = WeightsNUTS2$invw,
  minsize = 12, alpha = 0.05))
print(tr)
plot(tr, tp_args = list(which = 1))

## query coefficients
coef(tr, model = "tree")
coef(tr, model = "rho")
coef(tr, model = "all")
system.time({
ev <- eigenw(WeightsNUTS2$invw)
tr1 <- lagsarlmtree(ggdpcap ~ gdpcap0 + shgfcf + shsh + shsm |
    gdpcap0 + accessrail + accessroad + capital + regboarder + regcoast + regobj1 + cee + piigs,
  data = GrowthNUTS2, listw = WeightsNUTS2$invw, method = "eigen",
  control = list(pre_eig = ev), minsize = 12, alpha = 0.05)
})
coef(tr1, model = "rho")
}

\keyword{tree}
