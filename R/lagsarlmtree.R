utils::globalVariables(c(".tree", ".rhowy"))

lagsarlmtree <- function(formula, data, listw = NULL, method = "eigen",
  zero.policy = NULL, interval = NULL, control = list(),
  rhowystart = NULL, abstol = 0.001, maxit = 100, 
  dfsplit = TRUE, verbose = FALSE, plot = FALSE, ...)
{
  ## remember call
  cl <- match.call()

  ## process listw
  if(is.null(listw)) stop("a 'listw' object has to be provided")
  
  ## formula processing (full, tree, regression)
  ff <- Formula::as.Formula(formula)
  tf <- formula(ff, lhs = 1L, rhs = c(1L, 2L))
  rf <- formula(ff, lhs = 1L, rhs = 1L)
  ## intercept in lmtree?
  intercept <- attr(terms(rf), "intercept") > 0L
  if(!intercept) stop("linear models without intercept currently not supported")
  ## without tree
  rf0 <- rf
  ## with tree
  rf <- update(rf, . ~ 0 + .tree / .)

  ## initialization
  iteration <- 0L
  if (is.null(rhowystart)) {
    rhowystart <- lagsarlm(rf0, data = data, listw = listw, method = method,
      zero.policy = zero.policy, interval = interval, control = control)
    rhowystart <- rhowystart$rho * lag(listw, rhowystart$y)
  }
  data$.rhowy <- rhowystart  

  continue <- TRUE
  oldloglik <- -Inf

  ## iterate between lagsarlm and lmtree estimation
  while (continue) {
    iteration <- iteration + 1L

    ## lmtree
    tree <- lmtree(tf, data = data, offset = .rhowy, dfsplit = FALSE, ...)
    if(plot) plot(tree)
    data$.tree <- factor(predict(tree, newdata = data, type = "node"))

    ## estimate full lm model but force all coefficients from the
    ## .tree (and the overall intercept) to zero for the prediction
    lagsarlm <- if(length(levels(data$.tree)) > 1L) {
      lagsarlm(rf, data = data, listw = listw, method = method,
      zero.policy = zero.policy, interval = interval, control = control)
    } else {
      lagsarlm(rf0, data = data, listw = listw, method = method,
      zero.policy = zero.policy, interval = interval, control = control)
    }
    data$.rhowy <- lagsarlm$rho * lag(listw, lagsarlm$y)

    ## iteration information
    newloglik <- logLik(lagsarlm)
    continue <- (newloglik - oldloglik > abstol) & (iteration < maxit) 
    oldloglik <- newloglik
    if(verbose) print(newloglik)
  }
  
  ## collect results
  result <- list(
    formula = formula,
    call = cl,
    tree = tree,
    lagsarlm = lagsarlm,
    data = data,
    nobs = nrow(data),
    loglik = as.numeric(newloglik),
    df = attr(newloglik, "df"),
    dfsplit = dfsplit,
    iterations = iteration, 
    maxit = maxit,
    rhowystart = rhowystart, 
    abstol = abstol,
    listw = listw,
    mob.control = list(...)
  )
  class(result) <- "lagsarlmtree"
  return(result)
}

coef.lagsarlmtree <- function(object, model = c("all", "tree", "rho"), ...)
{
  model <- match.arg(model, c("all", "tree", "rho"))
  if(model == "all")  return(coef(object$lagsarlm, ...))
  if(model == "tree") return(coef(object$tree, ...))
  return(object$lagsarlm$rho)
}

plot.lagsarlmtree <- function(x, ...) {
  plot(x$tree, ...)
}

sctest.lagsarlmtree <- function(x, ...) {
  sctest.modelparty(x$tree, ...)
}

logLik.lagsarlmtree <- function(object, dfsplit = NULL, ...)
{
  if(is.null(dfsplit)) dfsplit <- object$dfsplit
  dfsplit <- as.integer(dfsplit) * (length(object$tree) - length(nodeids(object$tree, terminal = TRUE)))
  structure(object$loglik, df = object$df + dfsplit, nobs = object$nobs, class = "logLik")
}

print.lagsarlmtree <- function(x, title = "Spatial lag model tree", ...)
{
  print(x$tree, title = title, ...)
  cat("\nRho (from lagsarlm model):\n")
  print(coef(x, model = "rho"))
  invisible(x)
}

predict.lagsarlmtree <- function(object, newdata = NULL, type = "response", ...) { 
  if(is.null(newdata)) {
    newdata <- object$data
  }
  if(type == "node") {
    predict(object$tree, newdata = newdata, type = "node")
  } else {
    newdata$.tree <- predict(object$tree, newdata = newdata, type = "node")
    newdata$.tree <- factor(newdata$.tree,
    			    labels = levels(object$data$.tree))
    predict(object$lagsarlm, newdata = newdata, type = type, ...)
  }
}

impacts.lagsarlmtree <- function(obj, ...) {
  impacts(obj$lagsarlm, ...)
}


