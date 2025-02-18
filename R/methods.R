#' @export
predict.mermboost <- function(object,
                              newdata = NULL,
                              RE = TRUE,
                              type = c("link", "response", "class"),
                              which = NULL,
                              aggregate = c("sum", "cumsum", "none"),
                              ...) {
  
  ret <- object$predict(newdata = newdata,
                        RE = RE,
                        type = type,
                        which = which,
                        aggregate = aggregate,
                        ...)
  
  return(ret)
}




#' @export
predict.glmermboost <- function(object,
                                newdata = NULL,
                                RE = TRUE,
                                type = c("link", "response", "class"),
                                which = NULL,
                                aggregate = c("sum", "cumsum", "none"),
                                ...) {
  
  ret <- object$predict(newdata = newdata,
                        RE = RE,
                        type = type,
                        which = which,
                        aggregate = aggregate,
                        ...)
  
  return(ret)
}


#' @export
mstop.mer_cv <- function(object, ...) {
  return(object$m_opt)
}

#' @export
plot.mer_cv <- function(x, ...) {

  avg_risks <- x$avg_risks
  RISKMAT <- x$cv_risks
  iter <- length(avg_risks)

  plot(avg_risks, type = "l",
       ylim = c(min(RISKMAT), max(RISKMAT)),
       xlab = "Number of Iterations",
       ylab = "out-of-bag Risk", 
       ...)

  for (i in 1:ncol(RISKMAT)) {
    lines(1:iter, y = RISKMAT[,i], col = "gray")
  }
  lines(x = 1:length(avg_risks), y = avg_risks)
  abline(v = x$m_opt, lty = "dashed", col = "gray")
}


#' @export
ranef.glmermboost <- function(object, iteration = mstop(object), ...) {
  ranef.mermboost(object, iteration)
}
#' @export
ranef.mermboost <- function(object, iteration = mstop(object), ...) {
  object$nuisance()[[iteration]]$rf
}
#' @export
VarCorr.glmermboost <- function(x, sigma=1, iteration = mstop(x), ...) {
  VarCorr.mermboost(x, sigma, iteration)
}
#' @export
VarCorr.mermboost <- function(x, sigma=1, iteration = mstop(x), ...) {
  if(sigma!=1) {stop("Did not implement multiplier for the standard deviations.")}
  x$nuisance()[[iteration]]$Q
}
