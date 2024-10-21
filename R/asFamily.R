as.Family <- function(object, ...)
  UseMethod("as.Family")


.copy_variables <- function(call, envir, ran_parts) {
  nm <- names(call)
  nm <- nm[nm != ""]
  vars <- do.call("c", sapply(call[nm], all.vars))
  ret <- new.env()
  for(n in vars) {
    #if (n %in% names(envir))
    if (n %in% ls(envir))
      assign(n, get(n, envir), ret)
  }

  # Ensure ran_parts$family is correctly assigned
  if ("family" %in% names(ran_parts)) {
    assign("family", ran_parts$family, ret)
  }

  ret
}

#' @export
as.Family.merMod <- function(object, ran_parts) {
  q <- ran_parts$q
  id <- ran_parts$id
  Z0 <- ran_parts$Z0
  C <- ran_parts$Cm
  Z <- ran_parts$Z
  
  # object is an intercept object
  model <- object
  
  rf_orig <- lme4::ranef(model)[[id]]
  rf <- matrix(C %*% as.vector(t(rf_orig)), ncol = q, byrow = T)
  
  ff <- lme4::fixef(model)
  theta <- lme4::getME(model, "theta")
  Q <- lme4::VarCorr(model)
  name_of_fe <- "fixef"
  cl <- object@call
  ffm <- formula(object)
  env <- .copy_variables(cl, environment(ffm), ran_parts)
  environment(model@call$formula) <- env
  environment(attr(model@frame, "formula")) <- env

  N <- nrow(mf <- model.frame(object))
  Y <- model.response(mf)
  ffm <- as.formula(paste(deparse(ffm[[2]]), "~ 1 - 1"))
  cl$formula <- ffm

  fun <- gsub("^lme4::", "", deparse(cl[[1L]]))
  stopifnot(fun %in% c("lmer", "glmer"))
  fun <- gsub("er$", "", fun)
  
  cl[[1L]] <- as.name(fun)
  environment(cl$formula) <- env
  m0 <- eval(cl, envir = env)
  offset_. <- weights_. <- start_. <- lp_. <- NA
  fam <- family(object)


  risk <- function(y, f, w = 1) {
    ### newdata removes offset from predict
    if (any(w == 0)) {
      lp <- predict(model, newdata = mf, re.form = NA) + c(f) # ohne RE
    } else {
      lp <- predict(model, newdata = mf, re.form = NA) + rf_lp + c(f)
    }

    if (inherits(object, "lmerMod"))
      return(sum(w * (Y - lp)^2))

    mu <- fam$linkinv(lp)
    dev <- sum(fam$dev.resids(m0$y, mu, w))
    aic <- fam$aic(y = m0$y, n = 1, mu = mu, wt = w, dev = dev)
    loss_value <- aic / 2

    ### see as.Family.glm@risk
    return(loss_value)

  }

  ngradient <- function(y, f, w = 1) {
    if (length(f) == 1) f <- rep(f, N)
    if (length(w) == 1) w <- rep(w, N)

    ### this is not sufficient, data used to fit "object" needs to be
    ### the same as mboost(..., data); na.action can and will mess this
    ### up!!!
    stopifnot(NROW(f) == N)

    assign("offset_.", f, env)
    assign("weights_.", w, env)
    assign("start_.", lme4::getME(model, c("theta", "fixef")), env)
    if (is.na(start_.)) {start_. <- NULL}

    model <<- tryCatch(eval(update(model, verbose = 0, offset = offset_.,
                                   weights = weights_., start = start_.), env),
                       error = function(i) eval(model, env),
                       message = function(i) eval(model, env))

    rf_orig <- lme4::ranef(model)[[id]]
    ff <<- lme4::fixef(model) #- mean(f)

    theta <<- lme4::getME(model, "theta")
    Q <<- lme4::VarCorr(model)

    rf <<- matrix(C %*% as.vector(t(rf_orig)), ncol = q, byrow = T)

    rf_lp <<- as.matrix(Z %*% as.vector(t(rf)))
    #lp <- ff + rf_lp
    lp <- predict(model, re.form = NA) + rf_lp
    assign("lp_.", lp, env)

    # Update the nullmodel
    tmp <- eval(update(m0, offset = lp_., weights = weights_., evaluate = FALSE), env)

    if (logLik(tmp) < logLik(m0))
      warning("risk increase; decrease stepsize nu")
    m0 <<- tmp
    ### note: inherits(m0, "lm") will also be true for glms
    if (class(m0)[1] == "lm") {
      return(m0$residuals) # here comes the deviance
    }

    # when class = "lm" then nothing happens here
    ### for glms fake an intercept = 0 model
    class(tmp) <- c("merFamily", class(tmp))
    tmp$N <- N
    tmp$coefficients <- 0

    ### residual wrt a constant for _all_ observations
    #sandwich::estfun(tmp)
    residuals(tmp, type = "response")
  }

  Family(ngradient = ngradient,
         risk = risk,
         check_y = function(y) rep(TRUE, N),
         offset = function(y, w) 0,
         nuisance = function() return(list(ff = ff, rf = rf, theta = theta, Q = Q)),
         name = "glmer",
         response = function(f) fam$linkinv(f)
         )
}


#
# as.Family.lm <- function(object, ...) {
#
#   model <- object
#   cf <- coef(object)
#   X <- model.matrix(object)
#   Y <- model.response(model.frame(object))
#   N <- nrow(model.frame(object))
#
#   risk <- function(y, f, w = 1)
#     sum(w * (Y - (c(f) + X %*% cf))^2)
#
#   ngradient <- function(y, f, w = 1) {
#     if (length(f) == 1) f <- rep(f, N)
#     if (length(w) == 1) w <- rep(w, N)
#
#     ### this is not sufficient, data used to fit "object" needs to be
#     ### the same as mboost(..., data); na.action can and will mess this
#     ### up!!!
#     stopifnot(NROW(f) == N)
#
#     model <<- lm.wfit(x = X, y = Y, w = w, offset = c(f))
#     cf <<- coef(model)
#
#     ### residual wrt a constant for _all_ observations
#     model$residuals
#   }
#
#   Family(ngradient = ngradient, risk = risk,
#          check_y = function(y) rep(TRUE, N),
#          offset = function(y, w) 0,
#          nuisance = function() return(list(coefficients = cf)),
#          name = "lm",
#          response = function(f) f)
# }
#
# as.Family.glm <- function(object, ...) {
#
#   model <- object
#   cf <- coef(object)
#   X <- model.matrix(object)
#   stopifnot("(Intercept)" %in% colnames(X))
#   Y <- model.response(model.frame(object))
#   N <- nrow(model.frame(object))
#
#   fam <- family(object)
#   p <- object$rank
#   if (fam$family %in% c("gaussian", "Gamma", "inverse.gaussian"))
#     p <- p + 1
#
#   risk <- function(y, f, w = 1) {
#     mu <- fam$linkinv(c(f) + X %*% cf)
#     dev <- sum(fam$dev.resids(object$y, mu, w))
#     aic <- fam$aic(y = object$y, n = 1, mu = mu, wt = w)
#     ### see stats:::glm.wfit and stats:::logLik.glm
#     return(aic / 2)
#   }
#
#   ngradient <- function(y, f, w = 1) {
#     if (length(f) == 1) f <- rep(f, N)
#     if (length(w) == 1) w <- rep(w, N)
#
#     ### this is not sufficient, data used to fit "object" needs to be
#     ### the same as mboost(..., data); na.action can and will mess this
#     ### up!!!
#     stopifnot(NROW(f) == N)
#
#     tmp <- glm.fit(x = X, y = Y, weights = w, offset = c(f),
#                    family = fam)
#     #### see stats:::glm.wfit and stats:::logLik.glm
#     if ((p - tmp$aic / 2) < logLik(model))
#       warning("risk increase; decrease stepsize nu")
#
#     model[names(tmp)] <<- tmp
#     cf <<- coef(model)
#
#     ### residual wrt a constant for _all_ observations
#     sandwich::estfun(model)[,"(Intercept)"]
#   }
#
#   Family(ngradient = ngradient, risk = risk,
#          check_y = function(y) rep(TRUE, N),
#          offset = function(y, w) 0,
#          nuisance = function() return(list(coefficients = cf)),
#          name = "glm",
#          response = function(f) fam$linkinv(f))
# }



#' @export
model.matrix.merFamily <- function(object, ...)
  matrix(1, nrow = object$N, ncol = 1)






#
# as.Family.coxph <- function(object, ...) {
#
#   model <- object
#   cf <- coef(model)
#   ffm <- formula(object)
#   cl <- object$call
#   env <- .copy_variables(cl, environment(ffm))
#   environment(model$formula) <- env
#   N <- nrow(model.frame(object))
#
#   fm <- . ~ . + offset(offset_.)
#   m0 <- update(model, formula = fm,
#                subset = subset_., weights = weights_., evaluate = FALSE)
#   environment(m0$formula) <- env
#
#   offset_. <- weights_. <- subset_. <- cf_. <- NA
#
#   risk <- function(y, f, w = 1) {
#
#     assign("offset_.", f, env)
#     assign("weights_.", w, env)
#     assign("subset_.", which(w > 0), env)
#     assign("cf_.", cf, env)
#
#     ### do not re-estimate parameters in cf but rather evaluate
#     ### partial likelihood with cf and f
#     tmp <- eval(update(model, formula = fm,
#                        subset = subset_.,
#                        weights = weights_.,
#                        init = cf_.,
#                        control = coxph.control(iter.max = 0),
#                        evaluate = FALSE), env)
#     -logLik(tmp)
#   }
#
#   ngradient <- function(y, f, w = 1) {
#     if (length(f) == 1) f <- rep(f, N)
#     if (length(w) == 1) w <- rep(w, N)
#
#     ### this is not sufficient, data used to fit "object" needs to be
#     ### the same as mboost(..., data); na.action can and will mess this
#     ### up!!!
#     stopifnot(NROW(f) == N)
#
#     assign("offset_.", f, env)
#     assign("weights_.", w, env)
#     assign("subset_.", which(w > 0), env)
#
#     tmp <- eval(m0, env)
#     if (logLik(tmp) < logLik(model))
#       warning("risk increase; decrease stepsize nu")
#     model <<- tmp
#     cf <<- coef(model)
#
#     ### residual wrt a constant for _all_ observations
#     ret <- numeric(N)
#     ret[w > 0] <- residuals(model, type = "martingale")
#     ret
#   }
#
#   Family(ngradient = ngradient, risk = risk,
#          check_y = function(y) rep(TRUE, N),
#          offset = function(y, w) 0,
#          nuisance = function() return(list(coefficients = cf)),
#          name = "coxph",
#          ### conditional survivor function???
#          response = function(f) NA)
# }





