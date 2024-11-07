#' @export
glmermboost <- function(formula, data = list(),
                        weights = NULL,
                        offset = NULL, family = gaussian,
                        na.action = na.omit, contrasts.arg = NULL,
                        center = TRUE, control = boost_control(),
                        oobweights = NULL, ...) {

  data <- as.data.frame(data)
  
  if (!class(family) == "boost_family") {
    # Standardize the family input; if condition since accepting NBinomial
    family <- switch(
      class(family)[1],
      "character" = get(family),
      "function" = family(),
      "family" = family,
      stop("Invalid family argument")
    )
    
    if (inherits(family, "family") && family$family == "Gamma") {
      if (family$link == "inverse") {
        family <- Gamma(link = "log")
        message("Turned into log-link since inverse-link has problems concerning risk and gradient calculations.")
      }
    }
  }
  
  randoms <- organise_random(formula = formula, data = data, model = "glmm")
  data <- randoms$data
  ran_parts <- randoms$ran_parts
  ran_parts$family <- eval(family)

  # organise fixed effects
  formula_fix <- lme4::nobars(formula) # extract fix formula
  vars <- all.vars(formula_fix[[3]])
  
  X <- create_X(data, vars)

  # build correction matrix C
  Cm <- build_C(data = data, X = X, ran_parts = ran_parts)
  ran_parts$Cm <- Cm

  # create family for mboost
  fam <- build_family(ran_parts = ran_parts, weights = weights)

  # Finally, mboost
  mb_f <- as.formula(paste(formula_fix[[2]]," ~ ",
                           paste(deparse(formula_fix[[3]]), collapse = ""),
                           " - 1"))

  model <- glmboost(mb_f, data = data, weights = weights,
                    offset = offset, family = fam,
                    na.action = na.action, contrasts.arg = contrasts.arg,
                    center = center, control = control,
                    oobweights = oobweights, ...)

  #model$X <- extract(model, what = "design")
  coef_names <- sub("[A-Z]+$", "", names(coef(model)))
  
  
  # intercept calculation due to centering
  if (center == TRUE) {
    common_names <- intersect(colnames(model$model.frame()), names(model$center))
    
    newmf <- model$model.frame()
    # add centers
    for (name in common_names) {
      if (model$center[name] > 0) {
        newmf[, name] <- model$model.frame()[, name] +  model$center[name]
      } else {
        newmf[, name] <- model$model.frame()[, name] - abs(model$center[name])  # Subtracting the absolute value
      }
    }
    
    model$intercept <- model$nuisance()[[model$mstop()]]$ff - 
      mean(predict.glmboost(model, newdata = newmf))
  }
  
  model$Z <- ran_parts$Z
  model$lme4_family <- function() ran_parts$family
  model$formula <- formula
  model$ran_formula <- ran_parts$mm_f
  model$fix_formula <- formula_fix
  model$id <- ran_parts$id
  model$data <- data
  class(model) <- c("glmermboost", "glmboost", "mboost")
  
  model$predict <- local({
    mod <- model  # Store the model object within the local environment
    
    function(newdata = NULL,
             RE = TRUE,
             type = c("link", "response", "class"),
             which = NULL,
             aggregate = c("sum", "cumsum", "none"),
             ...) {
    
      ret <- internalpredict.glmermboost(object = mod,
                                         newdata = newdata,
                                         RE = RE,
                                         type = type,
                                         which = which,
                                         aggregate = aggregate,
                                         ...)
      
      return(ret)
    }
  })
  
  model$fitted <- function() model$predict(type = "response")
  
  
  return(model)
}


#' @export
mermboost <- function(formula, data = list(),
                      na.action = na.omit, weights = NULL,
                      offset = NULL, family = gaussian, control = boost_control(),
                      oobweights = NULL, baselearner = c("bbs", "bols", "btree", "bss", "bns"),
                      ...) {
  
  data <- as.data.frame(data)
  
  if (!class(family) == "boost_family") {
    # Standardize the family input; if condition since accepting NBinomial
    family <- switch(
      class(family)[1],
      "character" = get(family),
      "function" = family(),
      "family" = family,
      stop("Invalid family argument")
    )
    
    if (inherits(family, "family") && family$family == "Gamma") {
      if (family$link == "inverse") {
        family <- Gamma(link = "log")
        message("Turned into log-link since inverse-link has problems concerning risk and gradient calculations.")
      }
    }
  }

  randoms <- organise_random(formula = formula, data = data, model = "gamm")
  data <- randoms$data
  ran_parts <- randoms$ran_parts
  ran_parts$family <- eval(family)

  # organise fixed effects
  formula_fix <- lme4::nobars(formula) # extract fix formula
  vars <- all.vars(formula_fix[[3]])
  X <- create_X(data, vars)
  
  # build correction matrix C
  Cm <- build_C(data = data, X = X, ran_parts = ran_parts)
  ran_parts$Cm <- Cm

  # create family for mboost
  fam <- build_family(ran_parts = ran_parts, weights = weights)

  # Finally, mboost
  mb_f <- as.formula(paste(formula_fix[[2]]," ~ ",
                           paste(deparse(formula_fix[[3]]), collapse = "")))
  model <- mboost(mb_f, data = data, family = fam,
                  na.action = na.action, weights = weights,
                  offset = offset,
                  control = control,
                  oobweights = oobweights, baselearner = baselearner,
                  ...)

  model$Z <- ran_parts$Z
  model$lme4_family <- function() family
  model$formula <- formula
  model$ran_formula <- ran_parts$mm_f
  model$fix_formula <- formula_fix
  model$id <- ran_parts$id
  model$int_cond <- ran_parts$int_cond
  model$data <- data

  class(model) <- c("mermboost", "mboost")
  
  model$predict <- local({
    mod <- model  # Store the model object within the local environment
    
    function(newdata = NULL,
             RE = TRUE,
             type = c("link", "response", "class"),
             which = NULL,
             aggregate = c("sum", "cumsum", "none"),
             ...) {
      
      ret <- internalpredict.mermboost(object = mod,
                                         newdata = newdata,
                                         RE = RE,
                                         type = type,
                                         which = which,
                                         aggregate = aggregate,
                                         ...)
      
      return(ret)
    }
  })
  
  model$fitted <- function() model$predict(type = "response")

  return(model)
}



#' @export
mer_cvrisk <- function(object, folds, no_of_folds, cores = 1) {
  # check if both is given and stop
  if (!missing(folds)) {
    if (!missing(no_of_folds)) {
      stop("'folds' and 'no_of_folds' both specified")
    }
  }

  N <- nrow(object$data)

  if (!missing(no_of_folds)) {
    k <- no_of_folds
    fold_share <- 1 / k

    un_id <- unique(object$data[object$id])
    no_of_ind <- nrow(un_id)
    id_ind <- list()
    for (i in 1:no_of_ind) {
      id_ind[[i]] <- which(object$data[object$id] == i)
    }

    folds <- matrix(1, nrow = N, ncol = k)
    test_ind <- list()
    test_ind_size <- floor(no_of_ind * fold_share)
    if (test_ind_size == 0) {stop("You might want too many folds.")}
    for (i in 1:k) {
      test_ind <- sample(1:no_of_ind, test_ind_size)
      zero_ind <- unlist(id_ind[test_ind])
      folds[zero_ind,i] = 0
    }
  }

  k <- ncol(folds)
  # folds are given now
  iter <- object$mstop()
  lme4fam <- object$lme4_family
  ran_y <- rep(c(0,1), times = 5)
  ran_ids <- rep(1:2, each = 5)


  suppressMessages(suppressWarnings(
    if (identical(gaussian, lme4fam)) {
      fam <- family(ran_mod <- lme4::lmer(ran_y ~ 1 - 1 + (1 | ran_ids)))
    } else {
      # for negative binomial
      fam <- tryCatch({family(ran_mod <- lme4::glmer(ran_y ~ 1- 1 + (1 | ran_ids),
                                                     family = lme4fam))},
                      error = function(i) family(
                        ran_mod <- lme4::glmer.nb(ran_y ~ 1 - 1 + (1 | ran_ids))))
    }
  ))



  folds_l <- list()
  for (i in 1:ncol(folds)) {
    folds_l[[i]] <- folds[,i]
  }

  if (class(object)[1] == "glmermboost") {
    Y <- model.response(object$model.frame())
    if (is.factor(Y) && 0 %in% Y) {
      Y <- as.numeric(Y) - 1
    }

    cross_val <- function(giv_weights) {
      train_ind <- ifelse(giv_weights == 1, T, F)
      test_ind <- ifelse(giv_weights == 0, T, F)

      train <- object$data[train_ind,]
      test <- object$data[test_ind,]

      #########
      train_mod <- glmermboost(formula = object$formula,
                               data = train,
                               offset = object$offset,
                               family = object$lme4_family,
                               control = object$control)
      eta <- predict(train_mod,
                    newdata = test,
                    type = "link",
                    aggregate = "cumsum",
                    RE = F)
      aic <- c()

      for (j in 1:iter) {
        mu <- fam$linkinv(eta[,j])
        dev <- sum(fam$dev.resids(y = Y[test_ind], mu = mu, wt = 1))
        aic[j] <- fam$aic(y = Y[test_ind], n = 1, mu = mu, wt = 1, dev = dev)
      }
      return(aic / 2)
    }



    aics <- parallel::mclapply(folds_l, cross_val, mc.cores = cores)
    RISKMAT <- matrix(NA, nrow = iter, ncol = k)
    for (i in 1:length(aics)) {
      RISKMAT[,i] <- aics[[i]]
    }

  } else if (class(object)[1] == "mermboost") {
    Y <- as.matrix(object$data[all.vars(object$formula[[2]])])

    cross_val <- function(giv_weights) {
      train_ind <- ifelse(giv_weights == 1, T, F)
      test_ind <- ifelse(giv_weights == 0, T, F)

      train <- object$data[train_ind,]
      test <- object$data[test_ind,]

      train_mod <- mermboost(object$formula, 
                             data = train,
                             offset = object$offset,
                             family = object$lme4_family,
                             control = object$control)

      eta <- predict.mermboost(train_mod, 
                               newdata = test,
                               type = "link", 
                               aggregate = "cumsum",
                               RE = F)

      aic <- c()
      for (j in 1:iter) {
        mu <- fam$linkinv(eta[,j])
        dev <- sum(fam$dev.resids(y = Y[test_ind], mu = mu, wt = 1))
        aic[j] <- fam$aic(y = Y[test_ind], n = 1, mu = mu, wt = 1, dev = dev)
      }

      return(aic / 2)
    }

    aics <- parallel::mclapply(folds_l, cross_val, mc.cores = cores)
    RISKMAT <- matrix(NA, nrow = iter, ncol = k)
    for (i in 1:length(aics)) {
      RISKMAT[,i] <- aics[[i]]
    }

  }

  test_sizes <- colSums(1 - folds)

  RISKMAT <- sweep(RISKMAT, 2, test_sizes, "/")

  avg_risks <- rowMeans(RISKMAT)
  m_opt <- which.min(avg_risks)

  res <- list(cv_risks = RISKMAT,
              avg_risks = avg_risks,
              m_opt = m_opt,
              folds = folds)

  class(res) <- "mer_cv"

  return(res)
}


#' @export
find_ccc <- function(df, id_char) {
  splitted_list <- split(df, df[id_char])
  ccc_logical <- list()
  for (i in 1:length(splitted_list)) {
    ccc_logical[[i]] <- apply(splitted_list[[i]], 2,
                              function(x) length(unique(x)) == 1)
  }
  ccc_logical <- do.call(rbind, ccc_logical)

  ccc_row <- apply(ccc_logical, 2, function(x) all(x))

  id_ind <- which(names(ccc_row) == id_char)

  ccc_row <- ccc_row[-id_ind]
  return(ccc_row)
}


