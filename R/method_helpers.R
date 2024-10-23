internalpredict.glmermboost <- function(object,
                                newdata = NULL,
                                RE = TRUE,
                                type = c("link", "response", "class"),
                                which = NULL,
                                aggregate = c("sum", "cumsum", "none"),
                                ...) {
  
  type <- match.arg(type)
  aggregate <- match.arg(aggregate)
  
  
  # Determine the type of prediction based on conditions
  if (aggregate == "none") {
    message("No intercept and no random effects included.")
    type_pred <- type
  } else if (!is.null(which)) {
    type_pred <- type
  } else {
    type_pred <- "link"  # force to link when mermboost gets active
  }
  
  # Single call to predict.glmboost
  merm_preds <- predict.glmboost(object = object,
                                 newdata = newdata,
                                 type = type_pred,
                                 which = which,
                                 aggregate = aggregate,
                                 ...)
  
  # maybe stop and return mboost only
  if (aggregate == "none" || !is.null(which)) {
    return(merm_preds)
  }
  
  if (!is.null(newdata)) {
    if (RE) {
      if (!object$id %in% names(newdata)) {stop("Missing correct cluster-id.")}
      re_cond <- any(!(unlist(unique(newdata[object$id])) %in% object$data[object$id]))
      if (re_cond) {
        message("Not all new cluster-ids can be found in original data.
                Therefore, the argument for random effects changed to 'RE = F'.")
        RE <- FALSE
      }
      # we need a construction of a new Z somehow
      slope_covariates <- setdiff(all.vars(object$ran_formula[[3]]), object$id)
      Z0 <- newdata[slope_covariates]
      # did we use a random intercept?
      int_cond <- any(grepl("(Intercept)", names(object$nuisance()[[mstop(object)]]$theta)))
      if (int_cond) {Z0 <- cbind(1, Z0)}
      
      
      
      id_num <- newdata[object$id][,1]
      id_num_orig <- object$data[object$id][,1]
      un_id <- unique(id_num)
      un_id_orig <- unique(id_num_orig)
      RE_reorder_ind <- match(un_id, un_id_orig)
      n <- length(un_id)
      
      Z <- list()
      for (i in 1:n) {
        Z[[i]] <- as.matrix(Z0[id_num == un_id[i],])
      }
      Z <- Matrix::bdiag(Z)
      
      relevant_Z <- Z
    }
  } else if (is.null(newdata)) {
    # keep everything as it is, when no new data
    relevant_Z <- object$Z
    RE_reorder_ind <- 1:nrow(ranef(object))
  }
  
  # sum is just the last iteration
  if(aggregate == "sum") {
    merm_preds <- merm_preds + object$nuisance()[[mstop(object)]]$ff
    if (RE) {
      # create a gamma matrix, qn x 1
      GAMMA <- as.vector(t(ranef(object)[RE_reorder_ind,]))
      RE_pred <- matrix(relevant_Z %*% GAMMA)
      merm_preds <- merm_preds + RE_pred
    }
  } else if (aggregate == "cumsum") {
    if (class(merm_preds)[1] == "list") {
      stop("Tell the author the cumsum is a list again")
    }
    # for glmboost it is a matrix with dimensions observations x iterations
    # Loop through each column (iteration)
    for (i in 1:ncol(merm_preds)) {
      merm_preds[, i] <- merm_preds[, i] + object$nuisance()[[i]]$ff
      if (RE) {
        # Extract the random effects vector for the current iteration
        GAMMA <- as.vector(t(ranef(object, iteration = i)[RE_reorder_ind,]))
        
        # Compute RE_pred for the current iteration
        RE_pred <- matrix(relevant_Z %*% GAMMA)
        # Add RE_pred to the corresponding column of merm_preds
        merm_preds[, i] <- merm_preds[, i] + RE_pred
      }
    }
  }
  
  
  if (type %in% c("response", "class")) {
    my_merm_preds <- merm_preds
    my_fam <- object$lme4_family()
    # for aggregate sum this works
    merm_preds <- object$lme4_family()$linkinv(merm_preds)
    
    if (type == "class") {
      if (object$lme4_family$family == "binomial") {
        merm_preds <- round(merm_preds)
      } else {
        stop("Class-type predictions are only implemented for binomial family.")
      }
    }
  }
  
  return(merm_preds)
}


internalpredict.mermboost <- function(object,
                                      newdata = NULL,
                                      RE = TRUE,
                                      type = c("link", "response", "class"),
                                      which = NULL,
                                      aggregate = c("sum", "cumsum", "none"),
                                      ...) {
  
  type <- match.arg(type)
  aggregate <- match.arg(aggregate)
  
  
  # Determine the type of prediction based on conditions
  if (aggregate == "none") {
    message("No intercept and no random effects included.")
    type_pred <- type
  } else if (!is.null(which)) {
    type_pred <- type
  } else {
    type_pred <- "link"  # force to link when mermboost gets active
  }
  
  # Single call to predict.glmboost
  merm_preds <- predict.mboost(object = object,
                               newdata = newdata,
                               type = type_pred,
                               which = which,
                               aggregate = aggregate,
                               ...)
  
  # maybe stop and return mboost only
  if (aggregate == "none" || !is.null(which)) {
    return(merm_preds)
  }
  
  
  if (!is.null(newdata)) {
    if (RE) {
      if (!object$id %in% names(newdata)) {stop("Missing correct cluster-id.")}
      re_cond <- any(!(unlist(unique(newdata[object$id])) %in% object$data[object$id]))
      if (re_cond) {
        message("Not all new cluster-ids can be found in original data.
                Therefore, the argument for random effects changed to 'RE = F'.")
        RE <- FALSE
      }
      # we need a construction of a new Z somehow
      slope_covariates <- setdiff(all.vars(object$ran_formula[[3]]), object$id)
      Z0 <- newdata[slope_covariates]
      # did we use a random intercept?
      int_cond <- any(grepl("(Intercept)", names(object$nuisance()[[mstop(object)]]$theta)))
      if (int_cond) {Z0 <- cbind(1, Z0)}
      
      id_num <- newdata[object$id][,1]
      id_num_orig <- object$data[object$id][,1]
      un_id <- unique(id_num)
      un_id_orig <- unique(id_num_orig)
      RE_reorder_ind <- match(un_id, un_id_orig)
      n <- length(un_id)
      
      Z <- list()
      for (i in 1:n) {
        Z[[i]] <- as.matrix(Z0[id_num == un_id[i],])
      }
      Z <- Matrix::bdiag(Z)
      
      relevant_Z <- Z
    }
  } else if (is.null(newdata)) {
    # keep everything as it is, when no new data
    relevant_Z <- object$Z
    RE_reorder_ind <- 1:nrow(ranef(object))
  }
  
  if(RE) {
    if(aggregate == "sum") {
      # create a gamma matrix, qn x 1
      GAMMA <- as.vector(t(ranef(object)[RE_reorder_ind,]))
      RE_pred <- matrix(relevant_Z %*% GAMMA)
      merm_preds <- merm_preds + RE_pred
    } else if (aggregate == "cumsum") {
      if (class(merm_preds)[1] == "list") {
        stop("Tell the author the cumsum is a list again")
      }
      for (i in 1:object$mstop()) {
        # Extract the random effects vector for the current iteration
        GAMMA <- as.vector(t(ranef(object, iteration = i)[RE_reorder_ind,]))
        
        # Compute RE_pred for the current iteration
        RE_pred <- matrix(relevant_Z %*% GAMMA)
        # Add RE_pred to the corresponding column of merm_preds
        merm_preds[, i] <- merm_preds[, i] + RE_pred
      }
    }
  }
  
  
  if (type %in% c("response", "class")) {
    my_merm_preds <- merm_preds
    my_fam <- object$lme4_family()
    # for aggregate sum this works
    merm_preds <- object$lme4_family()$linkinv(merm_preds)
    
    if (type == "class") {
      if (object$lme4_family$family == "binomial") {
        merm_preds <- round(merm_preds)
      } else {
        stop("Class-type predictions are only implemented for binomial family.")
      }
    }
  }
  
  return(merm_preds)
}
