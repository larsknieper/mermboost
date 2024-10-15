check_data <- function(data, id, formula) {
  # Iterate through the columns and convert logical columns
  data[] <- lapply(data, function(x) {
    if (is.logical(x)) {
      return(as.integer(x))  # Convert TRUE to 1 and FALSE to 0
      message("Turned a logical into numeric.")
    } else {
      return(x)  # Leave non-logical columns unchanged
    }
  })

  rel_vars <- all.vars(formula[[3]])[which(all.vars(formula[[3]]) != id)]
  data_rel <- data[rel_vars]
  # Check if all columns, except the specified one, are numeric
  numeric_test <- all(sapply(data_rel[, !names(data_rel) %in% id], is.numeric))
  if (!numeric_test) {
    stop("Please make sure your covariates are numeric.")
  }

  return(data)
}


organise_random <- function(formula, data, model) {

  formula_random <- lme4::findbars(formula) # extract random formula

  # random parts list
  ran_parts <- list()
  ran_parts$id <- all.vars(formula_random[[1]][[3]]) # id variable

  data[ran_parts$id] <- as.factor(unlist(data[ran_parts$id]))
  if(class(unlist(data[ran_parts$id])) != "factor") {
    stop("Please turn the id into a factor.")
    }
  ran_parts$rs <- all.vars(formula_random[[1]][[2]]) # Random Slopes
  ran_parts$q <- length(ran_parts$rs) # no. of random effects
  ran_parts$int_cond <- grepl("1+", deparse(formula_random[[1]])) # is there a random intercept
  if (ran_parts$int_cond) {
    ran_parts$q <- ran_parts$q + 1
  }

  # data organisation by id
  data[ran_parts$id] <- as.numeric(unlist(data[ran_parts$id])) # add id to data frame
  if (all(data[order(data[,ran_parts$id]), ] != data)) {
    message("Data has been sorted by the id variable.")
    data <- data[order(data[,ran_parts$id]), ]
  }

  if (ran_parts$int_cond) {
    ran_part <- paste(c("1 ", ran_parts$rs), collapse = " + ")
  } else {
    ran_part <- paste(rs, collapse = " + ")
  }

  if (model == "glmm") {
    mm_f <- as.formula(paste(deparse(formula[[2]]),
                             "~ 1 + (", paste(ran_part),
                             " |", ran_parts$id, ")"))
  } else if (model == "gamm") {
    mm_f <- as.formula(paste(deparse(formula[[2]]),
                             "~ 1 - 1 + (", paste(ran_part),
                             " |", ran_parts$id, ")"))
  }

  ran_parts$mm_f <- mm_f

  # build Z0
  Z0 <- build_Z0(data = data, ran_parts = ran_parts)
  # build real Z
  ran_parts$Z <- build_Z(Z0 = Z0, data = data, id = ran_parts$id)

  # matrix for mixed model estimation
  Z0[,all.vars(formula[[2]])] <- data[, all.vars(formula[[2]])] # y
  Z0[,ran_parts$id] <- as.factor(data[,ran_parts$id])
  ran_parts$Z0 <- Z0

  return(list(data = data, ran_parts = ran_parts))

}


build_family <- function(ran_parts, weights) {
  Z0 <- ran_parts$Z0
  family <- ran_parts$family
  if (is.function(family)) {
    family <- family()
  }

  if (class(family) == "boost_family") {
    if (family@name == "Negative Negative-Binomial Likelihood") {
      gm <- lme4::glmer.nb(ran_parts$mm_f, data = Z0, weights = weights)
      # Overwrite ran_parts$family in the parent environment of 'main'
      assign("ran_parts",
             modifyList(ran_parts,
                        list(family =  family(gm))),
             envir = parent.frame())
    } else if (class(family) == "boost_family_glm") {
      stop("Please use a R-stats family.")
    }
  } else {
    if (inherits(family, "family") && family$family == "gaussian") {
      gm <- lme4::lmer(ran_parts$mm_f, data = Z0, weights = weights)
    } else {
      gm <- lme4::glmer(ran_parts$mm_f, data = Z0, family = family, weights = weights)
    }
  }

  fam <- as.Family.merMod(object = gm, ran_parts = ran_parts)
  return(fam)
}



build_C <- function(data, X, ran_parts) {
  id <- ran_parts$id
  rs <- ran_parts$rs
  q <- ran_parts$q
  int_cond <- ran_parts$int_cond
  N <- nrow(data)
  p <- ncol(X)

  id_split <- as.numeric(factor(as.matrix(data[id]),
                                levels = unique(as.matrix(data[id]))))
  n <- length(unique(id_split))

  ### extract cluster-constant covariates
  first <- rep(FALSE, N)
  for(i in 1:N){
    first[which.max(id_split==i)] <- TRUE
  }
  Xcor <- list()
  if (int_cond) {
    ccc = rep(FALSE, p)
    for(r in 1:p){
      if(Reduce('+', lapply(split(X[,r], id_split),
                            function(x)(length(unique(x))))) == n){ccc[r] <- TRUE}
    }

    Xcc <- as.matrix(X[first,ccc])

    # Check if any column consists of a single unique value/if intercept is a baselearner
    single_value_check <- any(apply(Xcc, 2, function(col) length(unique(col)) == 1))
    if (!single_value_check) {
      Xcc <- cbind(1, X[first,ccc])
    }

    Xcor[[1]] <- Xcc %*% solve(crossprod(Xcc)) %*% t(Xcc)
  }

  if(q > length(Xcor)) {
    x = matrix(rep(1, n), n, 1)
    for(s in (length(Xcor) + 1):q) {
      Xcor[[s]] <- x %*% solve(crossprod(x)) %*% t(x)
    }
  }

  p1 <- rep(seq(1, q * n, q), q) + rep(0:(q - 1), each = n)
  p2 <- rep(seq(1, q * n, n), n) + rep(0:(n - 1), each = q)
  P1 <- Matrix::sparseMatrix(seq_along(p1), p1)
  P2 <- Matrix::sparseMatrix(seq_along(p2), p2)

  Xcor <- Matrix::bdiag(Xcor)
  Cm <- P2 %*% (diag(n * q) - Xcor) %*% P1

  return(Cm)
}





build_Z0 <- function(data, ran_parts) {
  if (ran_parts$int_cond) {
    ran_parts$q <- ran_parts$q + 1
    Z0 <- as.data.frame(cbind(1, data[,ran_parts$rs]))
    names(Z0) <- c("int", ran_parts$rs)
  } else {
    Z0 <- data[,ran_parts$rs]
    names(Z0) <- ran_parts$rs
  }

  return(Z0)
}




build_Z <- function(Z0, data, id) {
  id_num <- data[id][,1]
  un_id <- unique(id_num)
  n <- length(un_id)

  Z <- list()
  for (i in 1:n) {
    Z[[i]] <- as.matrix(Z0[id_num == un_id[i],])
  }
  #Z <<- Z
  Z <- Matrix::bdiag(Z)

  return(Z)
}






