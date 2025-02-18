\name{glmermboost}
\alias{glmermboost}
\alias{glmermboost.formula}

\title{
Component-wise Gradient Boosting for Generalised Mixed Models
}
\description{
Gradient boosting for optimizing negative log-likelihoods as loss functions where component-wise
  linear models are utilized as base-learners and an estimation of random components is guaranteed
  via a maximum likelihood approach.
}
\usage{
glmermboost(formula, data = list(), weights = NULL,
          offset = NULL, family = gaussian,
          na.action = na.omit, contrasts.arg = NULL,
          center = TRUE, control = boost_control(), oobweights = NULL, ...)
}

\arguments{
  \item{formula}{ a symbolic description of the model to be fit in the
  lme4-format including random effects. }
  \item{data}{ a data frame containing the variables in the model. }
  \item{weights}{ an optional vector of weights to be used in the fitting
    process. }
  \item{offset}{ a numeric vector to be used as offset (optional).}
    \item{family}{!! This is in contrast to usual mboost -
    "only" a \code{\link{family}} object is possible - except for \code{NBinomial()}.}
  \item{na.action}{a function which indicates what should happen when the data
    contain \code{NA}s.}
  \item{contrasts.arg}{a list, whose entries are contrasts suitable for input
    to the \code{contrasts} replacement function and whose names are
    the names of columns of \code{data} containing factors.
    See \code{\link{model.matrix.default}}.}
  \item{center}{logical indicating of the predictor variables are centered before fitting.}
  \item{control}{ a list of parameters controlling the algorithm. For
    more details see \code{\link{boost_control}}. }
  \item{oobweights}{ an additional vector of out-of-bag weights, which is
    used for the out-of-bag risk (i.e., if \code{boost_control(risk =
      "oobag")}).}
  \item{\dots}{ additional arguments passed to \code{\link{mboost_fit}}; currently none.}
}

\details{
  The warning "model with centered covariates does not contain intercept" is correctly given -
  the intercept is estimated via the mixed model.

  A (generalized) linear mixed model is fitted using a boosting algorithm based on component-wise
  univariate linear models. Additionally, a mixed model gets estimated in every iteration and added
  to the current fit. The fit, i.e., the regression coefficients and random effects, can be
  interpreted in the usual way. This particualy methodology is described in
  Knieper et al. (2025).
}
\value{
  The description of \code{\link{glmboost}} holds while some methods are newly implemented
  like \code{\link{predict}}, \code{\link{plot}} and \code{\link{mstop}}. Only the former one requires
  an furhter argument. Additionally, methods  \code{\link{VarCorr}}, \code{\link{ranef}},
  \code{\link{VarCorr}} are implemented specifically.
}
\references{
Lars Knieper, Torsten Hothorn, Elisabeth Bergherr and Colin Griesbach (2025),
  mermboost....
  \emph{Journal }, \bold{volume}, pages.
}


\seealso{
  See \code{\link{mermboost}} for the same approach using additive models.

  See \code{\link{mer_cvrisk}} for a cluster-sensitive cross-validation.
}
\examples{
data(Orthodont)

# are there cluster-constant covariates?
find_ccc(Orthodont, "Subject")

# fit initial model
mod <- glmermboost(distance ~ age + Sex + (1 |Subject),
                   data = Orthodont, family = gaussian,
                   control = boost_control(mstop = 100))

# let mermboost do the cluster-sensitive cross-validation for you
norm_cv <- mer_cvrisk(mod, no_of_folds = 10)
opt_m <- mstop(norm_cv)

# fit model with optimal stopping iteration
mod_opt <- glmermboost(distance ~ age + Sex + (1 |Subject),
                   data = Orthodont, family = gaussian,
                   control = boost_control(mstop = opt_m))

# use the model as known from mboost
# in additional, there are some methods knwon from lme4
ranef(mod_opt)
VarCorr(mod_opt)


#######################

set.seed(123)

# Parameters
n_groups <- 10        # Number of groups
n_per_group <- 50     # Number of observations per group
beta_fixed <- c(0.5, -0.3, 0.7)  # Fixed effects for intercept, covariate1, covariate2
sigma_random <- 1     # Random effect standard deviation

# Simulate random effects (group-specific)
group_effects <- rnorm(n_groups, mean = 0, sd = sigma_random)

# Simulate covariates
covariate1 <- rnorm(n_groups * n_per_group)
covariate2 <- rnorm(n_groups * n_per_group)

# Simulate data
group <- rep(1:n_groups, each = n_per_group)
random_effect <- group_effects[group]

# Linear predictor including fixed effects and random effects
linear_predictor <- beta_fixed[1] + beta_fixed[2] * covariate1 +
                  beta_fixed[3] * covariate2 + random_effect
prob <- plogis(linear_predictor)  # Convert to probabilities

# Simulate binomial outcomes
y <- rbinom(n_groups * n_per_group, size = 1, prob = prob)

# Combine into a data frame
sim_data <- data.frame(group = group, y = y,
                       covariate1 = covariate1,
                       covariate2 = covariate2)
sim_data$group <- as.factor(sim_data$group)



mod3 <- glmermboost(y ~ covariate1 + covariate2 + (1 | group),
                    data = sim_data, family = binomial())
bin_cv <- mer_cvrisk(mod3, no_of_folds = 10)
mstop(bin_cv)
}


\keyword{mixed models}
\keyword{regression}
