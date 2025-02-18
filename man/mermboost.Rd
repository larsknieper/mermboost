\name{mermboost}
\alias{mermboost}
\title{ Gradient Boosting for Additive Mixed Models }
\description{
  Gradient boosting for optimizing negative log-likelihoods as loss functions,
  where component-wise arbitrary base-learners, e.g., smoothing procedures,
  are utilized as additive  base-learners. In addition, every iteration estimates
  random component via a maximum likelihood approach using the current fit.
}

\usage{
mermboost(formula, data = list(), na.action = na.omit, weights = NULL,
       offset = NULL, family = gaussian, control = boost_control(),
       oobweights = NULL, baselearner = c("bbs", "bols", "btree", "bss", "bns"),
       ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{ a symbolic description of the model to be fit in the
  lme4-format including random effects.  }

  \item{data}{ a data frame containing the variables in the model. }
  \item{na.action}{a function which indicates what should happen when
    the data contain \code{NA}s.}
    \item{weights}{ (optional) a numeric vector of weights to be used in
    the fitting process.}
  \item{offset}{ a numeric vector to be used as offset (optional).}
  \item{family}{!! This is in contrast to usual mboost -
  "only" a \code{\link{family}} object is possible - except for \code{NBinomial()}.}
  \item{control}{ a list of parameters controlling the algorithm. For
    more details see \code{\link{boost_control}}. }
  \item{oobweights}{ an additional vector of out-of-bag weights, which is
    used for the out-of-bag risk (i.e., if \code{boost_control(risk =
      "oobag")}). This argument is also used internally by
    \code{cvrisk}. }
  \item{baselearner}{ a character specifying the component-wise base
    learner to be used: \code{\link{bbs}} means P-splines with a
    B-spline basis (see Schmid and Hothorn 2008), \code{\link{bols}}
    linear models and \code{\link{btree}} boosts stumps.
    \code{bss} and \code{bns} are deprecated.
    Component-wise smoothing splines have been considered in Buehlmann
    and Yu (2003) and Schmid and Hothorn (2008) investigate P-splines
    with a B-spline basis. Kneib, Hothorn and Tutz (2009) also utilize
    P-splines with a B-spline basis, supplement them with their
    bivariate tensor product version to estimate interaction surfaces
    and spatial effects and also consider random effects base
    learners.}
  \item{\dots}{ additional arguments passed to \code{\link{mboost_fit}}; currently none.}
}

\details{

  A (generalized) additive mixed model is fitted using a boosting algorithm based on
  component-wise base-learners. Additionally, a mixed model gets estimated in every
  iteration and added to the current fit.

  The base-learners can either be specified via the \code{formula} object or via
  the \code{baselearner} argument. The latter argument is the default base-learner
  which is used for all variables in the formula, whithout explicit base-learner
  specification (i.e., if the base-learners are explicitly specified in \code{formula},
  the \code{baselearner} argument will be ignored for this variable).

  Of note, \code{"bss"} and \code{"bns"} are deprecated and only in the list for
  backward compatibility.

  Note that more base-learners (i.e., in addition to the ones provided
  via \code{baselearner}) can be specified in \code{formula}. See
  \code{\link{baselearners}} for details.
}
\value{
  The description of \code{\link{mboost}} holds while some methods are newly implemented
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
  See \code{\link{glmermboost}} for the same approach using additive models.

  See \code{\link{mer_cvrisk}} for a cluster-sensitive cross-validation.
}
\examples{
data("Orthodont")

# are there cluster-constant covariates?
find_ccc(Orthodont, "Subject")

mod <- mermboost(distance ~ bbs(age, knots = 4) + bols(Sex) + (1 |Subject),
                 data = Orthodont, family = gaussian,
                 control = boost_control(mstop = 100))

# let mermboost do the cluster-sensitive cross-validation for you
norm_cv <- mer_cvrisk(mod, no_of_folds = 10)
opt_m <- mstop(norm_cv)

# fit model with optimal stopping iteration
mod_opt <- mermboost(distance ~ bbs(age, knots = 4) + bols(Sex) + (1 |Subject),
                 data = Orthodont, family = gaussian,
                 control = boost_control(mstop = opt_m))

# use the model as known from mboost
# in additional, there are some methods knwon from lme4
ranef(mod_opt)
VarCorr(mod_opt)
}

