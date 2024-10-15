\name{methods}
\alias{mermboost_methods}

\alias{predict.mermboost}
\alias{predict.glmermboost}

\alias{ranef.glmermboost}
\alias{ranef.mermboost}

\alias{VarCorr.glmermboost}
\alias{VarCorr.mermboost}

\alias{mstop.mer_cv}
\alias{plot.mer_cv}


\title{ Methods for Gradient Boosting for Mixed Models Objects }
\description{ Methods for models fitted by mixed model boosting algorithms. }
\usage{
\method{predict}{mermboost}(object, newdata = NULL, RE = TRUE,
        type = c("link", "response", "class"), which = NULL,
        aggregate = c("sum", "cumsum", "none"), ...)
\method{predict}{glmermboost}(object, newdata = NULL, RE = TRUE,
        type = c("link", "response", "class"), which = NULL,
        aggregate = c("sum", "cumsum", "none"), ...)

\method{ranef}{mermboost} (object, iteration = mstop(object))
\method{ranef}{glmermboost} (object, iteration = mstop(object))

\method{VarCorr}{mermboost}(object, iteration = mstop(object))
\method{VarCorr}{glmermboost}(object, iteration = mstop(object))

\method{mstop}{mer_cv}(cv_object)
\method{plot}{mer_cv}(cv_object, main = NULL)
}

\arguments{
  \item{object}{ objects of class \code{glmermboost}, \code{mermboost} }
  \item{newdata}{ optionally, a data frame in which to look for variables with
          which to predict. In case the model was fitted using the \code{matrix}
          interface to \code{\link{glmermboost}}, \code{newdata} must be a \code{matrix}
          as well (an error is given otherwise).
          If \code{RE = TRUE} but not the same cluster-identifier is found in the newdata object,
          it gets set to FALSE, \code{RE = FALSE}. }
  \item{RE}{ a logical values (\code{TRUE/FALSE}) indicating whether to include random effects. }
  \item{which}{ a subset of base-learners to take into account for computing
                predictions or coefficients. If \code{which} is given
                (as an integer vector or characters corresponding
                 to base-learners) a list or matrix is returned.
                 This ignores the random effects. }
  \item{type}{ the type of prediction required.  The default is on the scale
          of the predictors; the alternative \code{"response"} is on
          the scale of the response variable.  Thus for a
          binomial model the default predictions are of log-odds
          (probabilities on logit scale) and \code{type = "response"} gives
          the predicted probabilities.  The \code{"class"} option returns
          predicted classes for binomial data.}
  \item{aggregate}{ a character specifying how to aggregate predictions
                    or coefficients of single base-learners. The default
                    returns the prediction or coefficient for the final number of
                    boosting iterations. \code{"cumsum"} returns a
                    matrix (one row per base-learner) with the
                    cumulative coefficients for all iterations
                    simultaneously (in columns). \code{"none"} returns a
                    list of matrices where the \eqn{j}th columns of the
                    respective matrix contains the predictions
                    of the base-learner of the \eqn{j}th boosting
                    iteration (and zero if the base-learner is not
                    selected in this iteration). Therefore, no random effects
                    are considered.}
  \item{iteration}{ an integer input that specifies from which iteration the
                    random component is to be drawn. }
  \item{cv_object}{ an object resulting from \code{\link{mer_cvrisk}}. }
  \item{main}{ user specific plot header. }
  \item{\dots}{ additional arguments passed to callies. }
}

\details{
  The methods should correspond the the \code{mboost} and \code{lme4} functions.
}


\references{
Lars Knieper, Torsten Hothorn, Elisabeth Bergherr and Colin Griesbach (2025),
  mermboost....
  \emph{Journal }, \bold{volume}, pages.
}



\seealso{
\code{\link{mstop.mer_cv}} and \code{\link{plot.mer_cv}}
}
\examples{
data(Orthodont)

mod <- glmermboost(distance ~ age + Sex + (1 |Subject),
                   data = Orthodont, family = gaussian,
                   control = boost_control(mstop = 50))

norm_cv <- mer_cvrisk(mod, no_of_folds = 10)
opt_m <- mstop(norm_cv)
plot(norm_cv)

any(predict(mod, RE = FALSE) == predict(mod, RE = TRUE))
all(predict(mod, RE = FALSE) ==
    predict.glmboost(mod) + mod$nuisance()[[mstop(mod)]]$ff
    )
ranef(mod)
VarCorr(mod, iteration = 10)
}



