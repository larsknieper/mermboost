\name{mer_cvrisk}
\alias{mer_cvrisk}
\alias{cvrisk.mermboost}

\title{ Cluster-sensitive Cross-Validation }
\description{
  Cross-validated estimation of the empirical risk for hyper-parameter selection.
  Folds are created cluster-sensitive, hence splitting data into train and tests sets
  considers the cluster-structure.
}
\usage{
mer_cvrisk(object, folds, no_of_folds, cores = 1)
}
\arguments{
  \item{object}{ an object of class \code{mermboost}.}
  \item{folds}{ a weight matrix with number of rows equal to the number
                of observations. The number of columns corresponds to
                the number of cross-validation runs. Can be computed
                using function \code{cv}.}
  \item{no_of_folds}{ creates the folds itself by taking the cluster structure into account. }
  \item{cores}{
    is passed on to \code{\link[parallel]{mclapply}} for parallel computing.
  }
}

\details{

  The number of boosting iterations is a hyper-parameter of the
  boosting algorithms implemented in this package. Honest,
  i.e., cross-validated, estimates of the empirical risk
  for different stopping parameters \code{mstop} are computed by
  this function which can be utilized to choose an appropriate
  number of boosting iterations to be applied.

  This function uses the cluster-identifier held in the \code{mermboost} object
  to split the data into cluster-sensitive folds if the corresponding argument
  \code{no_of_folds} is given.
  As this might lead to imbalanced splits the 1/0 matrix of folds can be given manually
  via the folds argument.

}

\value{
  An object of class \code{mer_cv}, containing the k-folds as a matrix,
  the corresponding estimates of the empirical risks, their average
  and the results optimal stopping iteration.
  \code{plot} and \code{mstop} methods are available.
}

\references{
Lars Knieper, Torsten Hothorn, Elisabeth Bergherr and Colin Griesbach (2025),
  mermboost....
  \emph{Journal }, \bold{volume}, pages.
}


\examples{
data(Orthodont)

mod <- mermboost(distance ~ bbs(age, knots = 4) + bols(Sex) + (1 |Subject),
                 data = Orthodont, family = gaussian,
                 control = boost_control(mstop = 100))

# let mermboost do the cluster-sensitive cross-validation for you
norm_cv <- mer_cvrisk(mod, no_of_folds = 10)
opt_m <- mstop(norm_cv)
plot(norm_cv)
}

