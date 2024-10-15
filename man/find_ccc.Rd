\name{find_ccc}
\alias{find_ccc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find cluster-constant covariates in a data set.
}
\description{
This function gives out logical indicators whether a variable is cluster-constant.
}
\usage{
find_ccc(df, id_char)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{ some data frame}
  \item{id_char}{ a character which is the column name for the cluster identifier. }
}
\details{
For uncorrected boosting of a mixed model the estimates for cluster constant covariates
might be biased as part of their effect is held in the random effects.
This bias is corrected by the underlying mermboost package.
}

\references{
Lars Knieper, Torsten Hothorn, Elisabeth Bergherr and Colin Griesbach (2025),
  mermboost....
  \emph{Journal }, \bold{volume}, pages.
}


\seealso{
\code{\link{glmermboost}} and \code{\link{mermboost}}
}

\examples{

data(Orthodont)
find_ccc(Orthodont, "Subject")
}


