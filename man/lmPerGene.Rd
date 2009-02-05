\name{lmPerGene}
\alias{lmPerGene}
\title{ Fit linear model for each gene }
\description{
  For each gene, \code{lmPerGene} fits the same, user-specified linear model.
  It returns the estimates of the model parameters and their variances
  for each fitted model. The function uses matrix algebra so it is much
  faster than repeated calls to \code{lm}.
}
\usage{
lmPerGene(eSet, formula, na.rm=TRUE)
}
\arguments{
  \item{eSet}{ An \code{ExpressionSet} object.}
  \item{formula}{ an object of class \code{\link{formula}} (or one that can be coerced to
          that class), specifying only the right-hand side starting with
      the '~' symbol. The LHS is automatically set as the expression
      levels provided in  \code{eSet}. The names of all predictors must
    exist in the phenotypic data of \code{eSet}.}
  \item{na.rm}{ Whether to remove missing observations. }
}
\details{
 This function efficiently computes the least squares fit of a linear
 regression to a set of gene expression values. We assume that there
 are \code{G} genes, on \code{n} samples, and that there are \code{p} variables in
 the regression equation.  So the result is that \code{G} different regressions
 are computed, and various summary statistics are returned.

  Since the independent variables are the same in each model fitting,
  instead of repeatedly fitting linear model for each gene,
  \code{lmPerGene} accelarates the fitting process by calculating the
  hat matrix \eqn{X(X'X)^(-1)X'} first.  Then matrix multiplication, and
  \code{solve} are to compute estimates of the model parameters.

  Leaving the formula blank (the default) will calculate an 
  intercept-only model, useful for generic pattern and outlier identification.
}
\value{
  A list with components:
  \item{eS}{The \code{ExpressionSet} used in the model fitting.}
  \item{x}{The design matrix of the coded predictor variables.}
  \item{Hmat}{The Hat matrix.}
  \item{coefficients}{A matrix of dimension \code{p} by \code{G}
   containing the estimated model parameters.}
  \item{sigmaSqr}{A vector of length $G$ containing the mean square error
    for that model, the sum of the residuals squared divided by \code{n - p}.}
  \item{coef.var}{A matrix of dimension \code{p} by \code{G} containing
   the estimated variances for the model parameters, for each regression.}
  \item{tstat}{A matrix of the same dimension as \code{coefficients},
  containing the $t$-statistics for each model estimate. This is simply
  \code{coefficients} divided by the square root of  \code{coef.var},
  and is provided for convenience.} 
}
\author{ Robert Gentleman, Assaf Oron }
\seealso{\code{\link{getResidPerGene}} to extract row-by-row residuals; \code{\link{gsealmPerm}} for
  code that utilizes  'lmPerGene' for gene-set-enrichment analysis (GSEA); and  \code{\link{CooksDPerGene}} for diagnostic functions on
  an object produced by 'lmPerGene'. Applying a by-gene regression in
  the manner performed here is a special case of a more generic
  linear-model framework available in the \code{\link[GlobalAncova]{GlobalAncova}}
  package; our assumption here is equivalent to a diagonal covariance structure
  between genes, with unequal variances.}

\examples{
data(sample.ExpressionSet)
layout(1)
lm1 = lmPerGene(sample.ExpressionSet,~sex)
qqnorm(lm1$coefficients[2,]/sqrt(lm1$coef.var[2,]),main="Sample Dataset: Sex Effect by Gene",ylab="Individual Gene t-statistic",xlab="Normal Quantile")
abline(0,1,col=2)
lm2 = lmPerGene(sample.ExpressionSet,~type+sex)
qqnorm(lm2$coefficients[2,]/sqrt(lm2$coef.var[2,]),main="Sample Dataset: Case vs. Control Effect by Gene, Adjusted for Sex",ylab="Individual Gene t-statistic",xlab="Normal Quantile")
abline(0,1,col=2)
}
\keyword{methods}
