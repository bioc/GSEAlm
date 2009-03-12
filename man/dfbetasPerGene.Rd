\name{dfbetasPerGene}
\alias{dffitsPerGene}
\alias{dfbetasPerGene}
\alias{CooksDPerGene}
\alias{Leverage}
\title{Linear-Model Deletion Diagnostics for Gene Expression (or similar)
  Data Structures}
\description{
This is an extension of standard linear-model diagnostics for use with
gene-expression datasets, in which the same model was run simultaneously
on each row of a response matrix.
}

\usage{
 dfbetasPerGene(lmobj)

 CooksDPerGene(lmobj)

 dffitsPerGene(lmobj)

 Leverage(lmobj)
}

\arguments{
  \item{lmobj}{ An object produced by \code{\link{lmPerGene}}. }
}
\details{ Deletion diagnostics gauge the influence of each observation
  upon model fit, by calculating values after removal of the observation
  and comparing to the complete-data version.

  DFFITS\eqn{_i} measures the distance on the response scale, between fitted
  values with and without observation y\_i, at point i. The distance
  is normalized by the regression standard error and the point's
  leverage (see below).

  Cook's D\eqn{_i} is the square of the distance, in parameter space, between
  parameter estimates witn and without observation y\_i, normalized and
  rescaled by standard errors and by a factor depending upon leverage.

  DFBETAS\eqn{_{i,j}} breaks the square root of Cook's D into its Euclidean
  components for each parameter j - but uses a somewhat different
  scaling function from Cook's D.

  The leverage is the diagonal of the "hat matrix"
  \eqn{X'(X'X)^{-1}X'}. This measure provides the relative weight of
  observation y\_i in the fitted value y-hat\_i. Typically observations
  with extreme X values (or belonging to smaller groups if model variables are
  categorical) will have high leverage.

  All these functions exist for standard regression, see
  \code{\link{influence.measures}}.

  The functions described here are extensions for the case in which the
  response is a matrix, and the same linear model is run on each row separately.

  For more details, see the references below.

  All functions are implemented in
  matrix form, which means they run quite fast.

}
\value{
 \code{\link{dfbetasPerGene}} A G x n x p array, where G, n are the
 number of rows and columns in the input's expression matrix, respectively,
 and p the number of parameters in the linear model (including intercept)

 \code{\link{CooksDPerGene}} A G x n matrix.

 \code{\link{dffitsPerGene}} A G x n matrix.

 \code{\link{Leverage}} A vector of length n, corresponding to the
 diagonal of the "hat matrix".
}

\references{

     Belsley, D. A., Kuh, E. and Welsch, R. E. (1980) Regression
     Diagnostics. New York: Wiley.

     Cook, R. D. and Weisberg, S. (1982) Residuals and Influence in
     Regression. London: Chapman and Hall.

     Williams, D. A. (1987) Generalized linear model diagnostics using
     the deviance and single case deletions. Applied Statistics *36*,
     181-191.

     Fox, J. (1997) Applied Regression, Linear Models, and Related
     Methods. Sage.

     LaMotte, L. R. (1999) Collapsibility hypotheses and diagnostic bounds in
     regression analysis. Metrika 50, 109-119.

     Jensen, D.R. (2001) Properties of selected subset diagnostics in
     regression. Statistics and Probability Letters 51, 377-388.


}

\author{Robert Gentleman, Assaf Oron}

\note{ The commonly-cited reference alert thresholds for diagnostic measures such as Cook's $D$
  and DFBETAS, found in older references, appear to be out of date. See LaMotte (1999) and Jensen
  (2001) for a more recent discussion. Our suggested practice is to
  inspect any samples or values that are visibly separate from the pack. }

\seealso{ \code{\link{influence.measures}} for the analogous simple
  regression diagnostic functions }
\examples{
data(sample.ExpressionSet)
layout(1)
lm1 = lmPerGene( sample.ExpressionSet,~score+type)
CD = CooksDPerGene(lm1)
### How does the distribution of mean Cook's distances across samples look?

boxplot(log2(CD) ~ col(CD),names=colnames(CD),ylab="Log Cook's
Distance",xlab="Sample")
### There are a few gross individual-observation outliers (which is why we plot on the log
### scale), but otherwise no single sample pops out as problematic. Here's
### one commonly-used alert level for problems:
lines(c(-5,30),rep(log2(2/sqrt(26)),2),col=2)


DFB = dfbetasPerGene(lm1)

### Looking for simultaneous two-effect outliers - 500 genes times 26
### samples makes 13000 data points on this plot

plot(DFB[,,2],DFB[,,3],main="DFBETAS for Score and Type (all genes)",xlab="Score Effect
Offset (normalized units)",ylab="Type Effect Offset (normalized units)",pch='+',cex=.5)
lines(c(-100,100),rep(0,2),col=2)
lines(rep(0,2),c(-100,100),col=2)

DFF = dffitsPerGene(lm1)
summary(apply(DFF,2,mean))

Lev = Leverage(lm1)
table(Lev)
### should have only two unique values because this is a dichotomous one-factor model
}
