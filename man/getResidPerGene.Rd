\name{getResidPerGene}
\alias{getResidPerGene}

\title{Row-by-Row Linear-Model Residuals for Gene Expression (or similar)
  Data Structures}

\description{
  This produces residuals of an identical linear model applied to each
  row of a gene expression matrix (or similar dataset). Computation
  speed is achieved via straightforward matrix algebra. Most
  commonly-used residual types are available.}

\usage{
getResidPerGene(lmobj, type = "extStudent")}
\arguments{
  \item{lmobj}{ An object produced by function \code{\link{lmPerGene}}. }
  \item{type}{ A string indicating the type of residual requeseted
  (defaults to externally-Studentized). }
}

\details{ Types of residuals now available:
  \describe{
    \item{"response"}{Response residuals, observed minus fitted}
    \item{"normalized"}{Response residuals divided by the estimated residual
      S.E.}
    \item{"intStudent"}{Internally Studentized residuals, often referred
      to as "Standardized"}
    \item{default}{Externally Studentized residuals, which can be used
    directly for outlier identification}
  }
}

\value{
  Returns a instance of \code{ExpressionSet} where the expression matrix
  contains the residuals. The \code{phenoData} are inherited from
  \code{lmobj$eS}.
}


\author{ Robert Gentleman, Assaf Oron }
\seealso{\code{\link{lmPerGene}}, \code{\link{resplot}},\code{\link{dfbetasPerGene}},\code{\link{influence.measures}} }
\examples{
data(sample.ExpressionSet)
lm1 = lmPerGene(sample.ExpressionSet,~sex)
r1 = getResidPerGene(lm1)
### now a boxplot of all residuals by sample
resplot(resmat=exprs(r1),fac=sample.ExpressionSet$sex)
### This plot is not very informative because of some gross outliers;
### try this instead
resplot(resmat=exprs(r1),fac=sample.ExpressionSet$sex,lims=c(-5,5))

}
\keyword{ methods }
