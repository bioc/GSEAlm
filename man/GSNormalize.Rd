\name{GSNormalize}
\alias{GSNormalize}
\alias{identity}
\alias{one}
\title{Aggregating and calculating expression statistics by Gene Set}
\description{
 Provides an interface for producing aggregate gene-set statistics, for gene-set-enrichment analysis (GSEA). 
 The function is best suited for mean or rescaled-mean GSEA approaches, but is hopefully generic enough 
 to enable other approaches as well.
}
\usage{
GSNormalize(dataset, incidence, gseaFun = crossprod, fun1 = "/", fun2 = sqrt, removeShift=FALSE, removeStat=mean, ...)
identity(x)
one(x)
}
\arguments{
  \item{dataset}{ a numeric matrix, typically of some gene-level statistics }
  \item{incidence}{ 0/1 incidence matrix indicating genes' membership in
  gene-sets}
  \item{gseaFun}{function name for the type of aggregation to take place, defaults to 'crossprod'. See 'Details' }
  \item{fun1}{function name for normalization, defaults to "/". See 'Details' }
  \item{fun2}{function name for scaling, defaults to 'sqrt'. See 'Details'}
  \item{removeShift}{logical: should normalization begin with a column-wise removal of the mean shift?}
  \item{removeStat}{(if above is TRUE) the column-wise statistic to be swept out of 'dataset'.}
  \item{...}{Additional arguments optionally passed on to 'gseaFun'.}
  \item{x}{any numerical value}
}
\details{
In gene-set-enrichment analysis (GSEA), the core step is
aggregating (or calculating) gene-set-level statistics from gene-set statistics. This
utility achieves the feat. It is tailored specifically for rescaled-sums of the type
suggested by Jiang and Gentleman (2007), but is designed as a generic template that should other GSEA approaches. 
In such cases, at this moment users should provide their own version of 'gseaFun'.

The default will generate sums of gene-level values divided by the square-root of the
gene-set size (in other words, gene-set means multiplied by the
square-root of gene-set size). The arithmetic works like this:

gene-set stat = gseaFun(t(incidence),dataset),...) 'fun1' fun2(gene-set size).

In case there is a known (or suspected) overall baseline shift (i.e.,
the mass of gene-level stats is not centered around zero) it may be
scientifically more meaningful to look for gene-set deviating from this
baseline rather than from zero. In this case, you can set 'removeShift=TRUE'.

Also provided are the 'identity' function (identity = function(x) x), so that
leaving 'gseaFun' and 'fun1' at their default and setting 'fun2 = identity' will generate gene-set
means -- and the 'one' function to neutralize the effect of both 'fun1' and
'fun2' (see note below).


}
\value{
'GSNormalize' returns a matrix with the same number of rows as 'incidence' and the same number
of columns as 'dataset' (if 'dataset' is a vector, the output will be a
vector as well).
The respective row and column names will carry through from 'dataset'
and 'incidence' to the output.

'identity' simply returns x.
'one' returns the number 1.
}
\references{
 Z. Jiang and R. Gentleman, "Extensions to Gene Set Enrichment
Analysis",Bioinformatics (23),306-313, 2007.}
\author{Assaf Oron }


\note{ If you want to create your own GSEA function for 'gseaFun', note that it should receive the transposed
incidence matrix as its first argument, and the gene-level stats as its second argument. 
In other words, both should have genes as rows. also, you can easily neutralize the effect of 'fun1', 'fun2' by setting
"fun2 = one".
}
\seealso{ \code{\link{gsealmPerm}}, which relies heavily on this
  function. The function  \code{\link[Category]{applyByCategory} from
    the \code{\link[Category]{Category}} package has similar functionality and is
    preferable when the applied function is
    complicated. \code{\link{GSNormalize}} is better optimized for
    matrix operations.}
\examples{
data(sample.ExpressionSet)
lm1 = lmPerGene(sample.ExpressionSet,~sex+type)

### Generating random pseudo-gene-sets
fauxGS=matrix(sample(c(0,1),size=50000,replace=TRUE,prob=c(.9,.1)),nrow=100)

### "tau-stats" for gene-SET-level type effect, adjusting for sex
fauxEffects=GSNormalize(lm1$coefficients[3,]/sqrt(lm1$coef.var[3,]),incidence=fauxGS)

qqnorm(fauxEffects)
### diagonal line represents zero-shift null; note that it doesn't fit
abline(0,1,col=2)
### a better option may be to run a diagonal through the middle of the
### data (nonzero-shift null, i.e. type may have an effect but it is the
### same for all gene-sets); note that if any outlier shows, it is a purely random one!

abline(median(fauxEffects),1,col=4)

#### Now try with baseline-shift removal

fauxEffects=GSNormalize(lm1$coefficients[3,]/sqrt(lm1$coef.var[3,]),incidence=fauxGS,removeShift=TRUE)

qqnorm(fauxEffects)
abline(0,1,col=2)


}

\keyword{ methods }
