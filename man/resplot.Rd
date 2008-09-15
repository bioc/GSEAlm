\name{resplot}
\alias{resplot}
\alias{restrip}
\alias{mnDiffPlot}

\title{ Simple Graphical Summaries for Gene Set Enrichment Analysis
  (GSEA)}
\description{
Diagnostic plots for GSEA. 'resplot' and 'restrip' group residuals (or
  expression levels) from a specific gene-set by sample. 'mnDiffPlot'
  shows mean expression differences for a dichotomous phenotype, by
  gene, for a specific gene set.
}
\usage{
resplot(GSname = "All", resmat, incidence = dumminc(resmat), fac,
     atomic = "Gene", core.text = "Residuals by Sample",
     yname = "Standardized Residual", xname = "Sample ID",
     ID = colnames(resmat), lims = 0, gnames = levels(factor(fac)),
     prefix = "", horiz = FALSE, colour=5,pch='+',...)

restrip(GSname = "All", resmat, incidence = dumminc(resmat), fac,
   atomic = "Gene", core.text = "Residuals by Sample",
  yname = "Standardized Residual", xname = "Sample ID", ID = colnames(resmat),
  gnames = levels(factor(fac)), prefix = "", colour=c(2:4,6), resort=TRUE,
  horiz = FALSE, resort.fun=num.positive, pch='+', ...)

mnDiffPlot(GSname = "All", exprmat, incidence = dumminc(exprmat), fac,
  atomic = "Gene", core.text = paste("Mean Expression Difference by",atomic),
  yname="Log Expression Ratio", xname="Log Expression",
  gnames = levels(factor(fac)), prefix = "", fitline=FALSE, varsize=FALSE,
  reverse=FALSE, ...)}
\arguments{
  \item{GSname}{Gene-set Name. See "Details".}
  \item{resmat,exprmat}{Numerical matrix with the values to be
  plotted. See "Details". }
  \item{incidence}{ Gene-set 0/1 membership matrix }
  \item{fac}{The phenotypical variable to plot by. Must be discrete. For
  'mnDiffPlot', must be dichotomous.}
  \item{atomic}{ string identifying the meaning of rows in the data
    matrix. Defaults to "Gene". }
  \item{core.text,gnames,prefix,xname,yname}{strings controlling the
    text of main and axis captions }
  \item{ID}{Group names associated with the data matrix columns }
  \item{lims}{plotting limits for the response axis}
  \item{horiz}{logical: whether the boxplots or strips should be
    horizontal (defaults to FALSE)}
  \item{colour}{color of boxplot filling ('resplot') or symbols
    ('restrip')}
  \item{pch}{the plotting symbol}
  \item{resort} {('restrip' only) whether to sort groups for better
    visibility}
  \item{resort.fun}{('restrip' only) what function to sort groups
  by. Ignored unless 'resort==TRUE'. See \code{\link{stripchart}} documentation for
  more details}

\item{fitline}{ ('mnDiffPlot' only) logical: whether a loess fit should be
    plotted}
  \item{varsize}{ ('mnDiffPlot' only) logical: whether symbol sizes should
    be proportional to (t-test style) standard errors}
   \item{reverse}{ ('mnDiffPlot' only)logical: whether the factor's order
     should be reversed so that the second level is on the x-axis rather
     than the first one}
  \item{\dots}{ Additional graphical parameters passed on to the generic
    plotting functions.}
}
\details{
These functions provide simple graphical summaries for processed
gene-expression data, or other similar datasets for which matrix form is
useful. They are tailored predominantly for GSEA, but are useful in
general as well.

'resplot' calls \code{\link{boxplot}} and 'restrip' calls
\code{\link{stripchart}}; both summarize *all* data points from those
rows in 'resmat' which are members in the gene-set specified by
'GSname'. The summary is by column. For each level of 'fac' there will
be a separate pane.

'mnDiffPlot' calls \code{\link{plot}}; it plots the mean differences, by
row, between columns belonging to the two groups specified by 'fac', as
a function of the mean values for the first group alone. Each row
translates to a single point on the graph. Again, the summary is only
for rows indicated by 'GSname'.

For gene-set selective plots to properly work, the incidence matrix
needs to have non-empty row names, and 'GSname' must match one of them.

If both 'GSname' and 'incidence' are left blank, automatic
utilities are called which help generate a summary of the entire matrix,
by column.

All functions plot a reference line signalling zero. 'mnDiffPlot' also
optionally plots a loess fit for expression differences (if 'fitline=TRUE').

}

\author{Assaf Oron }
\note{ One can use 'resplot'/'restrip' to plot raw expression values rather
than residuals; it all depends on what's in the data matrix.


}
\seealso{ \code{\link{boxplot}},\code{\link{plot}},\code{\link{stripchart}},\code{\link{par}},\code{\link[annotate]{GOmnplot}} }
\examples{

data(sample.ExpressionSet)
lm1 = lmPerGene(sample.ExpressionSet,~sex)
r1 = getResidPerGene(lm1)
### now a boxplot of all residuals by sample
resplot(resmat=exprs(r1),fac=sample.ExpressionSet$sex)
### This plot is not very informative because of some gross outliers;
### try this instead
resplot(resmat=exprs(r1),fac=sample.ExpressionSet$sex,lims=c(-5,5))

### stripchart for first 10 genes
restrip(resmat=exprs(r1)[1:10,],fac=sample.ExpressionSet$type,prefix="Not")

### note the wild trajectory of the loess fit:
mnDiffPlot(exprmat=exprs(sample.ExpressionSet),fac=sample.ExpressionSet$type,xname="Raw Expression",yname="Expression Difference",fitline=TRUE)
}

\keyword{ methods }
\keyword{ hplot }
