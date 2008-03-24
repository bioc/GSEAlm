#### This file has mostly diagnostic functions

### get residuals
### Cook's D
### DFBETAS, DFFITS
### Leverage (hat-mat diagonal)

### All have been vector-matrixized now for quick implementation

dffitsPerGene <- function(lmobj) {
  obj <- lmobj
  n <- ncol(exprs(obj$eS))
  p <- ncol(obj$x)
  e <- tcrossprod(diag(n) - obj$Hmat,exprs(obj$eS))
  SSE <- sapply(featureNames(obj$eS), function(i) { exprs(obj$eS)[i,] %*% e[,i] } )
  hii <- diag(obj$Hmat)

  se <- sapply(colnames(e),function(i) {
    sqrt((SSE[i] * (1-hii)- e[,i] * e[,i])/ (n-p-1))
  } )
  t <- e/se
  dffits <- apply(t,2,function(u) {
    u * sqrt(hii / (1-hii))
  } )
  dffits=t(dffits)

  return(dffits)
}


CooksDPerGene <- function(lmobj) {
     rs = getResidPerGene(lmobj,type="intStudent")
### residuals have to be internally Studentized for the shortcut formula to be correct
     p = ncol(lmobj$x)

     D = exprs(rs)^2
     sweep(D, 2, diag(lmobj$Hmat)/(p*(1-diag(lmobj$Hmat))), "*")
 }


dfbetasPerGene<- function(lmobj) {
### Shortcut formula from Jensen and Ramirez
  obj <- lmobj
  xx <- crossprod(obj$x)
  xxinv <- solve(xx)
  n <- ncol(exprs(obj$eS))
  p <- ncol(obj$x)

  Db <- array(0,dim=c(nrow(exprs(obj$eS)),n,p),dimnames=list(featureNames(obj$eS),sampleNames(obj$eS),rownames(obj$coefficients)))

  tees=exprs(getResidPerGene(lmobj))
  oneMinusH=1-Leverage(lmobj)
  xxx=tcrossprod(xxinv,obj$x)

  for (k in 1:p) {

      Db[,,k] = t(t(tees)*xxx[k,]/sqrt(xxinv[k,k]*oneMinusH))
  }

Db
}

getResidPerGene <-  function(lmobj, type="extStudent") {
### Residuals from a per-gene linear model
### Uses matrix algebra to make a faster calculation than doing them one by one
### Accepted types: response, normalized, internally Studentized (a.k.a. standardized) and externally Studentized (default)

    nSamp = ncol(lmobj$eS)
  dMat = diag(nSamp) - lmobj$Hmat
  e = crossprod(t(exprs(lmobj$eS)),dMat)


  p <- ncol(lmobj$x)

### The code (perhaps primitively) proceeds by successive adjustments needed to change the residual type, with a possible exit at each step:

if (type=="response")   return(new("ExpressionSet",exprs=e,phenoData=phenoData(lmobj$eS)))

### The following is a simple normalization by gene-specific res.std.err.:

  sigma = sqrt(rowSums(e^2)/(nSamp-p))

  e=e/sigma

if (type=="normalized")   return(new("ExpressionSet",exprs=e,phenoData=phenoData(lmobj$eS)))

###  Internal Studentization
stud.e = sweep(e, 2, sqrt(1-diag(lmobj$Hmat)), "/")

if (type=="intStudent" || type=="standardized")   return(new("ExpressionSet",exprs=stud.e,phenoData=phenoData(lmobj$eS)))

#### External Studentization (the default)

stud.e=stud.e*sqrt((nSamp-p-1)/(nSamp-p-stud.e^2))

  return(new("ExpressionSet",exprs=stud.e,phenoData=phenoData(lmobj$eS)))
}

Leverage <- function(lmobj) {
  Hmat <- lmobj$Hmat
  ans <- diag(Hmat)
  return(ans)
}
