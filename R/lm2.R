#### This file actually has the lmpergene function
### also:

### pvals from permutation matrix and qRequire (copied over
### from Category package where they are hidden functions)

### GSEA permutation test for one main effect in a multiple regression setting
### (generalization of gseattperm)



#in this function we  compute the linear regressions
#on a pergene basis
#note the NA handling is still not great, as we don't yet
# deal with NA's in the ExpressionSet

lmPerGene <- function(eSet,formula="", na.rm = TRUE) {

  if (formula=="") {

### Intercept-Only Model (default)
      nSamp = ncol(eSet)
      x=matrix(rep(1,nSamp),nrow=nSamp,ncol=1)
      colnames(x) = "Intercept"
      eS = eSet


  } else {
    xvarnames = all.vars(formula)
  badvars = !(xvarnames %in% varLabels(eSet))
  if( any(badvars) )
      stop("variable", xvarnames[badvars],
          "are not variable names in the supplied ExpressionSet")

  nvar <- length(xvarnames)


##drop any with missing values in some covariate
  if( na.rm ) {
      na = rep(FALSE, ncol(eSet))
      for (i in xvarnames)
         na = na | is.na(eSet[[i]])
      eS = eSet[, !na]
  } else
      eS = eSet

    nSamp=ncol(eS)
  xvar <- pData(eS)[, xvarnames, drop=FALSE]

    x = model.matrix(formula, data = xvar)
}
  xx <- crossprod(x)
  xxinv <- solve(xx)
  Hmat <- x %*% xxinv %*% t(x)
  dMat = diag(nSamp) - Hmat
  k = ncol(x)

  xy = exprs(eS) %*% x
  res = exprs(eS) %*% dMat
  beta = solve(xx, t(xy))
  varr = rowSums( exprs(eS) * res)/(nSamp - k)

  varbeta = matrix(diag(xxinv), ncol=k, nrow=nrow(eS),
     byrow=TRUE)
  varbeta = apply(varbeta, 2, function(x) x* varr)
  colnames(varbeta) = colnames(x)

  return(list(eS=eS, x = x, Hmat = Hmat, formula=formula,coefficients = beta,
              sigmaSqr = var, coef.var = t(varbeta)))
}

##### GSEA inference for main effect using multiple regression and permutation
##### Inference reported only for main effect, which MUST BE
##### THE FIRST VARIABLE IN THE FORMULA
##### This is because label permutations are done within each block defined by
##### level combinations of the other variables

##### This is an extension of "gseattperm"

gsealmPerm=function (eSet,formula="",mat,nperm,na.rm=TRUE,...) {

### For the most part we rely on 'lmPerGene' for formula validation, NA removal, etc. etc.

nSamp=ncol(eSet)

if (formula=="") {
        nvar=0
} else {
    xvarnames = all.vars(formula)
    nvar <- length(xvarnames)
}

### The observed t-values for the main effect

obsRaw=lmPerGene(eSet=eSet,formula=formula,na.rm=na.rm)

if (nvar>0) {
    observedStats= GSNormalize(obsRaw$coefficients[2,]/sqrt(obsRaw$coef.var[2,]),incidence=mat,...)
} else {
    observedStats= GSNormalize(t(obsRaw$coefficients),incidence=mat,fun2=identity,...)
}

### Permutation loop; we do the intercept-only case separately below

perm.eset=eSet

i <- 1L
if (nvar>0) {
    permMat <- matrix(0, nrow = nrow(eSet), ncol = nperm)

    while (i < (nperm + 1)) {

#### The crux (with nvar>=2): label permutation is done *within each covariate-combination level group separately*

        if (nvar>=2) {
            splitter=pData(eSet)[,xvarnames[2]]
            if (nvar>2) splitter=as.list(pData(eSet)[,xvarnames[2:nvar]])

            label.perm=unsplit(lapply(split(1:nSamp,splitter),sample),splitter)

### Now we only permute the labels of variable 1
            pData(perm.eset)[,xvarnames[1]]<-pData(eSet)[label.perm,xvarnames[1]]
        } else if (nvar==1) {
            pData(perm.eset)[,xvarnames[1]]<-pData(eSet)[sample(1:nSamp),xvarnames[1]]
        }
        temp.results<-lmPerGene(eSet=perm.eset,formula=formula,na.rm=na.rm)
# (na.rm=FALSE since we already dealt with na's)

### record t-score for permuted variable
        permMat[, i] <- temp.results$coefficients[2,]/sqrt(temp.results$coef.var[2,])

        i <- i + 1L
  }
  permMat <- GSNormalize(permMat,incidence=mat,...)
  rownames(permMat)=rownames(mat)

} else if (nvar==0) {

### Intercept only - row permutation and using raw expression means, no need for repeated calls to 'lmPerGene'

    permMat <- matrix(0, nrow = nrow(mat), ncol = nperm)
    rownames(permMat)=rownames(mat)

    for (i in 1:nperm)  permMat[,i]=GSNormalize(t(obsRaw$coefficients),incidence=mat[,sample(1:ncol(mat))],fun2=identity,...)
}
  return(pvalFromPermMat(observedStats, permMat))
}


#### unexported functions

pvalFromPermMat <- function(obs, perms) {
    N <- ncol(perms)
    pvals <- matrix(as.double(NA), nr=nrow(perms), ncol=2)
    dimnames(pvals) <- list(rownames(perms), c("Lower", "Upper"))

    tempObs <- rep(obs, ncol(perms))
    dim(tempObs) <- dim(perms)
    pvals[ , 1] <- rowSums(perms <= tempObs)/N
    pvals[ , 2] <- rowSums(perms >= tempObs)/N
    pvals
}
