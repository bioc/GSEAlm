### This file has utilities related to incidence matrices


### This function does the normalized aggregation over gene-sets;
### It is now hopefully generic enough to enable any GSEA flavor

GSNormalize<-function(dataset,incidence,gseaFun=crossprod,fun1="/",fun2=sqrt,removeShift=FALSE,removeStat=mean,...) {

    dataset=as.matrix(dataset)

### Removal of column-wise mean shift

    if (removeShift) {

        colStats=apply(dataset,2,removeStat)
        dataset=sweep(dataset,2,STATS=colStats)
    }

    if (ncol(incidence) != nrow(dataset)) stop ("GSNormalize: non-conforming matrices")


    outm=gseaFun(t(incidence),dataset,...)


    rownames(outm)=rownames(incidence)
    colnames(outm)=colnames(dataset)

    normby = fun2(rowSums(incidence))

    outm = sweep (outm,1,normby,FUN=fun1)

}



identity<-function(x) x
one <- function(x) 1
