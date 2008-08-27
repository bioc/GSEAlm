
### plotting functions for GSEA diagnostics





####################3
# plotting functions: no error check here
####################3

### Gene-set indiviudal residual Boxplots by sample, by factor level -
### for 2-level factors 
### Each point is a single residual from a single gene and sample

resplot=function(GSname="All", resmat,incidence=dumminc(resmat), fac,
atomic="Gene",core.text="Residuals by Sample", 
yname="Standardized Residual", xname="Sample ID", ID=colnames(resmat),
lims=0,gnames=levels(fac), prefix="", horiz=FALSE,
colour=5,pch='+',...)  {

setsize=sum(incidence[GSname,]>0)
atomics=paste(atomic,"s",sep="")
fac=factor(fac)
k=nlevels(fac)

lengths=table(fac)
tnames=levels(fac)

layout(1:k)
if (horiz) layout(t(1:k))

if (length(lims)==1) lims=range(resmat[incidence[GSname,]>0,])
lass=ifelse(horiz,1,3)
if(horiz) {
            timp=xname
            xname=yname
            yname=timp }

for (a in 1:k) {
    title1=paste(prefix,GSname,gnames[a],core.text, "(",setsize,atomics,")")
    boxplot(split(resmat[incidence[GSname,]>0,fac==tnames[a]],col(resmat[incidence[GSname,]>0,fac==tnames[a]])),main=title1,cex=.7,cex.axis=.8,names=ID[fac==tnames[a]],las=lass,pch=pch,range=1,ylim=lims,xlab=xname,ylab=yname,boxwex=0.8*min(1,log(lengths[a])/log(mean(lengths))),horizontal=horiz,col=colour,...)
    if (horiz) { lines(rep(0,2),c(-5,length(ID)),col=2)
        } else lines(c(-5,length(ID)),rep(0,2),col=2)
}
}

### Same, but stripcharts for smaller-size gene sets

restrip=function(GSname="All", resmat, incidence=dumminc(resmat), fac,
atomic="Gene", core.text="Residuals by Sample",
yname="Standardized Residual", xname="Sample ID", ID=colnames(resmat),
gnames=levels(fac), prefix="", colour=c(2:4,6), resort=TRUE,
horiz=FALSE, resort.fun=num.positive, pch='+',...) {

vert=!horiz    
setsize=sum(incidence[GSname,]>0)
myraw=resmat[incidence[GSname,]>0,]
atomics=paste(atomic,"s",sep="")


if (resort) {
   colsort=sort(apply(myraw,2,resort.fun),index=T)
   myraw=myraw[,colsort$ix]
   fac=fac[colsort$ix]
   ID=ID[colsort$ix]
}

fac=factor(fac)
k=nlevels(fac)
lengths=table(fac)
tnames=levels(fac)

layout(1:k)
if (!vert) layout(t(1:k))

#if (lims==0) lims=range(resmat[incidence[GSname,]>0,])
lass=ifelse(!vert,1,3)
if(!vert) {
            timp=xname
            xname=yname
            yname=timp }

for (a in 1:k) {
    title1=paste(prefix,GSname,gnames[a],core.text, "(",setsize,atomics,")")

    stripchart(split(myraw[,fac==tnames[a]],col(myraw[,fac==tnames[a]])),cex=.8,cex.axis=.8,group.names=ID[fac==tnames[a]],xlab=xname,ylab=yname,las=lass,pch=pch,vert=vert,col=colour,...)
    title(title1)
    if (vert) { lines(c(-5,length(ID)),rep(0,2))
     } else lines(rep(0,2),c(-5,length(ID)))
}

}

############# An upgrade (or downgrade?) to the MNplot

mnDiffPlot=function(GSname="All", exprmat, incidence=dumminc(exprmat),
fac, atomic="Gene", core.text=paste("Mean Expression Difference by",
                    atomic),
yname="Log Expression Ratio", xname="Log Expression",
gnames=levels(factor(fac)), prefix="", fitline=FALSE, varsize=FALSE,
reverse=FALSE, ...) {

layout(1)
fac=factor(fac)
lengths=table(fac)
tnames=levels(fac)
setsize=sum(incidence[GSname,]>0)
atomics=paste(atomic,"s",sep="")

if (reverse) {
    lengths=rev(lengths)
    gnames=rev(gnames)
    tnames=rev(tnames)
}
title0=paste(prefix,GSname,core.text, "(",setsize,atomics,")")
title2=paste(yname,"(",gnames[2],"/",gnames[1],")")

xvals=apply(exprmat[incidence[GSname,]>0,fac==tnames[1]],1,mean)
yvals=apply(exprmat[incidence[GSname,]>0,fac==tnames[2]],1,mean)-xvals

if (varsize) {
    xerrs=apply(exprmat[incidence[GSname,]>0,fac==tnames[1]],1,var)/lengths[1]
    yerrs=apply(exprmat[incidence[GSname,]>0,fac==tnames[2]],1,var)/lengths[2]+xerrs
}
sizes=1
if(varsize) sizes=sqrt(yerrs/mean(yerrs))
ylims=range(yvals,na.rm=TRUE)
if(ylims[1]>0) ylims[1]=0
if(ylims[2]<0) ylims[2]=0

plot(xvals,yvals,main=title0,xlab=paste (gnames[1],":",xname),ylab=title2,cex=sizes,ylim=ylims,...)
abline(0,0,col=2)
if (fitline) lines(sort(xvals),predict(loess(yvals~xvals),newdata=sort(xvals)),col=4)
}


####################### unexported utilities

### dummy incidence matrix for getting all resids

dumminc=function(resmat) {

myout=matrix(1,nrow=1,ncol=nrow(resmat))
rownames(myout)="All"
myout}

num.positive=function(x) sum(x>0)
num.extreme=function(x,lo,hi) sum(x<lo | x>hi)
num.low=function(x,lo) sum(x<lo)
num.high=function(x,hi) sum(x>hi)

