## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(fig.width = 4.5, fig.height = 3.5)

## ---- message = FALSE----------------------------------------------------
library(deconvSeq)

## ------------------------------------------------------------------------
data("data_celltypes_rnaseq") 

## ------------------------------------------------------------------------
set.seed(1234)
b0 = getb0.rnaseq(dge.celltypes, design.rnaseq, ncpm.min=1, nsamp.min=4)

## ------------------------------------------------------------------------
nb0 = 50
resultx1 = getx1.rnaseq(nb0,b0,dge.celltypes)

## ------------------------------------------------------------------------
x2 = as.data.frame(design.rnaseq,row.names=sample.id.rnaseq)
cr = getcorr(resultx1$x1,x2)
plot(cr, ylim=c(0,1), ylab="Correlation")

## ----results="hide"------------------------------------------------------
data("data_tissue_rnaseq") 
dge.tissue = getdge(cnts.tissue,design=NULL, ncpm.min=1, nsamp.min=4)

## ----results="hide"------------------------------------------------------
nb0=50
resultx1.tissue = getx1.rnaseq(nb0,b0, dge.tissue)

## ------------------------------------------------------------------------
x1 = cbind(lymph=resultx1.tissue$x1$tissuefBcell+resultx1.tissue$x1$tissuefTcell, 
 mono = resultx1.tissue$x1$tissuefMonocytes, gran = resultx1.tissue$x1$tissuefGranulocytes)
x2 = as.matrix(cbc.rnaseq/100)
cr = getcorr(x1,x2)
plot(cr, ylim=c(0,1), ylab="Correlation")

## ------------------------------------------------------------------------
data("data_celltypes_rrbs") 

## ----results="hide"------------------------------------------------------
set.seed(1234)
b0 = getb0.biseq(methmat, design.rrbs)

## ------------------------------------------------------------------------
nb0=250
resultx1 = getx1.biseq(nb0,b0,methmat,sample.id.rrbs,celltypes.rrbs)

## ------------------------------------------------------------------------
x2 = as.data.frame(design.rrbs,row.names=sample.id.rrbs)
cr=getcorr(resultx1$x1,x2)
plot(cr, ylim=c(0,1), ylab="Correlation")

## ------------------------------------------------------------------------
data("data_tissue_rrbs")

## ------------------------------------------------------------------------
nb0=250
resultx1.tissue = getx1.biseq(nb0,b0,methmat.tissue,sample.id.tissue,celltypes.rrbs)

## ------------------------------------------------------------------------
x1 = cbind(lymph=resultx1.tissue$x1[,1]+resultx1.tissue$x1[,2], mono = resultx1.tissue$x1[,3], gran = resultx1.tissue$x1[,4])
x2 = as.matrix(cbc.rrbs/100)
cr = getcorr(x1,x2)
plot(cr, ylim=c(0,1), ylab="Correlation")

