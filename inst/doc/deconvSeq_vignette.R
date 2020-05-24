## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(fig.width = 4.5, fig.height = 3.5)

## ---- message = FALSE---------------------------------------------------------
library(deconvSeq)

## -----------------------------------------------------------------------------
file1 = system.file("extdata","sample1_genecounts.txt", package="deconvSeq")
file2 = system.file("extdata","sample2_genecounts.txt", package="deconvSeq")
countmat = getrnamat(filnames=c(file1,file2),sample.id=c("sample1","sample2"))

## -----------------------------------------------------------------------------
data("data_celltypes_rnaseq") 

## -----------------------------------------------------------------------------

celltypes=c("Bcells","Tcells","Monocytes","Neutrophils")
#cell types of single cell type RNAseq or methylation files
file.celltypes = c(rep("Monocytes",4),rep("Tcells",4), rep("Bcells",3), rep("Neutrophils",13))
design=data.frame(Bcells=ifelse(file.celltypes=="Bcells",1,0),Tcells=ifelse(file.celltypes=="Tcells",1,0),Monocytes=ifelse(file.celltypes=="Monocytes",1,0),Neutrophils=ifelse(file.celltypes=="Neutrophils",1,0))


## -----------------------------------------------------------------------------

celltypes=c("Bcells","Tcells","Monocytes","Neutrophils")
file.celltypes = factor(c(rep("Monocytes",4),rep("Tcells",4), rep("Bcells",3), rep("Neutrophils",13)), levels = celltypes)
design <- model.matrix(~-1+file.celltypes) 


## ----results="hide"-----------------------------------------------------------
dge.celltypes = getdge(cnts.celltypes, design.rnaseq, ncpm.min=1, nsamp.min=4)

## -----------------------------------------------------------------------------
set.seed(1234)
b0 = getb0.rnaseq(dge.celltypes, design.rnaseq, ncpm.min=1, nsamp.min=4, sigg=NULL)

## -----------------------------------------------------------------------------
resultx1 = getx1.rnaseq(NB0="top_bonferroni",b0,dge.celltypes)

## -----------------------------------------------------------------------------
resultx1 = getx1.rnaseq(NB0=50,b0,dge.celltypes)

## -----------------------------------------------------------------------------
x2 = as.data.frame(design.rnaseq,row.names=sample.id.rnaseq)
cr = getcorr(resultx1$x1,x2)
plot(cr, ylim=c(0,1), ylab="Correlation")

## ----results="hide"-----------------------------------------------------------
data("data_tissue_rnaseq") 
dge.tissue = getdge(cnts.tissue,design=NULL, ncpm.min=1,nsamp.min=4)

## ----results="hide"-----------------------------------------------------------
resultx1.tissue = getx1.rnaseq(NB0=50,b0,dge.tissue)

## ----warning=FALSE, message=FALSE---------------------------------------------
x1 = cbind(lymph=resultx1.tissue$x1$tissuefBcell+resultx1.tissue$x1$tissuefTcell, 
 mono = resultx1.tissue$x1$tissuefMonocytes, gran = resultx1.tissue$x1$tissuefGranulocytes)
x2 = as.matrix(cbc.rnaseq/100)
cr = getcorr(x1,x2)
plot(cr, ylim=c(0,1), ylab="Correlation")

## -----------------------------------------------------------------------------
data("data_scrnaseq") 

## ----warning=FALSE, results='hide',message=FALSE------------------------------
cnts.sc = prep_scrnaseq(cnts.scrnaseq, genenametype = "hgnc_symbol",cellcycle=NULL,count.threshold=0.05)

## -----------------------------------------------------------------------------
cnts.sc.G1 = getcellcycle(cnts.sc,"G1")

## -----------------------------------------------------------------------------
cnts.sc.G1.train = cnts.sc.G1[,c(which(substr(colnames(cnts.sc.G1),3,6)=="Tcon")[1:250],which(substr(colnames(cnts.sc.G1),3,6)=="Treg")[1:150])]
cnts.sc.G1.valid = cnts.sc.G1[,-which(colnames(cnts.sc.G1) %in% colnames(cnts.sc.G1.train))]
tissue.sc = substr(colnames(cnts.sc.G1.train),3,6)
names(tissue.sc) = colnames(cnts.sc.G1.train)
sample.id.sc = colnames(cnts.sc.G1.train)
design.sc = model.matrix(~-1+as.factor(tissue.sc))
colnames(design.sc) = levels(as.factor(tissue.sc))
rownames(design.sc) = names(tissue.sc)
design.sc = design.sc[colnames(cnts.sc.G1.train),]

## ----results="hide"-----------------------------------------------------------
dge.sc = getdge(cnts.sc.G1.train,design.sc,ncpm.min=1, nsamp.min=4, method="bin.loess")
b0.sc = getb0.rnaseq(dge.sc, design.sc, ncpm.min=1, nsamp.min=4)

## ----results="hide"-----------------------------------------------------------
tissue_s.sc = substr(colnames(cnts.sc.G1.valid),3,6)
names(tissue_s.sc) = colnames(cnts.sc.G1.valid)
sample.id_s.sc = colnames(cnts.sc.G1.valid)
design_s.sc = model.matrix(~-1+as.factor(tissue_s.sc))
colnames(design_s.sc) = levels(as.factor(tissue_s.sc))
rownames(design_s.sc) = names(tissue_s.sc)
design_s.sc = design_s.sc[colnames(cnts.sc.G1.valid),]
dge_s.sc = getdge(cnts.sc.G1.valid, design_s.sc, ncpm.min=1, nsamp.min=4, method="bin.loess")
resultx1_s.sc = getx1.rnaseq(NB0=1500,b0.sc, dge_s.sc)

## ----warning=FALSE, message=FALSE---------------------------------------------
x2 = as.data.frame(design_s.sc,row.names=sample.id_s.sc)
sc = getcorr(resultx1_s.sc$x1,x2)
getmeancorr(sc)

## ----results="hide"-----------------------------------------------------------
#scRNAseq data
singlecelldata = cnts.sc.G1.train 
#known single cell types of the scRNAseq data
celltypes.sc = tissue.sc 
#tissue data with unknown cell types
tissuedata = cnts.sc.G1.valid 
#obtain design matrix from scRNAseq data 
design.singlecell = model.matrix(~-1+as.factor(celltypes.sc))
colnames(design.singlecell) = levels(as.factor(celltypes.sc))
rownames(design.singlecell) = names(celltypes.sc)
#obtain projection matrix
dge.singlecell = getdge(singlecelldata,design.singlecell,ncpm.min=1, nsamp.min=4, method="bin.loess")
b0.singlecell = getb0.rnaseq(dge.singlecell, design.singlecell, ncpm.min=1, nsamp.min=4)
#obtain cell type proportions in tissue
dge_tissue.sc = getdge(tissuedata, NULL, ncpm.min=1, nsamp.min=4, method="bin.loess")
resultx1_tissue.sc = getx1.rnaseq(NB0=1500,b0.singlecell, dge_tissue.sc)

## ---- message=FALSE, results=FALSE--------------------------------------------
#file1 = system.file("extdata","sample1_methratio.txt", package="deconvSeq")
#file2 = system.file("extdata","sample2_methratio.txt", package="deconvSeq")
#methmat = getmethmat(filnames=c(file1,file2), sample.id=c("sample1","sample2"), filtype="bsmap")

## -----------------------------------------------------------------------------
data("data_celltypes_rrbs") 

## -----------------------------------------------------------------------------

celltypes=c("Bcells","Tcells","Monocytes","Neutrophils")
#cell types of single cell type RNAseq or methylation files
file.celltypes = c(rep("Monocytes",4),rep("Tcells",4), rep("Bcells",3), rep("Neutrophils",13))
design=data.frame(Bcells=ifelse(file.celltypes=="Bcells",1,0),Tcells=ifelse(file.celltypes=="Tcells",1,0),Monocytes=ifelse(file.celltypes=="Monocytes",1,0),Neutrophils=ifelse(file.celltypes=="Neutrophils",1,0))


## -----------------------------------------------------------------------------

celltypes=c("Bcells","Tcells","Monocytes","Neutrophils")
file.celltypes = factor(c(rep("Monocytes",4),rep("Tcells",4), rep("Bcells",3), rep("Neutrophils",13)), levels = celltypes)
design <- model.matrix(~-1+file.celltypes) 


## ----results="hide"-----------------------------------------------------------
set.seed(1234)
b0 = getb0.biseq(methmat, design.rrbs, sigg=NULL)

## -----------------------------------------------------------------------------
resultx1 = getx1.biseq(NB0="top_bonferroni",b0,methmat,sample.id.rrbs,celltypes.rrbs)

## -----------------------------------------------------------------------------
resultx1 = getx1.biseq(NB0=250,b0,methmat,sample.id.rrbs,celltypes.rrbs)

## -----------------------------------------------------------------------------
x2 = as.data.frame(design.rrbs,row.names=sample.id.rrbs)
cr=getcorr(resultx1$x1,x2)
plot(cr, ylim=c(0,1), ylab="Correlation")

## -----------------------------------------------------------------------------
data("data_tissue_rrbs")

## -----------------------------------------------------------------------------
resultx1.tissue = getx1.biseq(NB0=250,b0,methmat.tissue,sample.id.tissue,celltypes.rrbs)

## ----warning=FALSE, message=FALSE---------------------------------------------
x1 = cbind(lymph=resultx1.tissue$x1[,1]+resultx1.tissue$x1[,2], mono = resultx1.tissue$x1[,3], gran = resultx1.tissue$x1[,4])
x2 = as.matrix(cbc.rrbs/100)
cr = getcorr(x1,x2)
plot(cr, ylim=c(0,1), ylab="Correlation")

