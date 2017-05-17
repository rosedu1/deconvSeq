##deconvolution of RNAseq
#
#library(edgeR)
#library(statmod) #for glmnb.fit
#library(nloptr) #for nloptr, nonlinear optimization with constraints
#
#
#set.seed(1234)
#
##======================================================
##INPUT FILES AND PARAMETERS  
##======================================================
#
##individual cell types
#setwd("/proj/regeps/regep00/studies/aneurysms/analyses/rerdu/rnaseq_leucegene_grch38/")
#fildir=paste0("SRR102",2945:3003)
#filnames=sapply(1:length(fildir), function(x) paste0("/proj/regeps/regep00/studies/aneurysms/analyses/rerdu/rnaseq_leucegene_grch38/results/",fildir[x],"/htseq_out/",fildir[x],"_accepted_hits.sortedname.gene.counts"))
#

#' generate ordered sample identifiers
#' @export
sample.idVec = function() paste0("SRR102",c(2945:2986,2990:3003))

#' generate vector of tissue type labels (ordered)
#' @export
tissueVec = function()
   c(rep("Bcells",14),rep("Granulocytes",14),rep("Monocytes",14),rep("Tcells",14))

#' generate vector of cell types present
#' @export
celltypeVec = function() c("Bcells","Tcells","Monocytes","Granulocytes")


#whole blood, putting in more samples than just the 5 samples with cbc seems to help with the stability of the predictions
#setwd("/proj/regeps/regep00/studies/aneurysms/analyses/rerdu/rnaseq_grch38_blood_cell")
#fildir_s=c("RD1_0314R1D00004_A160712690","RD2_0512R1D00008_A160712691", "RD3_0611R1D00002_A160712692","RD4_0811R1D00011_A160712693","RD5_0812R1D00010_A160712694", "RD6_0912R1D00013_A160712695","RD7_0913R1D00010_A160712696","RD8_1012R1D00003_A160712697","SRR1076511","SRR1313514","SRR1314113","SRR1314916","SRR1376619","SRR1433860","SRR1447523","SRR1453456","SRR1469614","SRR1476300","SRR1489926")
#filnames_s=c()
#filnames_s = c(sapply(1:8, function(x) paste0("/proj/regeps/regep00/studies/aneurysms/analyses/rerdu/rnaseq_grch38_blood_cell/results/",fildir_s[x],"/htseq_out/",fildir_s[x],"_accepted_hits.sortedname.gene.counts")), sapply(9:19, function(x) paste0("/proj/regeps/regep00/studies/aneurysms/analyses/rerdu/rnaseq_gtex_blood_grch38/results/",fildir_s[x],"/htseq_out/",fildir_s[x],"_accepted_hits.sortedname.gene.counts")) )
#sample.id_s = fildir_s
#tissuef_s=factor(c(rep("aneu_blood",8), rep("control_blood",11)), levels=c("control_blood", "aneu_blood"))
#design_s = NULL

#' generate data frame with CBC
#' @export
cbcDF = function() data.frame(sample.id=c("RD1_0314R1D00004_A160712690","RD5_0812R1D00010_A160712694", "RD6_0912R1D00013_A160712695","RD7_0913R1D00010_A160712696","RD8_1012R1D00003_A160712697"), lymph=c(21.4,24.9,0,9.3,6.1), mono=c(6.1,6.5,5,7.6,6.4), gran=c(72.3,68.6,94,83.1,87.5))



#will keep genes that achieve ncpm count per million in at least nsamp samples (used ncpm=1 and nsamp=4)
#ncpm = 1
#nsamp = 4  




#==================================================
#INITIALIZE VARIABLES
#==================================================
#tissuei = which(tissue %in% celltypes)
#tissuef = factor(tissue[tissuei], levels=celltypes)
#design <- model.matrix(~-1+tissuef) 


#==================================================
#FUNCTIONS
#==================================================
#' obtain EdgeR DGE object from count matrix
#' @param countmat a matrix of counts with rows corresponding to genes and columns to samples
#' @param design output of model.matrix
#' @param ncpm.min filtering threshold for cpm
#' @param nsamp.min minimum number of samples that must have cpm exceeding ncpm.min for a gene to be retained
#' @return dge DGE object
#' @examples
#' data("leucegene_btmg_b0_counts")
#' okids = sample.idVec()
#' cnts = leucegene_btmg_b0_counts[, okids]
#' tissf = factor(tissueVec(), levels=celltypeVec())
#' des = model.matrix(~-1+tissf)
#' dge = getdge(cnts, des, ncpm.min=1, nsamp.min=4)
#' @export
getdge <- function(countmat, design, ncpm.min, nsamp.min){
	

	dge=DGEList(countmat,genes=rownames(countmat))

	#filter
	keep <- rowSums(cpm(dge)>ncpm.min)>=nsamp.min
	dge <- dge[keep,]

	#recompute library sizes
	dge$samples$lib.size <- colSums(dge$counts)

	#normalize. Default method is TMM
	dge <- calcNormFactors(dge)

	#estimate dispersion
	#first estimate common dispersion
	dge <- estimateGLMCommonDisp(dge, design, verbose=TRUE)

	#estimate gene specific dispersions
	dge <- estimateGLMTrendedDisp(dge,design)
	dge <- estimateGLMTagwiseDisp(dge,design)

	return(dge)	
}




#EdgeR analysis, obtain B0, the projection matrix
#' compute B0, projection matrix, given counts
#' @param countmat a matrix of counts with rows corresponding to genes and columns to samples
#' @param design output of model.matrix
#' @param ncpm.min filtering threshold for cpm
#' @param nsamp.min minimum number of samples that must have cpm exceeding ncpm.min for a gene to be retained
#' @return b0 projection matrix
#' @return pfstat f-stat p-values for each gene
#' @return dge DGE object
#' @examples
#' data("leucegene_btmg_b0_counts")
#' okids = sample.idVec()
#' cnts = leucegene_btmg_b0_counts[, okids]
#' tissf = factor(tissueVec(), levels=celltypeVec())
#' des = model.matrix(~-1+tissf)
#' b0 = getb0(cnts, des, ncpm.min=1, nsamp.min=4)
#' summary(b0)
#' @export
getb0 <- function(countmat, design, ncpm.min=1, nsamp.min=4){

	dge <- getdge(countmat, design, ncpm.min, nsamp.min)
	
	#offset = log(lib.size *norm.factors) using TMM algorithm from Robinson and Oshlack 2010.
	offset = getOffset(dge)
	#tagwise dispersion, ie genewise dispersion
	dispersion = dge$tagwise.dispersion

	#number of genes
	ngenes = nrow(dge$counts)

	#number of cell types
	ntypes = ncol(design)  
	nsamples = nrow(design)
	b0 = matrix(0,nrow=ngenes,ncol=ntypes) #initialize b0
	design0 <- matrix(1,nrow=nrow(design),1) #null design matrix
	fstat = c()
	pfstat = c()

	for(k in 1:ngenes){
		y0= dge$counts[k,]
		w0 = design
		#full model. Glmnb.fit from statmod package, uses Levenberg-Marquardt damping algorithm just like edgeR
		glm.full = glmnb.fit(w0, y0, dispersion=dispersion[k], offset=offset, start.method="mean", tol=1e-06, maxit=100, trace=FALSE)
		#null model 
		glm.0 = glmnb.fit(design0, y0, dispersion=dispersion[k], offset=offset, start.method="mean", tol=1e-06, maxit=100, trace=FALSE)

		#fstatistic = F_(p2-p1)_(N-p2)= [(D_null-D_full)/(# cell types -1)]/[D_full/(# samples - # cell types)]         
		fstatk = ((glm.0$deviance-glm.full$deviance)/(ntypes-1))/(glm.full$deviance/(nsamples-ntypes))
	
		#pvalue = 1-pf(fstat,df1,df2) = pf(fstat,df1,df2,lower.tail=FALSE)	
		pk = 1-pf(fstatk,ntypes-1,nsamples-ntypes)
	
		#coef are b0k
		b0k = glm.full$coefficients
	
		fstat = c(fstat,fstatk)
		pfstat = c(pfstat,pk)
		b0[k,] = b0k
	}

	colnames(b0)=names(glm.full$coefficients)
	rownames(b0) = dge$genes$genes

	#adjust for multiple testing
	pfstat.adjust = p.adjust(pfstat,method="BH")

	list(b0,pfstat.adjust,nsamples, dge, offset, dispersion)

	return(list(b0=b0,pfstat.adjust=pfstat.adjust, dge=dge))
}


#x1,x1b are mixture matrices where rows are samples, columns are cell types
#row names and column names for x1 and x2 should match
#' compute correlation between estimated and actual composition
#' @param x1 composition matrix where rows are samples and columns are cell types
#' @param x2 second composition matrix where rows and columns match x1
#' @return cor.type correlation by cell type
#' @return cor.sample correlation by sample
#' @examples
#' data("leucegene_btmg_b0_counts")
#' okids = sample.idVec()
#' cnts = leucegene_btmg_b0_counts[, okids]
#' tissf = factor(tissueVec(), levels=celltypeVec())
#' des = model.matrix(~-1+tissf)
#' b0 = getb0(cnts, des, ncpm.min=1, nsamp.min=4)
#' nb0 = 3000
#' dge = getdge(cnts, des, ncpm.min=1, nsamp.min=4)
#' resultx1 = getx1(nb0,b0,dge)
#' x2 = as.data.frame(des,row.names=okids)
#' result.corr = getcorr(resultx1$x1,x2)
#' @export
getcorr <- function(x1, x2){
	
	#correlation of each cell type over all samples between predicted and actual cell mixture
	cor.type = mapply(cor,x1[,1:ncol(x1)],x2[,1:ncol(x2)])
	names(cor.type)=colnames(x1)
	
	#correlation of each sample over all cell types between predicted and actual cell mixture
	cor.sample = mapply(cor,split(as.matrix(x1),row(x1)),split(as.matrix(x2),row(x2)))
	
	return(list(cor.type=cor.type, cor.sample=cor.sample))
}



#MAXITER is maximum iterations before calculation is stopped and considered nonconvergent
#use this to get x1 mixture data for tissue using cell type projection matrix b0
#' compute mixture data given projection matrix
#' @param NB0 number of genes to be retained, ordered by f-stat p-value
#' @param resultb0 output of \code{\link{getb0}}
#' @param dge output of \code{\link{getdge}}
#' @param MAXITER integer number of iterations allowed
#' @return x1 cell mixture of sample
#' @return converged convergence
#' @export
getx1<-function(NB0, resultb0, dge, MAXITER=1000){
	
	b0 = resultb0$b0
	pfstat = resultb0$pfstat.adjust

	b0r = b0[order(pfstat)[1:NB0],]
	ntypes = ncol(b0r)
	
	#for default, start with equal proportions of all cell types but can give initial x0 starting mixture
	x0=rep(1/ntypes, ntypes)
	
	offset_s = getOffset(dge)
	dispersion_s = dge$tagwise.dispersion #genewise dispersion
	nsamples_s = nrow(dge$samples)

	#initialize mixture matrix
	x1 = matrix(0,nrow=nsamples_s,ncol=ntypes)
	one_vec = rep(1,ntypes)
	gx <- function(x) {one_vec %*% x-1} #inequality constraint
	converged = c() #convergence if iter < maxiter

	for(i in 1:nsamples_s){
		#look at the top NB0 genes only
		y1 = dge$counts[which(rownames(dge$counts) %in% rownames(b0r)),i] 
		#obj is negative log likelihood, x is the mixture vector
		obj <- function(x){-sum(dnbinom(y1,size=dispersion_s, mu=exp(b0r %*% x + offset_s[i])))}
		#use constraint that gx<=0 or sum_i x_i <=1
		#constrained lower bound to 0, upper bound to 1 for each x_i and sum_i of x_i <=1, leaving room for other cellular components that may not have been accounted for
		opt_result = nloptr(x0=x0, eval_f=obj,lb=rep(0,ntypes), ub=rep(1,ntypes), eval_g_ineq=gx, opts=list(algorithm = "NLOPT_LN_COBYLA", maxeval=MAXITER))
	
		x1[i,]=opt_result$solution #predicted coefficients
		converged = c(converged,ifelse(opt_result$iterations < MAXITER,1,0))
	}

	colnames(x1)=colnames(b0r)
	rownames(x1)=rownames(dge$samples)
	x1=as.data.frame(x1)
	
	return(list(x1=x1,converged=converged))
}








notnow = function() {

#==================================================
#MAIN
#==================================================

#-------------------------------------------
#project individual cell types onto cell types to get mixture fractions using the top NB0 genes
#input data files from filnames are count data output from htseq
#output projection matrix is in resultb0: /proj/regeps/regep00/studies/aneurysms/analyses/rerdu/rnaseq_leucegene_grch38/b0_leucegene_btmg.rda
#-------------------------------------------

#=================================================
#get projection matrix b0
#=================================================
#given count data for individual cell types, get projection matrix b0 and pvalues (as calculated from F statistics) for that gene
resultb0 = getb0(filnames, sample.id, design, ncpm, nsamp)


#=================================================
#validate individual cell type by calculating cell mixture
#=================================================

#nb0 is the top nb0 number of genes to be used in projection matrix for individual cell types
nb0list = seq(25,500,25)

#calculate the cell type mixture x1 for each sample given projection matrix in resultb0. The samples are individual cell types so we can validate the prediction with the actual cell types using the design matrix
#result's components are x1, converged
result=c()
result = sapply(nb0list, function(x) getx1(x, resultb0, dge_s=resultb0$dge, MAXITER=1000))



#get correlations
#--------------------
cor.result=c()
cor.result = sapply(1:length(nb0list), function(i) getcorr(result[1,i]$x1,as.data.frame(design,row.names=sample.id)))

#correlation for each predicted cell type with actual cell type over all samples
cortype = cbind(nb0list,matrix(unlist(cor.result[1,]),ncol=length(celltypes),byrow=TRUE))
colnames(cortype)=c("nb0",celltypes)

#correlation of the predicted cell mixture for each sample with the actual cell mixture
corsample = cbind(nb0list,matrix(unlist(cor.result[2,]), ncol=length(sample.id), byrow=TRUE))
colnames(corsample)=c("nb0",sample.id)


#plots
#--------------------
par(oma=c(2,3,1,1))
matplot(cortype[,"nb0"], cortype[,c(2:5)],type="b", xlab="Number of genes", ylab="Correlation",col=c("green","blue","red","purple"), pch=19, cex=2, cex.axis=2, cex.lab=2)
legend("bottomright",c(celltypes), col=c("green","blue","red","purple"), cex=2, pch=19)
#dev.copy2pdf(file="/proj/regeps/regep00/studies/aneurysms/analyses/rerdu/rnaseq_leucegene_grch38/cortype_leucegene_btmg_nb025to500_ESGNonly.pdf")

par(oma=c(2,3,1,1))
matplot(corsample[,"nb0"], corsample[,c(2:ncol(corsample))],type="l", xlab="Number of genes", ylab="Correlation")
#dev.copy2pdf(file="/proj/regeps/regep00/studies/aneurysms/analyses/rerdu/rnaseq_leucegene_grch38/corsample_leucegene_btmg_nb025to500_ESGNonly.pdf")



#=================================================
#get cell mixture for whole tissue and validate
#=================================================

#-------------------------------------------
#project whole tissue onto cell types, given projection matrix
#input resultb0 from above section: /proj/regeps/regep00/studies/aneurysms/analyses/rerdu/rnaseq_leucegene_grch38/decon_b0_rnaseq.rda
#input tissue datafiles in filnames_s
#output /proj/regeps/regep00/studies/aneurysms/analyses/rerdu/rnaseq_leucegene_grch38/decon_b0_rnaseq_wholeblood.rda
#-------------------------------------------


#get dge data type from count data for whole tissue
dge_s = getdge(filnames_s, sample.id_s, design_s, ncpm, nsamp)

#number of genes to be used in projection matrix for whole tissue
nb0list_s=c(3000)

#calculate cell type mixture x1 for each sample
#result.tissue's components are x1 and converged.
result.tissue=c()
result.tissue = sapply(nb0list_s, function(x) getx1(x, resultb0, dge_s, MAXITER=1000))



#calculate correlations of cell mixture prediction with known differential cell counts
#---------------------------------------
cor.sample.tissue=c()
for(j in 1:length(nb0list_s)){
	mix = result.tissue[1,j]$x1
	mix2 = cbind(mix[,"tissuefBcells"]+mix[,"tissuefTcells"], mix[,"tissuefMonocytes"], mix[,"tissuefGranulocytes"])
	colnames(mix2)=c("lymph","mono","gran")
	rownames(mix2)=rownames(mix)
	cor.sample.tissue = rbind(cor.sample.tissue,sapply(1:nrow(cbc), function(k) { cor(as.numeric(mix2[as.character(cbc$sample.id[k]),]),as.numeric(cbc[k,2:4]) ) }))
}
cor.sample.tissue = cbind(nb0list_s, cor.sample.tissue)
colnames(cor.sample.tissue)=c("nb0",as.character(cbc$sample.id))


cor.sample.tissue


#nb0list_s=1000, no non-ENSG
#cor.sample.tissue
#      nb0 RD1_0314R1D00004_A160712690 RD5_0812R1D00010_A160712694
#[1,] 1000                   0.9753363                   0.9647361
#     RD6_0912R1D00013_A160712695 RD7_0913R1D00010_A160712696
#[1,]                   0.9988821                   0.9998056
#     RD8_1012R1D00003_A160712697
#[1,]                  0.08661279


#nb0list_s=2000m no non-ENSG
#cor.sample.tissue
#      nb0 RD1_0314R1D00004_A160712690 RD5_0812R1D00010_A160712694
#[1,] 2000                   0.9753363                   0.9999972
#     RD6_0912R1D00013_A160712695 RD7_0913R1D00010_A160712696
#[1,]                   0.9757208                   0.9998056
#     RD8_1012R1D00003_A160712697
#[1,]                   0.6849306


#nb0list_s=3000 no non-ENSG
#cor.sample.tissue
#      nb0 RD1_0314R1D00004_A160712690 RD5_0812R1D00010_A160712694
#[1,] 3000                   0.9753363                   0.9687035
#     RD6_0912R1D00013_A160712695 RD7_0913R1D00010_A160712696
#[1,]                   0.9988821                   0.9998056
#     RD8_1012R1D00003_A160712697
#[1,]                   0.9909303



}





