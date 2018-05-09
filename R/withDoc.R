##======================================================
##INPUT FILES AND PARAMETERS  for RNAseq
##======================================================
#' generate ordered sample identifiers for RNAseq data
#' @export
sample.idVec = function() paste0("SRR102",c(2945:2986,2990:3003))

#' generate vector of tissue type labels (ordered) for RNAseq data
#' @export
tissueVec = function()
   c(rep("Bcells",14),rep("Granulocytes",14),rep("Monocytes",14),rep("Tcells",14))

#' generate vector of cell types present
#' @export
celltypeVec = function() c("Bcells","Tcells","Monocytes","Granulocytes")

#' generate ordered sample identifiers for RNAseq data for whole blood
#' @export
sample.id.tissue.Vec = function() c(paste0("sample",1:8),"SRR1076511","SRR1313514","SRR1314113","SRR1314916","SRR1376619","SRR1433860","SRR1447523","SRR1453456","SRR1469614","SRR1476300","SRR1489926")

#' generate data frame with CBC
#' @export
cbcDF = function() data.frame(sample.id=c("sample1","sample5", "sample6","sample7","sample8"), lymph=c(21.4,24.9,0,9.3,6.1), mono=c(6.1,6.5,5,7.6,6.4), gran=c(72.3,68.6,94,83.1,87.5))






##======================================================
##INPUT FILES AND PARAMETERS  for RRBS
##======================================================

#' generate ordered sample identifiers for RRBS data
#' @export
sample.idVec2 = function() c("SRR1104848","SRR1104855","SRR1104856","SRR1104857", "SRR1104838","SRR1104839","SRR1104841","SRR1104842","SRR1104852", "SRR1104853","SRR1104854","SRR2960993","SRR1508407","SRR1508408","SRR1508409","SRR1508410", "SRR1508411","SRR1508412","SRR1508413","SRR1508414","SRR1508415", "SRR1508416","SRR1508417","SRR1508418")

#' generate vector of tissue type labels (ordered) for RRBS data
#' @export
tissueVec2 = function() c(rep("Monocytes",4),rep("Tcells",4), rep("Bcells",3), rep("Neutrophils",13))
   
#' generate vector of cell types present for RRBS data
#' @export
celltypeVec2 = function() c("Bcells","Tcells","Monocytes","Neutrophils")


#' generate ordered sample identifiers for RNAseq data for whole blood
#' @export
sample.id.tissue.Vec2 = function() c("sample1","sample5", "sample6","sample7","sample8")






#==================================================
#FUNCTIONS
#==================================================
#' obtain EdgeR DGE object from count matrix
#' @param countmat a matrix of counts with rows corresponding to genes and columns to samples
#' @param design output of model.matrix
#' @param ncpm.min filtering threshold for cpm
#' @param nsamp.min minimum number of samples that must have cpm exceeding ncpm.min for a gene to be retained
#' @return dge DGE object
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



#make one big methylation matrix
getmethmat <- function(cpgmethfil){
	nfil = length(cpgmethfil)
	methmat=c()
	
	for(i in 1:nfil){
		load(cpgmethfil[i])
		methmati = getData(methi)
		methmat = rbind(methmat, methmati)	
	}
	
	return(methmat)
}




#EdgeR analysis, obtain B0, the projection matrix
#' compute B0, projection matrix, given counts
#' @import statmod
#' @param countmat a matrix of counts with rows corresponding to genes and columns to samples
#' @param design output of model.matrix
#' @param ncpm.min filtering threshold for cpm
#' @param nsamp.min minimum number of samples that must have cpm exceeding ncpm.min for a gene to be retained
#' @return b0 projection matrix
#' @return pfstat f-stat p-values for each gene
#' @return dge DGE object
#' @export
getb0.rnaseq <- function(countmat, design, ncpm.min=1, nsamp.min=4){

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
#' #RNA sequencing example
#'
#' set.seed(1234)
#' data("celltypes_cnts")
#' okids = sample.idVec()
#' cnts = celltypes_cnts[, okids]
#' tissf = factor(tissueVec(), levels=celltypeVec())
#' des = model.matrix(~-1+tissf)
#' b0 = getb0.rnaseq(cnts, des, ncpm.min=1, nsamp.min=4)
#' nb0 = 50
#' dge = getdge(cnts, des, ncpm.min=1, nsamp.min=4)
#' resultx1 = getx1.rnaseq(nb0,b0,dge)
#' x2 = as.data.frame(des,row.names=okids)
#' getcorr(resultx1$x1,x2)
#'
#' #Validate projection matrix on whole tissue
#' data("wholeblood_cnts")
#' okids.tissue = sample.id.tissue.Vec()
#' cnts.tissue = wholeblood_cnts[,okids.tissue]
#' nb0=3000
#' dge.tissue = getdge(cnts.tissue,design=NULL, ncpm.min=1, nsamp.min=4)
#' resultx1 = getx1.rnaseq(nb0,b0, dge.tissue)
#' cbc = cbcDF()
#' x1 = cbind(lymph=resultx1$x1[as.character(cbc$sample.id),1]+resultx1$x1[as.character(cbc$sample.id),2], mono = resultx1$x1[as.character(cbc$sample.id),3], gran = resultx1$x1[as.character(cbc$sample.id),4])
#' x2=as.matrix(cbc[,2:4])/100
#' getcorr(x1,x2)

#'
#' #Bisulfite sequencing example
#'
#' set.seed(1234)
#' data("methmat")
#' data("methmat_geo_small")
#' methmat = methmat_geo
#' okids = sample.idVec2()
#' tissf = factor(tissueVec2(), levels=celltypeVec2())
#' des = model.matrix(~-1+tissf)
#' b0 = getb0.biseq(methmat, des)
#' nb0=1000
#' resultx1 = getx1.biseq(nb0,b0,methmat,okids,celltypeVec2())
#' x2 = as.data.frame(des,row.names=okids)
#' getcorr(resultx1$x1,x2)
#'
#' #Validate projection matrix on whole tissue
#' data("methmat_wholeblood_small")
#' methmat.tissue = methmat_wholeblood
#' okids.tissue = sample.id.tissue.Vec2()
#' nb0=1000
#' resultx1 = getx1.biseq(nb0,b0,methmat.tissue,okids.tissue,celltypeVec2())
#' cbc = cbcDF()
#' x1 = cbind(lymph=resultx1$x1[as.character(cbc$sample.id),1]+resultx1$x1[as.character(cbc$sample.id),2], mono = resultx1$x1[as.character(cbc$sample.id),3], gran = resultx1$x1[as.character(cbc$sample.id),4])
#' x2=as.matrix(cbc[,2:4])/100
#' getcorr(x1,x2)
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
#use this to get x1 mixture data for tissue using cell type projection matrix b0 (RNAseq)
#' compute mixture data given projection matrix (RNAseq)
#' @import nloptr
#' @param NB0 number of genes to be retained, ordered by f-stat p-value
#' @param resultb0 output of \code{\link{getb0.rnaseq}}
#' @param dge output of \code{\link{getdge}}
#' @param MAXITER integer number of iterations allowed
#' @return x1 cell mixture of sample
#' @return converged convergence
#' @export
getx1.rnaseq<-function(NB0, resultb0, dge, MAXITER=1000){
	
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



#' compute B0, projection matrix, given methylation counts
#' @param methmat a matrix of counts with rows corresponding to methylation sites and columns. Columns are: chr start   end strand coverage1 numCs1 numTs1 coverage2 numCs2 numTs2 ....
#' @param design output of model.matrix
#' @param celltypes vector of cell types 
#' @return b0 projection matrix. Coefficients are beta.
#' @export
getb0.biseq <- function(methmat, design, mixedef=FALSE){
	formula = paste0(c(colnames(design),0),collapse="+")
	b0 = my_diffMethFromDesign_matrix(methmat,design,formula,cur_treatment=NULL,mixedef,fstat=TRUE)	
	return(b0)		
}



#' compute mixture data given projection matrix (bisulfite sequencing)
#' @param NB0 number of genes to be retained, ordered by f-stat p-value
#' @param b0 output of \code{\link{getb0.biseq}}
#' @param methmat a matrix of counts with rows corresponding to methylation sites and columns. Columns are: chr start   end strand coverage1 numCs1 numTs1 coverage2 numCs2 numTs2 ....
#' @param sample.id vector of sample IDs
#' @param celltypes vector of cell types
#' @param MAXITER integer number of iterations allowed
#' @return x1 cell mixture of sample
#' @return converged convergence
#' @export
getx1.biseq <- function(NB0,b0,methmat,sample.id,celltypes,MAXITER=10000){
	b0r = b0[order(b0$pfstat)[1:NB0],]
	ntypes=length(grep("beta",colnames(b0)))
	nsamples=length(sample.id)

	#initialize mixture matrix
	x1 = matrix(0,nrow=nsamples,ncol=ntypes)
	one_vec = rep(1,ntypes)
	x0 = rep(1/ntypes, ntypes) #for initial mixture, use equal amounts
	gx <- function(x) {one_vec %*% x-1} #inequality constraint
	converged = c() #convergence if iter < maxiter
	miter = c()


	#get only the NB0 genes from methmat
	meth.small = merge(x=b0r,y=methmat,by=c("chr","start","end","strand"))
	meth.cn = colnames(meth.small)

	#get names of columns of the coefficients, ie beta
	type.coef = colnames(b0r)[grep("beta",colnames(b0r))]

	for(i in 1:nsamples){
		#look at the top NB0 genes only
	
		y1 = meth.small[,c(1:4,which(meth.cn %in% paste0("numCs",i)),which(meth.cn %in% paste0("numTs",i)))]
		colnames(y1)=c("chr","start","end","strand","Cs","Ts")
		
		#obj is negative log likelihood, x is the mixture vector
		b0r.beta = as.matrix(meth.small[,type.coef])
		obj1 <- function(x){-sum(dbinom(y1$Cs,size=(y1$Cs+y1$Ts), prob=1/(1+exp(-(b0r.beta %*% x))) ))}
		#use constraint that gx<=0 or sum_i x_i <=1
		#constrained lower bound to 0, upper bound to 1 for each x_i and sum_i of x_i <=1, leaving room for other cellular components that may not have been accounted for
		opt_result = nloptr(x0=x0, eval_f=obj1,lb=rep(0,ntypes), ub=rep(1,ntypes), eval_g_ineq=gx, opts=list(algorithm = "NLOPT_LN_COBYLA", maxeval=MAXITER))
	
		x1[i,]=opt_result$solution #predicted coefficients
		converged = c(converged,ifelse(opt_result$iterations < MAXITER,1,0))
		miter = c(miter,opt_result$iterations)
	}

	colnames(x1)=colnames(b0r.beta)
	rownames(x1)=sample.id
	x1=as.data.frame(x1)
	
	return(list(x1=x1,converged=converged))
}

#==========================












