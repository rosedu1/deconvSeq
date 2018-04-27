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
getdge <- function(countmat, design, ncpm.min=1, nsamp.min=4){
	

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
#' compute b0, projection matrix, given counts
#' @param dge DGE object
#' @param design output of model.matrix
#' @param ncpm.min filtering threshold for cpm
#' @param nsamp.min minimum number of samples that must have cpm exceeding ncpm.min for a gene to be retained
#' @return b0 projection matrix
#' @return p0 pvalues matrix
#' @return dge DGE object
#' @return converged convergence status. 1=converged.
#' @return pfstat F-stat p-values for each gene
#' @export
getb0.rnaseq <- function(dge, design, ncpm.min=1, nsamp.min=4){

	#dge <- getdge(countmat, design, ncpm.min, nsamp.min)
	
	#offset = log(lib.size *norm.factors) using TMM algorithm from Robinson and Oshlack 2010.
	offset = getOffset(dge)
	offset = offset-mean(offset)
	#tagwise dispersion, ie genewise dispersion
	dispersion = dge$tagwise.dispersion

	#number of genes
	ngenes = nrow(dge$counts)

	#number of cell types
	ntypes = ncol(design)  
	nsamples = nrow(design)
	b0 = matrix(0,nrow=ngenes,ncol=ntypes) #initialize b0
	design0 <- matrix(1,nrow=nrow(design),1) #null design matrix
	fstat = rep(NA,ngenes)
	pfstat = rep(NA, ngenes)
	names(pfstat) = rownames(dge$counts)
	p0=matrix(0,nrow=ngenes,ncol=ntypes) #initialize p0, matrix of pvalues for each cell type
	converged = c()

	for(k in 1:ngenes){
		y0= dge$counts[k,]
		w0 = design
		p0k=c()

		#use glm.nb from statmod 
		thisdat = as.data.frame(cbind(y0,w0))
		thisformula = as.formula(paste("y0 ~ ", paste(colnames(w0), collapse="+"),"+offset(offset)+0"))

		#full model
		glm.full=tryCatch(glm.nb(thisformula,data=thisdat, link="log", init.theta=1/dispersion[k]),warning=function(w) {NA}, error=function(e) {NA} )
		
		#null model
		formula0 = as.formula("y0 ~ 1")
		glm.0=tryCatch(glm.nb(formula0,data=thisdat, link="log", init.theta=1/dispersion[k]),warning=function(w) {NA}, error=function(e) {NA} )

		if(!is.na(unlist(glm.full)[[1]])){
			test = summary(glm.full)
			p0k = test$coefficients[,4]
			b0k = test$coefficients[,1]
			converged = c(converged, as.numeric(glm.full$converged)) #1 means converged for glm.nb
	
			b0[k,] = b0k
			p0[k,] = p0k
			
			if(!is.na(unlist(glm.0)[[1]])){
			
				#fstatistic = F_(p2-p1)_(N-p2)= [(D_null-D_full)/(# cell types -1)]/[D_full/(# samples - # cell types)]         
				fstatk = ((glm.0$deviance-glm.full$deviance)/(ntypes-1))/(glm.full$deviance/(nsamples-ntypes))
	
				#pvalue = 1-pf(fstat,df1,df2) = pf(fstat,df1,df2,lower.tail=FALSE)	
				pk = 1-pf(fstatk,ntypes-1,nsamples-ntypes)
						
				fstat[k] = fstatk
				pfstat[k] = pk

			} else {
				fstat[k]=NA
				pfstat[k]=NA
			}
			
			
		} else {
			converged = c(converged, NA) 
			b0[k,] = NA
			p0[k,] = NA

		}
	}

	colnames(b0)=names(glm.full$coefficients)
	rownames(b0) = dge$genes$genes
	
	colnames(p0)=names(glm.full$coefficients)
	rownames(p0) = dge$genes$genes
	

	return(list(b0=b0,p0=p0, dge=dge, converged=converged,pfstat=pfstat))
}





#MAXITER is maximum iterations before calculation is stopped and considered nonconvergent
#use this to get x1 mixture data for tissue using cell type projection matrix b0 (RNAseq)
#' compute mixture data given projection matrix (RNAseq)
#' @param NB0 number of genes to be retained, ordered by f-stat p-value
#' @param resultb0 output of \code{\link{getb0.rnaseq}}
#' @param dge output of \code{\link{getdge}}
#' @param MAXITER integer number of iterations allowed
#' @param x0 initial cell type composition for fitting
#' @return x1 cell mixture of sample
#' @return converged convergence. 1=converged
#' @export
getx1.rnaseq<-function(NB0, resultb0, dge_s, MAXITER=1000, x0=NULL){
	
	#keep only the genes that are also in dge_s
	b0r.1 = resultb0$b0
	pfstat.1 = resultb0$pfstat
	
	b0r = b0r.1[which(rownames(b0r.1) %in% dge_s$genes[,1]),]
	pfstat = pfstat.1[which(rownames(b0r.1) %in% dge_s$genes[,1])]
	b0r = b0r[order(pfstat)[1:NB0],]
	ntypes = ncol(b0r)
	
	#for default, start with equal proportions of all cell types but can give initial x0 starting mixture
	if(is.null(x0)) x0=rep(1/ntypes, ntypes)
	
	offset_s = getOffset(dge_s) #length of nsamples
	offset_s = offset_s-mean(offset_s)

	dispersion_s.index=c()
	for(i in 1:nrow(b0r)){
		dispersion_s.index=c(dispersion_s.index, which(dge_s$genes[,1] %in% rownames(b0r)[i]))
	}
	dispersion_s = dge_s$tagwise.dispersion[dispersion_s.index] #genewise dispersion
	
	nsamples_s = nrow(dge_s$samples)

	#initialize mixture matrix
	x1 = matrix(0,nrow=nsamples_s,ncol=ntypes)
	one_vec = rep(1,ntypes)
	gx <- function(x) {one_vec %*% x-1} #inequality constraint
	converged = c() #convergence if iter < maxiter

	for(i in 1:nsamples_s){
		#look at the top NB0 genes only
		#y1 = dge_s$counts[which(rownames(dge_s$counts) %in% rownames(b0r)),i] 
		y1 = dge_s$counts[rownames(b0r),i]
		#obj is negative log likelihood, x is the mixture vector
		offset_s.i = offset_s[i]
		obj <- function(x){-sum(dnbinom(y1,size=dispersion_s, mu=exp(b0r %*% x + offset_s.i), log=TRUE))}
		opt_result=tryCatch(solnp(pars=x0,fun=obj,eqfun=gx,eqB=0,LB=rep(0,ntypes),UB=rep(1,ntypes), control=list(outer.iter=MAXITER, inner.iter=MAXITER, trace=0)),warning=function(w) {NA}, error=function(e) {NA} )

		if(!is.na(unlist(opt_result)[[1]])){	
			x1[i,]=opt_result$pars
			converged = c(converged, opt_result$convergence) #for solnp converged = 0 means converged
		} else {
			x1[i,]=rep(NA,ntypes)
			converged = c(converged, NA)
		}
	}

	colnames(x1)=colnames(b0r)
	rownames(x1)=rownames(dge_s$samples)
	x1=as.data.frame(x1)
	
	return(list(x1=x1,converged=converged))
}



#' compute b0, projection matrix, given methylation counts
#' @param methmat a matrix of counts with rows corresponding to methylation sites and columns. Columns are: chr start   end strand coverage1 numCs1 numTs1 coverage2 numCs2 numTs2 ....
#' @param design design matrix, output of model.matrix
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
#' @param x0 initial cell type composition for fitting
#' @return x1 cell mixture of sample
#' @return converged convergence
#' @export
getx1.biseq <- function(NB0,b0,methmat,sample.id,celltypes,MAXITER=10000,x0=NULL){
	b0r = b0[order(b0$pfstat)[1:NB0],]
	ntypes=length(grep("beta",colnames(b0)))
	nsamples=length(sample.id)

	#initialize mixture matrix
	x1 = matrix(0,nrow=nsamples,ncol=ntypes)
	one_vec = rep(1,ntypes)
	if(is.null(x0)) x0=rep(1/ntypes, ntypes) #for initial mixture, use equal amounts
	gx <- function(x) {one_vec %*% x-1} #equality constraint
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
		obj1 <- function(x){-sum(dbinom(y1$Cs,size=(y1$Cs+y1$Ts), prob=1/(1+exp(-(b0r.beta %*% x))),log=TRUE ))}
		#use constraint that gx=0 or sum_i x_i =1
		#constrained lower bound to 0, upper bound to 1 for each x_i and sum_i of x_i =1		
		opt_result=solnp(pars=x0,fun=obj1,eqfun=gx,eqB=0,LB=rep(0,ntypes),UB=rep(1,ntypes), control=list(outer.iter=MAXITER, inner.iter=MAXITER, trace=0))

		x1[i,]=opt_result$pars
		converged = c(converged, opt_result$convergence) #for solnp converged = 0 means converged
	
	}

	colnames(x1)=colnames(b0r.beta)
	rownames(x1)=sample.id
	x1=as.data.frame(x1)
	
	return(list(x1=x1,converged=converged))
}





#x1,x1b are mixture matrices where rows are samples, columns are cell types
#row names and column names for x1 and x2 should match
#' compute correlation between estimated and actual composition
#' @param x1 composition matrix where rows are samples and columns are cell types
#' @param x2 second composition matrix where rows and columns match x1
#' @return correlation by sample
#' @examples
#' #RNA sequencing example
#'
#' set.seed(1234)
#' data("data_celltypes_rnaseq") #contains DGE object dge.celltypes, count matrix cnts.celltypes, and design matrix des for cell types design.rnaseq, sample.id.rnaseq
#' dge.celltypes = getdge(cnts.celltypes, design.rnaseq, ncpm.min=1, nsamp.min=4)
#' b0 = getb0.rnaseq(dge.celltypes, design.rnaseq, ncpm.min=1, nsamp.min=4)
#' nb0 = 50
#' resultx1 = getx1.rnaseq(nb0,b0,dge.celltypes)
#' x2 = as.data.frame(design.rnaseq,row.names=sample.id.rnaseq)
#' getcorr(resultx1$x1,x2)
#'
#' #Validate projection matrix on whole tissue
#' data("data_tissue_rnaseq") #contains DGE object dge.tissue, count matrix cnts.tissue, and CBC cbc.rnaseq
#' nb0=50
#' dge.tissue = getdge(cnts.tissue,design=NULL, ncpm.min=1, nsamp.min=4)
#' resultx1.tissue = getx1.rnaseq(nb0,b0, dge.tissue)
#' x1 = cbind(lymph=resultx1.tissue$x1$tissuefBcell+resultx1.tissue$x1$tissuefTcell, mono = resultx1.tissue$x1$tissuefMonocytes, gran = resultx1.tissue$x1$tissuefGranulocytes)
#' x2 = as.matrix(cbc.rnaseq/100)
#' getcorr(x1,x2)

#'
#' #Bisulfite sequencing example
#'
#' set.seed(1234)
#' data("data_celltypes_rrbs") #contains methylation matrix methmat, design.rrbs, celltypes.rrbs, sample.id.rrbs
#' b0 = getb0.biseq(methmat, design.rrbs)
#' nb0=250
#' resultx1 = getx1.biseq(nb0,b0,methmat,sample.id.rrbs,celltypes.rrbs)
#' x2 = as.data.frame(design.rrbs,row.names=sample.id.rrbs)
#' getcorr(resultx1$x1,x2)
#'
#' #Validate projection matrix on whole tissue
#' data("data_tissue_rrbs") #contains methylation matrix methmat.tissue, sample.id.tissue, actual differential cell count cbc.rrbs
#' nb0=250
#' resultx1.tissue = getx1.biseq(nb0,b0,methmat.tissue,sample.id.tissue,celltypes.rrbs)
#' x1 = cbind(lymph=resultx1.tissue$x1[,1]+resultx1.tissue$x1[,2], mono = resultx1.tissue$x1[,3], gran = resultx1.tissue$x1[,4])
#' x2 = as.matrix(cbc.rrbs/100)
#' getcorr(x1,x2)
#' @export
getcorr <- function(x1, x2){
	
	#correlation of each sample over all cell types between predicted and actual cell mixture
	cor.sample = mapply(cor,split(as.matrix(x1),row(x1)),split(as.matrix(x2),row(x2)))
	
	return(cor.sample)
}








