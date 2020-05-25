#' @rawNamespace import(SingleCellExperiment, except = c(cpm,weights))
#' @rawNamespace import(SummarizedExperiment, except = c(merge,as.data.frame,start,end))
#' @rawNamespace import(scater, except = c(filter,plotMDS))
#'@import edgeR stats methylKit scran 
#'@importFrom MASS glm.nb
#'@importFrom lme4 glmer
#'@importFrom Rsolnp solnp
#'@importFrom Rmisc CI
#'@importFrom psych fisherz paired.r fisherz2r
#'@importFrom biomaRt getBM useMart useDataset
#'@importFrom utils read.delim
NULL


#====================
# EXAMPLE DATA
#====================


#' This is the example cell type data for RNA sequencing
#'
#' \itemize{
#'	\item dge.celltypes DGE object from \(EdgeR\) of RNAseq results
#'	\item cnts.celltypes count matrix
#'	\item design.rnaseq design matrix 
#'	\item sample.id.rnaseq sample IDs
#' }
#' @name data_celltypes_rnaseq
#' @docType data
#' @keywords datasets
#' @usage data("data_celltypes_rnaseq")
NULL

#' This is the example tissue data for RNA sequencing
#'
#' \itemize{
#'	\item dge.tissue DGE object (from EdgeR) of RNAseq results 
#'	\item cnts.tissue count matrix
#'	\item cbc.rnaseq known cell type mixture for tissue
#' }
#' @name data_tissue_rnaseq
#' @docType data
#' @keywords datasets
#' @usage data("data_tissue_rnaseq")
NULL

#' This is the example data for scRNA sequencing
#'
#' \itemize{
#'	\item cnts.scrnaseq count matrix for full scRNAseq dataset
#' }
#' @name data_scrnaseq
#' @docType data
#' @keywords datasets
#' @usage data("data_scrnaseq")
NULL


#' This is the example cell type data for bisulfite sequencing
#'
#' \itemize{
#'	\item methmat methylation matrix 
#'	\item celltypes.rrbs vector of kinds of cell types
#'	\item design.rrbs design matrix 
#'	\item sample.id.rrbs sample IDs 
#' }
#' @name data_celltypes_rrbs
#' @docType data
#' @keywords datasets
#' @usage data("data_celltypes_rrbs")
NULL

#' This is the example tissue data for bisulfite sequencing
#'
#' \itemize{
#'	\item methmat.tissue methylation matrix 
#'	\item sample.id.tissue sample IDs 
#'	\item cbc.rrbs known cell type mixture for tissue
#' }
#' @name data_tissue_rrbs
#' @docType data
#' @keywords datasets
#' @usage data("data_tissue_rrbs")
NULL

#' DGE object (from EdgeR) of RNAseq results for cell types
"dge.celltypes"

#' Count matrix for cell types (RNAseq)
"cnts.celltypes"

#' Design matrix for cell types (RNAseq)
"design.rnaseq"

#' Sample IDs for cell types (RNAseq)
"sample.id.rnaseq"

#' DGE object (from EdgeR) of RNAseq results for tissue
"dge.tissue"

#' Count matrix for tissue (RNAseq)
"cnts.tissue"

#' Known cell type mixture for tissue (RNAseq)
"cbc.rnaseq"

#' Count matrix for single cell RNAseq (scRNAseq)
"cnts.scrnaseq"

#' methylation matrix for cell types (bisulfite sequencing)
"methmat"

#' vector of kinds of cell types (bisulfite sequencing)
"celltypes.rrbs"

#' Design matrix for cell types (bisulfite sequencing)
"design.rrbs"

#' Sample IDs for cell types (bisulfite sequencing)
"sample.id.rrbs"

#' methylation matrix for tissue (bisulfite sequencing)
"methmat.tissue"

#' Sample IDs for tissue (bisulfite sequencing)
"sample.id.tissue"

#' Known cell type mixture for tissue (bisulfite sequencing)
"cbc.rrbs"

#' Count matrix for single cell RNAseq
"cnts.scrnaseq"


#==================================================
#FUNCTIONS
#==================================================

#' Combine count matrix for individual samples
#'
#' Takes count matrix for individual samples (e.g. output of HTSeq output). Format: first column is Ensembl gene id, second column is count, no header. Removes genes with all 0 reads and rows that are not Ensembl gene IDs. Output is a count matrix where rows are genes and columns are samples.
#' Sample input files:  sample_genecounts1.txt, sample_genecounts2.txt
#' @param filnames filenames of individual samples
#' @param sample.id sample IDs 
#' @return count matrix with genes in rows. First column is gene name and 
#' @examples
#' file1 = system.file("extdata","sample1_genecounts.txt", package="deconvSeq")
#' file2 = system.file("extdata","sample2_genecounts.txt", package="deconvSeq")
#' countmat = getrnamat(filnames=c(file1,file2),sample.id=c("sample1","sample2"))
#' @export
getrnamat <- function(filnames, sample.id){
	x_t=list()
	for(i in 1:length(filnames)){
		x_t[[i]] = read.delim(filnames[i], header=FALSE)
			colnames(x_t[[i]])=c("gene",sample.id[i])
	}

	x2_t=x_t[[1]]
	for(i in 2:length(filnames)){
		x2_t = merge(x2_t,x_t[[i]],by.x="gene",by.y="gene")	
	}
	rownames(x2_t) = x2_t$gene
		
	#delete rows that are all 0. First column is gene name
	z=rowSums(x2_t[2:ncol(x2_t)])
	z2 = which(z==0) 
	if(length(z2!=0)) x2_t = x2_t[-z2,] 
	
	#delete rows that are not ENSG*
	y = which(substr(rownames(x2_t),1,1)=="_")
	if(length(y!=0)) x2_t = x2_t[-y,]

	x4_t = x2_t[,-1]

	return(x4_t)	
}





#' Obtain EdgeR DGE object from count matrix
#' @param countmat a matrix of counts with rows corresponding to 
#'   genes and columns to samples
#' @param design output of model.matrix
#' @param ncpm.min filtering threshold for cpm
#' @param nsamp.min minimum number of samples that must have cpm 
#' @param method method used to esimate trended dispersion in EdgeR: "auto" (default), "bin.spline", "bin.loess", "power", "spline"
#'   exceeding ncpm.min for a gene to be retained
#' @return dge DGE object
#' @export
getdge <- function(countmat, design, ncpm.min=1, nsamp.min=4, method="auto"){
	

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
	dge <- estimateGLMTrendedDisp(dge,design, method=method)
	dge <- estimateGLMTagwiseDisp(dge,design)

	return(dge)	
}




#EdgeR analysis, obtain B0, the projection matrix
#' compute b0, projection matrix, given counts
#' @param dge DGE object
#' @param design output of model.matrix
#' @param ncpm.min filtering threshold for cpm
#' @param nsamp.min minimum number of samples that must have cpm 
#'   exceeding ncpm.min for a gene to be retained
#' @param sigg vector of predetermined significant genes
#' @return b0 projection matrix
#' @return p0 pvalues matrix
#' @return dge DGE object
#' @return converged convergence status. 1=converged.
#' @return pfstat F-stat p-values for each gene
#' @export
getb0.rnaseq <- function(dge, design, ncpm.min=1, nsamp.min=4, sigg=NULL){

	#dge <- getdge(countmat, design, ncpm.min, nsamp.min)
	
	#offset = log(lib.size *norm.factors) using TMM algorithm from Robinson and Oshlack 2010.
	offset = getOffset(dge)
	offset = offset-mean(offset)

	if(!is.null(sigg)){
		dge$counts = dge$counts[sigg,]
		dge$genes = dge$genes[sigg,]
		dge$trended.dispersion = dge$trended.dispersion[match(sigg,rownames(dge$counts))]
		dge$tagwise.dispersion = dge$tagwise.dispersion[match(sigg,rownames(dge$counts))]
	}

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

	colnb0 =0 
	
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

			if(colnb0[[1]]==0) colnb0=names(glm.full$coefficients)
			
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

	#colnames(b0)=names(glm.full$coefficients)
	colnames(b0)=colnb0
	rownames(b0) = rownames(dge$counts)
	
	colnames(p0)=colnb0
	rownames(p0) =  rownames(dge$counts)
	
	
	return(list(b0=b0,p0=p0, dge=dge, converged=converged,pfstat=pfstat))
}





#MAXITER is maximum iterations before calculation is stopped and considered nonconvergent
#use this to get x1 mixture data for tissue using cell type projection matrix b0 (RNAseq)
#' compute mixture data given projection matrix (RNAseq)
#' @param NB0 number of genes to be retained, ordered by f-stat p-value. Other options: "all" uses all genes, "top_bonferroni" uses genes with adjusted p-values <0.05 after bonferroni correction, "top_fdr" uses genes with adjusted p-values <0.05 after FDR correction. Default is "top-bonferroni"
#' @param resultb0 output of \code{\link{getb0.rnaseq}}
#' @param dge_s output of \code{\link{getdge}}
#' @param MAXITER integer number of iterations allowed
#' @param x0 initial cell type composition for fitting
#' @return x1 cell mixture of sample
#' @return converged convergence. 1=converged
#' @export
getx1.rnaseq<-function(NB0="top_bonferroni", resultb0, dge_s, MAXITER=1000, x0=NULL){
	
	if(NB0=="all") {
		NB0 = length(resultb0$pfstat)
	} else if (NB0=="top_bonferroni"){
		NB0 = length(which(sort(p.adjust(resultb0$pfstat,"bonferroni"))<0.05))
	} else if (NB0=="top_fdr"){
		NB0 = length(which(sort(p.adjust(resultb0$pfstat,"fdr"))<0.05))		
	}
	
	#keep only the genes that are also in dge_s
	b0r.1 = resultb0$b0
	pfstat.1 = resultb0$pfstat
	
	b0r = b0r.1[which(rownames(b0r.1) %in% dge_s$genes[,1]),]
	pfstat = pfstat.1[which(rownames(b0r.1) %in% dge_s$genes[,1])]
	b0r = b0r[order(pfstat)[1:NB0],]
	ntypes = ncol(b0r)
	
	#for default, start with equal proportions of all cell types 
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




#' Makes methylation count matrix 
#'
#' Makes methylation count matrix using output from either BSMAP or Bismark. 
#'
#' Removes sites with 0 counts in all samples.
#'
#' Input file from BSMAP has columns: chromosome, position, strand, context, ratio, eff_CT_count, C_count, CT_count, rev_G_count, rev_GA_count, CI_lower, CI_upper.
#'
#' Input file from Bismark coverage for CG only has columns "chrBase","chr","base","strand","coverage","freqC","freqT"
#'
#' Sample input files from BSMAP: sample1_methratio.txt, sample1_mehtratio.txt
#'
#' For BSMAP input files, extracts CpGs and output in filname.CpG in the current working directory.
#' Adds "Chr" to chromosome name and output in filname.CpG_chr.
#' @param filnames input filenames
#' @param sample.id sample IDs
#' @param filtype "bsmap", "bismark". Default is "bsmap"
#' @return methylation count matrix where rows are CpG sites and columns are: chromosome, start, end, strand, number of Ts+Cs for sample 1, number of Cs for sample 1, number of Ts for sample 1, ....
#' @examples
#'  #file1 = system.file("extdata","sample1_methratio.txt", package="deconvSeq")
#'  #file2 = system.file("extdata","sample2_methratio.txt", package="deconvSeq")
#'  #methmat = getmethmat(filnames=c(file1,file2), sample.id=c("sample1","sample2"), filtype="bsmap")
#' @export
getmethmat <- function(filnames, sample.id, filtype="bsmap"){
	
	if(filtype=="bsmap"){
		WD = getwd()
		for(i in 1:length(filnames)){
			system(paste0("awk '($4~/^CG/)' ",filnames[i]," > ",WD,"/", basename(filnames[i]),".CpG"))
		}
	
		awkfil = system.file("extdata","convert_chr_name_grch38.awk", package="deconvSeq")
		nfil=c()
		for(i in 1:length(filnames)){
			system(paste0("awk -f ", awkfil,"  ",basename(filnames[i]),".CpG"," > ", WD,"/" ,basename(filnames[i]),".CpG_chr"))
			nfil=c(nfil,paste0(basename(filnames[i]),".CpG_chr"))
		}

		#parameters for assembly and context are for methRead and do not affect final count matrix
		obj=methRead(as.list(nfil), sample.id=as.list(sample.id), assembly="grch38",header=FALSE, context="CpG",resolution="base", pipeline=list(fraction=TRUE,chr.col=1,start.col=2,end.col=2,coverage.col=6,strand.col=3, freqC.col=5), treatment=rep(0,length(filnames)))
		
	} else if (filtype=="bismark"){
		
		obj=methRead(as.list(filnames), sample.id=as.list(sample.id), assembly="grch38",header=TRUE, context="CpG",resolution="base", treatment=rep(0,length(filnames)),pipeline="bismarkCoverage")

		
	} else {
		
		cat("Wrong file type. Filetypes are 'bsmap' or 'bismark'\n")
		return()
	}

	meth=unite(obj,min.per.group=NULL)
	meth.data = getData(meth)
	
	#remove rows in which all Ts are 0 or all Cs are 0 
	#meth.filter = meth.data[-which(rowSums(meth.data[,meth@numTs.index])==0 | rowSums(meth.data[,meth@numCs.index])==0),]
	if(length(which(rowSums(meth.data[,meth@numCs.index])==0))==0 & length(which(rowSums(meth.data[,meth@numTs.index])==0))){
		meth.filter = meth.data[-which(rowSums(meth.data[,meth@numTs.index])==0 | rowSums(meth.data[,meth@numCs.index])==0),]
	} else {
		meth.filter = meth.data
	}

	if(filtype=="bsmap") cat("files *.CpG and *.CpG_chr are written to current working directory\n")
	
	return(meth.filter)
}



#' compute b0, projection matrix, given methylation counts
#' @param methmat a matrix of counts with rows corresponding to methylation sites and columns. Columns are: chr start   end strand coverage1 numCs1 numTs1 coverage2 numCs2 numTs2 ....
#' @param design design matrix, output of model.matrix
#' @param sigg predetermined signature CpG sites. Format sites as chromosome 
#'    name, chromosome location, strand: chrN_position(+/-). For example, chr1_906825-
#' @return b0 projection matrix. Coefficients are beta.
#' @export
getb0.biseq <- function(methmat, design, sigg=NULL){
	if (!is.null(sigg)) {
		methmat = methmat[sigg,]
	}
	
	formula = paste0(c(colnames(design),0),collapse="+")
	b0 = my_diffMethFromDesign_matrix(methmat,design,formula)	
	return(b0)		
}



#' compute mixture data given projection matrix (bisulfite sequencing)
#' @param NB0 number of genes to be retained, ordered by f-stat p-value. Other options: "all" uses all genes, "top_bonferroni" uses genes with adjusted p-values <0.05 after bonferroni correction, "top_fdr" uses genes with adjusted p-values <0.05 after FDR correction. Default is "top-bonferroni"
#' @param b0 output of \code{\link{getb0.biseq}}
#' @param methmat a matrix of counts with rows corresponding to methylation 
#'   sites and columns. Columns are: chr start   end strand coverage1 numCs1 
#'   numTs1 coverage2 numCs2 numTs2 ....
#' @param sample.id vector of sample IDs
#' @param celltypes vector of cell types
#' @param MAXITER integer number of iterations allowed
#' @param x0 initial cell type composition for fitting
#' @return x1 cell mixture of sample
#' @return converged convergence
#' @export
getx1.biseq <- function(NB0="top_bonferroni",b0,methmat,sample.id,celltypes,MAXITER=10000,x0=NULL){
	
	if(NB0=="all") {
		NB0 = length(b0$pfstat)
	} else if (NB0=="top_bonferroni"){
		NB0 = length(which(sort(p.adjust(b0$pfstat,"bonferroni"))<0.05))
	} else if (NB0=="top_fdr"){
		NB0 = length(which(sort(p.adjust(b0$pfstat,"fdr"))<0.05))		
	}

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






#' Quality control of scRNAseq
#'
#' do quality control of scRNAseq based on library size, feature counts, chrM, cell cycle phase, and count threshold
#' @param scrna_mat count matrix  for single cell RNAseq
#' @param genenametype nomenclature for genes in count matrix: "hgnc_symbol" or "ensembl_id"
#' @param cellcycle filter for specific cell cycle phase: "G1", "G2M", or "S". Default is NULL, for no cell cycle phase filtering.
#' @param count.threshold remove genes where average count is less than count.threshold. Default is NULL, for no count threshold filtering
#' @return count matrix after quality control
#' @export
prep_scrnaseq <- function(scrna_mat, genenametype = "hgnc_symbol",cellcycle=NULL,count.threshold=NULL){
	if(genenametype != "hgnc_symbol" & genenametype !="ensembl_id" ){
		cat("Error: genenametype must be either 'hgnc_symbol' or 'ensembl_id'\n")
		return()
	}
	
	#get chromosomal location
	ensembl=useMart("ensembl")
	ensembl=useDataset("hsapiens_gene_ensembl",mart=ensembl)
	#note that getBM output is in sorted order
	
	if(genenametype=="hgnc_symbol"){
		genelocation = getBM(attributes=c("hgnc_symbol", "chromosome_name", "start_position", "end_position","strand","ensembl_gene_id"), filters="hgnc_symbol", values=list(rownames(scrna_mat)), mart=ensembl)
		#some genes do not have ensembl IDs while others have more than one ensemble ids. Remove duplicated gene symbols 
		kg = which(duplicated(genelocation$hgnc_symbol))
		genelocation = genelocation[-kg,]
		rownames(genelocation)=genelocation$hgnc_symbol
		#reorder genelocation to match rownames(scrna_mat)
		gl = genelocation[match(rownames(scrna_mat),rownames(genelocation)),]		
		rownames(gl)=rownames(scrna_mat)
		
		#convert rownames of scrna_mat to ENSEMBL IDs, leave as symbol if there is no matching ensembl_id
		scrna_mat.en = scrna_mat
		rownames(scrna_mat.en)=gl$ensembl_gene_id
		rownames(scrna_mat.en)[which(is.na(rownames(scrna_mat.en)))] = rownames(scrna_mat)[which(is.na(rownames(scrna_mat.en)))] 
		#for(i in 1:nrow(scrna_mat.en)) if(is.na(rownames(scrna_mat.en)[i])) rownames(scrna_mat.en)[i]=rownames(scrna_mat)[i]
		#for 2 genes that map to same ensembl ID, concatenate symbol and ensembl id
		rownames(scrna_mat.en)[which(duplicated(rownames(scrna_mat.en)))] = paste0(gl$ensembl_gene_id[which(duplicated(rownames(scrna_mat.en)))],".",gl$hgnc_symbol[which(duplicated(rownames(scrna_mat.en)))])



	} else if(genenametype=="ensembl_id"){
		genelocation = getBM(attributes=c("hgnc_symbol", "chromosome_name", "start_position", "end_position","strand","ensembl_gene_id"), filters="ensembl_gene_id", values=list(rownames(scrna_mat)), mart=ensembl)
		kg = which(duplicated(genelocation$ensembl_gene_id))
		genelocation = genelocation[-kg,]
		rownames(genelocation)=genelocation$ensembl_gene_id
		gl = genelocation[rownames(scrna_mat),]
		rownames(gl)=rownames(scrna_mat)
		
		scrna_mat.en= scrna_mat

	} else {
		cat('genenametype should be "hgnc_symbol" or "ensembl_id"\n')
		return()
	}

	sce <- SingleCellExperiment(list(counts=scrna_mat.en))  #counts must be matrix

	#Quality Control on cells
	#-------------------------------
	#library(scater)
		
	mito <- which(gl$chromosome_name=="M")
	#update defunct functions calculateQCMetrics to perCellQCMetrics
	#sce <- calculateQCMetrics(sce, feature_controls=list(Mt=mito))
	#libsize.drop <- isOutlier(sce$total_counts, nmads=3, type="lower",log=TRUE)
	#feature.drop <- isOutlier(sce$total_features_by_counts, nmads=3, type="lower",log=TRUE,)
	#mt.drop = isOutlier(sce$pct_counts_Mt, nmads=3, type="higher")
	
	sce.percell = perCellQCMetrics(sce,subsets=list(Mt=mito))
	libsize.drop <- isOutlier(sce.percell[,"sum"], nmads=3, type="lower",log=TRUE)
	feature.drop <- isOutlier(sce.percell[,"detected"], nmads=3, type="lower",log=TRUE,)
	mt.drop = isOutlier(sce.percell[,"subsets_Mt_percent"], nmads=3, type="higher")

	keep <- !(libsize.drop | feature.drop | mt.drop)
	#data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop), ByMt=sum(mt.drop), Remaining=sum(keep))
	sce$PassQC <- keep
	sce <- sce[,keep]

	#classification of cell cycle phase. Consider looking only at subsets of cells in the same phase
	#library(scran)
	if(!is.null(cellcycle)){
		set.seed(1234)
		hs.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))
		#mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran")) #for mouse
		assignments <- cyclone(sce, hs.pairs, gene.names=rownames(sce))
		sce$phases <- assignments$phases
		sce = sce[,which(sce$phases==cellcycle)]

	}

	if(!is.null(count.threshold)){
		#replaced defunct calcAverage with calculateAverage
		#ave.counts <- calcAverage(sce, use_size_factors=FALSE)
		ave.counts <- calculateAverage(sce, size_factors=NULL)
		#remove genes with average counts less than count.threshold
		counts.keep <- ave.counts >= count.threshold
		sce <- sce[counts.keep,]
	}

	#remove genes that are not expressed in any cell, this is also done as part of getdge
	num.cells <- nexprs(sce, byrow=TRUE)
	to.keep <- num.cells > 0
	sce <- sce[to.keep,]

	exprs_mat <- assay(sce, i = "counts")

	return(exprs_mat)
}


#cellcycle is G1, G2M or S
#species is human or mouse
#' Cell cycle filter
#'
#' filter for cells with particular cell cycle phase: "G1", "G2M", or "S"
#' @param scrna_mat count matrix  for single cell RNAseq with ensembl gene IDs
#' @param cellcycle filter for specific cell cycle phase: "G1", "G2M", or "S". Default is NULL, for no cell cycle phase filtering.
#' @param species "human" or "mouse"
#' @return count matrix containing only cells with the selected cell cycle phase
#' @export
getcellcycle <- function(scrna_mat, cellcycle="G1",species="human"){
	set.seed(1234)
	sce <- SingleCellExperiment(list(counts=scrna_mat)) 
		if(species=="human"){
				cc.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))
		} else if (species=="mouse"){
				cc.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran")) #for mouse
		} else {
			cat("Error: Species must be 'human' or 'mouse'\n")
			return()
		}
		assignments <- cyclone(sce, cc.pairs, gene.names=rownames(sce))
		sce$phases <- assignments$phases
		sce = sce[,which(sce$phases==cellcycle)]
		exprs_mat.cellcycle = assay(sce,i="counts")[,which(sce$phases==cellcycle)]

	return(exprs_mat.cellcycle)	
}


#' Mean Correlation
#'
#' get mean correlation by first doing Fisher transform to z, averaging z, then reverse transform back
#' do quality control of scRNAseq based on library size, feature counts, chrM, cell cycle phase, and count threshold
#' @param rho correlations
#' @return mean correlation and 95% CI
#' @export
getmeancorr <- function(rho){
	rho = rho[which(!is.na(rho))]
	rho = ifelse(abs(rho)>0.999999,sign(rho)*0.999999,rho) #to avoid the infinity problem
	z = fisherz(rho)
	mz = mean(z, na.rm=TRUE)
	sdz = sd(z, na.rm=TRUE)
	mrho = fisherz2r(mz)
	#sdrho = fisherz2r(sdz)
	cirho = fisherz2r(CI(fisherz(rho),ci=0.95)) #upper CI, mean, lower CI

	return(cirho)
}




#x1,x1b are mixture matrices where rows are samples, columns are cell types
#row names and column names for x1 and x2 should match
#' compute correlation between estimated and actual composition for each sample
#' @param x1 composition matrix where rows are samples and columns are cell types
#' @param x2 second composition matrix where rows and columns match x1
#' @return correlation by sample
#' @export
getcorr <- function(x1, x2){
	
	#correlation of each sample over all cell types between predicted and actual cell mixture
	cor.sample = mapply(cor,split(as.matrix(x1),row(x1)),split(as.matrix(x2),row(x2)))
	
	return(cor.sample)
}







