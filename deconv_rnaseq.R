#deconvolution of RNAseq

library(edgeR)
library(statmod) #for glmnb.fit
library(nloptr) #for nloptr, nonlinear optimization with constraints


set.seed(1234)

#======================================================
#INPUT FILES AND PARAMETERS  
#======================================================

#individual cell types
setwd("/proj/regeps/regep00/studies/aneurysms/analyses/rerdu/rnaseq_leucegene_grch38/")
fildir=paste0("SRR102",2945:3003)
filnames=sapply(1:length(fildir), function(x) paste0("/proj/regeps/regep00/studies/aneurysms/analyses/rerdu/rnaseq_leucegene_grch38/results/",fildir[x],"/htseq_out/",fildir[x],"_accepted_hits.sortedname.gene.counts"))
sample.id=paste0("SRR102",2945:3003)
tissue=c(rep("Bcells",14),rep("Granulocytes",14),rep("Monocytes",14),rep("CD34",3),rep("Tcells",14))
celltypes = c("Bcells","Tcells","Monocytes","Granulocytes")

#whole blood
setwd("/proj/regeps/regep00/studies/aneurysms/analyses/rerdu/rnaseq_grch38_blood_cell")
fildir_s=c("RD1_0314R1D00004_A160712690","RD5_0812R1D00010_A160712694", "RD6_0912R1D00013_A160712695","RD7_0913R1D00010_A160712696","RD8_1012R1D00003_A160712697")
filnames_s=c()
filnames_s = c(sapply(1:length(fildir_s), function(x) paste0("/proj/regeps/regep00/studies/aneurysms/analyses/rerdu/rnaseq_grch38_blood_cell/results/", fildir_s[x], "/htseq_out/", fildir_s[x], "_accepted_hits.sortedname.gene.counts")))
sample.id_s = fildir_s

#differential cell counts for whole blood	
cbc = data.frame(sample.id=c("RD1_0314R1D00004_A160712690","RD5_0812R1D00010_A160712694", "RD6_0912R1D00013_A160712695","RD7_0913R1D00010_A160712696","RD8_1012R1D00003_A160712697"), lymph=c(21.4,24.9,0,9.3,6.1), mono=c(6.1,6.5,5,7.6,6.4), gran=c(72.3,68.6,94,83.1,87.5))



#will keep genes that achieve ncpm count per million in at least nsamp samples (used ncpm=1 and nsamp=4)
ncpm = 1
nsamp = 4  

#number of genes to be used in projection matrix for individual cell types
#nb0list = seq(10,200,10)
nb0list=c(10,100)
#number of genes to be used in projection matrix for whole tissue
#nb0list_s = seq(10,200,10)
nb0list_s=c(10,100)



#==================================================
#INITIALIZE VARIABLES
#==================================================
tissuei = which(tissue %in% celltypes)
tissuef = factor(tissue[tissuei], levels=celltypes)
design <- model.matrix(~-1+tissuef) 
design_s <- NULL



#==================================================
#FUNCTIONS
#==================================================
#EdgeR analysis, obtain B0
getb0 <- function(filnames, sample.id, design, ncpm, nsamp){

	x=list()
	for(i in 1:length(filnames)){
		x[[i]] = read.delim(filnames[i], header=FALSE)
		colnames(x[[i]])=c("gene",sample.id[i])
	}


	x2=x[[1]]
	for(i in 2:length(filnames)){
		x2 = merge(x2,x[[i]],by.x="gene",by.y="gene")	
	}
	rownames(x2)=x2$gene

	#delete rows that are all 0. First column is gene name
	z=rowSums(x2[2:ncol(x2)])
	z2 = which(z==0) 
	x2 = x2[-z2,] 

	x4 = x2[,-1]
	x4 = x4[,tissuei]

	dge=DGEList(x4,genes=rownames(x4))

	#filter
	keep <- rowSums(cpm(dge)>ncpm)>=nsamp
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







#use the most significant NB0 genes for projection matrix and outputs cell mixture x1, convergence, correlation by cell type, correlation by sample
#MAXITER is maximum iterations before calculation is stopped and considered nonconvergent
getx1<-function(NB0, b0, design, dge_s, pfstat, MAXITER=1000){
	b0r = b0[order(pfstat)[1:NB0],]
	
	offset_s = getOffset(dge_s)
	dispersion_s = dge_s$tagwise.dispersion #genewise dispersion
	nsamples_s = nrow(dge_s$samples)

	#initialize mixture matrix
	ntypes = ncol(b0)
	x1 = matrix(0,nrow=nsamples_s,ncol=ntypes)
	one_vec = rep(1,ntypes)
	x0 = rep(1/ntypes, ntypes) #for initial mixture, use equal amounts
	gx <- function(x) {one_vec %*% x-1} #inequality constraint
	converged = c() #convergence if iter < maxiter



	for(i in 1:nsamples_s){
		#look at the top NB0 genes only
		y1 = dge_s$counts[which(rownames(dge_s$counts) %in% rownames(b0r)),i] 
		#obj is negative log likelihood, x is the mixture vector
		obj <- function(x){-sum(dnbinom(y1,size=dispersion_s, mu=exp(b0r %*% x + offset_s[i])))}
		#use constraint that gx<=0 or sum_i x_i <=1
		#constrained lower bound to 0, upper bound to 1 for each x_i and sum_i of x_i <=1, leaving room for other cellular components that may not have been accounted for
		opt_result = nloptr(x0=x0, eval_f=obj,lb=rep(0,ntypes), ub=rep(1,ntypes), eval_g_ineq=gx, opts=list(algorithm = "NLOPT_LN_COBYLA", maxeval=MAXITER))
	
		x1[i,]=opt_result$solution #predicted coefficients
		converged = c(converged,ifelse(opt_result$iterations < MAXITER,1,0))
	}

	colnames(x1)=colnames(b0r)
	rownames(x1)=rownames(dge_s$samples)
	x1=as.data.frame(x1)
	des = as.data.frame(design)
	
	#correlation of each cell type over all samples between predicted and actual cell mixture
	cor.type = mapply(cor,x1[,1:ncol(x1)],des[,1:ncol(x1)])
	names(cor.type)=celltypes
	#correlation of each sample over all cell types between predicted and actual cell mixture
	cor.sample = mapply(cor,split(as.matrix(x1),row(x1)),split(as.matrix(des),row(des)))

	return(list(x1=x1,converged=converged,cor.type=cor.type,cor.sample=cor.sample))
}	





#given fastq files, get dge
getdge <- function(filnames_s, sample.id_s, design_s, ncpm, nsamp){
	
	x_s=list()
	for(i in 1:length(filnames_s)){
		x_s[[i]] = read.delim(filnames_s[i], header=FALSE)
		colnames(x_s[[i]])=c("gene",sample.id_s[i])
	}


	x2_s=x_s[[1]]
	for(i in 2:length(filnames_s)){
		x2_s = merge(x2_s,x_s[[i]],by.x="gene",by.y="gene")	
	}
	rownames(x2_s) = x2_s$gene

	x4_s = x2_s[,-1]

	dge_s=DGEList(x4_s,genes=rownames(x4_s))

	#filter
	keep_s <- rowSums(cpm(dge_s)>ncpm)>=nsamp
	dge_s <- dge_s[keep_s,]

	#recompute library sizes
	dge_s$samples$lib.size <- colSums(dge_s$counts)

	#normalize. Default method is TMM
	dge_s <- calcNormFactors(dge_s)
	offset_s = getOffset(dge_s)

	#estimate dispersion
	#first estimate common dispersion
	dge_s <- estimateGLMCommonDisp(dge_s, design_s, verbose=TRUE)

	#estimate gene specific dispersions
	dge_s <- estimateGLMTrendedDisp(dge_s,design_s)
	dge_s <- estimateGLMTagwiseDisp(dge_s,design_s)

	return(dge_s)	
}






#use this to get x1 mixture data for tissue using cell type b0
getx1_tissue<-function(NB0, b0, dge_s, MAXITER=1000){
	b0r = b0[order(pfstat)[1:NB0],]
	ntypes = ncol(b0r)
	
	#for default, start with equal proportions of all cell types but can give initial x0 starting mixture
	x0=rep(1/ntypes, ntypes)
	
	offset_s = getOffset(dge_s)
	dispersion_s = dge_s$tagwise.dispersion #genewise dispersion
	nsamples_s = nrow(dge_s$samples)

	#initialize mixture matrix
	x1 = matrix(0,nrow=nsamples_s,ncol=ntypes)
	one_vec = rep(1,ntypes)
	gx <- function(x) {one_vec %*% x-1} #inequality constraint
	converged = c() #convergence if iter < maxiter

	for(i in 1:nsamples_s){
		#look at the top NB0 genes only
		y1 = dge_s$counts[which(rownames(dge_s$counts) %in% rownames(b0r)),i] 
		#obj is negative log likelihood, x is the mixture vector
		obj <- function(x){-sum(dnbinom(y1,size=dispersion_s, mu=exp(b0r %*% x + offset_s[i])))}
		#use constraint that gx<=0 or sum_i x_i <=1
		#constrained lower bound to 0, upper bound to 1 for each x_i and sum_i of x_i <=1, leaving room for other cellular components that may not have been accounted for
		opt_result = nloptr(x0=x0, eval_f=obj,lb=rep(0,ntypes), ub=rep(1,ntypes), eval_g_ineq=gx, opts=list(algorithm = "NLOPT_LN_COBYLA", maxeval=MAXITER))
	
		x1[i,]=opt_result$solution #predicted coefficients
		converged = c(converged,ifelse(opt_result$iterations < MAXITER,1,0))
	}

	colnames(x1)=colnames(b0r)
	rownames(x1)=rownames(dge_s$samples)
	x1=as.data.frame(x1)
	
	return(list(x1=x1,converged=converged))
}









#==================================================
#MAIN
#==================================================

#-------------------------------------------
#project individual cell types onto cell types to get mixture fractions using the top NB0 genes
#-------------------------------------------
resultb0 = getb0(filnames, sample.id, design, ncpm, nsamp)

dge = resultb0$dge
b0 = resultb0$b0
pfstat = resultb0$pfstat.adjust


#may need to change nb0list examined to get optimal one
result=c()
result = sapply(nb0list, function(x) getx1(x, b0, design, dge, pfstat, MAXITER=1000))


#plot results
cortype = cbind(nb0list, matrix(unlist(result[3,]),ncol=4, byrow=TRUE))
colnames(cortype)=c("nb0","Bcells","Tcells","Monocytes","Granulocytes")
par(oma=c(2,3,1,1))
matplot(cortype[,"nb0"], cortype[,c(2:5)],type="b", xlab="Number of genes", ylab="Correlation",col=c("green","blue","red","purple"), pch=19, cex=2, cex.axis=2, cex.lab=2)
legend("bottomright",c("Bcells","Tcells","Monocytes","Granulocytes"), col=c("green","blue","red","purple"), cex=2, pch=19)

corsample = cbind(nb0list,matrix(unlist(result[4,]), ncol=length(result[4][[1]]), byrow=TRUE))
colnames(corsample)=c("nb0",sample.id[c(1:42,46:59)])
matplot(corsample[,"nb0"], corsample[,c(2:ncol(corsample))],type="l", xlab="Number of genes", ylab="Correlation")


#-------------------------------------------
#project whole tissue onto cell types, given projection matrix
#-------------------------------------------

dge_s = getdge(filnames_s, sample.id_s, design_s, ncpm, nsamp)

result.tissue=c()
result.tissue = sapply(nb0list_s, function(x) getx1_tissue(x, b0, dge_s, MAXITER=1000))


#calculate correlations
cor.sample.tissue=c()
for(j in 1:length(nb0list_s)){
	mix = result.tissue[1,j]$x1
	mix2 = cbind(mix[,"tissuefBcells"]+mix[,"tissuefTcells"], mix[,"tissuefMonocytes"], mix[,"tissuefGranulocytes"])
	colnames(mix2)=c("lymph","mono","gran")
	rownames(mix2)=rownames(mix)
	cor.sample.tissue = rbind(cor.sample.tissue,sapply(1:nrow(cbc), function(k) { cor(as.numeric(mix2[k,]),as.numeric(cbc[which(as.character(cbc$sample.id) %in% rownames(mix2)[k]),2:4]))}))
}
cor.sample.tissue = cbind(nb0list_s, cor.sample.tissue)
colnames(cor.sample.tissue)=c("nb0",as.character(cbc$sample.id))







