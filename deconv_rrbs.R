#deconvolution for RRBS data
#input are outputs from bsmap filtered for CpGs with the following format (*_CpG_chr.txt):
#    1) chromosome
#    2) coordinate (1-based)
#    3) strand
#    4) sequence context (CG|CHG|CHH)
#    5) methylation ratio, calculated as #C_counts / #eff_CT_counts
#    6) number of effective total C+T counts on this locus (#eff_CT_counts). This excludes the real Ts and only counts the converted Ts.
#            CT_SNP="no action", #eff_CT_counts = #CT_counts
#            CT_SNP="correct", #eff_CT_counts = #CT_counts * (#rev_G_counts / #rev_GA_counts)
#    7) number of total C counts on this locus (#C_counts) 
#    8) number of total C+T counts on this locus (#CT_counts)
#    9) number of total G counts on this locus of reverse strand (#rev_G_counts)
#    10) number of total G+A counts on this locus of reverse strand (#rev_GA_counts)
#    11) Neutrophilslower bound of 95% confidence interval of methylation ratio, calculated by Wilson score interval for binomial proportion.
#    12) upper bound of 95% confidence interval of methylation ratio, calculated by Wilson score interval for binomial proportion.





library(methylKit)
library(lme4)
library(blme)
library(nloptr) 
library(hash) #hash tables for color keys for plots
library(parallel)

source("/udd/rerdu/my_diffMethFromDesign6.R")






#===============================================
#INPUT FILES AND PARAMETERS
#===============================================


#note that RNAseq file used "Granulocytes" rather than "Neutrophils"
#individua cell types
basefil="geo_btmn"
setwd("/proj/regeps/regep00/studies/aneurysms/analyses/rerdu/rrbs_geo_grch38")
sample.id=c("SRR1104848","SRR1104855","SRR1104856","SRR1104857", "SRR1104838","SRR1104839","SRR1104841","SRR1104842","SRR1104852", "SRR1104853","SRR1104854","SRR2960993","SRR1508407","SRR1508408","SRR1508409","SRR1508410", "SRR1508411","SRR1508412","SRR1508413","SRR1508414","SRR1508415", "SRR1508416","SRR1508417","SRR1508418")
filnames = paste0("/proj/regeps/regep00/studies/aneurysms/analyses/rerdu/rrbs_geo_grch38/", sample.id, "_CpG_chr.txt")
treatmentc = c(rep("Monocytes",4),rep("Tcells",4), rep("Bcells",3), rep("Neutrophils",13))
treatment = as.numeric(as.character(factor(treatmentc, levels=c("Bcells","Tcells","Monocytes","Neutrophils"),labels=c(0,1,2,3))))
design=data.frame(Bcells=ifelse(treatmentc=="Bcells",1,0),Tcells=ifelse(treatmentc=="Tcells",1,0),Monocytes=ifelse(treatmentc=="Monocytes",1,0),Neutrophils=ifelse(treatmentc=="Neutrophils",1,0))
formula='Bcells+Tcells+Monocytes+Neutrophils+0'
celltypes=c("Bcells","Tcells","Monocytes","Neutrophils")

#whole blood
basefil_s="whole_blood"
sample.id_s=c("10009987","10003442","10003445","10013595", "10003462")
filnames_s = paste0("/proj/regeps/regep00/studies/aneurysms/analyses/rerdu/rrbs_grch38/", sample.id_s, "_CpG_chr.txt")	
treatment_s = rep(0,length(sample.id_s))
cbc = data.frame(sample.id=c("10013595","10003462", "10003442","10009987","10003445"), lymph=c(21.4,24.9,0,9.3,6.1), mono=c(6.1,6.5,5,7.6,6.4), gran=c(72.3,68.6,94,83.1,87.5))
rownames(cbc)=cbc$sample.id
cbc=cbc[as.character(sample.id_s),]



#parameters
assembly="grch38"
context="CpG"
#destranded=FALSE  #can use destranded = TRUE to get better coverage
mixedef = FALSE  
fstat = TRUE
ncores = 8  #number of cores to use
	
	
	



#===============================================
#FUNCTIONS
#===============================================

#takes bsmap output and divides into 20 files of meth objects for use with methylkit
#output is basefil_cpg_meth*.rda
makemeth <- function(filnames, sample.id, treatment, assembly, context){

	obj=read(as.list(filnames), sample.id=as.list(sample.id), assembly=assembly,header=FALSE, context=context,resolution="base", pipeline=list(fraction=TRUE,chr.col=1,start.col=2,end.col=2,coverage.col=6,strand.col=3, freqC.col=5), treatment=treatment)
	meth=unite(obj,min.per.group=NULL)
	meth.data = getData(meth)

	#remove rows in which all Ts are 0 or all Cs are 0 
	meth.filter = meth.data[-which(rowSums(meth.data[,meth@numTs.index])==0 | rowSums(meth.data[,meth@numCs.index])==0),]

	meth2 = new("methylBase",meth.filter, sample.ids=meth@sample.ids,assembly=meth@assembly, context=meth@context, treatment=meth@treatment,coverage.index=meth@coverage.index, numCs.index=meth@numCs.index, numTs.index=meth@numTs.index, destranded=meth@destranded, resolution=meth@resolution)
	meth2.data = getData(meth2)

	size = nrow(meth2.data)
	#***divide into 20 files for parallel processing
	nbin = 20 
	binsize = as.integer(size/nbin)

	datalist = list()
	for(i in 1:nbin){
		i1 = ((i-1)*binsize+1)
		i2 = i*binsize
		if(i < nbin){
			datalist[[i]] = meth2.data[i1:i2,]	
		} else {
			datalist[[i]] = meth2.data[i1:size,]	
		}
	}

	for(i in 1:nbin){
		methi = new("methylBase",datalist[[i]], sample.ids=meth@sample.ids,assembly=meth@assembly, context=meth@context, treatment=meth@treatment,coverage.index=meth@coverage.index, numCs.index=meth@numCs.index, numTs.index=meth@numTs.index, destranded=meth@destranded, resolution=meth@resolution)

		outfil = paste0(basefil,"_cpg_meth",i,".rda")
		cat(outfil,"\n")
		#save(methi,file=outfil)
	}
}






#input: *_cpg_meth*.rda
#output: *_cpg_meth*pv.csv
pmeth <- function(cpgmethfil, design, formula, ncores, mixedef,fstat){
	
	load(cpgmethfil)

	#split this up
	CHUNKSIZE=1000
	chunks = floor(dim(methi)[1]/CHUNKSIZE)
	reschunk = dim(methi)[1] %% CHUNKSIZE
	if(reschunk>0) chunks = chunks+1
	#initialize output file
	outfil = paste0(substring(cpgmethfil,1,nchar(cpgmethfil)-4),"pv_test.csv")

	for(i in 1:chunks){
		cat(i,"\n")
		if(i< chunks){
			j = (i-1)*CHUNKSIZE+1
			k = i*CHUNKSIZE
			methchunk = methi[j:k,]
		} else {
			j = (i-1)*CHUNKSIZE+1
			k = dim(methi)[1]
			methchunk = methi[j:k,]
		}
		mdiff_dat = my_diffMethFromDesign(methchunk,design,formula,ncores,mixedef,fstat)
		if(i==1){
			write.table(mdiff_dat, file=outfil, sep=",", quote=FALSE, row.names=FALSE, col.names=TRUE, append=FALSE)
		} else {
			write.table(mdiff_dat, file=outfil, sep=",", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
		}
	}	
}






#get cell type mixture for individual cell type data, and correlations
getx1 <- function(NB0,b0,meth,celltypes,MAXITER=10000){
	
	b0r = b0[order(b0$pfstat)[1:NB0],]
	ntypes=length(grep("beta",colnames(b0)))	
	meth.data = getData(meth)
	nsamples=length(meth@sample.ids)

	#initialize mixture matrix
	x1 = matrix(0,nrow=nsamples,ncol=ntypes)
	one_vec = rep(1,ntypes)
	x0 = rep(1/ntypes, ntypes) #for initial mixture, use equal amounts
	gx <- function(x) {one_vec %*% x-1} #inequality constraint
	converged = c() #convergence if iter < maxiter
	miter = c()

	#get only the NB0 genes from meth.data
	meth.small = merge(x=b0r,y=meth.data,by=c("chr","start","end","strand"))
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
	rownames(x1)=meth@sample.ids
	x1=as.data.frame(x1)
		
	cor.type=c()
	for(i in 1:length(celltypes)){
		cor.type = c(cor.type,cor(x1[,paste0("beta.",celltypes[i])], design[,celltypes[i]]))
	}
	names(cor.type)=celltypes
	
	cor.sample = sapply(1:nrow(x1), function(i){cor(as.numeric(x1[i,paste0("beta.",celltypes)]), as.numeric(design[i,celltypes]))}) 
	
	return(list(x1=x1,converged=converged,cor.type=cor.type,cor.sample=cor.sample))
}





#use this to get x1 for tissue using b0
getx1.tissue <- function(NB0,b0,meth_s,celltypes,MAXITER=10000){
	b0r = b0[order(b0$pfstat)[1:NB0],]
	ntypes=length(grep("beta",colnames(b0)))
	meth.data = getData(meth_s)
	nsamples=length(meth_s@sample.ids)

	#initialize mixture matrix
	x1 = matrix(0,nrow=nsamples,ncol=ntypes)
	one_vec = rep(1,ntypes)
	x0 = rep(1/ntypes, ntypes) #for initial mixture, use equal amounts
	gx <- function(x) {one_vec %*% x-1} #inequality constraint
	converged = c() #convergence if iter < maxiter
	miter = c()


	#get only the NB0 genes from meth.data
	meth.small = merge(x=b0r,y=meth.data,by=c("chr","start","end","strand"))
	meth.cn = colnames(meth.small)

	#get names of columns of the coefficients, ie beta
	type.coef = colnames(b0r)[grep("beta",colnames(b0r))]

	#save.image("ec_smc_fib.rda")

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
	rownames(x1)=meth_s@sample.ids
	x1=as.data.frame(x1)
	
	return(list(x1=x1,converged=converged))
}








#=================================================
#MAIN
#=================================================

makemeth(filnames, sample.id, treatment, assembly, context)

cpgmethfilnames = paste0(basefil,"_cpg_meth",1:20,".rda")

#get b0
mclapply(cpgmethfilnames, function(x) pmeth(x, design, formula, ncores, mixedef,fstat))

filnames_b0 = sapply(1:20, function(x) paste0(basefil,"_cpg_meth",x,"pv.csv"))
b0=c()
for(i in 1:20){
	idat = read.csv(filnames_b0[i],as.is=TRUE)
	b0 = rbind(b0,idat)	
}




#-----------------------------------
#validate individual cell type
#-----------------------------------

obj=read(as.list(filnames), sample.id=as.list(sample.id), assembly=assembly,header=FALSE, context=context,resolution="base", pipeline=list(fraction=TRUE,chr.col=1,start.col=2,end.col=2,coverage.col=6,strand.col=3, freqC.col=5), treatment=treatment)
meth=unite(obj,min.per.group=NULL)

#get cell type mixture for nb0list number of genes
nb0list = seq(100,2000,100)
result=c()
result = sapply(nb0list, function(x) getx1(x,b0,meth,celltypes=celltypes))

#make plots
cortype = cbind(nb0list, matrix(unlist(result[3,]),ncol=4, byrow=TRUE))
colnames(cortype)=c("nb0","Bcells","Tcells","Monocytes","Granulocytes")
	
par(oma=c(2,3,1,1))
matplot(cortype[,"nb0"], cortype[,c(2:5)],type="b", xlab="Number of CpG sites", ylab="Correlation",col=c("green","blue","red","purple"), pch=19, cex=2, cex.axis=2, cex.lab=2)
legend("bottomright",c("Bcells","Tcells","Monocytes","Granulocytes"), col=c("green","blue","red","purple"), cex=2, pch=19)

corsample = cbind(nb0list,matrix(unlist(result[4,]), ncol=24, byrow=TRUE))
colnames(corsample)=c("nb0",sample.id)
par(oma=c(2,3,1,1))
matplot(corsample[,"nb0"], corsample[,c(2:25)],type="l", xlab="# Sites", ylab="Correlation",cex=2, cex.axis=2, cex.lab=2,lwd=2)




#-----------------------------------
#validate whole tissue
#-----------------------------------

obj_s=read(as.list(filnames_s), sample.id=as.list(sample.id_s), assembly=assembly,header=FALSE, context=context, resolution="base", pipeline=list(fraction=TRUE,chr.col=1,start.col=2,end.col=2,coverage.col=6,strand.col=3, freqC.col=5), treatment=treatment_s)
meth_s=unite(obj_s,min.per.group=NULL)
nsamples_s=length(meth_s@sample.ids)
ntypes = length(celltypes)

#get cell type mixture for nb0list number of genes
nb0list_s = seq(100,2000,100)
result.tissue=c()
result.tissue = sapply(nb0list_s, function(x) getx1.tissue(x,b0,meth_s,celltypes=celltypes))
	
#calculate correlations
cor.sample.tissue=c()
for(j in 1:length(nb0list_s)){
	mix = result.tissue[1,j]$x1[as.character(cbc$sample.id),]
	mix2 = cbind(mix[,"beta.Bcells"]+mix[,"beta.Tcells"], mix[,"beta.Monocytes"], mix[,"beta.Neutrophils"])
	colnames(mix2)=c("lymph","mono","gran")
	rownames(mix2)=rownames(mix)
	cor.sample.tissue = rbind(cor.sample.tissue,sapply(1:nrow(cbc), function(k) { cor(as.numeric(mix2[k,]),as.numeric(cbc[which(as.character(cbc$sample.id) %in% rownames(mix2)[k]),2:4]))}))
}
cor.sample.tissue = cbind(nb0list_s, cor.sample.tissue)
colnames(cor.sample.tissue)=c("nb0",as.character(cbc$sample.id))

















