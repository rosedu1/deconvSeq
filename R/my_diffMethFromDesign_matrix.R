#modified from methylKit
f1<-function(cutoff,rawp){sum(rawp<cutoff)/length(rawp)}

my_diffMethFromDesign_matrix<-function(cur_methmat, cur_design, cur_formula, cur_treatment=NULL, mixedef=FALSE, fstat=TRUE){

	#has this format:  chr start   end strand coverage1 numCs1 numTs1 coverage2 numCs2 numTs2	
	numCs.index = grep("Cs",colnames(cur_methmat))
	numTs.index = grep("Ts",colnames(cur_methmat))
	c_data = cur_methmat[,numCs.index]
	t_data = cur_methmat[,numTs.index]
	num_base=nrow(cur_methmat)


	if(!is.null(cur_treatment)){
			#get the indices (column number) for numCs and numTs in each set
		#note that the diff meth only compares treatments 0 and 1
		set1.Cs=numCs.index[cur_treatment==1]
		set2.Cs=numCs.index[cur_treatment==0]
		set1.Ts=numTs.index[cur_treatment==1]
		set2.Ts=numTs.index[cur_treatment==0]

    	# calculate mean methylation change
		if(length(set1.Cs) > 1){ #if more than one subject in treatment group
			mom.meth1=100*rowMeans(cur_methmat[,set1.Cs]/cur_methmat[,set1.Cs-1],na.rm=TRUE) # get means of means, column set.Cs-1 is the total coverage column (C+T)
			pm.meth1=100*rowSums(cur_methmat[,set1.Cs])/rowSums(cur_methmat[,set1.Cs-1],na.rm=TRUE) # get weighted means, this is total Cs for all samples/total coverage for all samples
		}else{
			mom.meth1    = 100*(cur_methmat[,set1.Cs]/cur_methmat[,set1.Cs-1]) # get % methylation
			pm.meth1     = mom.meth1
		}

		if(length(set2.Cs)>1){
			mom.meth2=100*rowMeans(cur_methmat[,set2.Cs]/cur_methmat[,set2.Cs-1],na.rm=TRUE)
			pm.meth2=100*rowSums(cur_methmat[,set2.Cs])/rowSums(cur_methmat[,set2.Cs-1],na.rm=TRUE) # get weighted means
		}else{
			mom.meth2    = 100*(cur_methmat[,set2.Cs]/cur_methmat[,set2.Cs-1]) # get % methylation
			pm.meth2     = mom.meth2
		}
		pm.mean.diff=pm.meth1-pm.meth2
		mom.mean.diff=mom.meth1-mom.meth2
			
	} else {
		pm.mean.diff = NA
		mom.mean.diff = NA
	}
                          


	## Make a function that given an index, returns the methylation difference,
	##coefficient value (for the last term in the formula)
	## and the wald test p-value grabbed from the glm object at that index
    #***return all coefficients and p values for each covariate
###############
	cur_function<-function(x){
		
		if(x %% 100==0) cat(x,"\n")
		# x is an index
		# it is assumed that the appropriate data has been loaded into
		# c_data, t_data, cur_design, and cur_formula
		cur_data=data.frame(cur_design, t(c_data[ x, ]), t(t_data[ x, ]))
		colnames(cur_data)=c(colnames(cur_design), 'Cs', 'Ts')
	
		#initialize error markers
		botherror=0
		botherror.0=0
	
	
		#obj=glm(as.formula(paste0("cbind(Cs, Ts)~", cur_formula)), data=cur_data, family=binomial(link=logit) )
		#obj=glmer(as.formula(paste0("cbind(Cs, Ts)~", cur_formula)), data=cur_data, family=binomial(link=logit))
		#obj=bglmer(as.formula(paste0("cbind(Cs, Ts)~", cur_formula)), data=cur_data, family=binomial(link=logit),fixef.prior = normal)
		
		#use glmer unless there is a warning or error, then use bglmer, and if there are still convergence problems then output NA
	
	
		if(mixedef==TRUE){
		#use mixed effects model
        #if bglmer results in an error or warning, return 1
			#bf <- tryCatch(bglmer(as.formula(paste0("cbind(Cs, Ts)~", cur_formula)), data=cur_data, family=binomial(link=logit),fixef.prior = normal), warning=function(w) {1}, error=function(e) {1})
            #if glmer results in an error or a warning, use bf for bglmer instead
			obj <- tryCatch(glmer(as.formula(paste0("cbind(Cs, Ts)~", cur_formula)), data=cur_data, family=binomial(link=logit)),error=function(e) {1}, warning=function(w) {1})
		
			#botherror= ifelse(class(obj) == "glmerMod" | class(obj)=="bglmerMod",0,1)
			botherror= ifelse(class(obj) == "glmerMod",0,1)
            
            if((fstat==TRUE) && (botherror==0)){
                #bf.0 <- tryCatch(bglmer(as.formula(paste0("cbind(Cs, Ts)~", 1)), data=cur_data, family=binomial(link=logit),fixef.prior = normal), warning=function(w) {1}, error=function(e) {1})
                #if glmer results in an error or a warning, use bf for bglmer instead
                obj.0 <- tryCatch(glmer(as.formula(paste0("cbind(Cs, Ts)~", 1)), data=cur_data, family=binomial(link=logit)),error=function(e) {1}, warning=function(w) {1})
                botherror.0= ifelse(class(obj.0) == "glmerMod",0,1)
                if(botherror.0==0){ 
                	fstat = anova(obj,obj.0,test="LRT") #for likelihood ratio test
                	pfstat = fstat[2,"Pr(>Chisq)"] #pvalue of full model vs null model
                } else {
                	fstat = NA
                	pfstat = NA
                }
            }
            
		} else {
		#use regular glm	
		#only modified this for cell type deconvolution
			#obj <- tryCatch( glm(as.formula(paste0("cbind(Cs, Ts)~", cur_formula)), data=cur_data, family=binomial(link=logit)), error=function(e) {1}, warning=function(w) {1})
			obj <- tryCatch( glm(as.formula(paste0("cbind(Cs, Ts)~", cur_formula)), data=cur_data, family=binomial(link=logit)), error=function(e) {1})
			botherror = ifelse(class(obj)[1] == "glm",0,1)
            
            if((fstat==TRUE) && (botherror==0)){
                #obj.0 <- tryCatch( glm(as.formula(paste0("cbind(Cs, Ts)~", 1)), data=cur_data, family=binomial(link=logit)), error=function(e) {1}, warning=function(w) {1})
                obj.0 <- tryCatch( glm(as.formula(paste0("cbind(Cs, Ts)~", 1)), data=cur_data, family=binomial(link=logit)), error=function(e) {1})
                botherror.0 = ifelse(class(obj.0)[1] == "glm",0,1)
                if(botherror.0==0){
                	lrt = anova(obj,obj.0,test="LRT") #for likelihood ratio test
                	pfstat = lrt[2,"Pr(>Chi)"] #pvalue of full model vs null model
                }
                
                #also do submodels with each cell type
                pfstat.type=c()
                ntypes=names(obj$coefficients)
                for (i in 1:length(ntypes)){
                		obj.type = tryCatch( glm(as.formula(paste0("cbind(Cs, Ts)~", ntypes[i])), data=cur_data, family=binomial(link=logit)), error=function(e) {1})
                		botherror.type = ifelse(class(obj.type)[1] == "glm",0,1)	
                		if(botherror.type==0){
                			lrt.type = 	anova(obj.type,obj.0,test="LRT") #for likelihood ratio test
                			pfstat.type=c(pfstat.type,lrt.type[2,"Pr(>Chi)"])
                		} else {
                			pfstat.type=c(pfstat.type,NA)
                		}
                }
                names(pfstat.type) = ntypes
                
            }
		}
		
		if(botherror==0 && botherror.0==0){
            
            #for mixed effects model, only fixed effects results are obtained
            beta <- obj$coefficients #includes coeff that are NAs
            p.value <- summary(obj)$coefficients[,'Pr(>|z|)'] #note that summary removes NAs
            se <- summary(obj)$coefficients[,'Std. Error']
            
            #put NAs back so that beat, p.value, and se are the same length
            if (length(p.value) != length(beta)){
            		test = rep(NA,length(beta))
            	names(test) = names(beta)
            	test[names(p.value)]=p.value
            	p.value = test
            }
                        
            if (length(se) != length(beta)){
            		test = rep(NA,length(beta))
            	names(test) = names(beta)
            	test[names(se)]=se
            	se = test
            }
            
            
            
            #for case of the one SMC sample having 0 Cs and Ts, there are no coefficient for SMC
              
            if (fstat==TRUE){
            		return(c(beta,se,p.value,pfstat.type,pfstat))
            } else {
            		return(c(beta,se,p.value))
            }
            
            
            #return( c(beta_last,p.value) )
		} else {
             return( NA )
            
		}
	}
	
###########################

	res= t(sapply(1:num_base, cur_function))
	
	
	#get length of res[[i]]
	if(class(res)=="list"){
		#if there is an NA (botherror=1 or bothererror.0=1) as res then this becomes a list
		lres=0
		i=1
		while(lres==0){
			if(!is.na(res[[i]][1])) lres=length(res[[i]])
			i=i+1
		}
	
		na.ind = which(unlist(lapply(res,function(x)is.na(x)[1])))
		for(j in na.ind) res[[j]]=rep(NA,lres)
		res2 = matrix(unlist(res),ncol=lres,byrow=TRUE)
		colnames(res2)=names(res[[i-1]])	
			
	} else {
		if(class(res[1,1])=="list"){
			#if no (single) NA then res is a 1xn matrix where each element is a list, ie res[1,1] is a list
			#take care of when there are fewer than the full set of coefficients and make that row NA
			
			reslen = sapply(res, length)
			reslentab = table(reslen)
			if(length(reslentab)>1){
				#1/2/17 the real length should be the one with the longest length, not necessarily the one with the highest frequency
				#realntype =  which(reslentab %in% max(reslentab))
				realntype =  which(names(reslentab) %in% max(names(reslentab)))
				
				bresi = which(reslen %in% as.numeric(names(reslentab)[which(!names(reslentab) %in% names(reslentab)[realntype])]) )
				
				for(k in 1:length(bresi)){
					res[[bresi[k]]] = rep(NA,as.numeric(names(reslentab)[realntype]))
				}
				
				res2 = matrix(unlist(res),ncol=as.numeric(names(reslentab)[realntype]),byrow=TRUE)
				#1/2/17 also fixed this so that index is not just one beyond the last bad index but the first good index
				realindex = which(reslen %in% as.numeric(names(reslentab)[realntype]))[1]
				colnames(res2) =  names(res[[realindex]])

			} else {
				
				res2 = matrix(unlist(res),ncol=as.numeric(names(reslentab)),byrow=TRUE)
				colnames(res2) = names(res[[1]])
				
			}
			
			
		} else {
			#for res as a matrix
			lres=ncol(res)
			res2 = res
		}
	}
	


if(fstat==TRUE){
	#number of covariate columns, last column is pvaue from full vs null model LRT
	ncov = ncol(res2)-1  
	result=data.frame(cur_methmat[,c('chr','start','end','strand')],beta = res2[ ,c(1:(ncov/4))],se=res2[,c((ncov/4+1):(2*ncov/4))],pvalue=res2[ , c((2*ncov/4+1):(3*ncov/4))], pfstat.type=res2[,(3*ncov/4+1):ncov], pfstat=res2[,ncov+1], meth.diff=pm.mean.diff,weighted.meth.diff=mom.mean.diff, stringsAsFactors=F)

} else {
	ncov = ncol(res2)
	result=data.frame(cur_methmat[,c('chr','start','end','strand')],beta = res2[ ,c(1:(ncov/3))],se=res2[,c((ncov/3+1):(2*ncov/3))],pvalue=res2[ , c((2*ncov/3+1):ncov)], meth.diff=pm.mean.diff, weighted.meth.diff=mom.mean.diff, stringsAsFactors=F)

}



return(result)

}


