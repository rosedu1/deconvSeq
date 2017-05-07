#modified diffMeth:uses glm and pvalues from glm, outputs mean methylation difference (mean(Ci/(Ci+Ti))) and weighted methylation difference (C in all samples)/(C+T in all samples) in methylDiff object as output

#example:
#design=data.frame(gender=c('M','M','M', 'F', 'F','F'), treatment=c(0,1,0,1,0,1))
#formula='gender+treatment'
#num.cores=4
# class(meth)
#[1] "methylBase"
#attr(,"package")
#[1] "methylKit"

#meth= unite(obj)
#myDiff=diffMethFromDesign(meth, design, formula, num.cores)
#design=data.frame(subject=c(1,1,2,2), aneu=c(1,0,1,0))
#formula='(1|subject)+aneu'
#**replaced glm with mixed linear model from lme4: glmer
#glmer(cbind(Cs,Ts)~(1|subject)+aneu, data=cur_data, family=binomial(link=logit))
#**subjects is the random effects.  Put the covariate of interest as the last one in the formula

#8/1/15 use glmer.  If separation occurs, try bglmer with normal prior for fixed effects, if this fails, then give NA
#Gelman et al. (2008) suggest that the input variables of a categorical regression are standardised and that the associated regression parameters are assumed independent in the prior. Gelman et al. (2008) recommend a scaled t-distribution with a single degree of freedom (scaled Cauchy) and a scale of 10 for the intercept and 2.5 for the regression parameters.
#Gelman, A. et al. (2008) The Annals of Appled Statistics 2 4 1360-1383
#the default in bglmer ie fixef.prior=null uses Gelman's prior with Fixef prior: normal(sd = c(10, 2.5), corr = 0, common.scale = FALSE)

#8/4/15
#output just the dataframe without qvalues so that we can combine them later.  For use in pmeth.R

#8/9/15
#add mixedef option so that we can use it without mixed effects and just use regular glm

#7/15/16
#add return of coefficients, standard errors and pvalues for all covariates not just the last covariate, so that we can also use this to do deconvolution. NOTE that meth.diff and weighted.meth.diff is still based on difference between treatment groups 0 and 1 (independent of design matrix)
#add fstat option to get pvalue for full model vs null model (to be used for deconvolution)
#fixed a typo in x in my_diffMethFromDesign2.R, where beta=res[,2] gives the wrong column, should have been res[,1]

#7/23/16 have trouble differentiating similar tissues. May  not be a problem when looking at different cell types (ie good for picking out blood but bad when trying to identify STA or aneurysm). Idea: also look at LRT of each cell type alone and pick out those as well for b0. ie add pfstat for model with each cell type alone

#7/27/16 return results even if there is a warning, only halt for errors

#12/22/16 add repeated measures. Cur_design is only used to make covariates part of cur_data, so just make it data rather than a design matrix
#12/26/16 bglmer is way too slow for mixed models, about 4.6s per calculation as opposed to glmer. Remove bglmer.

#1/2/17 my_diffMethFromDesign5.R fix res part since there may be more NAs than actual results for the repeated measures calculations

library(blme)
f1<-function(cutoff,rawp){sum(rawp<cutoff)/length(rawp)}

my_diffMethFromDesign<-function(cur_meth, cur_design, cur_formula, num.cores, mixedef=TRUE, fstat=TRUE){

	## Assumes you have a methylBase object of interest called 'cur_meth'
	## Pulls out the T and C data here so you don't need to keep doing it later
	#***make sure last column in design is the variable of interest


	c_data=getData(cur_meth)[ , cur_meth@numCs.index]
	t_data=getData(cur_meth)[ , cur_meth@numTs.index]
	num_base=dim(cur_meth)[[1]]
    #last_design_col=as.integer(cur_design[ , dim(cur_design)[[2]]])



	#****added by RD
	#get the indices (column number) for numCs and numTs in each set
	set1.Cs=cur_meth@numCs.index[cur_meth@treatment==1]
	set2.Cs=cur_meth@numCs.index[cur_meth@treatment==0]
	set1.Ts=cur_meth@numTs.index[cur_meth@treatment==1]
	set2.Ts=cur_meth@numTs.index[cur_meth@treatment==0]
	subst=S3Part(cur_meth) #same as getData(cur_meth)
     # calculate mean methylation change
	if(length(set1.Cs) > 1){ #if more than one subject in treatment group
		mom.meth1=100*rowMeans(subst[,set1.Cs]/subst[,set1.Cs-1],na.rm=TRUE) # get means of means, column set.Cs-1 is the total coverage column (C+T)
		pm.meth1=100*rowSums(subst[,set1.Cs])/rowSums(subst[,set1.Cs-1],na.rm=TRUE) # get weighted means, this is total Cs for all samples/total coverage for all samples
	}else{
		mom.meth1    = 100*(subst[,set1.Cs]/subst[,set1.Cs-1]) # get % methylation
		pm.meth1     = mom.meth1
	}

	if(length(set2.Cs)>1){
		mom.meth2=100*rowMeans(subst[,set2.Cs]/subst[,set2.Cs-1],na.rm=TRUE)
		pm.meth2=100*rowSums(subst[,set2.Cs])/rowSums(subst[,set2.Cs-1],na.rm=TRUE) # get weighted means
	}else{
		mom.meth2    = 100*(subst[,set2.Cs]/subst[,set2.Cs-1]) # get % methylation
		pm.meth2     = mom.meth2
	}
	pm.mean.diff=pm.meth1-pm.meth2
	mom.mean.diff=mom.meth1-mom.meth2
                          


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
			#coef_numb=length(obj$coefficients)length
			#deviance <- obj$null.deviance - obj$deviance
			#****use output of glmer, last row is covariate of interest
			coef_numb=nrow(summary(obj)$coefficients)

			#beta0<-<- obj$coefficients[1]

			## I assume the last term in the formula has the coefficient of interest,
			##and grab the wald test value
			## associated with the glm object

			#beta_last <- obj$coefficients[coef_numb,1]
			#****
            #beta_last <- summary(obj)$coefficients[coef_numb,1]
            #p.value<- summary(obj)$coefficients[ coef_numb, 'Pr(>|z|)']
			#p.value <- 1-pchisq(deviance,df=1)
			#meth_diff=100*(ilogit(sum(obj$coefficients)-beta_last)-ilogit(sum(obj$coefficients)))
			#fitted gives predicted value using original data and original scale, while predict gives response using any data and using scale before inverse link function is applied
			#obj$fitted.values = fitted.values(obj) = fitted(obj) = ilogit(predict(obj))
			#obj$linear.predictors = predict(obj)
			#***I will use the old meth diff calculations based on number Cs and Ts.  This new version gives beta_last*C in all samples/(C+T in all samples)
			#meth_diff=100*mean(obj$fitted.values-ilogit(obj$linear.predictors-beta_last*last_design_col))
			#return( c(meth_diff, beta_last,p.value) )
			#***
            
            #for mixed effects model, only fixed effects results are obtained
            beta <- summary(obj)$coefficients[,'Estimate']
            p.value <- summary(obj)$coefficients[,'Pr(>|z|)']
            se <- summary(obj)$coefficients[,'Std. Error']
            
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

	#res=t(simplify2array(parallel::mclapply(1:100, cur_function, mc.cores=num.cores)))

	## Actually do the most of the work, and get the data frame into shape in
	##the same step
	## Obviously you can up the number of cores
	#***get rid of parallel processig
	#res=t(simplify2array(parallel::mclapply(1:num_base, cur_function, mc.cores=num.cores)))
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
#		k=0
#		j=1
#		while (k ==0){
#			if(j %in% na.ind) {j=j+1}
#			else {k==1}
#		}	
#		colnames(res2)=names(res[[j]])	
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
	
	
	## Get the q values
	#cat('Qvales')
	#pvals  = fix.q.values.glm(res[ , 3],slim=TRUE)
	#***
	#remove the NA rows
	#resi = which(!is.na(res[,2]))
	#pvals  = fix.q.values.glm(res[resi , 2],slim=TRUE)


	## Make the final data frame and return it instead of methylDiff

#	x=data.frame(getData(cur_meth)[,c('id','chr','start','end','strand')],pvalue=pvals[ , 1], qvalue=pvals[,2],meth.diff=res[ , 1],stringsAsFactors=F)
#x=data.frame(getData(cur_meth)[,c('chr','start','end','strand')],pvalue=pvals[ , 1], qvalue=pvals[,2],meth.diff=res[ , 1],stringsAsFactors=F)
#***
#x=data.frame(getData(cur_meth)[resi,c('chr','start','end','strand')],pvalue=pvals[ , 1], qvalue=pvals[,2],meth.diff=pm.mean.diff[resi],weighted.meth.diff=mom.mean.diff[resi], stringsAsFactors=F)

if(fstat==TRUE){
	#number of covariate columns, last column is pvaue from full vs null model LRT
	ncov = ncol(res2)-1  
	result=data.frame(getData(cur_meth)[,c('chr','start','end','strand')],beta = res2[ ,c(1:(ncov/4))],se=res2[,c((ncov/4+1):(2*ncov/4))],pvalue=res2[ , c((2*ncov/4+1):(3*ncov/4))], pfstat.type=res2[,(3*ncov/4+1):ncov], pfstat=res2[,ncov+1], meth.diff=pm.mean.diff,weighted.meth.diff=mom.mean.diff, stringsAsFactors=F)

} else {
	ncov = ncol(res2)
	result=data.frame(getData(cur_meth)[,c('chr','start','end','strand')],beta = res2[ ,c(1:(ncov/3))],se=res2[,c((ncov/3+1):(2*ncov/3))],pvalue=res2[ , c((2*ncov/3+1):ncov)], meth.diff=pm.mean.diff,weighted.meth.diff=mom.mean.diff, stringsAsFactors=F)

}



return(result)

	# make a dataframe and return it

	## Make the new methylDiff object

	#output=new("methylDiff",x,sample.ids=cur_meth@sample.ids, assembly=cur_meth@assembly, context=cur_meth@context, treatment=cur_meth@treatment, destranded=cur_meth@destranded, resolution=cur_meth@resolution)
	#output
}


