# includes functions: compAllocVtoCompAlloc, ls.mvrm, label.count, plot.generic, plot.mvrm   
compAllocVtoCompAlloc <- function(G,p,compAllocV){
    compAlloc<-array(-19,p*(p-1)/2)     
    move <- 0
    for (k in 0:(G-1)){
        temp <- 0
        for (i in 1:(p-1)){
            for (j in (i+1):p){
                if (compAllocV[i] == compAllocV[j] && compAllocV[i] == k){ 
                    compAlloc[(i-1)*p+j-(i)*(i+1)/2] <- move
                    temp <- temp+1
		        }
            }
        }
        if (temp > 0) move<-move+1
    }      
    move2<-0
    for (k in 0:(G-2)){
        for (l in (k+1):(G-1)){
            temp <- 0
            for (i in 1:(p-1)){
                for (j in (i+1):p){                 
			        if (((compAllocV[i] == k && compAllocV[j] == l) || (compAllocV[i] == l && compAllocV[j] == k))){ 
                            compAlloc[(i-1)*p+j-(i)*(i+1)/2] <- move+move2
                            temp<-temp+1		            
	                }
                }
            }
            if (temp > 0) move2<-move2+1
        }
    }
    return(compAlloc)
}

ls.mvrm <- function(x){    
    file1 <- paste(x$DIR,"BNSP.nu.ls.txt",sep="") 
    file2 <- paste(x$DIR,"BNSP.eta.ls.txt",sep="")
    file3 <- paste(x$DIR,"BNSP.nmembers.ls.txt",sep="") 
    file4 <- paste(x$DIR,"BNSP.clusters.txt",sep="")    
    nu<-mvrm2mcmc(x,"nu")
    eta<-mvrm2mcmc(x,"eta")
    nmembers<-mvrm2mcmc(x,"nmembers")
    nparams<-x$LGc+(x$LGc+1)+1 #gamma + eta + nmembers
    nu <- nu[,c(mapply(seq,seq(1,x$LGc),length.out=x$H,by=x$LGc))]
    eta <- eta[,c(mapply(seq,seq(1,x$LGc+1),length.out=x$H,by=x$LGc+1))]
    mat<-cbind(nu,eta,nmembers)
    mcmc.pars <- array(data = mat, dim = c(x$nSamples, x$H, nparams))
    compAlloc<-mvrm2mcmc(x,"compAlloc")
    compAllocV<-mvrm2mcmc(x,"compAllocV")
    ifelse(x$G == 1, z<-compAlloc+1, z<-compAllocV+1)
    probs<-mvrm2mcmc(x,"probs")
    ifelse(x$G ==1, DIM<-x$H, DIM<-x$G)
    pr <- array(data = probs, dim = c(x$d, x$nSamples, DIM))
    #pr <- array(data = probs, dim = c(x$p, x$nSamples, x$H))
    p <- aperm(pr, c(2,1,3))
    if (x$G == 1){ 
        ls <- label.switching(method = "STEPHENS", p = p, mcmc = mcmc.pars, z = z, K = x$H)
        mcmc.pars.per<-permute.mcmc(mcmc.pars, ls$permutations$STEPHENS)$output
    }
    if (x$G > 1){   
	    ls <- label.switching(method = "STEPHENS", p = p, z = z, K = x$G)
        #Permute component allocations for the variables
        compAllocVPer<-compAllocV
        for (sw in 1:x$nSamples)
            for (k in 0:(x$G-1)) 
                compAllocVPer[sw,compAllocV[sw,]==ls$permutations[[1]][sw,k+1]-1]<-k
        #Permute component allocations for the correlations
        compAllocPer<-matrix(0,ncol=x$d,nrow=x$nSamples)
        for (sw in 1:x$nSamples){
            if (sum(compAllocVPer[sw,]==compAllocV[sw,])==x$p) compAllocPer[sw,]<-compAlloc[sw,] 
            else compAllocPer[sw,]<-compAllocVtoCompAlloc(x$G,x$p,compAllocVPer[sw,])
        }
        #Find clustering
        tab <- table(apply(compAllocPer+1, 1, paste, collapse=", "))
        ls$clusters <- as.numeric(unlist(strsplit(names(which.max(tab)), ', ')))
        
        #From permutation matrix of variables to permuation matrix of the correlations
        fpm<-matrix(-99,nrow=x$nSamples,ncol=x$H)
        for (sw in 1:x$nSamples){
            b<-compAlloc[sw,]
            a<-compAllocPer[sw,]
            fpv<-array(-99,x$H)
            if ((max(a,b)+2)<=x$H) fpv[(max(a,b)+2):x$H]<-c((max(a,b)+2):x$H)
            for (k in 0:max(a,b))
                fpv[k+1]<-a[min(which(b==k))]+1
            fpm[sw,] <- fpv
        }
        #With fpm, permute parameters
        mcmc.pars.per<-permute.mcmc(mcmc.pars, fpm)$output               
	}
    nu<-NULL
    if (x$LGc > 0){ 
        nu<-mcmc.pars.per[,,1:x$LGc]
        nu <-matrix(nu,nrow=x$nSamples)    
        nu <- nu[,c(mapply(seq,seq(1,x$H),length.out=x$LGc,by=x$H))]
	}
    eta<-mcmc.pars.per[,,(1+x$LGc):(nparams-1)]
    eta <-matrix(eta,nrow=x$nSamples)    
    eta <- eta[,c(mapply(seq,seq(1,x$H),length.out=(x$LGc+1),by=x$H))]    
    nmembers<-mcmc.pars.per[,,nparams]
    write.table(nu,file=file1,row.names = FALSE,col.names = FALSE)
    write.table(eta,file=file2,row.names = FALSE,col.names = FALSE)
    write.table(nmembers,file=file3,row.names = FALSE,col.names = FALSE)
    write.table(ls$clusters,file=file4,row.names = FALSE,col.names = FALSE, quote=FALSE)   
}

label.count <- function(term,labels,counts,model){
    if (is.character(term)){        
        int.label1<-int.label2<-grep(term,labels,fixed=TRUE)
        k<-1
        while (!length(int.label1)==1 && (nchar(term)>k)){
            int.label1<-grep(substr(term,1,nchar(term)-k),labels,fixed=TRUE)
            k<-k+1
	    }	    
        k<-1
        while (!length(int.label2)==1 && (nchar(term)>k)){
            int.label2<-grep(substr(term,k,nchar(term)),labels,fixed=TRUE)
            k<-k+1
	    }	    	    
	    if ((!length(int.label1) == 1) || (!length(int.label2) == 1) || !(int.label1==int.label2))
	        stop(cat(cat("`term` in", model, "model should be an integer between 1 and", length(labels),"or one of: "), cat(labels,sep=", ",fill=TRUE)),call. = FALSE)
	    int.label<-int.label1
    }
    if (is.numeric(term)) int.label<-term
    if ((!length(int.label) == 1) || (int.label > length(labels)))  
        stop(cat(cat("`term` in", model, "model should be an integer between 1 and", length(labels),"or one of: "), cat(labels,sep=", ",fill=TRUE)),call. = FALSE)
    count<-counts[int.label]
    return(c(int.label,count))
}	

plot.generic <- function(mvrmObj, MEAN, STDEV, CORM, CORS, DEP, SCALE,                          
                         response, intercept, grid, 
                         centre, quantiles, static, contour, centreEffects, plotOptions, 
                         int.label, count, vars, label, is.D, formula.term, which.Spec, assign,
                         knots, storeMeanVector, storeIndicator, data, plotEmptyCluster){
    p<-mvrmObj$p    
    for (i in 1:length(formula.term)){
        if (!grepl("knots",formula.term[i]) && which.Spec[i] > 0){
    	    formula.term[i]<-substr(formula.term[i],1,nchar(formula.term[i])-1)
            formula.term[i]<-paste(formula.term[i],", knots=knots[[",which.Spec[i],"]])",sep="")
	    }	
    }
    V<-NULL
    for (k in 1:length(int.label)) V<-c(V,which(assign %in% int.label[k]))
	if (STDEV) V <- V - 1
    if ((intercept || sum(is.D)) && MEAN) V<-c(1,V)    
    indeces<-V
    if (STDEV) indeces<-c(1,indeces+1)
    if (count==1){
	    if (!is.D){
	        min1<-min(with(data,eval(as.name(vars[1]))))
		    max1<-max(with(data,eval(as.name(vars[1]))))
		    newR1<-seq(min1,max1,length.out=grid)
            newData<-data.frame(newR1)
            colnames(newData)<-vars
            if (sum(!is.na(formula.term))){
                dm<-DM(formula=reformulate(formula.term),data=newData,n=NROW(newData),knots=knots,
                       meanVector=storeMeanVector[indeces],indicator=storeIndicator[indeces],centre=TRUE)
                Ds<-dm$X
                Dynamic<-dm$Dynamic
                amplitude<-dm$DSP[1]
                startSin<-dm$DSP[2]
                harmonics<-dm$DSP[3]
                #print(dim(Ds))
                #print(startSin)
		    }            
            if (!sum(!is.na(formula.term))) Ds<-matrix(1,1,1)
            if (((!1%in%V)||STDEV) && sum(!is.na(formula.term))) Ds<-Ds[,-1]
        }
        if (is.D){
            newData<-data.frame(sort(unique(with(data,eval(as.name(vars[1]))))))
            colnames(newData)<-vars
            Ds<-DM(formula=reformulate(formula.term),data=newData,n=NROW(newData),knots=NULL,
                   meanVector=storeMeanVector[indeces],indicator=storeIndicator[indeces],centre=TRUE)$X
            if (STDEV) Ds<-Ds[,-1]        
	    }
    }
    if (count==2){
	    if (sum(is.D)){
	        dv<-which(is.D==1)
	        cv<-which(is.D==0)
	        min1<-min(with(data,eval(as.name(vars[cv]))))
			max1<-max(with(data,eval(as.name(vars[cv]))))
            newR1<-seq(min1,max1,length.out=grid)
            newR2<-unique(with(data,eval(as.name(vars[dv]))))
            newData<-expand.grid(newR1,newR2)
            colnames(newData)<-vars[c(cv,dv)]
            #which.Spec <- which.Spec[which.Spec > 0] 
            #whichKnots <- which.Spec
            #Dstar<-knots[[whichKnots]]
            #Ds<-DM(formula=reformulate(formula.term),data=newData,n=NROW(newData),knots=Dstar,
            #       meanVector=storeMeanVector[indeces],indicator=storeIndicator[indeces],centre=TRUE)$X
            Ds<-DM(formula=reformulate(formula.term),data=newData,n=NROW(newData),knots=knots,
                   meanVector=storeMeanVector[indeces],indicator=storeIndicator[indeces],centre=TRUE)$X
            if (STDEV) Ds<-Ds[,-1]
	    }
        if (sum(is.D)==0){
		    min1<-min(with(data,eval(as.name(vars[1]))))
			min2<-min(with(data,eval(as.name(vars[2]))))
			max1<-max(with(data,eval(as.name(vars[1]))))
			max2<-max(with(data,eval(as.name(vars[2]))))
            newR1<-seq(min1,max1,length.out=grid)
            newR2<-seq(min2,max2,length.out=grid)
            newData<-expand.grid(newR1,newR2)
            colnames(newData)<-vars
            #whichKnots <- which.Spec
            #Dstar<-knots[[whichKnots]]
            #Ds<-DM(formula=reformulate(formula.term),data=newData,n=NROW(newData),knots=Dstar,
            #       meanVector=storeMeanVector[indeces],indicator=storeIndicator[indeces],centre=TRUE)$X
            Ds<-DM(formula=reformulate(formula.term),data=newData,n=NROW(newData),knots=knots,
                   meanVector=storeMeanVector[indeces],indicator=storeIndicator[indeces],centre=TRUE)$X
            if (STDEV) Ds<-Ds[,-1]
		}
    }  
    if (MEAN){   
		dim2<-1
		if (CORM) dim2<-mvrmObj$H 
		#print("jbf")
		#print(dim(Ds))
        fit<-array(0,dim=c(mvrmObj$nSamples,dim2,NROW(Ds)))
        #fit<-matrix(0,nrow=mvrmObj$nSamples,ncol=NROW(Ds))
        #meanReg<-array(0,dim=c(x$nSamples,x$H,x$LUT)) 
	    if (CORM==0 && DEP==0){
			etaFN <- file.path(paste(mvrmObj$DIR,"BNSP.beta.txt",sep=""))
			how.many <- p*(mvrmObj$LG+1)
			multiplier <- mvrmObj$LG+1
			#print(dim(Dynamic))
			if (!is.null(Dynamic)){
				harmonicsFN <- file.path(paste(mvrmObj$DIR,"BNSP.Hbeta.txt",sep=""))
			    how.many2 <- 2*harmonics 
			} 
		}else if (CORM==1){
		    if (mvrmObj$H==1) etaFN <- file.path(paste(mvrmObj$DIR,"BNSP.eta.txt",sep=""))
		    if (mvrmObj$H > 1){
			    etaFN <- file.path(paste(mvrmObj$DIR,"BNSP.eta.ls.txt",sep=""))
                clusFN <- paste(mvrmObj$DIR,"BNSP.clusters.txt",sep="")
                if (!file.exists(etaFN) || !file.exists(clusFN)) 
                    quiet(ls.mvrm(mvrmObj))                  
                clus<-array(unlist(read.table(clusFN)))
                tabs<-table(clus)
		        labs<-array(0,mvrmObj$H)
                labs[as.numeric(names(tabs))]<-tabs 
                V <- rep(V,mvrmObj$H) + rep(seq(0,max(V)*(mvrmObj$H-1),by=max(V)),each=length(V))
			}
		    how.many <- (mvrmObj$LGc+1)*mvrmObj$H
		    multiplier<-0
		}else{
			etaFN <- file.path(paste(mvrmObj$DIR,"BNSP.psi.txt",sep=""))
			how.many <- p*p*mvrmObj$LK
			multiplier <- mvrmObj$LK		
		}			
        eFile<-file(etaFN,open="r")
        if (!is.null(Dynamic)) 
            eFile2<-file(harmonicsFN,open="r")
	    for (i in 1:mvrmObj$nSamples){
            eta<-scan(eFile,what=numeric(),n=how.many,quiet=TRUE)
            if (!is.null(Dynamic)){
			    etaHar<-scan(eFile2,what=numeric(),n=how.many2,quiet=TRUE)						            
			    #ifelse(sum(Ds[,1] == 1) == dim(Ds)[1], addIntercept <- 1, addIntercept <- 0)
			    #print(addIntercept)
			    #Ds[,-c(1:startSin)] <- 0
			    Ds[,c((startSin+1):(startSin+amplitude+1))] <- 0
			    for (l in 1:how.many2)
			        Ds[,c((startSin+1):(startSin+amplitude+1))] <- Ds[,c((startSin+1):(startSin+amplitude+1))] + etaHar[l] * Dynamic[,(1+(amplitude+1)*(l-1)):((amplitude+1)*l)]
			    #if (addIntercept) Ds <- cbind(1,Ds)
			    #print(Ds)
			    #print(dim(Ds))
			}
			#print(dim(Ds))
			
			#print("a")
			#print(dim(Ds))
			#print(length(c(eta[V+(response-1)*multiplier])))
			
            fit[i,,] <- matrix(c(eta[V+(response-1)*multiplier]),byrow=TRUE,ncol=dim(Ds)[2]) %*% t(as.matrix(Ds))            
            
            #print("b")
            #if (i==1) {print(fit[i,,]);print(matrix(c(eta[V+(response-1)*multiplier]),nrow=1));print(t(as.matrix(Ds)))}
            #fit[i,,]<-as.matrix(Ds)%*%matrix(c(eta[V+(response-1)*multiplier]))
            #meanReg[i,,] <- eta[i,,]%*%t(uXc)            
            if (CORM && mvrmObj$FT==1) fit[i,,] <- tanh(fit[i,,])                        
            if (centreEffects) fit[i,,]<-fit[i,,]-mean(fit[i,,])
	    }
	    close(eFile)
	    if (!is.null(Dynamic)) close(eFile2)
    }
	if (STDEV){	
	    fit<-array(0,dim=c(mvrmObj$nSamples,1,NROW(Ds)))
	    if (CORS==0 && SCALE==0){
	        alphaFN <- file.path(paste(mvrmObj$DIR,"BNSP.alpha.txt",sep=""))
            sigma2FN <- file.path(paste(mvrmObj$DIR,"BNSP.sigma2.txt",sep=""))
            how.many1 <- p*mvrmObj$LD
            how.many2 <- p
		}
		if (CORS==1){
			alphaFN <- file.path(paste(mvrmObj$DIR,"BNSP.omega.txt",sep=""))
            sigma2FN <- file.path(paste(mvrmObj$DIR,"BNSP.sigma2R.txt",sep=""))		
            how.many1 <- mvrmObj$LDc
            how.many2 <- 1
		}
		if (SCALE==1){
		    alphaFN <- file.path(paste(mvrmObj$DIR,"BNSP.psi.txt",sep=""))
            sigma2FN <- file.path(paste(mvrmObj$DIR,"BNSP.phi2.txt",sep=""))
            how.many1 <- p*mvrmObj$LC
            how.many2 <- p
		}
        aFile<-file(alphaFN,open="r")	    
        s2File<-file(sigma2FN,open="r")
        for (i in 1:mvrmObj$nSamples){
            alpha<-scan(aFile,what=numeric(),n=how.many1,quiet=TRUE)
            
		    #ifelse(intercept || !is.na(formula.term), 
		    #       s2<-scan(s2File,what=numeric(),n=how.many2,quiet=TRUE)[response], s2<-1)            
            #if (!is.na(formula.term) && intercept)
            #    fit[i,,]<-s2*exp(as.matrix(Ds)%*%matrix(c(alpha[V+(response-1)*(how.many1/p)])))
            #if (!is.na(formula.term) && !intercept)
            #    fit[i,,]<-exp(as.matrix(Ds)%*%matrix(c(alpha[V+(response-1)*(how.many1/p)])))
            #if (is.na(formula.term)) fit[i,,]<-s2
            
            ifelse(intercept, s2<-scan(s2File,what=numeric(),n=how.many2,quiet=TRUE)[response], s2<-1)
            fit[i,,]<-s2*exp(as.matrix(Ds)%*%matrix(c(alpha[V+(response-1)*(how.many1/p)])))
            
            if (!SCALE) fit[i,,]<-sqrt(fit[i,,])
            if (centreEffects) fit[i,,]<-fit[i,,]/mean(fit[i,,])
	    }
	    close(aFile)
	    close(s2File)  
	}
	if (count==1 && !is.D){
	    newX<-unlist(newData)
	    centreM<-drop(apply(fit,c(2,3),centre))
		if (sum(!is.na(formula.term)))	        
	        dataM<-data.frame(rep(newX,dim(fit)[2]),c(t(centreM)))	            			
	    if (!sum(!is.na(formula.term)))
	        dataM<-data.frame(expand.grid(newX,centreM))	            	
	    nms <- c(vars,"centreM")
	    if (!is.null(quantiles)){
			nms<-c(nms,"QLM","QUM")
			QM<-apply(fit,c(2,3),quantile,probs=quantiles,na.rm=TRUE)
			if (sum(!is.na(formula.term)))	        
			    dataM<-data.frame(dataM,c(t(QM[1,,])),c(t(QM[2,,])))			
			if (!sum(!is.na(formula.term)))
			    dataM<-data.frame(dataM,rep(QM[1,,],each=length(newX)),rep(QM[2,,],each=length(newX)))
		}	    
	    colnames(dataM) <- nms	    
	    if (CORM==1 && mvrmObj$H > 1){
            dataM<-data.frame(dataM,group=factor(rep(1:mvrmObj$H,each=grid)),size=rep(labs,each=grid))		       
            dataM<-dataM[order(-dataM$size, dataM[,1]),]            
            labs2<-labs
            if (!plotEmptyCluster) {
				dataM<-dataM[dataM$size > 0,]
				labs2<-labs[labs>0]
			}            
        }        	    	    	    	    	    
	    #dataRug<-data.frame(b=with(data,eval(as.name(vars))))	    
	    plotElM<-c(geom_line(aes_string(x=as.name(vars),y=centreM),col=4,alpha=0.5),
                   #geom_rug(mapping=aes(x=dataRug$b),data=dataRug,alpha=0.3),
                   list(ylab(label)))        
        if (!is.null(quantiles))
            plotElM<-c(geom_ribbon(data=dataM,aes_string(x=vars, ymin="QLM", ymax="QUM"),alpha=0.2),plotElM)			    	                    
        if (CORM==1 && mvrmObj$H > 1){
	        plotElM<-c(geom_line(aes_string(x=as.name(vars),y="centreM",group="group",linetype="group"),col=4,alpha=0.5),
                       #geom_rug(mapping=aes(x=dataRug$b),data=dataRug,alpha=0.3),
                       list(ylab(label)))        
            if (!is.null(quantiles))
                plotElM<-c(geom_ribbon(data=dataM,aes_string(x=vars,ymin="QLM",ymax="QUM",group="group"),alpha=0.2),plotElM)			    	            
        }           
        ggM<-ggplot(data=dataM)
	    plotM<-ggM + plotElM + plotOptions        
        if (CORM==1 && mvrmObj$H > 1){
			cns<-noquote(unlist(lapply(seq(1,p-1),function (k) paste(rep(k,p-k),seq(k+1,p),sep=""))))
			cgsa<-noquote(paste(lapply(1:length(labs2), function(k)cns[clus==k]))) 		
            plotM <- ggM + plotElM + guides(fill=FALSE) +                 
                     scale_linetype_manual(name = "correlation groups: ",                      
                     labels = cgsa, 
                     values=seq(1:length(labs2))) + 
                     theme(legend.position = "bottom") +
                     plotOptions 
		}
	    	    
	}     
	if (count==1 && is.D){	
		lvs<-levels(with(data,eval(as.name(vars[1]))))
		if (is.null(lvs)) lvs<-sort(unique(with(data,eval(as.name(vars[1])))))
		                          #sort(unique(with(data,eval(as.name(vars[1])))))
	    df<-data.frame(x=rep(lvs,each=mvrmObj$nSamples),y=c(fit))
        plotElM<-c(geom_boxplot(),list(xlab(label),ylab("")))
        #ggM<-ggplot(data=df,aes(x=factor(df$x),y=df$y))
        ggM<-ggplot(df,aes(factor(df[,1]),df[,2]))
#	    ggM<-ggplot(data=dataM,aes_string(x=vars[1],y=vars[2],z=centreM))    

	    plotM<-ggM + plotElM + plotOptions
	}     
	if (count==2 && sum(is.D)==1){	    
	    centreM<-drop(apply(fit,c(2,3),centre))	
		disc.var<-vars[which(is.D==1)]
		cont.var<-vars[which(is.D==0)]
		dataM<-data.frame(newData,centreM)
		nms <- c(colnames(newData),"centreM")
		if (!is.null(quantiles)){ 			
			nms<-c(nms,"QLM","QUM")
			QM<-apply(fit,c(2,3),quantile,probs=quantiles,na.rm=TRUE)
			dataM<-data.frame(dataM,c(t(QM[1,,])),c(t(QM[2,,])))					    		    
		}
		colnames(dataM) <- nms
		
		#DG<-data.frame(b=with(data,eval(as.name(cont.var))),c=with(data,eval(as.name(disc.var))))
		
		#plotElM<-c(geom_line(aes_string(x=as.name(vars),y="centreM",group="group",linetype="group"),col=4,alpha=0.5),
        #               geom_rug(mapping=aes(x=dataRug$b),data=dataRug,alpha=0.3),
        #               list(ylab(label)))        
        
        #if (!is.null(quantiles))
        #    plotElM<-c(geom_ribbon(data=dataM,aes_string(x=vars,ymin="QLM",ymax="QUM",group="group"),alpha=0.2),plotElM
		
		plotElM<-c(geom_line(aes_string(x=cont.var,y=centreM,
		           group=disc.var,linetype=disc.var,col=disc.var),alpha=0.5),
                   #geom_rug(mapping=aes(x=DG$b,group=c,linetype=c),data=DG,alpha=0.3),
                   list(ylab(label[which(is.D==0)])))
        
        if (!is.null(quantiles))
            plotElM<-c(geom_ribbon(data=dataM,aes_string(x=cont.var,ymin="QLM",ymax="QUM",
                                   group=disc.var,col=disc.var),alpha=0.2),plotElM)
		
		ggM<-ggplot(data=dataM)
	    plotM<-ggM + plotElM + plotOptions
    }
    if (count==2 && sum(is.D)==0){
	    centreM<-drop(apply(fit,c(2,3),centre))
	    if (contour==1){	        
	        dataM<-data.frame(newData,centreM)
		    nms <- c(colnames(newData),"centreM")		    
		    colnames(dataM) <- nms
		    ggM<-ggplot(data=dataM,aes_string(x=vars[1],y=vars[2],z=centreM))    
	        level<-1
	        plotM<-ggM + geom_contour(data=dataM,aes(colour = stat(level))) + plotOptions
	    }
		if (static && contour==0){
		    defaultList<-list(x=as.numeric(newR1),y=as.numeric(newR2),z=matrix(centreM,length(newR1),length(newR2)),colvar=matrix(centreM,length(newR1),length(newR2)))
            along="xy";
            space=0.6;
            if (MEAN) optionalList<-list(xlab=vars[1],ylab=vars[2],zlab=label,along=along,space=space,add=FALSE,bty="g",main=paste("mean of", mvrmObj$varsY[[response]])) 
            if (STDEV) optionalList<-list(xlab=vars[1],ylab=vars[2],zlab=label,along=along,space=space,add=FALSE,bty="g",main=paste("stdev of", mvrmObj$varsY[[response]]))
            allOptions<-c(defaultList,plotOptions,optionalList)
            do.call("ribbon3D",allOptions[!duplicated(names(allOptions))])
		}
        if (!static && contour==0){
            a<-as.matrix(cbind(newData,centreM))
            if (is.null(plotOptions$col)) col=rainbow(16,2/3)
            else col<-plotOptions$col
            plotCentreM <- centreM
            if (min(centreM)<=0) plotCentreM <- centreM + abs(min(centreM)) + 1
            ra<-ceiling(length(col)*plotCentreM/max(plotCentreM))
            defaultList<-list(x=a,col=col[ra])
            if (MEAN) optionalList<-list(size=0.4,bg=1,axisLabels=c(vars[1],label,vars[2]),main=paste("mean of", mvrmObj$varsY[[response]]))
            if (STDEV) optionalList<-list(size=0.4,bg=1,axisLabels=c(vars[1],label,vars[2]),main=paste("stdev of", mvrmObj$varsY[[response]]))
            allOptions<-c(defaultList,plotOptions,optionalList)
            do.call("scatterplot3js",allOptions[!duplicated(names(allOptions))])
		}
    }  
    if (contour==1)
		return(plotM)
}

plot.mvrm <- function(x, model, term, response, response2, 
                      intercept = TRUE, grid = 30, centre = mean, 
                      quantiles = c(0.1, 0.9), contour = TRUE, static = TRUE, 
                      centreEffects = FALSE, plotOptions = list(), nrow, ask = FALSE,
                      plotEmptyCluster = FALSE, combine = FALSE, ...)
{
    oldpar <- NULL
    on.exit(par(oldpar))
    #
	if (!is.function(centre)) stop("centre must be a function, usually mean or median")
	if (length(quantiles)==1) quantiles <- c(quantiles,1-quantiles)
    if (length(quantiles) > 2) stop("up to two quantiles")
    if (!is.null(quantiles)) quantiles <- unique(sort(quantiles))
    grid<-round(grid)
    if (missing(response)) response <- c(1:x$p)
    if (missing(response2)) response2 <- c(1:x$p)
    if (max(response) > x$p) stop("argument response exceeds the number of responses")
    #
    if (missing(model)) model<-"all"
    MEAN <- 0; STDEV <- 0; CORM <- 0; CORS <- 0; DEP <- 0; SCALE <- 0
    if ((model=="all" || model=="mean") && (x$NG > 0)) MEAN <- 1
    if ((model=="all" || model=="stdev") && (x$ND > 0)) STDEV <- 1
    if ((model=="all" || model=="corm") && x$p > 1 && x$LUT > 1) CORM <- 1
    if ((model=="all" || model=="cors") && x$p > 1 && x$LUT > 1) CORS <- 1
    if ((model=="all" || model=="dep") && (x$NK > 0)) DEP <- 1
    if ((model=="all" || model=="scale") && (x$NC > 0)) SCALE <- 1
    #
    if (MEAN==0 && STDEV==0 && CORM==0 && CORS==0 && DEP==0 && SCALE==0) stop("no terms to plot");    
    termM <- NULL; termSD <- NULL; termMC<-NULL; termSDC<-NULL; termD<-NULL; termSC<-NULL 
    countM <- NULL; countSD <- NULL; countMC <- NULL; countSDC <- NULL; countD <- NULL; countSC <- NULL
    int.labelM <- NULL; int.labelSD <- NULL; int.labelD <- NULL; int.labelSC <- NULL
    if (missing(term)){
		if (MEAN) {termM<-1:x$NG; countM<-x$countX; int.labelM <- 1:x$NG}
		if (STDEV) {termSD<-1:x$ND; countSD<-x$countZ; int.labelSD <- 1:x$ND}
		if (DEP) {termD<-1:(x$NK-1); countD<-x$countC[termD]; int.labelD <- termD}
		if (SCALE) {termSC<-1:x$NC; countSC<-x$countW; int.labelSC <- 1:x$NC}
	}else{
		if (MEAN){
			termM<-term
			for (i in termM){
				#print(i)
	  		    lcX <- label.count(i,x$labelsX,x$countX,"mean")
			    int.labelM[i]<-lcX[1]
			    countM[i]<-lcX[2]
			}	
		}
		if (STDEV){
			termSD<-term
			for (i in termSD){
	  		    lcZ <- label.count(i,x$labelsZ,x$countZ,"stdev")
			    int.labelSD[i]<-lcZ[1]
			    countSD[i]<-lcZ[2]
			}			
		}
		if (DEP){
			termD<-term
			for (i in termD){
	  		    lcC <- label.count(i,x$labelsC,x$countC,"autoregressive")
			    int.labelD[i]<-lcC[1]
			    countD[i]<-lcC[2]
			}			
		}
		if (SCALE){
			termSC<-term
			for (i in termSC){
	  		    lcW <- label.count(i,x$labelsW,x$countW,"scale")
			    int.labelSC[i]<-lcW[1]
			    countSC[i]<-lcW[2]
			}			
		}
	}	
	if (CORM) {termMC<-1; countMC<-x$countXc[termMC]} 
	if (CORS) {termSDC<-1; countSDC<-x$countZc[termSDC]}
    #
    if ((contour == 0) && 
        (sum(c(countM,countSD,countMC,countSDC,countD,countSC)>1) > 0) && 
        (length(c(countM,countSD,countMC,countSDC,countD,countSC)) > 1)){
        warning("for 3-d plots, specify one model, one term and one response at a time")
        contour <- 1
	}
	#
	if (combine && length(term) == 1) combine <- FALSE
	if (combine) termM <- termSD <- DEP <- SCALE <- 1
    #
    my_plots <- list()
    count <- 1
    for (r in response){
		if (MEAN)
        for (i in termM){
			#lcX <- label.count(i,x$labelsX,x$countX,"mean")
            #int.labelX <- lcX[1]
            #countX <- lcX[2]
            if (!combine){
				int.labelX <- int.labelM[i]                
                countX <- countM[i]
                varsArg <- x$varsX[[int.labelX]]
                is.DArg <- x$is.Dx[[int.labelX]]
                which.SpecArg <- x$which.SpecX[[int.labelX]]
		    }else{
			    int.labelX <- term
			    varsArg <- unique(unlist(x$varsX[int.labelX]))
			    countX <- length(varsArg)
			    is.DArg <- unlist(x$is.Dx[int.labelX])
			    if (!sum(is.DArg)) is.DArg <- 0
			    which.SpecArg <- unlist(x$which.SpecX[int.labelX])			    
			    if ((length(unique(varsArg)) > 1) && (sum(!is.DArg) > 1)) stop("can't combine more than 2 continuous terms; consider interaction terms")	
			    #print(c(int.labelX))
			    #print(c(countX))
			    #print(c(varsArg))
			    #print(c(is.DArg))
			    #print(c(which.SpecArg))
			}            
            if (countX==1) contour <- 1
            plotOptions2 <- list(ggtitle(paste("mean of", x$varsY[r])),plotOptions)            
            my_plots[[count]] <- plot.generic(mvrmObj=x, MEAN=1, STDEV=0, CORM=0, CORS=0, DEP=0, SCALE=0,                                               
                                              response=r, intercept=intercept, grid=grid,
                                              centre=centre, quantiles=quantiles, static=static, 
                                              contour=contour,
                                              centreEffects=centreEffects, plotOptions=plotOptions2,                                              
                                              int.label=int.labelX, count=countX, vars=varsArg,                                               
                                              label=x$labelsX[int.labelX], is.D=is.DArg, 
                                              formula.term=x$formula.termsX[int.labelX],                                               
                                              which.Spec=which.SpecArg, assign=x$assignX, 
                                              knots=x$Xknots, storeMeanVector=x$storeMeanVectorX,                                              
                                              storeIndicator=x$storeIndicatorX, data=x$data,
                                              plotEmptyCluster=plotEmptyCluster)
            if ((ask==TRUE) && (length(my_plots)>0)) print(my_plots[[count]])
            if (count==1)
                oldpar <- c(oldpar, par(ask=ask))
            count <- count + 1
        }
        if (STDEV)
        for (i in termSD){
			#lcZ <- label.count(i,x$labelsZ,x$countZ,"stdev")
            #int.labelZ <- lcZ[1]
            #countZ <- lcZ[2]
            if (!combine){
				int.labelZ <- int.labelSD[i]                
                countZ <- countSD[i]
                varsArg <- x$varsZ[[int.labelZ]]
                is.DArg <- x$is.Dz[[int.labelZ]]
                which.SpecArg <- x$which.SpecZ[[int.labelZ]]
		    }else{
			    int.labelZ <- term
			    countZ <- 2
			    varsArg <- unlist(x$varsZ[int.labelZ])
			    is.DArg <- unlist(x$is.Dz[int.labelZ])
			    which.SpecArg <- unlist(x$which.SpecZ[int.labelZ])			    
			}
            if (combine && (!sum(is.DArg)==1)) stop("can only combine a continuous with a discrete term")            
            if (countZ==1) contour <- 1            
            plotOptions2 <- list(ggtitle(paste("stdev of", x$varsY[[r]])),plotOptions)
            if (x$LUT > 1) plotOptions2 <- list(ggtitle(paste("innov stdev of", x$varsY[[r]])),plotOptions)            
            my_plots[[count]] <- plot.generic(x, MEAN=0, STDEV=1, CORM=0, CORS=0, DEP=0, SCALE=0,
                                              response=r, intercept=intercept, grid=grid,
                                              centre=centre, quantiles=quantiles, static=static, 
                                              contour=contour,
                                              centreEffects=centreEffects, plotOptions=plotOptions2,
                                              int.label=int.labelZ, count=countZ, vars=varsArg, 
                                              label=x$labelsZ[int.labelZ], is.D=is.DArg, 
                                              formula.term=x$formula.termsZ[int.labelZ], 
                                              which.Spec=which.SpecArg, assign=x$assignZ, 
                                              knots=x$Zknots, storeMeanVector=x$storeMeanVectorZ,
                                              storeIndicator=x$storeIndicatorZ, data=x$data,
                                              plotEmptyCluster=plotEmptyCluster)
            if ((ask==TRUE) && (length(my_plots)>0)) print(my_plots[[count]])
            if (count==1)
                oldpar <- c(oldpar, par(ask=ask))
            count <- count + 1
        }
        if (SCALE)
        for (i in termSC){
			#lcW <- label.count(i,x$labelsW,x$countW,"scale")
            #int.labelW <- lcW[1]
            #countW <- lcW[2]            
            if (!combine){
				int.labelW <- int.labelSC[i] 
                countW <- countSC[i]
                varsArg <- x$varsW[[int.labelW]]
                is.DArg <- x$is.Dw[[int.labelW]]
                which.SpecArg <- x$which.SpecW[[int.labelW]]
		    }else{
			    int.labelW <- term
			    countW <- 2
			    varsArg <- unlist(x$varsW[int.labelW])
			    is.DArg <- unlist(x$is.Dw[int.labelW])
			    which.SpecArg <- unlist(x$which.SpecW[int.labelW])			    
			}
            if (combine && (!sum(is.DArg)==1)) stop("can only combine a continuous with a discrete term")                       
            if (countW==1) contour <- 1
            plotOptions2 <- list(ggtitle(paste("scale of", x$varsY[[r]])),plotOptions)
            my_plots[[count]] <- plot.generic(x, MEAN=0, STDEV=1, CORM=0, CORS=0, DEP=0, SCALE = 1,
                                              response=r, intercept=intercept, grid=grid,
                                              centre=centre, quantiles=quantiles, static=static, 
                                              contour=contour,
                                              centreEffects=centreEffects, plotOptions=plotOptions2,
                                              int.label=int.labelW, count=countW, vars=varsArg, 
                                              label=x$labelsW[int.labelW], is.D=is.DArg, 
                                              formula.term=x$formula.termsW[int.labelW], 
                                              which.Spec=which.SpecArg, assign=x$assignW, 
                                              knots=x$Wknots, storeMeanVector=x$storeMeanVectorW,
                                              storeIndicator=x$storeIndicatorW, data=x$data,
                                              plotEmptyCluster=plotEmptyCluster)
            if ((ask==TRUE) && (length(my_plots)>0)) print(my_plots[[count]])
            if (count==1)
                oldpar <- c(oldpar, par(ask=ask))
            count <- count + 1
        }
	}	
	if (DEP)
    	for (r1 in response){
	    	for (r2 in response2){		    
				Pair<-(r1-1)*x$p + r2
                for (i in termD){
			        #lcC <- label.count(i,x$labelsC,x$countC,"autoregressive")
                    #int.labelC <- lcC[1]
                    #countC <- lcC[2]
                    if (!combine){
				        int.labelC <- int.labelD[i]                
                        countC <- countD[i]
                        varsArg <- x$varsC[[int.labelC]]
                        is.DArg <- x$is.Dc[[int.labelC]]
                        which.SpecArg <- x$which.SpecC[[int.labelC]]
		            }else{
			            int.labelC <- term
			            countC <- 2
			            varsArg <- unlist(x$varsC[int.labelC])
			            is.DArg <- unlist(x$is.Dc[int.labelC])
			            which.SpecArg <- unlist(x$which.SpecC[int.labelC])			    
			        }
                    if (combine && (!sum(is.DArg)==1)) stop("can only combine a continuous with a discrete term")                                 
                    if (countC==1) contour <- 1
                    plotOptions2 <- list(ggtitle(paste("autoreg", x$varsY[r1], x$varsY[r2])), plotOptions)                                
                    my_plots[[count]] <- plot.generic(mvrmObj=x, MEAN=1, STDEV=0, CORM=0, CORS=0, DEP=1, SCALE=0, 
                                                      response=Pair, intercept=intercept, grid=grid,
                                                      centre=centre, quantiles=quantiles, static=static, 
                                                      contour=contour,
                                                      centreEffects=centreEffects, plotOptions=plotOptions2,
                                                      int.label=int.labelC, count=countC, vars=varsArg, 
                                                      label=x$labelsC[int.labelC], is.D=is.DArg, 
                                                      formula.term=x$formula.termsC[int.labelC], 
                                                      which.Spec=which.SpecArg, assign=x$assignC, 
                                                      knots=x$Cknots, storeMeanVector=x$storeMeanVectorC,
                                                      storeIndicator=x$storeIndicatorC, data=x$lag,
                                                      plotEmptyCluster=plotEmptyCluster)
                    if ((ask==TRUE) && (length(my_plots)>0)) print(my_plots[[count]])
                    if (count==1)
                        oldpar <- c(oldpar, par(ask=ask))
                    count <- count + 1
		        }
            }
	    }
    if (CORM){
        for (i in termMC){	
			contour <- 1		
            plotOptions2 <- list(ggtitle(paste("mean cor")),plotOptions)            
            which.Spec <- 0
            if (length(x$which.SpecXc)) which.Spec <- x$which.SpecXc[[1]] 
            is.D <- 0
            if (length(x$is.Dxc)) is.D <- x$is.Dxc[[1]]
            vars <- x$varTime
            if (length(x$varsXc)) vars <- x$varsXc[[1]]
            my_plots[[count]] <- plot.generic(mvrmObj=x, MEAN=1, STDEV=0, CORM=1, CORS=0, DEP=0, SCALE=0,
                                              response=1, intercept=intercept, 
                                              grid=grid, centre=centre, quantiles=quantiles, static=static, contour=contour,
                                              centreEffects=centreEffects, plotOptions=plotOptions2,
                                              int.label=1, count=1, vars=vars, 
                                              label=x$labelsXc[1], is.D=is.D, 
                                              formula.term=x$formula.termsXc[1], 
                                              which.Spec=which.Spec, assign=x$assignXc,
                                              knots=x$Xcknots, storeMeanVector=x$storeMeanVectorXc,
                                              storeIndicator=x$storeIndicatorXc, data=x$data,
                                              plotEmptyCluster=plotEmptyCluster)                        
            if ((ask==TRUE) && (length(my_plots)>0)) print(my_plots[[count]])
            if (count==1)
                oldpar <- c(oldpar, par(ask=ask))
            count <- count + 1
        }
	}	
	if (CORS){
        for (i in termSDC){		
			contour <- 1	
            plotOptions2 <- list(ggtitle(paste("stdev cor")),plotOptions)                        
            which.Spec <- 0
            if (length(x$which.SpecZc)) which.Spec <- x$which.SpecZc[[1]] 
            is.D <- 0
            if (length(x$is.Dzc)) is.D <- x$is.Dzc[[1]]
            vars <- x$varTime
            if (length(x$varsZc)) vars <- x$varsXc[[1]]
            my_plots[[count]] <- plot.generic(mvrmObj=x, MEAN=0, STDEV=1, CORM=0, CORS=1, DEP=0, SCALE=0, 
                                              response=1, intercept=intercept, 
                                              grid=grid, centre=centre, quantiles=quantiles, static=static, contour=contour,
                                              centreEffects=centreEffects, plotOptions=plotOptions2,
                                              int.label=1, count=1, vars=vars, 
                                              label=x$labelsZc[1], is.D=is.D, 
                                              formula.term=x$formula.termsZc[1], 
                                              which.Spec=which.Spec, assign=x$assignZc,
                                              knots=x$Zcknots, storeMeanVector=x$storeMeanVectorZc,
                                              storeIndicator=x$storeIndicatorZc, data=x$data,
                                              plotEmptyCluster=plotEmptyCluster)
            if ((ask==TRUE) && (length(my_plots)>0)) print(my_plots[[count]])
            if (count==1)
                oldpar <- c(oldpar, par(ask=ask))
            count <- count + 1
        }
	}	
	if (missing(nrow)){
	    if (CORM==0 && CORS==0 && DEP==0){nrow <- min(x$p,length(my_plots))}else{nrow <- ceiling(sqrt(count-1))}
	}		
	if ((ask==FALSE) && (length(my_plots)>0)) quiet(print(grid.arrange(grobs = my_plots, nrow = nrow)))
}
