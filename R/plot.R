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

plot.generic <- function(mvrmObj, MEAN, STDEV, CORM, CORS, DEP, term, response, intercept, grid, 
                         centre, quantiles, static, contour, centreEffects, plotOptions, 
                         int.label, count, vars, label, is.D, formula.term, which.Spec, assign,
                         knots, storeMeanVector, storeIndicator, data, plotEmptyCluster){
    p<-mvrmObj$p
    if (!grepl("knots",formula.term) && which.Spec>0){			
        z<-length(formula.term) #use this to fix both sm and s smooths
    	formula.term[z]<-substr(formula.term[z],1,nchar(formula.term[z])-1)
        formula.term[z]<-paste(formula.term[z],", knots=knots)")
	}	
    V<-which(assign %in% int.label)
	if (STDEV) V <- V - 1
    if ((intercept || sum(is.D)) & MEAN) V<-c(1,V)
    indeces<-V
    if (STDEV) indeces<-c(1,indeces+1)
    if (count==1){
	    if (!is.D){
	        min1<-min(with(data,eval(as.name(vars[1]))))
		    max1<-max(with(data,eval(as.name(vars[1]))))
		    newR1<-seq(min1,max1,length.out=grid)
            newData<-data.frame(newR1)
            colnames(newData)<-vars
            whichKnots <- which.Spec
            if (length(knots)>0 && whichKnots>0) {Dstar<-data.frame(knots=knots[[whichKnots]])}else{Dstar<-NULL}            
            if (!is.na(formula.term))
                Ds<-DM(formula=reformulate(formula.term),data=newData,n=NROW(newData),knots=Dstar,
                       meanVector=storeMeanVector[indeces],indicator=storeIndicator[indeces])$X
            if (is.na(formula.term)) Ds<-matrix(1,1,1)
            if (((!1%in%V)||STDEV) && !is.na(formula.term)) Ds<-Ds[,-1]
        }
        if (is.D){
            newData<-data.frame(sort(unique(with(data,eval(as.name(vars[1]))))))
            colnames(newData)<-vars
            Ds<-DM(formula=reformulate(formula.term),data=newData,n=NROW(newData),knots=NULL,
                   meanVector=storeMeanVector[indeces],indicator=storeIndicator[indeces])$X
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
            whichKnots <- which.Spec
            Dstar<-knots[[whichKnots]]
            Ds<-DM(formula=reformulate(formula.term),data=newData,n=NROW(newData),knots=Dstar,
                   meanVector=storeMeanVector[indeces],indicator=storeIndicator[indeces])$X
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
            whichKnots <- which.Spec
            Dstar<-knots[[whichKnots]]
            Ds<-DM(formula=reformulate(formula.term),data=newData,n=NROW(newData),knots=Dstar,
                   meanVector=storeMeanVector[indeces],indicator=storeIndicator[indeces])$X
            if (STDEV) Ds<-Ds[,-1]
		}
    }  
    if (MEAN){   
		dim2<-1
		if (CORM) dim2<-mvrmObj$H 
        fit<-array(0,dim=c(mvrmObj$nSamples,dim2,NROW(Ds)))
        #fit<-matrix(0,nrow=mvrmObj$nSamples,ncol=NROW(Ds))
        #meanReg<-array(0,dim=c(x$nSamples,x$H,x$LUT)) 
	    if (CORM==0 && DEP==0){
			etaFN <- file.path(paste(mvrmObj$DIR,"BNSP.beta.txt",sep=""))
			how.many <- p*(mvrmObj$LG+1)
			multiplier <- mvrmObj$LG+1
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
	    for (i in 1:mvrmObj$nSamples){
            eta<-scan(eFile,what=numeric(),n=how.many,quiet=TRUE)
            fit[i,,] <- matrix(c(eta[V+(response-1)*multiplier]),byrow=TRUE,ncol=dim(Ds)[2]) %*% t(as.matrix(Ds))            
            #if (i==1) {print(fit[i,,]);print(matrix(c(eta[V+(response-1)*multiplier]),nrow=1));print(t(as.matrix(Ds)))}
            #fit[i,,]<-as.matrix(Ds)%*%matrix(c(eta[V+(response-1)*multiplier]))
            #meanReg[i,,] <- eta[i,,]%*%t(uXc)            
            if (CORM && mvrmObj$out[[60]]==1) fit[i,,] <- tanh(fit[i,,])                        
            if (centreEffects) fit[i,,]<-fit[i,,]-mean(fit[i,,])
	    }
	    close(eFile)
    }
	if (STDEV){	
	    fit<-array(0,dim=c(mvrmObj$nSamples,1,NROW(Ds)))
	    if (CORS==0){
	        alphaFN <- file.path(paste(mvrmObj$DIR,"BNSP.alpha.txt",sep=""))
            sigma2FN <- file.path(paste(mvrmObj$DIR,"BNSP.sigma2.txt",sep=""))
            how.many1 <- p*mvrmObj$LD
            how.many2 <- p
		}else{
			alphaFN <- file.path(paste(mvrmObj$DIR,"BNSP.omega.txt",sep=""))
            sigma2FN <- file.path(paste(mvrmObj$DIR,"BNSP.sigma2R.txt",sep=""))		
            how.many1 <- mvrmObj$LDc
            how.many2 <- 1
		}
        aFile<-file(alphaFN,open="r")	    
        s2File<-file(sigma2FN,open="r")
        for (i in 1:mvrmObj$nSamples){
            alpha<-scan(aFile,what=numeric(),n=how.many1,quiet=TRUE)
		    ifelse(intercept || !is.na(formula.term), s2<-scan(s2File,what=numeric(),n=how.many2,quiet=TRUE)[response], s2<-1)            
            if (!is.na(formula.term))
                fit[i,,]<-sqrt(s2*exp(as.matrix(Ds)%*%matrix(c(alpha[V+(response-1)*mvrmObj$LD]))))            
            if (is.na(formula.term))
                fit[i,,]<-sqrt(s2)
            if (centreEffects) fit[i,,]<-fit[i,,]/mean(fit[i,,])
	    }
	    close(aFile)
	    close(s2File)  
	}
	if (count==1 && !is.D){
	    newX<-unlist(newData)
	    centreM<-drop(apply(fit,c(2,3),centre))
		if (!is.na(formula.term))	        
	        dataM<-data.frame(rep(newX,dim(fit)[2]),c(t(centreM)))	            			
	    if (is.na(formula.term))
	        dataM<-data.frame(expand.grid(newX,centreM))	            	
	    nms <- c(vars,"centreM")
	    if (!is.null(quantiles)){
			nms<-c(nms,"QLM","QUM")
			QM<-apply(fit,c(2,3),quantile,probs=quantiles,na.rm=TRUE)
			if (!is.na(formula.term))	        
			    dataM<-data.frame(dataM,c(t(QM[1,,])),c(t(QM[2,,])))			
			if (is.na(formula.term))
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
        ggM<-ggplot(data=df,aes(x=factor(df$x),y=df$y))
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
		           group=disc.var,linetype=disc.var),col=4,alpha=0.5),
                   #geom_rug(mapping=aes(x=DG$b,group=c,linetype=c),data=DG,alpha=0.3),
                   list(ylab(label)))
        
        if (!is.null(quantiles))
            plotElM<-c(geom_ribbon(data=dataM,aes_string(x=cont.var,ymin="QLM",ymax="QUM",
                                   group=disc.var),alpha=0.2),plotElM)
		
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
                      plotEmptyCluster = FALSE, ...)
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
    if (max(response) > x$p) stop("argument response exceeds the number of responses");
    
    if (missing(model)) model<-"all"
    MEAN <- 0; STDEV <- 0; CORM <- 0; CORS <- 0; DEP <- 0;
    if ((model=="all" || model=="mean") && (x$NG > 0)) MEAN <- 1
    if ((model=="all" || model=="stdev") && (x$ND > 0)) STDEV <- 1
#    if ((model=="all" || model=="corm") && (x$NGc > 0)) CORM <- 1
#    if ((model=="all" || model=="cors") && (x$NDc > 0)) CORS <- 1
    if ((model=="all" || model=="corm") && x$p > 1) CORM <- 1
    if ((model=="all" || model=="cors") && x$p > 1) CORS <- 1
    if ((model=="all" || model=="dep") && (x$NK > 0)) DEP <- 1
    
    if (MEAN==0 && STDEV==0 && CORM==0 && CORS==0 && DEP==0) stop("no terms to plot");    
    termM <- NULL; termSD <- NULL; termMC<-1; termSDC<-1; termD<-NULL  
    if (missing(term)) {termM<-1:x$NG; termSD<-1:x$ND; termD<-1:(x$NK-1)}
    if (!missing(term)) {termM <- termSD <- termD <- term}
    
    #
    if (((MEAN+STDEV+length(termM)+length(termSD)+length(response))>3) && contour == 0)
        warning("for 3-d plots, specify one model, one term and one response at a time")
    if (((MEAN+STDEV+length(termM)+length(termSD)+length(response))>3)) contour <- 1
    #
    my_plots <- list()
    count <- 1
    for (r in response){
		if (MEAN)
        for (i in termM){
			lcX <- label.count(i,x$labelsX,x$countX,"mean")
            int.labelX <- lcX[1]
            countX <- lcX[2]
            plotOptions2 <- list(ggtitle(paste("mean of", x$varsY[r])),plotOptions)            
            my_plots[[count]] <- plot.generic(mvrmObj=x, MEAN=1, STDEV=0, CORM=0, CORS=0, DEP=0,
                                              term=i, response=r, intercept=intercept, grid=grid,
                                              centre=centre, quantiles=quantiles, static=static, contour=contour,
                                              centreEffects=centreEffects, plotOptions=plotOptions2,
                                              int.label=int.labelX, count=countX, vars=x$varsX[[int.labelX]], 
                                              label=x$labelsX[int.labelX], is.D=x$is.Dx[[int.labelX]], 
                                              formula.term=x$formula.termsX[int.labelX], 
                                              which.Spec=x$which.SpecX[[int.labelX]], assign=x$assignX, 
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
			lcZ <- label.count(i,x$labelsZ,x$countZ,"stdev")
            int.labelZ <- lcZ[1]
            countZ <- lcZ[2]
            plotOptions2 <- list(ggtitle(paste("stdev of", x$varsY[[r]])),plotOptions)
            if (x$LUT > 1)plotOptions2 <- list(ggtitle(paste("innov stdev of", x$varsY[[r]])),plotOptions)
            my_plots[[count]] <- plot.generic(x, MEAN=0, STDEV=1, CORM=0, CORS=0, DEP=0,
                                              term=i, response=r, intercept=intercept, grid=grid,
                                              centre=centre, quantiles=quantiles, static=static, contour=contour,
                                              centreEffects=centreEffects, plotOptions=plotOptions2,
                                              int.label=int.labelZ, count=countZ, vars=x$varsZ[[int.labelZ]], 
                                              label=x$labelsZ[int.labelZ], is.D=x$is.Dz[[int.labelZ]], 
                                              formula.term=x$formula.termsZ[int.labelZ], 
                                              which.Spec=x$which.SpecZ[[int.labelZ]], assign=x$assignZ, 
                                              knots=x$Zknots, storeMeanVector=x$storeMeanVectorZ,
                                              storeIndicator=x$storeIndicatorZ, data=x$data,
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
			        lcC <- label.count(i,x$labelsC,x$countC,"autoregressive")
                    int.labelC <- lcC[1]
                    countC <- lcC[2]
                    plotOptions2 <- list(ggtitle(paste("autoreg", x$varsY[r1], x$varsY[r2])), plotOptions)            
                    my_plots[[count]] <- plot.generic(mvrmObj=x, MEAN=1, STDEV=0, CORM=0, CORS=0, DEP=1, 
                                                      term=i, response=Pair, intercept=intercept, grid=grid,
                                                      centre=centre, quantiles=quantiles, static=static, contour=contour,
                                                      centreEffects=centreEffects, plotOptions=plotOptions2,
                                                      int.label=int.labelC, count=countC, vars=x$varsC[[int.labelC]], 
                                                      label=x$labelsC[int.labelC], is.D=x$is.Dc[[int.labelC]], 
                                                      formula.term=x$formula.termsC[int.labelC], 
                                                      which.Spec=x$which.SpecC[[int.labelC]], assign=x$assignC, 
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
            plotOptions2 <- list(ggtitle(paste("mean cor")),plotOptions)            
            which.Spec <- 0
            if (length(x$which.SpecXc)) which.Spec <- x$which.SpecXc[[1]] 
            is.D <- 0
            if (length(x$is.Dxc)) is.D <- x$is.Dxc[[1]]
            vars <- x$varTime
            if (length(x$varsXc)) vars <- x$varsXc[[1]]
            my_plots[[count]] <- plot.generic(mvrmObj=x, MEAN=1, STDEV=0, CORM=1, CORS=0, DEP=0,
                                              term=-99, response=1, intercept=intercept, 
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
            plotOptions2 <- list(ggtitle(paste("stdev cor")),plotOptions)                        
            which.Spec <- 0
            if (length(x$which.SpecZc)) which.Spec <- x$which.SpecZc[[1]] 
            is.D <- 0
            if (length(x$is.Dzc)) is.D <- x$is.Dzc[[1]]
            vars <- x$varTime
            if (length(x$varsZc)) vars <- x$varsXc[[1]]
            my_plots[[count]] <- plot.generic(mvrmObj=x, MEAN=0, STDEV=1, CORM=0, CORS=1, DEP=0, 
                                              term=-99, response=1, intercept=intercept, 
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
	    if (CORM==0 && CORS==0 && DEP==0){nrow <- length(x$p)}else{nrow <- ceiling(sqrt(count-1))}
	}	     	
	if ((ask==FALSE) && (length(my_plots)>0)) quiet(print(grid.arrange(grobs = my_plots, nrow = nrow)))
}
