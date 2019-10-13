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

plot.generic <- function(mvrmObj, model, term, response, intercept, grid, centre, quantiles, static, contour,
                         centreEffects, plotOptions, 
                         int.label,count,vars,label,is.D,formula.term,which.Spec,assign,
                         knots,storeMeanVector,storeIndicator,data,MEAN,STDEV){
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
            
            Ds<-DM(formula=reformulate(formula.term),data=newData,n=NROW(newData),knots=Dstar,
                   meanVector=storeMeanVector[indeces],indicator=storeIndicator[indeces])$X
            if ((!1%in%V)||STDEV) Ds<-Ds[,-1]
        }
        if (is.D){
            newData<-data.frame(unique(with(data,eval(as.name(vars[1])))))
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
        fit<-matrix(0,nrow=mvrmObj$nSamples,ncol=NROW(Ds))
	    etaFN <- file.path(paste(mvrmObj$DIR,"BNSP.beta.txt",sep=""))
        eFile<-file(etaFN,open="r")
	    for (i in 1:mvrmObj$nSamples){
            eta<-scan(eFile,what=numeric(),n=p*(mvrmObj$LG+1),quiet=TRUE)
            fit[i,]<-as.matrix(Ds)%*%matrix(c(eta[V+(response-1)*(mvrmObj$LG+1)]))
            if (centreEffects) fit[i,]<-fit[i,]-mean(fit[i,])
	    }
	    close(eFile)
    }
	if (STDEV){	
	    fit<-matrix(0,nrow=mvrmObj$nSamples,ncol=NROW(Ds))
	    alphaFN <- file.path(paste(mvrmObj$DIR,"BNSP.alpha.txt",sep=""))
        aFile<-file(alphaFN,open="r")
	    sigma2FN <- file.path(paste(mvrmObj$DIR,"BNSP.sigma2.txt",sep=""))
        s2File<-file(sigma2FN,open="r")
        for (i in 1:mvrmObj$nSamples){
            alpha<-scan(aFile,what=numeric(),n=p*mvrmObj$LD,quiet=TRUE)
		    if (intercept) s2<-scan(s2File,what=numeric(),n=p,quiet=TRUE)[response]
		    else s2=1
            fit[i,]<-sqrt(s2*exp(as.matrix(Ds)%*%matrix(c(alpha[V+(response-1)*mvrmObj$LD]))))
            if (centreEffects) fit[i,]<-fit[i,]/mean(fit[i,])
	    }
	    close(aFile)
	    close(s2File)  
	}
	if (count==1 && !is.D){
	    centreM<-apply(fit,2,centre)
	    QM<-NULL
	    if (!is.null(quantiles)) QM<-drop(t(apply(fit,2,quantile,probs=quantiles,na.rm=TRUE)))
	    newX<-newData
	    dataM<-data.frame(cbind(newX,centreM,QM))
	    nms <- c(vars,"centreM")
	    if (!is.null(quantiles)) nms<-c(nms,"QLM","QUM")
	    colnames(dataM) <- nms
	    dataRug<-data.frame(b=with(data,eval(as.name(vars))))
	    plotElM<-c(geom_line(aes_string(x=as.name(vars),y=centreM),col=4,alpha=0.5),
                   geom_rug(mapping=aes(x=dataRug$b),data=dataRug,alpha=0.3),
                   list(ylab(label)))
        if (!is.null(quantiles))
            plotElM<-c(geom_ribbon(data=dataM,aes_string(x=vars, ymin="QLM", ymax="QUM"),alpha=0.2),plotElM)
		ggM<-ggplot(data=dataM)
	    plotM<-ggM + plotElM + plotOptions
	}     
	if (count==1 && is.D){	
		lvs<-levels(with(data,eval(as.name(vars[1]))))
		if (is.null(lvs)) lvs<-unique(with(data,eval(as.name(vars[1]))))
	    df<-data.frame(x=rep(lvs,each=mvrmObj$nSamples),y=c(fit))
        plotElM<-c(geom_boxplot(),list(xlab(label),ylab("")))
        ggM<-ggplot(data=df,aes(x=factor(df$x),y=df$y))
	    plotM<-ggM + plotElM + plotOptions
	}     
	if (count==2 && sum(is.D)==1){
	    centreM<-apply(fit,2,centre)
		QM<-NULL
		if (!is.null(quantiles)) QM<-drop(t(apply(fit,2,quantile,probs=quantiles,na.rm=TRUE)))
		disc.var<-vars[which(is.D==1)]
		cont.var<-vars[which(is.D==0)]
		dataM<-data.frame(newData,centreM,QM)
		nms <- c(colnames(newData),"centreM")
		if (!is.null(quantiles)) nms<-c(nms,"QLM","QUM")
		colnames(dataM) <- nms
		DG<-data.frame(b=with(data,eval(as.name(cont.var))),c=with(data,eval(as.name(disc.var))))
		plotElM<-c(geom_line(aes_string(x=cont.var,y=centreM,
		           group=disc.var,colour=disc.var,linetype=disc.var),alpha=0.80),
                   geom_rug(mapping=aes(x=DG$b,group=c,colour=c,linetype=c),data=DG,alpha=0.3),
                   list(ylab(label)))
        if (!is.null(quantiles))
            plotElM<-c(geom_ribbon(data=dataM,aes_string(x=cont.var,ymin="QLM",ymax="QUM",
                                   group=disc.var,fill=as.name(disc.var)),alpha=0.2),plotElM)
		ggM<-ggplot(data=dataM)
	    plotM<-ggM + plotElM + plotOptions
    }
    if (count==2 && sum(is.D)==0){
	    centreM<-apply(fit,2,centre)
	    if (contour==1){	        
	        dataM<-data.frame(newData,centreM)
		    nms <- c(colnames(newData),"centreM")		    
		    colnames(dataM) <- nms
		    ggM<-ggplot(data=dataM,aes_string(x=vars[1],y=vars[2],z=centreM))    
	        plotM<-ggM + geom_contour(data=dataM,aes(colour = stat(centreM))) + plotOptions
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


plot.mvrm <- function(x, model, term, response, intercept = TRUE, grid = 30,
                      centre = mean, quantiles = c(0.1, 0.9), contour = TRUE, static = TRUE, 
                      centreEffects = FALSE, plotOptions = list(), nrow, ask = FALSE,...)
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
    if (max(response) > x$p) stop("argument response exceeds the number of responses");
    if (missing(model)) model<-"both"
    MEAN <- 0; STDEV <- 0
    if ((model=="both" || model=="mean") && (x$NG > 0)) MEAN <- 1
    if ((model=="both" || model=="stdev") && (x$ND > 0)) STDEV <- 1
    if (MEAN==0 && STDEV==0) stop("no terms to plot; only intercepts in the model");    
    termM <- NULL; termSD <- NULL; 
    if (missing(term)) {termM<-1:x$NG; termSD<-1:x$ND}
    if (!missing(term)) {termM <- termSD <- term}
    if (missing(nrow)) nrow <- length(response)
    #
    if (((MEAN+STDEV+length(termM)+length(termSD)+length(response))>3) && contour == 0)
        warning("for 3-d plots, specify one model, one term and one response at a time")
    if (((MEAN+STDEV+length(termM)+length(termSD)+length(response))>3)) contour <- 1
    #
    my_plots <- list()
    count <- 1
    for (r in response) {
		if (MEAN)
        for (i in termM){
			lcX <- label.count(i,x$labelsX,x$countX,"mean")
            int.labelX <- lcX[1]
            countX <- lcX[2]
            plotOptions2 <- list(ggtitle(paste("mean of", x$varsY[r])),plotOptions)            
            my_plots[[count]] <- plot.generic(mvrmObj=x, model="mean", term=i, response=r, intercept=intercept, grid=grid,
                                              centre=centre, quantiles=quantiles, static=static, contour=contour,
                                              centreEffects=centreEffects, plotOptions=plotOptions2,
                                              int.label=int.labelX, count=countX, vars=x$varsX[[int.labelX]], 
                                              label=x$labelsX[int.labelX], is.D=x$is.Dx[[int.labelX]], 
                                              formula.term=x$formula.termsX[int.labelX], 
                                              which.Spec=x$which.SpecX[[int.labelX]], assign=x$assignX, 
                                              knots=x$Xknots, storeMeanVector=x$storeMeanVectorX,
                                              storeIndicator=x$storeIndicatorX, data=x$data, MEAN=1, STDEV=0)
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
            my_plots[[count]] <- plot.generic(x, model="stdev", term=i, response=r, intercept=intercept, grid=grid,
                                              centre=centre, quantiles=quantiles, static=static, contour=contour,
                                              centreEffects=centreEffects, plotOptions=plotOptions2,
                                              int.label=int.labelZ, count=countZ, vars=x$varsZ[[int.labelZ]], 
                                              label=x$labelsZ[int.labelZ], is.D=x$is.Dz[[int.labelZ]], 
                                              formula.term=x$formula.termsZ[int.labelZ], 
                                              which.Spec=x$which.SpecZ[[int.labelZ]], assign=x$assignZ, 
                                              knots=x$Zknots, storeMeanVector=x$storeMeanVectorZ,
                                              storeIndicator=x$storeIndicatorZ, data=x$data, MEAN=0, STDEV=1)
            if ((ask==TRUE) && (length(my_plots)>0)) print(my_plots[[count]])
            if (count==1)
                oldpar <- c(oldpar, par(ask=ask))
            count <- count + 1
        }
	}
	if ((ask==FALSE) && (length(my_plots)>0)) quiet(print(grid.arrange(grobs = my_plots, nrow = nrow)))
}
