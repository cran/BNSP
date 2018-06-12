sm<-function (...,nknots=10,knots=NULL,bs="rd"){    
    pf <- parent.frame()    
    vars<-as.list(substitute(list(...)))[-1]
    d<-length(vars)
    if (d > 2) stop("Up to bivariate covariates supported")
    term<-NULL
    for (i in 1:d){ 
	    term[i]<-deparse(vars[[i]],backtick=TRUE,width.cutoff=500)
	    term[i] <- attr(terms(reformulate(term[i])),"term.labels")
    }
    if (length(unique(term))!=d) stop("Repeated variables are not permitted")        
    RNK <- round(nknots) 
    if (RNK!=nknots) warning("Number of knots should be integer and has been rounded")
    nknots <- RNK
    label<-paste("sm(",term[1],sep="")
    if (d>1) for (i in 2:d) label<-paste(label,",",term[i],sep="")
    label<-paste(label,")",sep="")  
    is.D<-rep(0,d)
    x2<-NULL
    for (i in 1:d){ 
        mm<-model.matrix(~eval(vars[[i]],pf))
        lvs<-levels(eval(vars[[i]],pf))
        if (length(lvs)==0) {colnames(mm)[2]<-c(term[i]); is.D[i]<-0}
        if (length(lvs)>=2) {colnames(mm)<-paste(term[i],lvs,sep=""); is.D[i]<-1}
        x2<-cbind(x2,mm[,-1,drop=FALSE])
    }
    if (d==1){
        if (is.null(knots)){             
            knots<-seq(from = 0, to = 1, length = nknots + 2)[-c(1, nknots + 2)]
            knots<-seq(from = 0, to = 1, length = nknots)
            knots<-unique(quantile(x2,knots,type=1))                
        }        
        nknots<-length(knots)
        X<-x2        
        if (bs=="pl"){ 
            for (i in 1:nknots) X<-cbind(X,pmax(x2-knots[i],rep(0,length(x2))))
		}
        else if (bs=="rd"){
			for (i in 1:nknots){ 
			    D<-(x2-knots[i])^2
			    D2<-D*log(D)
			    D2[which(D==0)]<-0
			    X<-cbind(X,D2)
			}
		} 
		else stop("chosen bs not supported")    
    } 
    else if (d==2){
        if (!is.null(knots) && !is.matrix(knots)) stop("for bivariate smoothers knots must in matrix form")
        if (!is.null(knots)) nknots<-c(dim(knots))
        if (is.null(knots) && length(nknots)==1) nknots<-rep(nknots,2)        
        if (is.null(knots)){             
			knots1<-seq(from = 0, to = 1, length = nknots[1] + 2)[-c(1, nknots[1] + 2)]
            knots1<-seq(from = 0, to = 1, length = nknots[1])
            knots2<-seq(from = 0, to = 1, length = nknots[2] + 2)[-c(1, nknots[2] + 2)]
            knots2<-seq(from = 0, to = 1, length = nknots[2])
            if (dim(x2)[2] == 2){                
                knots1<-unique(quantile(x2[,1],knots1,type=1))            
                knots2<-unique(quantile(x2[,2],knots2,type=1)) 
                knots<-as.matrix(expand.grid(knots1,knots2))
			}else if (is.D[1] == 1){
			    knots1<-unique(x2[,-dim(x2)[2]])
			    knots2<-unique(quantile(x2[,dim(x2)[2],drop=FALSE],knots2,type=1))
			    knots <- cbind(knots1[rep(1:nrow(knots1), length(knots2)), ], rep(knots2, each = nrow(knots1)))
			}else if (is.D[2] == 1){
				knots1<-unique(quantile(x2[,1],knots1,type=1))
			    knots2<-unique(x2[,-1])
			    knots <- cbind(rep(knots1, nrow(knots2)), knots2[rep(1:nrow(knots2), length(knots1)), ])
			}
			colnames(knots)<-colnames(x2)
		}      
        X<-x2
        if (bs=="rd"){ 
            for (i in 1:NROW(knots)){ 
                D<-apply((x2-matrix(knots[i,],nrow=NROW(x2),ncol=NCOL(x2),byrow=TRUE))^2,1,sum)
                D2<-D*log(D)
                D2[which(D==0)]<-0
                X<-cbind(X,D2)
			}
		} else stop("for bivariate smoothers only rd supported")
    }
    X<-data.frame(X)
    colnames(X)[1:NCOL(x2)]<-colnames(x2)
    colnames(X)[-c(1:NCOL(x2))]<-paste(label,1:(NCOL(X)-NCOL(x2)),sep=".")
    XK<-list(X=X,knots=knots)
    return(XK)    
}

DM<-function(formula,data,mm,ns,knots=NULL,meanVector,indicator){
    Rknots <- list()
    specials <- c('sm')    
    trms<-terms.formula(formula,specials=specials)
    attr(trms,"intercept")<-1
    nFactors<-dim(attr(trms,"factors"))[2]
    if (attr(trms,"response") && mm){ 
        y<-with(data,eval(trms[[2]])) 
        n<-length(y)
	} else n<-ns    
    if (!is.null(nFactors)){
        whereSpecials <- unique(unlist(apply(attr(trms,"factors")[unlist(attr(trms,"specials")), ,drop = F] > 0,1,which)))
        if (length(whereSpecials) < nFactors){
            if (length(whereSpecials) > 0){ 
                trms2<-drop.terms(trms, dropx = whereSpecials)
                trms3<-drop.terms(trms, dropx = -whereSpecials)
		    }else{ 
                trms2<-trms
                trms3<-NULL
		    }
		    X<-model.matrix(trms2,data=data)
	    }else{ 
            X<-matrix(1,ncol=1,nrow=n)
            colnames(X)<-"(Intercept)"
			trms3<-delete.response(trms)
	    }
        if (!is.null(trms3)) 
            for (i in 1:dim(attr(trms3,"factors"))[2]){ 
                XK<-with(data,eval(trms3[i][[2]]))        
                X<-cbind(X,XK$X)
                Rknots[[i]] <- XK$knots
			}
    }else if (attr(trms,"response") && mm){ 
        X<-model.matrix(trms,data=data)    
	}else{ 
        X<-matrix(1,ncol=1,nrow=n)
        colnames(X)<-"(Intercept)"
	}
	#print("From dm")
	#print(head(X))	
	
	if (missing(meanVector)) meanVector<-apply(as.matrix(X),2,mean)
	if (missing(indicator)){ 
	    unique.values<-unlist(lapply(apply(X,2,unique),length))
	    indicator<-which(unique.values<=2 & unique.values<=n)
	}
			
	#main<-colnames(X)
	#print(main)
	#main<-main[main %in% colnames(data)]
	#print(main)
	#is.F<-vector()
	#if (length(main) > 0) for (i in 1:length(main)) is.F[i]<-is.factor(data[,main[i]])
	#fact<-main[is.F]
	#indicator<-sort(c(indicator,which(colnames(X)==fact)))
	#print(indicator)	
		
	X[,-indicator] <- X[,-indicator] - matrix(1,nrow=n)%*%matrix(meanVector[-indicator],nrow=1)     
    
    if (mm) {Ret<-as.matrix(cbind(y,X))}
    else {Ret<-as.matrix(X)}
    
    return(list(yX=Ret,Rknots=Rknots,meanVector=meanVector,indicator=indicator))
}

match.call.defaults <- function(...) {
    call <- evalq(match.call(expand.dots = FALSE), parent.frame(1))
    formals <- evalq(formals(), parent.frame(1))
    for(k in setdiff(names(formals), names(call)))
        call[k] <- list( formals[[k]] )
    match.call(sys.function(sys.parent()), call)
}

mvrm <- function(formula,data,sweeps,burn=0,thin=1,seed,StorageDir,
                 c.betaPrior="IG(0.5,0.5*n)",c.alphaPrior="IG(1.1,1.1)",
                 pi.muPrior="Beta(1,1)",pi.sigmaPrior="Beta(1,1)",sigmaPrior="HN(2)",...){
    #Fix thin etc
    if (thin <= 0) thin=1
    thin <- as.integer(thin)
    sweeps <- as.integer(sweeps)
    burn <- as.integer(burn)
    if (missing(sweeps)) stop("provide sweeps argument")
    nSamples <- length(seq(1,(sweeps-burn),by=thin))
    #Specials
    specials <- c('sm')
    # Match call
    call <- match.call(expand.dots = FALSE)
    call2 <- match.call.defaults()
    #Data
    if (missing(data)) stop("provide data argument") 
    data <- as.data.frame(data)        
    #Formula
    if (length(as.Formula(formula))[1] > 1) stop("more than one response provided") 
    if (length(as.Formula(formula))[2] > 2) stop("more than two regression models provided")
    if (length(as.Formula(formula))[2] == 1) formula <- as.Formula(formula, ~1)
    formula.save <- formula
    formula.v<-formula(as.Formula(formula),lhs=0,rhs=2)
    formula<-formula(as.Formula(formula),lhs=1,rhs=1)
    # Design matrices, response, indicators
    XYK<-DM(formula,data,1,5,NULL)
    XY<-XYK$yX
    Y<-XY[,1]    
    n<-length(Y)
    X<-XY[,-1]
    Xknots<-XYK$Rknots
    storeMeanVectorX<-XYK$meanVector
    storeIndicatorX<-XYK$indicator
    main<-colnames(X)[-c(1,grep(specials,colnames(X)))]
    mainApp<-sapply(main,grep,colnames(X)[-1],simplify=FALSE)
    gX<-vector()
    count<-0
    if (length(mainApp) >=1){
        gX[mainApp[[1]]]<-1
        count<-1
    }
    if (length(mainApp) >= 2){
        for (i in 2:(length(mainApp))){
            if (sum(mainApp[[i-1]] %in% mainApp[[i]])==0) 
                count<-count+1
            gX[mainApp[[i]]]<-count
        }
    }
    NG<-count
    LG<-NCOL(X)-1 
    vecLG<-table(gX)
    cusumVecLG<-c(0,cumsum(vecLG))
    #    
    ZK<-DM(formula.v,data,0,n,NULL)
    Z<-ZK$yX  
    attr(Z,"assign")<-NULL
    Zknots<-ZK$Rknots
    storeMeanVectorZ<-ZK$meanVector
    storeIndicatorZ<-ZK$indicator        
    #       
    main<-colnames(Z)[-c(1,grep(specials,colnames(Z)))]
    mainApp<-sapply(main,grep,colnames(Z)[-1],simplify=FALSE)
    gZ<-vector()
    count<-0
    if (length(mainApp) >=1){
        gZ[mainApp[[1]]]<-1
        count<-1
    } 
    if (length(mainApp) >= 2){
        for (i in 2:(length(mainApp))){
            if (sum(mainApp[[i-1]] %in% mainApp[[i]])==0) 
                count<-count+1
            gZ[mainApp[[i]]]<-count
        }
    }
    ND<-count
    LD<-NCOL(Z)-1 
    vecLD<-table(gZ)
    cusumVecLD<-c(0,cumsum(vecLD))
    MVLD<-max(cusumVecLD)
    #Prior for c.beta
    sp<-strsplit(c.betaPrior,"IG\\(")
    sp<-strsplit(sp[[1]][2],"\\)")
    sp<-strsplit(sp[[1]][1],",")
    cetaParams<-c(as.numeric(sp[[1]][1]),eval(parse(text=sp[[1]][2])))
    #Prior for c.alpha
    sp<-strsplit(c.alphaPrior,"IG\\(")
    sp<-strsplit(sp[[1]][2],"\\)")
    sp<-strsplit(sp[[1]][1],",")
    calphaParams<-as.numeric(sp[[1]])
    #Prior for pi.mu
    if (!length(pi.muPrior)==NG) if (!length(pi.muPrior)==1) stop("Beta prior for pi.mu must have 1 distribution or as many distributions as the number of smooth terms in the the mean model")
    pimu<-NULL
    for (k in 1:length(pi.muPrior)){
        sp<-strsplit(pi.muPrior[k],"Beta\\(")
        sp<-strsplit(sp[[1]][2],"\\)")
        sp<-strsplit(sp[[1]][1],",")
        pimu<-rbind(pimu,as.numeric(sp[[1]]))
    }
    pimu<-c(pimu)
    if (length(pi.muPrior)==1) pimu<-rep(pimu,each=NG)
    #Prior for pi.sigma
    if (!length(pi.sigmaPrior)==ND) if (!length(pi.sigmaPrior)==1) stop("Beta prior for pi.sigma must have 1 distribution or as many distributions as the number of smooth terms in the the variance model")
    pisigma<-NULL
    for (k in 1:length(pi.sigmaPrior)){
        sp<-strsplit(pi.sigmaPrior[k],"Beta\\(")
        sp<-strsplit(sp[[1]][2],"\\)")
        sp<-strsplit(sp[[1]][1],",")
        pisigma<-rbind(pisigma,as.numeric(sp[[1]]))
    }
    pisigma<-c(pisigma)
    if (length(pi.sigmaPrior)==1) pisigma<-rep(pisigma,each=ND)   
    #Prior for sigma2
    specials<-c("HN","IG")
    sp<-strsplit(sigmaPrior,"\\(")
    if (sp[[1]][1] %in% specials){ 
        if (match(sp[[1]][1],specials)==1) HN<-1
        if (match(sp[[1]][1],specials)==2) HN<-0
    } else stop("unrecognised prior for sigma2")
    sp<-strsplit(sp[[1]][2],"\\)")
    sp<-strsplit(sp[[1]][1],",")
    sigmaParams<-as.numeric(sp[[1]])
    #Seed
    if (missing(seed)) seed<-as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31)
    # Storage directory & files
    WF <- 1
    if (!missing(StorageDir)){
        StorageDir <- path.expand(StorageDir)
        ncharwd <- nchar(StorageDir)}
    if (!missing(StorageDir)) if (!(substr(StorageDir,ncharwd,ncharwd)=="/")) StorageDir <- paste(StorageDir,"/",sep="")
    if (!missing(StorageDir)) if (!file.exists(StorageDir)) dir.create(StorageDir)
    if (missing(StorageDir)) stop("provide a storage directory via argument StorageDir")
    if (missing(StorageDir)) {WF <- 0; StorageDir <- paste(getwd(),"/",sep="")}
    on.exit(if (WF==0) file.remove(
    paste(StorageDir,"BNSP.alpha.txt",sep=""),
    paste(StorageDir,"BNSP.calpha.txt",sep=""),
    paste(StorageDir,"BNSP.beta.txt",sep=""),
    paste(StorageDir,"BNSP.sigma2.txt",sep=""),
    paste(StorageDir,"BNSP.delta.txt",sep=""),
    paste(StorageDir,"BNSP.gamma.txt",sep=""),
    paste(StorageDir,"BNSP.cbeta.txt",sep="")))
    #Tuning Parameters
    #if (missing(f)) 
    f<-1
    #if (missing(g)) 
    g<-10
    #if (missing(h)) 
    h<-5
	if (!length(h)==ND) h<-rep(mean(h),ND)
    #Block size selection
    #if (missing(blockSizeProbG)){ 
    blockSizeProbG <- rep(0,LG) 
    blockSizeProbG[1:5]<-c(10,25,30,25,10)
	#}
    #if (missing(blockSizeProbD)){
	blockSizeProbD <- rep(0,LD)
    blockSizeProbD[1:5]<-c(10,25,30,25,10)
    #}
    maxBSG <- max(which(blockSizeProbG>0))
    maxBSD <- max(which(blockSizeProbD>0))
    #Call C
    out<-.C("mvrmC",
            as.integer(seed),as.character(StorageDir),as.integer(WF),
            as.integer(sweeps),as.integer(burn),as.integer(thin),
            as.double(Y),as.double(as.matrix(X)),as.double(as.matrix(Z)),as.integer(n),as.integer(LG),as.integer(LD),
            as.double(blockSizeProbG),as.integer(maxBSG),as.double(blockSizeProbD),as.integer(maxBSD), 
            as.double(f),as.double(g),as.double(h),
            as.integer(NG),as.integer(ND),as.integer(vecLG),as.integer(vecLD),
            as.integer(cusumVecLG),as.integer(cusumVecLD),as.integer(MVLD),
            as.double(cetaParams),as.double(calphaParams),as.double(pimu),as.double(pisigma),
            as.integer(HN),as.double(sigmaParams),as.double(0))
    #Output
    fit <- list(call=call,call2=call2,formula=formula.save,seed=seed,
                data=data,X=X,Xknots=Xknots,Z=Z,Zknots=Zknots,LG=LG,LD=LD,
                mcpar=c(as.integer(burn+1),as.integer(seq(from=burn+1,by=thin,length.out=nSamples)[nSamples]),as.integer(thin)),
                nSamples=nSamples,storeMeanVectorX=storeMeanVectorX,storeMeanVectorZ=storeMeanVectorZ,
                #f=out[[17]][1],g=out[[18]][1],h=out[[19]],
                DIR=StorageDir,deviance=out[[33]][1]/nSamples,nullDeviance=-2*logLik(lm(Y ~ 1)),
                storeIndicatorX=storeIndicatorX,storeIndicatorZ=storeIndicatorZ)
    class(fit) <- 'mvrm'
    return(fit)
}

mvrm2mcmc <- function(mvrmObj,labels){
    all.labels <- c("alpha","calpha","cbeta","delta","beta","gamma","sigma2")
    mtch<-match(labels,all.labels)
	R<-NULL
    if (any(mtch==1) && mvrmObj$LD > 0) R<-matrix(unlist(read.table(paste(mvrmObj$DIR,"BNSP.alpha.txt",sep=""))),ncol=mvrmObj$LD,dimnames=list(c(),colnames(mvrmObj$Z)[-1]))
    if (any(mtch==2)) R<-matrix(unlist(read.table(paste(mvrmObj$DIR,"BNSP.calpha.txt",sep=""))),ncol=1,dimnames=list(c(),c("calpha")))
    if (any(mtch==3)) R<-matrix(unlist(read.table(paste(mvrmObj$DIR,"BNSP.cbeta.txt",sep=""))),ncol=1,dimnames=list(c(),c("ceta")))
    if (any(mtch==4) && mvrmObj$LD > 0) R<-matrix(unlist(read.table(paste(mvrmObj$DIR,"BNSP.delta.txt",sep=""))),ncol=mvrmObj$LD,dimnames=list(c(),paste("delta",seq(1,mvrmObj$LD),sep="_")))
    if (any(mtch==5) && mvrmObj$LG > 0) R<-matrix(unlist(read.table(paste(mvrmObj$DIR,"BNSP.beta.txt",sep=""))),ncol=mvrmObj$LG+1,dimnames=list(c(),colnames(mvrmObj$X)))
    if (any(mtch==6) && mvrmObj$LG > 0) R<-matrix(unlist(read.table(paste(mvrmObj$DIR,"BNSP.gamma.txt",sep=""))),ncol=mvrmObj$LG,dimnames=list(c(),paste("gamma",seq(1,mvrmObj$LG),sep="_")))
    if (any(mtch==7)) R<-matrix(unlist(read.table(paste(mvrmObj$DIR,"BNSP.sigma2.txt",sep=""))),ncol=1,dimnames=list(c(),c("sigma2")))	
	if (!is.null(R)){
	    attr(R, "mcpar") <- mvrmObj$mcpar
        attr(R, "class") <- "mcmc"
	}
	return(R)
}

plot.mvrm <- function(x, model="mean", term, intercept=TRUE, grid=30, 
                      centre=mean, quantiles=c(0.1, 0.9), static=TRUE, 
                      centreEffects=FALSE, plotOptions=list(), ...){
    mvrmObj=x
    x=y=NULL
    specials <- c('sm')
	plotLabel<-term
	if (!is.function(centre)) stop("centre must be a function, usually mean or median")
    n<-NROW(mvrmObj$X)
    grid<-round(grid)
    mainM<-colnames(mvrmObj$X)[-c(1,grep(specials,colnames(mvrmObj$X)))]
    mainV<-colnames(mvrmObj$Z)[-c(1,grep(specials,colnames(mvrmObj$Z)))]
    main<-unique(c(mainM,mainV))
	A<-sub(specials,"",term)
	B<-sub("\\(","",A)
	C<-sub("\\)","",B)
	sl<-strsplit(C,",")
	sl<-sub(" ","",sl[[1]])
	match1<-unlist(lapply(sl,match,colnames(mvrmObj$data)))
    match1<-match1[!is.na(match1)]
	label<-colnames(mvrmObj$data)[match1]
	count<-length(label)
	label2<-main[unlist(lapply(sl,grep,main))]
	new.formula1<-formula(as.Formula(mvrmObj$formula),lhs=1,rhs=1)#mean formula for later use
	new.formula2<-formula(as.Formula(mvrmObj$formula),lhs=0,rhs=2)#variance formula for later use
	MEAN<-0
    if (model=="mean" || model=="both"){ #check if chosen label is in the mean model
        if (!sum(label2 %in% mainM)==length(label2)) stop("chosen term doesn't appear in the mean model")
        if (length(label2)==0) stop("chosen term doesn't appear in the mean model")
        TEST<-1
        trms<-terms.formula(new.formula1,specials=specials)
        termLabels<-attr(trms,"term.labels")
        for (i in 1:length(termLabels)){
            test<-termLabels[i]
            A<-sub(specials,"",test)
            B<-sub("\\(","",A)
            C<-sub("\\)","",B)
            sl<-strsplit(C,",")
            sl<-gsub(" ", "", sl[[1]], fixed = TRUE)
            variables.in.term<-sl[sl %in% colnames(mvrmObj$data)]
            if (min(variables.in.term==label)){ 
                TEST<-0
                idTermM<-i
			}
         }  
         if (TEST) stop("chosen term doesn't appear in mean model")
         MEAN<-1
    }
    STDEV<-0
    if (model=="stdev" || model=="both"){ #check if chosen label is in the variance model
        if (!sum(label2 %in% mainV)==length(label2)) stop("chosen term doesn't appear in the variance model")
        if (length(label2)==0) stop("chosen term doesn't appear in the variance model")
        TEST<-1
        trms<-terms.formula(new.formula2,specials=specials)
        termLabels<-attr(trms,"term.labels")
        idTermS<-(-99)
        for (i in 1:length(termLabels)){
            test<-termLabels[i]
            A<-sub(specials,"",test)
            B<-sub("\\(","",A)
            C<-sub("\\)","",B)
            sl<-strsplit(C,",")
            sl<-gsub(" ", "", sl[[1]], fixed = TRUE)
            variables.in.term<-sl[sl %in% colnames(mvrmObj$data)]
            if (min(variables.in.term==label)){ 
                TEST<-0
                idTermS<-i
			}
        }  
        if (TEST) stop("chosen term doesn't appear in variance model")
        STDEV<-1
    }
    is.D<-vector()
    for (i in 1:count) is.D[i]<-is.factor(mvrmObj$data[,label[i]])
    unique.values<-vector()
    for (i in 1:count) unique.values[i]<-length(unique(mvrmObj$data[,label[i]]))
	indicator<-which(unique.values<=2)
    is.D[indicator]<-1     	
    if (count==2 && static && model=="both") par(mfrow=c(1,2))
	ML<-VL<-NULL
    if (MEAN) ML<-sort(unique(unlist(lapply(label2,grep,colnames(mvrmObj$X)))))
    if (STDEV) VL<-sort(unique(unlist(lapply(label2,grep,colnames(mvrmObj$Z)))))
    if (length(c(ML,VL))==0) stop("no matching variables")
    #if (length(ML)==1 && MEAN==1) stop("to plot 1-dim parameter use function plot()")
    #if (length(VL)==1 && STDEV==1) stop("to plot 1-dim parameter use function plot()")
    if ((intercept || sum(is.D)) & MEAN) ML<-c(1,ML)
    if (length(quantiles)==1) quantiles <- c(quantiles,1-quantiles)
    if (length(quantiles) > 2) stop("up to two quantiles")
    if (!is.null(quantiles)) quantiles <- unique(sort(quantiles))
    if (MEAN){
        if (count==1){  
			DsM<-unique(mvrmObj$X[,ML,drop=FALSE])
            adjG<-min(grid,NROW(DsM))          
            if (!is.D) DsM<-DsM[order(DsM[,label2[1]]),,drop=FALSE]	        
	        DsM<-DsM[seq(1,NROW(DsM),length.out=adjG),,drop=FALSE]	
	        if (is.D) DsM<-DsM[do.call(order,as.data.frame(DsM[,NCOL(DsM):1])),,drop=FALSE]        
		}
	    if (count==2){    
	        if (sum(is.D)){
	            DsM<-unique(mvrmObj$X[,ML])
                adjG<-min(grid^2,NROW(DsM))            
                DsM<-DsM[do.call(order,as.data.frame(DsM[,label2])),]                       
                #if (!is.D) DsM<-DsM[order(DsM[,label2[1]]),]	        
	            DsM<-DsM[seq(1,NROW(DsM),length.out=adjG),]		        	        
			}
	        #if (is.D) DsM<-DsM[do.call(order,as.data.frame(DsM[,NCOL(DsM):1])),]        	        	        	       
	        #adjG<-grid; if (is.D[1]) adjG<-length(unique(mvrmObj$data[,label[1]]))	        	        
	        #newR1<-seq(min(mvrmObj$X[,label2[1]]),max(mvrmObj$X[,label2[1]]),length.out=adjG)
            #if (is.D[1]) newR1<-factor(newR1)
            #adjG<-grid; if (is.D[2]) adjG<-length(unique(mvrmObj$data[,label[2]]))
            #newR2<-seq(min(mvrmObj$X[,label2[2]]),max(mvrmObj$X[,label2[2]]),length.out=adjG)
            #if (is.D[2]) newR2<-factor(newR2)
            #newData<-data.frame(expand.grid(newR1,newR2))
            #colnames(newData)<-label2            
            #trms<-terms.formula(mvrmObj$formula,specials=specials)
            #trms3<-drop.terms(trms, dropx = -idTermM)
            #DsM<-DM(trms3,newData,0,NROW(newData)) 
            #DM(formula(paste("~",plotLabel)),newData,0,225)                        
            #DsM<-DM(mvrmObj$formula,newData,0,NROW(newData))             
            #DsM<-DsM[,ML]
            #DsM<-unique(mvrmObj$X[,ML])
            #adjG<-min(grid^2,NROW(DsM))            	        
	        #DsM<-DsM[seq(1,NROW(DsM),length.out=adjG),]		        	        
	        #if (is.D) DsM<-DsM[do.call(order,as.data.frame(DsM[,NCOL(DsM):1])),]        	        	        	       
	        #adjG<-grid; if (is.D[1]) adjG<-length(unique(mvrmObj$data[,label[1]]))	        	        
	        #newR1<-seq(min(mvrmObj$X[,label2[1]]),max(mvrmObj$X[,label2[1]]),length.out=adjG)	        	        
            #if (is.D[1]) newR1<-factor(newR1)
            #adjG<-grid; if (is.D[2]) adjG<-length(unique(mvrmObj$data[,label[2]]))
            #newR2<-seq(min(mvrmObj$X[,label2[2]]),max(mvrmObj$X[,label2[2]]),length.out=adjG)            
            #if (is.D[2]) newR2<-factor(newR2)
            if (sum(is.D)==0){
                newR1<-seq(min(mvrmObj$data[,label2[1]]),max(mvrmObj$data[,label2[1]]),length.out=grid)
                newR2<-seq(min(mvrmObj$data[,label2[2]]),max(mvrmObj$data[,label2[2]]),length.out=grid)
                newData<-expand.grid(newR1,newR2)
                colnames(newData)<-label2                        
                
                trms<-terms.formula(new.formula1,specials=specials)
                whereSpecials <- unique(unlist(apply(attr(trms,"factors")[unlist(attr(trms,"specials")), ,drop = F] > 0,1,which)))
                whichKnots <- which(idTermM == whereSpecials)
                #trms3<-drop.terms(trms, dropx = -idTermM)            
                #XK<-with(mvrmObj$data,eval(trms3[1][[2]]))
                #D<-matrix(XK$knots,ncol=2)
                #Dstar<-XK$knots
                #print(Dstar)
                Dstar<-mvrmObj$Xknots[[whichKnots]]
                
                #DsM<-DM(trms3,newData,0,NROW(newData)) 
                #print(formula(paste("~",sub(")",",knots=D)",plotLabel))))
                #print(formula(parse(paste("~",sub(")",",knots=D)",plotLabel)))))
                DsM<-DM(formula(paste("~",sub(")",",knots=knots)",plotLabel))),data=newData,0,NROW(newData),knots=Dstar)$yX                   
                #DsM<-DM(mvrmObj$formula,newData,0,NROW(newData))             
                #DsM<-DsM[,ML]    
		    }                    
        }
        fitM<-matrix(0,nrow=mvrmObj$nSamples,ncol=NROW(DsM))
		etaFN <- file.path(paste(mvrmObj$DIR,"BNSP.beta.txt",sep="")) 
        eFile<-file(etaFN,open="r")
		for (i in 1:mvrmObj$nSamples){
            eta<-scan(eFile,what=numeric(),n=mvrmObj$LG+1,quiet=T)
            fitM[i,]<-DsM%*%matrix(c(eta[ML]))
            if (centreEffects) fitM[i,]<-fitM[i,]-mean(fitM[i,])
		}
		close(eFile)
		if (count==1 && !is.D){# && length(ML) > 2){
		    centreM<-apply(fitM,2,centre)
		    QM<-NULL
		    if (!is.null(quantiles)) QM<-drop(t(apply(fitM,2,quantile,probs=quantiles,na.rm=TRUE)))
		    newX<-DsM[,label2]+mean(mvrmObj$data[,label])		
		    dataM<-data.frame(cbind(newX,centreM,QM))		
		    nms <- c(label2,"centreM")
		    if (!is.null(quantiles)) nms<-c(nms,"QLM","QUM")
		    colnames(dataM) <- nms
		    plotElM<-c(geom_line(aes_string(x=label2,y=centreM),col=4,alpha=0.5),
                       geom_rug(data=data.frame(mvrmObj$data),aes_string(x=label,y=NULL),alpha=0.3),list(ylab(plotLabel)))
            if (!is.null(quantiles))
                plotElM<-c(geom_ribbon(data=dataM,aes_string(x=label2, ymin="QLM", ymax="QUM"),alpha =0.2),plotElM)
		    ggM<-ggplot(data=dataM)
	        plotM<-ggM + plotElM + ggtitle("mean") + plotOptions
	    }
	    if (count==1 && is.D){
	        df<-data.frame(x=rep(levels(mvrmObj$data[,label]),each=mvrmObj$nSamples),y=c(fitM))	        
            plotElM<-c(geom_boxplot(),list(xlab(label),ylab("")))
            ggM<-ggplot(data=df,aes(x=x,y=y))
	        plotM<-ggM + plotElM + ggtitle("mean") + plotOptions
	    }
	    if (count==1 && !is.D && length(ML)<3 && 1==0){	     #   
	        df<-data.frame(x=rep(0,length(c(fitM))),y=c(fitM))
            plotElM<-c(geom_boxplot(),list(ylab(label),xlab("")))
            ggM<-ggplot(data=df,aes(x=x,y=y))
	        plotM<-ggM + plotElM + ggtitle("mean") + 
	        theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
	        plotOptions
	    }
	    if (count==2 && sum(is.D)==1){
		    centreM<-apply(fitM,2,centre)
		    QM<-NULL
		    if (!is.null(quantiles)) QM<-drop(t(apply(fitM,2,quantile,probs=quantiles,na.rm=TRUE)))		    
		    #if (!is.D[1]) newR1<-newR1+mean(mvrmObj$data[,label[1]])
		    #if (!is.D[2]) newR2<-newR2+mean(mvrmObj$data[,label[2]])		    
		    #dataM<-data.frame(cbind(expand.grid(newR1,newR2),centreM,QM))				    
		    #if (!is.D[1]) DsM[,label2[1]]<-DsM[,label2[1]]+mean(mvrmObj$data[,label[1]])
		    #if (!is.D[2]) DsM[,label2[2]]<-DsM[,label2[2]]+mean(mvrmObj$data[,label[2]])		    
		    #if (is.D[1]) DsM[,label2[1]]<-factor(DsM[,label2[1]])
		    #if (is.D[2]) DsM[,label2[2]]<-factor(DsM[,label2[2]])		    
		    lvs<-levels(mvrmObj$data[,label[which(is.D==1)]])
		    nms<-paste(label[which(is.D==1)],lvs[-1],sep="")
		    pos<-match(nms,colnames(DsM))
		    fc<-lvs[c(max.col(DsM[,pos,drop=FALSE],ties.method="first")+1*(!apply(DsM[,pos,drop=FALSE],1,sum)==0))]
		    dataM<-data.frame(fc,DsM[,label[which(is.D==0)]]+mean(mvrmObj$data[,label[which(is.D==0)]]),centreM,QM)
            #if (is.D[1]) dataM[,label2[1]]<-factor(dataM[,label2[1]])
		    #if (is.D[2]) dataM[,label2[2]]<-factor(dataM[,label2[2]])		   		    
		    nms <- c(label[which(is.D==1)],label[which(is.D==0)],"centreM")
		    if (!is.null(quantiles)) nms<-c(nms,"QLM","QUM")
		    colnames(dataM) <- nms
		    plotElM<-c(geom_line(aes_string(x=label[which(is.D==0)],y=centreM,group=label[which(is.D==1)],colour=label[which(is.D==1)],linetype=label[which(is.D==1)]),alpha=1),
                       geom_rug(data=data.frame(mvrmObj$data),aes_string(x=label[which(is.D==0)],y=NULL),alpha=0.3),list(ylab(plotLabel)))
            if (!is.null(quantiles))
                plotElM<-c(geom_ribbon(data=dataM,aes_string(x=label[which(is.D==0)], ymin="QLM", ymax="QUM",group=label[which(is.D==1)],fill=label[which(is.D==1)]),alpha=0.2),plotElM)
		    ggM<-ggplot(data=dataM)
	        plotM<-ggM + plotElM + ggtitle("mean") + plotOptions
	    }
	    if (count==2 && sum(is.D)==0){
		    centreM<-apply(fitM,2,centre)		    
		    #newR1<-seq(min(mvrmObj$X[,label2[1]]),max(mvrmObj$X[,label2[1]]),length.out=adjG)
		    #newR2<-seq(min(mvrmObj$X[,label2[2]]),max(mvrmObj$X[,label2[2]]),length.out=adjG)
		    newR1<-newR1+mean(mvrmObj$data[,label[1]])
		    newR2<-newR2+mean(mvrmObj$data[,label[2]])
		    if (static){		        
		        defaultList<-list(x=as.numeric(newR1),y=as.numeric(newR2),z=matrix(centreM,length(newR1),length(newR2)),colvar=matrix(centreM,length(newR1),length(newR2)))
                along="xy"; if (is.D[1]) along="y"; if (is.D[2]) along="x"
                space=0.6; if (sum(is.D)) space<-0.8                                                        
                optionalList<-list(xlab=label2[1],ylab=label2[2],zlab=plotLabel,along=along,space=space,add=FALSE,bty="g",main="mean")
                allOptions<-c(defaultList,plotOptions,optionalList)
                do.call("ribbon3D",allOptions[!duplicated(names(allOptions))])
                if (!is.null(quantiles) && 1==0){
                    defaultList<-list(x=mvrmObj$X[,label2[1]],y=mvrmObj$X[,label2[2]],z=matrix(centreM+1,grid,grid),colvar=matrix(centreM+1,grid,grid))                                                
                    optionalList<-list(along="xy",space=0.6,add=TRUE)
                    allOptions<-c(defaultList,plotOptions,optionalList)
                    do.call("ribbon3D",allOptions[!duplicated(names(allOptions))])
				}
		    }
            if (!static){             
                newData<-as.matrix(expand.grid(newR1,newR2))
                a<-cbind(newData,centreM)
                if (is.null(plotOptions$col)) col=rainbow(16,2/3)
                else col<-plotOptions$col
                plotCentreM <- centreM
                if (min(centreM)<=0) plotCentreM <- centreM + abs(min(centreM)) + 1
                ra<-ceiling(length(col)*plotCentreM/max(plotCentreM))
                defaultList<-list(x=a,col=col[ra])
                optionalList<-list(size=0.4,bg=1,axisLabels=c(label2[1],plotLabel,label2[2]),main="mean")
                allOptions<-c(defaultList,plotOptions,optionalList)
                plotM<-do.call("scatterplot3js",allOptions[!duplicated(names(allOptions))])
			}
		}
	}
	if (STDEV){
        if (count==1){
            DsV<-unique(mvrmObj$Z[,VL,drop=FALSE])
            adjG<-min(grid,NROW(DsV))          
            if (!is.D) DsV<-DsV[order(DsV[,label2[1]]),,drop=FALSE]	        
	        DsV<-DsV[seq(1,NROW(DsV),length.out=adjG),,drop=FALSE]
	        if (is.D) DsV<-DsV[do.call(order,as.data.frame(DsV[,NCOL(DsV):1])),,drop=FALSE]	       
		}
		if (count==2){
	        if (sum(is.D)){
	            DsV<-unique(mvrmObj$Z[,VL])
                adjG<-min(grid^2,NROW(DsV))  
                DsV<-DsV[do.call(order,as.data.frame(DsV[,label2])),]                       
                #if (!is.D) DsM<-DsM[order(DsM[,label2[1]]),]	        
	            DsV<-DsV[seq(1,NROW(DsV),length.out=adjG),]
			}		        
	        #adjG<-grid; if (is.D[1]) adjG<-length(unique(mvrmObj$data[,label[1]]))	
	        #newR1<-seq(min(mvrmObj$Z[,label2[1]]),max(mvrmObj$Z[,label2[1]]),length.out=adjG)
            #if (is.D[1]) newR1<-factor(newR1)
            #adjG<-grid; if (is.D[2]) adjG<-length(unique(mvrmObj$data[,label[2]]))
            #newR2<-seq(min(mvrmObj$Z[,label2[2]]),max(mvrmObj$Z[,label2[2]]),length.out=adjG)
            #if (is.D[2]) newR2<-factor(newR2)
            #newData<-data.frame(expand.grid(newR1,newR2))
            #colnames(newData)<-label2
            #trms<-terms.formula(mvrmObj$formula.v,specials=specials)
            #trms3<-drop.terms(trms, dropx = -idTermS)
            #DsV<-DM(trms3,newData,0,NROW(newData))    
            #DsV<-DM(mvrmObj$formula.v,newData,0,NROW(newData))
            #DsV<-DsV[,VL] 
            if (sum(is.D)==0){
                newR1<-seq(min(mvrmObj$data[,label2[1]]),max(mvrmObj$data[,label2[1]]),length.out=grid)
                newR2<-seq(min(mvrmObj$data[,label2[2]]),max(mvrmObj$data[,label2[2]]),length.out=grid)
                newData<-expand.grid(newR1,newR2)
                colnames(newData)<-label2                        
                #trms<-terms.formula(new.formula2,specials=specials)
                #trms3<-drop.terms(trms, dropx = -idTermS)          
                #XK<-with(mvrmObj$data,eval(trms3[1][[2]]))                
                #Dstar<-XK$knots
                trms<-terms.formula(new.formula2,specials=specials)
                whereSpecials <- unique(unlist(apply(attr(trms,"factors")[unlist(attr(trms,"specials")), ,drop = F] > 0,1,which)))
                whichKnots <- which(idTermS == whereSpecials)
                Dstar<-mvrmObj$Zknots[[whichKnots]]
                DsV<-DM(formula(paste("~",sub(")",",knots=knots)",plotLabel))),data=newData,0,NROW(newData),knots=Dstar)$yX    
                DsV<-DsV[,-1]
		    }                                
        }
        fitV<-matrix(0,nrow=mvrmObj$nSamples,ncol=NROW(DsV))
		alphaFN <- file.path(paste(mvrmObj$DIR,"BNSP.alpha.txt",sep="")) 
        aFile<-file(alphaFN,open="r")
		sigma2FN <- file.path(paste(mvrmObj$DIR,"BNSP.sigma2.txt",sep="")) 
        s2File<-file(sigma2FN,open="r")
        for (i in 1:mvrmObj$nSamples){
        	alpha<-scan(aFile,what=numeric(),n=mvrmObj$LD,quiet=T)
		    if (intercept) s2<-scan(s2File,what=numeric(),n=1,quiet=T)
		    else s2=1
            fitV[i,]<-sqrt(s2*exp(DsV%*%matrix(c(alpha[VL-1]))))
            if (centreEffects) fitV[i,]<-fitV[i,]/mean(fitV[i,])
		}
		close(aFile)
		close(s2File)		
		if (count==1 && !is.D){
		    centreV<-apply(fitV,2,centre)
		    QV<-NULL
		    if (!is.null(quantiles)) QV<-drop(t(apply(fitV,2,quantile,probs=quantiles,na.rm=TRUE)))
		    newV<-DsV[,label2]+mean(mvrmObj$data[,label])
		    dataV<-data.frame(newV,centreV,QV)
            nms<-c(label2,"centreV")
		    if (!is.null(quantiles)) nms<-c(nms,"QLV","QUV")
		    colnames(dataV) <- nms        
		    plotElV<-c(geom_line(aes_string(x=label2,y=centreV),col=4,alpha=0.5),
                       geom_rug(data=data.frame(mvrmObj$data),aes_string(x=label,y=NULL),alpha=0.3),list(ylab(plotLabel)))
            if (!is.null(quantiles))
                plotElV<-c(geom_ribbon(data=dataV,aes_string(x=label2, ymin="QLV", ymax="QUV"),alpha =0.2),plotElV)   
		    ggV<-ggplot(data=dataV)
		    plotV<-ggV + plotElV + ggtitle("st dev") + plotOptions
		}
		if (count==1 && is.D){
	        df<-data.frame(x=rep(levels(mvrmObj$data[,label]),each=mvrmObj$nSamples),y=c(fitV))
            plotElV<-c(geom_boxplot(),list(xlab(label),ylab("")))
            ggV<-ggplot(data=df,aes(x=x,y=y))
	        plotV<-ggV + plotElV + ggtitle("st dev") + plotOptions
	    }
	    if (count==2 && sum(is.D)==1){		    						    
		    centreV<-apply(fitV,2,centre)
		    QV<-NULL
		    if (!is.null(quantiles)) QV<-drop(t(apply(fitV,2,quantile,probs=quantiles,na.rm=TRUE)))		
		    #if (!is.D[1]) DsV[,label2[1]]<-DsV[,label2[1]]+mean(mvrmObj$data[,label[1]])
		    #if (!is.D[2]) DsV[,label2[2]]<-DsV[,label2[2]]+mean(mvrmObj$data[,label[2]])	
		    lvs<-levels(mvrmObj$data[,label[which(is.D==1)]])
		    nms<-paste(label[which(is.D==1)],lvs[-1],sep="")
		    pos<-match(nms,colnames(DsV))
		    fc<-lvs[c(max.col(DsV[,pos,drop=FALSE],ties.method="first")+1*(!apply(DsV[,pos,drop=FALSE],1,sum)==0))]		    		    
		    dataV<-data.frame(fc,DsV[,label[which(is.D==0)]]+mean(mvrmObj$data[,label[which(is.D==0)]]),centreV,QV)		    
            #if (is.D[1]) dataV[,label2[1]]<-factor(dataV[,label2[1]])
		    #if (is.D[2]) dataV[,label2[2]]<-factor(dataV[,label2[2]])			    
		    #nms <- c(label,"centreV")
		    nms <- c(label[which(is.D==1)],label[which(is.D==0)],"centreM")
		    if (!is.null(quantiles)) nms<-c(nms,"QLV","QUV")
		    colnames(dataV) <- nms
		    plotElV<-c(geom_line(aes_string(x=label[which(is.D==0)],y=centreV,group=label[which(is.D==1)],linetype=label[which(is.D==1)],colour=label[which(is.D==1)]),alpha=0.5),
                       geom_rug(data=data.frame(mvrmObj$data),aes_string(x=label[which(is.D==0)],y=NULL),alpha=0.3),list(ylab(plotLabel)))
            if (!is.null(quantiles))
                plotElV<-c(geom_ribbon(data=dataV,aes_string(x=label[which(is.D==0)], ymin="QLV", ymax="QUV",group=label[which(is.D==1)],fill=label[which(is.D==1)]),alpha=0.2),plotElV)
		    ggV<-ggplot(data=dataV)
	        plotV<-ggV + plotElV + ggtitle("st dev") + plotOptions
	    }
		if (count==2 && sum(is.D)==0){
		    centreV<-apply(fitV,2,centre)
		    #newR1<-seq(min(mvrmObj$Z[,label2[1]]),max(mvrmObj$Z[,label2[1]]),length.out=adjG)
			#newR1<-newR1+mean(mvrmObj$data[,label[1]])
		    #newR2<-seq(min(mvrmObj$Z[,label2[2]]),max(mvrmObj$Z[,label2[2]]),length.out=adjG)
		    #newR2<-newR2+mean(mvrmObj$data[,label[2]])
		    
		    newR1<-newR1+mean(mvrmObj$data[,label[1]])
		    newR2<-newR2+mean(mvrmObj$data[,label[2]])
		    
		    if (static){		        				
				defaultList<-list(x=as.numeric(newR1),y=as.numeric(newR2),z=matrix(centreV,length(newR1),length(newR2)),colvar=matrix(centreV,length(newR1),length(newR2)))
                along="xy"; if (is.D[1]) along="y"; if (is.D[2]) along="x"
                space=0.6; if (sum(is.D)) space<-0.8
                optionalList<-list(xlab=label2[1],ylab=label2[2],zlab=plotLabel,along=along,space=space,add=FALSE,bty="g",main="st dev")
                allOptions<-c(defaultList,plotOptions,optionalList)
                do.call("ribbon3D",allOptions[!duplicated(names(allOptions))])
                if (!is.null(quantiles) && 1==0){
                    defaultList<-list(x=mvrmObj$X[,label2[1]],y=mvrmObj$X[,label2[2]],z=matrix(centreM+1,grid,grid),colvar=matrix(centreM+1,grid,grid))                                                
                    optionalList<-list(along="xy",space=0.6,add=TRUE)
                    allOptions<-c(defaultList,plotOptions,optionalList)
                    do.call("ribbon3D",allOptions[!duplicated(names(allOptions))])
				}
		    }
            if (!static){                
                newData<-as.matrix(expand.grid(newR1,newR2)) 
                a<-cbind(newData,centreV)
                if (is.null(plotOptions$col)) col=rainbow(16,2/3)
                else col<-plotOptions$col
                ra<-ceiling(length(col)*centreV/max(centreV))
                defaultList<-list(x=a,col=col[ra])
                optionalList<-list(size=0.4,bg=1,axisLabels=c(label2[1],plotLabel,label2[2]),main="st dev")
                allOptions<-c(defaultList,plotOptions,optionalList)
                plotV<-do.call("scatterplot3js",allOptions[!duplicated(names(allOptions))])                
			}
		}
    }
	if (count==1 || sum(is.D)==1){
	    if (model=="mean") return(plotM)
	    else if (model=="stdev") return(plotV)
	    else if (model=="both") return(grid.arrange(plotM, plotV, ncol=2))
	}	
	if (count==2 && !static && sum(is.D)==0){
	    if (MEAN) print(plotM)
	    if (STDEV) print(plotV)
	}	
}

predict.mvrm <- function(object,newdata,interval=c("none","credible","prediction"),level=0.95,nSamples=100, ...){
    #newdata
    if (missing(newdata) || is.null(newdata)) newdata<-object$data
    newdata <- as.data.frame(newdata)
    specials <- c('sm') 
    if (colnames(newdata)[1]=="newdata") colnames(newdata)<-colnames(object$X)[-c(1,grep(specials,colnames(object$X)))]
    #Formula
    formula <- object$formula
    formula.v<-formula(as.Formula(formula),lhs=0,rhs=2)
    formula<-formula(as.Formula(formula),lhs=1,rhs=1)
    #Reformulate and get X design matrix
    trms<-terms.formula(formula,specials=specials)
    whereSpecials <- unique(unlist(apply(attr(trms,"factors")[unlist(attr(trms,"specials")), ,drop = F] > 0,1,which)))
    nFactors<-dim(attr(trms,"factors"))[2]
    if (length(whereSpecials) == nFactors) {formula2<-~1}else 
        {trms2<-drop.terms(trms, dropx = whereSpecials)
        formula2<-reformulate(attr(trms2, "term.labels"))}    
    for (k in 1:length(whereSpecials)){
        term<-attr(trms,"term.labels")[whereSpecials][k]
        A<-sub(specials,"",term)
        B<-sub("\\(","",A)
        C<-sub("\\)","",B)
        sl<-strsplit(C,",")
        sl<-sub(" ","",sl[[1]])
        match1<-unlist(lapply(sl,match,colnames(object$data)))
        match1<-match1[!is.na(match1)]
        label<-colnames(object$data)[match1]
        K<-paste("sm(",paste(label,collapse=","),",knots= knots[[",k,"]])") 
        terms.f2<-terms(formula2)
        term.labels.f2<-attr(terms.f2,"term.labels")
        terms.f2<-paste(term.labels.f2,collapse=" + ")
        formula2<-reformulate(paste(terms.f2,"+", K,collapse=" + "))
    }
    nd<-object$data[0,match(colnames(newdata),colnames(object$data)),drop=FALSE]
    nd[1:NROW(newdata),] <- newdata
    X<-DM(formula2,nd,0,NROW(newdata),object$Xknots,object$storeMeanVectorX,object$storeIndicatorX)$yX
	fitM<-matrix(0,nrow=object$nSamples,ncol=NROW(newdata))
    etaFN <- file.path(paste(object$DIR,"BNSP.beta.txt",sep="")) 
    eFile<-file(etaFN,open="r")
	for (i in 1:object$nSamples){
        eta<-scan(eFile,what=numeric(),n=object$LG+1,quiet=T)
        fitM[i,]<-c(X%*%matrix(c(eta)))
	}
	close(eFile)
	predictions<-fit<-apply(fitM,2,mean)
	interval <- match.arg(interval)
	if (interval=="credible"){
		QMb<-apply(fitM,2,quantile,probs=c((1-level)/2,1-(1-level)/2),na.rm=TRUE)
		QM<-matrix(QMb,nrow=NROW(newdata),byrow=T)
		colnames(QM) <- rownames(QMb)
		predictions<-cbind(fit,QM)	
		colnames(predictions) <- c("fit", "lwr", "upr")	
	}
	if (interval=="prediction"){
	    trms<-terms.formula(formula.v,specials=specials)
        nFactors<-dim(attr(trms,"factors"))[2]
        if (!is.null(nFactors)){ 
            whereSpecials <- unique(unlist(apply(attr(trms,"factors")[unlist(attr(trms,"specials")), ,drop = F] > 0,1,which)))}
        else whereSpecials <- NULL
        if (is.null(nFactors) || length(whereSpecials) == nFactors) {formula2<-~1}else 
            {trms2<-drop.terms(trms, dropx = whereSpecials)
            formula2<-reformulate(attr(trms2, "term.labels"))}        
        if (length(whereSpecials) > 0)
            for (k in 1:length(whereSpecials)){
                term<-attr(trms,"term.labels")[whereSpecials][k]
                A<-sub(specials,"",term)
                B<-sub("\\(","",A)
                C<-sub("\\)","",B) 
                sl<-strsplit(C,",")
                sl<-sub(" ","",sl[[1]])
                match1<-unlist(lapply(sl,match,colnames(object$data)))
                match1<-match1[!is.na(match1)]
                label<-colnames(object$data)[match1]
                K<-paste("sm(",paste(label,collapse=","),",knots= knots[[",k,"]])") 
                terms.f2<-terms(formula2)
                term.labels.f2<-attr(terms.f2,"term.labels")
                terms.f2<-paste(term.labels.f2,collapse=" + ")
                formula2<-reformulate(paste(terms.f2,"+", K,collapse=" + "))
            }
        Z<-DM(formula2,newdata,0,NROW(newdata),object$Zknots,object$storeMeanVectorZ,object$storeIndicatorZ)$yX
        Z<-Z[,-1]
        fitV<-matrix(0,nrow=object$nSamples,ncol=NROW(newdata))     
		alphaFN <- file.path(paste(object$DIR,"BNSP.alpha.txt",sep="")) 
        aFile<-file(alphaFN,open="r")
		sigma2FN <- file.path(paste(object$DIR,"BNSP.sigma2.txt",sep="")) 
        s2File<-file(sigma2FN,open="r")
        for (i in 1:object$nSamples){
        	alpha<-scan(aFile,what=numeric(),n=object$LD,quiet=T)
		    s2<-scan(s2File,what=numeric(),n=1,quiet=T)
            fitV[i,]<-sqrt(s2*exp(Z%*%matrix(c(alpha))))            
		}
		close(aFile)
		close(s2File)	
		sample.Pred <- rnorm(nSamples*object$nSamples*NROW(newdata),c(fitM),c(fitV))	
		sample.Pred <- matrix(sample.Pred,ncol=NROW(newdata)*object$nSamples,byrow=T)
		QM<-NULL
		for (k in 1:NROW(newdata))
		    QM<-rbind(QM,quantile(sample.Pred[,(1+object$nSamples*(k-1)):(object$nSamples*k)],probs=c((1-level)/2,1-(1-level)/2),na.rm=TRUE))
		predictions<-cbind(fit,QM)
		colnames(predictions) <- c("fit", "lwr", "upr")
	}
	return(data.frame(predictions))
}

print.mvrm <- function(x,  digits = 5, ...) {
  cat("\nCall:\n")
  print(x$call)
  cat("\n")
  cat(x$nSamples,"posterior samples\n")
  G<-as.data.frame(mvrm2mcmc(x,"gamma"))
  D<-as.data.frame(mvrm2mcmc(x,"delta"))
  colnames(G)<- colnames(x$X)[-1] 
  colnames(D)<- colnames(x$Z)[-1] 
  if (x$LG > 0){
      cat("\nMean model: marginal probabilities of term inclusion:\n")
      print(apply(G,2,mean), digits=5)
  }
  if (x$LD > 0){
      cat("\nVariance model: marginal probabilities of term inclusion:\n")
      print(apply(D,2,mean), digits=5)
  }
}	

summary.mvrm <- function(object, nModels = 5, digits = 5, ...) {
  cat("\nSpecified model for the mean and variance:\n")
  print(object$formula, showEnv = FALSE)
  Prior<-NULL
  for (k in 9:13) Prior<-c(Prior,paste(names(as.list(object$call2))[k],"=",object$call2[[k]]))
  cat("\nSpecified priors:\n")
  print(noquote(sub("Prior","",Prior)))
  cat("\nTotal posterior samples:",object$nSamples,"; burn-in:",object$mcpar[1]-1,"; thinning:",object$mcpar[3],"\n")
  cat("\nFiles stored in",object$DIR,"\n")
  deviance <- matrix(c(object$nullDeviance,object$deviance),ncol=1)
  rownames(deviance) <- c("Null deviance:", "Mean posterior deviance:")
  colnames(deviance) <- c("")
  print(deviance, digits = digits)
  cat("\nJoint mean/variance model posterior probabilities:\n")
  G<-as.data.frame(mvrm2mcmc(object,"gamma"))
  D<-as.data.frame(mvrm2mcmc(object,"delta"))
  colnames(G)<- sub(" ",".",paste("mean",colnames(object$X)[-1])) 
  if (object$LD > 0) colnames(D)<- sub(" ",".",paste("var",colnames(object$Z)[-1])) 
  if (object$LG > 0 && object$LD > 0) Joint<-cbind(G,D)
  if (object$LG > 0 && object$LD == 0) Joint<-G
  if (object$LG == 0 && object$LD > 0) Joint<-D
  g<-count(Joint)
  g<-g[order(g$freq,decreasing=TRUE),]
  rownames(g)<-seq(1,NROW(g))
  g$prob<-100*g$freq/sum(g$freq)
  g$cumulative<-cumsum(g$prob)
  print(g[1:nModels,])
  cat("Displaying", nModels, "models of the",NROW(g),"visited\n")
  cat(nModels,"models account for",sub(" ","",paste(g$cumulative[nModels],"%")),"of the posterior mass\n")
}
