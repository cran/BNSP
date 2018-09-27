s<-function(...,data,knots=NULL,absorb.cons=FALSE, scale.penalty=TRUE,n=nrow(data),dataX=NULL,
            null.space.penalty=FALSE,sparse.cons=0,diagonal.penalty=FALSE,apply.by=TRUE,modCon=0,
            k=-1,fx=FALSE,bs="tp",m=NA,by=NA,xt=NULL,id=NULL,sp=NULL,pc=NULL){
    a<-mgcv::s(...,k=k,fx=fx,bs=bs,m=m,by=by,xt=xt,id=id,sp=sp,pc=pc)
    a$by<-deparse(substitute(by), backtick = TRUE, width.cutoff = 500)
    b<- mgcv::smoothCon(a,data=data,knots=knots,absorb.cons=absorb.cons,
              scale.penalty=scale.penalty,n=n,dataX=dataX,
              null.space.penalty=null.space.penalty,sparse.cons=sparse.cons,
              diagonal.penalty=diagonal.penalty,apply.by=apply.by,modCon=modCon) 
    return(b)
}

te<-function(...,data,knots=NULL,absorb.cons=FALSE, scale.penalty=TRUE,n=nrow(data),dataX=NULL,
            null.space.penalty=FALSE,sparse.cons=0,diagonal.penalty=FALSE,apply.by=TRUE,modCon=0,
            k=NA,bs="cr",m=NA,d=NA,by=NA,fx=FALSE,np=TRUE,xt=NULL,id=NULL,sp=NULL,pc=NULL){
    a<-mgcv::te(...,k=k,bs=bs,m=m,d=d,by=by,fx=fx,np=np,xt=xt,id=id,sp=sp,pc=pc)
    a$by<-deparse(substitute(by), backtick = TRUE, width.cutoff = 500)
    b<- mgcv::smoothCon(a,data=data,knots=knots,absorb.cons=absorb.cons,
              scale.penalty=scale.penalty,n=n,dataX=dataX,
              null.space.penalty=null.space.penalty,sparse.cons=sparse.cons,
              diagonal.penalty=diagonal.penalty,apply.by=apply.by,modCon=modCon) 
    return(b)
}

ti<-function(...,data,knots=NULL,absorb.cons=FALSE, scale.penalty=TRUE,n=nrow(data),dataX=NULL,
            null.space.penalty=FALSE,sparse.cons=0,diagonal.penalty=FALSE,apply.by=TRUE,modCon=0,
            k=NA,bs="cr",m=NA,d=NA,by=NA,fx=FALSE,np=TRUE,xt=NULL,id=NULL,sp=NULL,mc=NULL,pc=NULL){
    a<-mgcv::ti(...,k=k,bs=bs,m=m,d=d,by=by,fx=fx,np=np,xt=xt,id=id,sp=sp,mc=mc,pc=pc)
    a$by<-deparse(substitute(by), backtick = TRUE, width.cutoff = 500)
    b<- mgcv::smoothCon(a,data=data,knots=knots,absorb.cons=absorb.cons,
              scale.penalty=scale.penalty,n=n,dataX=dataX,
              null.space.penalty=null.space.penalty,sparse.cons=sparse.cons,
              diagonal.penalty=diagonal.penalty,apply.by=apply.by,modCon=modCon) 
    return(b)
}

sm<-function(...,k=10,knots=NULL,bs="rd"){    
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
    nknots<-k
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
    X<-x2
    if (d==1 & nknots>0){
        if (is.null(knots)){             
            knots<-seq(from = 0, to = 1, length = nknots + 2)[-c(1, nknots + 2)]
            knots<-seq(from = 0, to = 1, length = nknots)
            knots<-unique(quantile(x2,knots))
            knots<-data.frame(knots=knots)
        } 
        if (!is.null(knots)){
			knots<-data.frame(knots=knots)
		}
        nknots<-dim(knots)[1]
        if (bs=="pl"){ 
            for (i in 1:nknots) X<-cbind(X,pmax(x2-knots[i],rep(0,length(x2))))
		}
        else if (bs=="rd"){
			for (i in 1:nknots){				 
			    D<-(x2-knots[i,])^2    
			    D2<-D*log(D)
			    D2[which(D==0)]<-0			    
			    X<-cbind(X,D2)
			}
		} 
		else stop("chosen bs not supported")    
    } 
    else if (d==2 & nknots>0){
        #if (!is.null(knots) && !is.matrix(knots)) stop("for bivariate smoothers knots must in matrix form")
        if (!is.null(knots)) nknots<-c(dim(knots))
        if (is.null(knots) && length(nknots)==1) nknots<-rep(nknots,2)        
        if (is.null(knots)){             
			knots1<-seq(from = 0, to = 1, length = nknots[1] + 2)[-c(1, nknots[1] + 2)]
            knots1<-seq(from = 0, to = 1, length = nknots[1])
            knots2<-seq(from = 0, to = 1, length = nknots[2] + 2)[-c(1, nknots[2] + 2)]
            knots2<-seq(from = 0, to = 1, length = nknots[2])
            if (dim(x2)[2] == 2){                
                knots1<-unique(quantile(x2[,1],knots1))
                knots2<-unique(quantile(x2[,2],knots2))
                knots<-as.matrix(expand.grid(knots1,knots2))
			}else if (is.D[1] == 1){
			    knots1<-unique(x2[,-dim(x2)[2]])
			    knots2<-unique(quantile(x2[,dim(x2)[2],drop=FALSE],knots2))
			    knots <- cbind(knots1[rep(1:nrow(knots1), length(knots2)), ], rep(knots2, each = nrow(knots1)))
			}else if (is.D[2] == 1){
				knots1<-unique(quantile(x2[,1],knots1))
			    knots2<-unique(x2[,-1])
			    knots <- cbind(rep(knots1, nrow(knots2)), knots2[rep(1:nrow(knots2), length(knots1)), ])
			}
			colnames(knots)<-colnames(x2)
		}      
        if (bs=="rd"){ 
            for (i in 1:NROW(knots)){ 
				#print(x2-matrix(knots[i,],nrow=NROW(x2),ncol=NCOL(x2),byrow=TRUE))
                #print("h")
                D<-apply((x2-matrix(as.numeric(knots[i,]),nrow=NROW(x2),ncol=NCOL(x2),byrow=TRUE))^2,1,sum)
                D2<-D*log(D)
                D2[which(D==0)]<-0
                X<-cbind(X,D2)
			}
		} else stop("for bivariate smoothers only rd supported")
    }
    X<-data.frame(X)
    colnames(X)[1:NCOL(x2)]<-colnames(x2)
    colnames(X)[-c(1:NCOL(x2))]<-paste(label,1:(NCOL(X)-NCOL(x2)),sep=".")
    XK<-list(X=X,knots=knots,count=d,vars=unlist(term),is.D=is.D,label=label)
    return(XK)    
}

DM<-function(formula,data,n,knots=NULL,predInd=FALSE,meanVector,indicator,mvrmObj){
    is.M <- NULL
    count <- NULL
    k<-0
    Rknots <- list()
    vars <- list()
    is.D <- list()
    which.Spec<-list()
    specials <- c('sm','s','te','ti')    
    trms<-terms.formula(formula,specials=specials)
    attr(trms,"intercept")<-1
    nFactors<-dim(attr(trms,"factors"))[2]
    if (attr(trms,"response")){ 
        y<-with(data,eval(attr(trms,"variables")[[2]]))
        n<-length(y)
	}else{ 
	    y <- NULL
	    n<-n
	} 
	labels <- colnames(attr(trms,"factors"))
	formula.terms<-attr(trms,"term.labels")
	Design.X<-matrix(1,ncol=1,nrow=n)
    colnames(Design.X)<-"(Intercept)"
    assign <- 0    
    if (!is.null(nFactors)){
        for (i in 1:nFactors){
            trms2<-drop.terms(trms, dropx = -i)            
            is.spec<-sum(unlist(lapply(attr(trms2,"specials"),sum))) > 0
            is.sm<-!is.null(attr(trms2,"specials")$sm)
            is.s<-!is.null(attr(trms2,"specials")$s) || !is.null(attr(trms2,"specials")$te) ||!is.null(attr(trms2,"specials")$ti)
            is.M.i<-0 
            if (!is.spec && (length(with(data,eval(attr(trms2,"variables"))[[1]]))/n > 1)) is.M.i <-1     
            is.M <- c(is.M, is.M.i)            
            if (!is.spec){                 
                Design.Xt<-model.matrix(trms2,data=data)             
                remove.cols<-which(c(apply(Design.Xt,2,sd)==0))
                if (length(remove.cols)>0) Design.Xt<-Design.Xt[,-remove.cols,drop=FALSE]
                count<-c(count,1)
                vars[[i]]<-as.character(attr(trms2,"variables")[[2]])
                is.D.i<-0
                lvs<-levels(with(data,eval(attr(trms2,"variables")))[[1]])                
                if (length(lvs)>=2) is.D.i<-1
                if (length(c(unique(with(data,eval(attr(trms2,"variables")))[[1]])))==2) is.D.i<-1                
                is.D[[i]]<-is.D.i
                which.Spec[[i]]<--99
			}
            if (is.sm){				
			    XK<-with(data,eval(trms2[1][[2]]))			    
                Design.Xt<-XK$X
                remove.cols<-which(c(apply(Design.Xt,2,sd)==0))
                Rknots[[k<-k+1]]<- XK$knots      
                count<-c(count,XK$count)
                vars[[i]]<-XK$vars 
                is.D[[i]]<-XK$is.D
                labels[i]<-XK$label 
                which.Spec[[i]]<-k  
			}
			if (is.s){							    
			    dataPr<-match(c("data"), names(trms2[1][[2]]), 0L)
			    trmWD<-trms2[1][[2]]			    
			    if (!predInd){ 
			        if (!dataPr){  
     			        trmWD<-deparse(trms2[1][[2]], backtick = TRUE, width.cutoff = 500)
     			        trmWD<-substr(trmWD,1,nchar(trmWD)-1)
     			        trmWD<-paste(trmWD,", data=data)")
				    }
     			    evalS<-eval(parse(text=trmWD))		     			    
     			    NL<-1
     			    if (!evalS[[1]]$by=="NA") NL<-max(1,length(levels(with(data,eval(as.name(evalS[[1]]$by))))))
     			    Design.Xt<-NULL     			    
     			    for (j in 1:NL) Design.Xt<-cbind(Design.Xt,evalS[[j]]$X)     			         			    
     			    vars[[i]]<-evalS[[1]]$term
     			    Rknots[[k<-k+1]]<-data.frame(evalS[[1]]$xp)
			        try(colnames(Rknots[[k]])<-vars[[i]],silent=TRUE)
			        count<-c(count,evalS[[1]]$dim)				    		    	    
			        is.D.i<-rep(0,evalS[[1]]$dim)		    
			        labels[i]<-evalS[[1]]$label
			        if (!evalS[[1]]$by=="NA"){
				        count[length(count)]<-count[length(count)]+1
				        vars[[i]]<-c(evalS[[1]]$term,evalS[[1]]$by)			    			       
			            is.D.i<-c(is.D.i,1)
			            labels[i]<-substr(labels[i],1,nchar(labels[i])-1)	
				    }      
			        is.D[[i]]<-is.D.i			    			        			        
			        which.Spec[[i]]<-k
			        is.M[i] <- 1 #convention to avoid centering objects from mgcv	  
				}				
				if (predInd){ 			        			        			        		        
			        if (!dataPr){
     			        trmWD<-deparse(trms2[1][[2]], backtick = TRUE, width.cutoff = 500)
     			        trmWD<-substr(trmWD,1,nchar(trmWD)-1)
     			        trmWD<-paste(trmWD,", data=mvrmObj$data)")
				    }
     			    evalS<-eval(parse(text=trmWD))
     			    NL<-1
     			    if (!evalS[[1]]$by=="NA") NL<-max(1,length(levels(with(data,eval(as.name(evalS[[1]]$by))))))
     			    Design.Xt<-NULL     			    
     			    for (j in 1:NL) Design.Xt<-cbind(Design.Xt,mgcv::PredictMat(evalS[[j]],data=data.frame(data)))  
				}			    			    			   			    
			    remove.cols<-which(c(apply(Design.Xt,2,sd)==0))
			    colnames(Design.Xt)<-paste(labels[i],1:NCOL(Design.Xt),sep=".")
			}
			Design.X<-cbind(Design.X,Design.Xt)			
			assign<-c(assign,rep(i,NCOL(Design.Xt)))
        }
	}
	if (missing(meanVector)) meanVector<-apply(as.matrix(Design.X),2,mean)
	if (missing(indicator)){ 
	    unique.values<-unlist(lapply(apply(Design.X,2,unique),length))
	    indicator<-c(unique.values<=2)	    
	}
	if (sum(is.M)){
	    mat.arg<-which(is.M==1)
	    for (i in 1:length(mat.arg))
		    indicator[which(assign==mat.arg[i])]<-TRUE    
	}
	Design.X[,!indicator]<-Design.X[,!indicator]-matrix(1,nrow=n)%*%matrix(meanVector[!indicator],nrow=1)       
	together<-list()
	e<-0
	if (!is.null(nFactors)){
        ptns<-which(lapply(vars,length) > 1) #candidates for together
        if (length(ptns) > 0)
		    for (i in 1:nFactors)
	            for (j in ptns)
	                if (sum(vars[[i]] %in% vars[[j]]) > 0 && (!i==j)) together[[e<-e+1]]<-c(i,j)    
	}   
    return(list(y=y,X=Design.X,assign=assign,Rknots=Rknots,meanVector=meanVector,
                indicator=indicator,labels=labels,count=count,vars=vars,is.D=is.D,
                which.Spec=which.Spec,formula.terms=formula.terms,together=together))
}

match.call.defaults <- function(...) {
    call <- evalq(match.call(expand.dots = FALSE), parent.frame(1))
    formals <- evalq(formals(), parent.frame(1))
    for(k in setdiff(names(formals), names(call)))
        call[k] <- list( formals[[k]] )
    match.call(sys.function(sys.parent()), call)
}

mvrm <- function(formula,data=list(),sweeps,burn=0,thin=1,seed,StorageDir,
                 c.betaPrior="IG(0.5,0.5*n)",c.alphaPrior="IG(1.1,1.1)",
                 pi.muPrior="Beta(1,1)",pi.sigmaPrior="Beta(1,1)",sigmaPrior="HN(2)",...){
    if (thin <= 0) thin=1
    thin <- as.integer(thin)
    sweeps <- as.integer(sweeps)
    burn <- as.integer(burn)
    if (missing(sweeps)) stop("provide sweeps argument")
    nSamples <-0
    if (sweeps > 0) nSamples <- length(seq(1,(sweeps-burn),by=thin))
    # Match call
    call <- match.call(expand.dots = FALSE)
    call2 <- match.call.defaults()
    #Data
    if (!is.list(data) && !is.data.frame(data)) 
        data <- as.data.frame(data)       
    #Formula
    if (length(as.Formula(formula))[1] > 1) stop("more than one response provided") 
    if (length(as.Formula(formula))[2] > 2) stop("more than two regression models provided")
    if (length(as.Formula(formula))[2] == 1) formula <- as.Formula(formula, ~1)
    formula.save <- formula
    formula.v<-formula(as.Formula(formula),lhs=0,rhs=2)
    formula<-formula(as.Formula(formula),lhs=1,rhs=1)
    # Design matrices, response, indicators    
    XYK<-DM(formula=formula,data=data)    
    Y<-XYK$y    
    n<-length(Y)
    X<-XYK$X
    Xknots<-XYK$Rknots
    storeMeanVectorX<-XYK$meanVector
    storeIndicatorX<-XYK$indicator
    LG<-NCOL(X)-1 
    vecLG<-table(XYK$assign)[-1]
    NG<-length(vecLG)
    cusumVecLG<-c(0,cumsum(vecLG))
    assignX<-XYK$assign
    labelsX<-XYK$labels
    countX<-XYK$count
    varsX<-XYK$vars
    is.Dx<-XYK$is.D
    which.SpecX<-XYK$which.Spec
    formula.termsX<-XYK$formula.terms
    togetherX<-XYK$together
    #
    ZK<-DM(formula=formula.v,data=data,n=n)
    Z<-ZK$X  
    Zknots<-ZK$Rknots
    storeMeanVectorZ<-ZK$meanVector
    storeIndicatorZ<-ZK$indicator       
    LD<-NCOL(Z)-1 
    vecLD<-table(ZK$assign)[-1]
    ND<-length(vecLD)
    cusumVecLD<-c(0,cumsum(vecLD))
    if (LD==0) MVLD <- 0
    if (LD > 0) MVLD<-max(vecLD)
    assignZ<-ZK$assign
    labelsZ<-ZK$labels
    countZ<-ZK$count
    varsZ<-ZK$vars
    is.Dz<-ZK$is.D
    which.SpecZ<-ZK$which.Spec
    formula.termsZ<-ZK$formula.terms
    togetherZ<-ZK$together
    #       
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
                data=data,Y=Y,X=X,Xknots=Xknots,Z=Z,Zknots=Zknots,LG=LG,LD=LD,
                mcpar=c(as.integer(burn+1),as.integer(seq(from=burn+1,by=thin,length.out=nSamples)[nSamples]),as.integer(thin)),
                nSamples=nSamples,storeMeanVectorX=storeMeanVectorX,storeMeanVectorZ=storeMeanVectorZ,
                f=out[[17]][1],g=out[[18]][1],h=out[[19]][1:ND],
                DIR=StorageDir,deviance=out[[33]][1]/nSamples,nullDeviance=-2*logLik(lm(Y ~ 1)),
                storeIndicatorX=storeIndicatorX,storeIndicatorZ=storeIndicatorZ,assignX=assignX,
                assignZ=assignZ,labelsX=labelsX,labelsZ=labelsZ,countX=countX,countZ=countZ,
                varsX=varsX,varsZ=varsZ,is.Dx=is.Dx,is.Dz=is.Dz,which.SpecX=which.SpecX,
                which.SpecZ=which.SpecZ,formula.termsX=formula.termsX,formula.termsZ=formula.termsZ,
                togetherX=togetherX,togetherZ=togetherZ)
    class(fit) <- 'mvrm'
    return(fit)
}

mvrm2mcmc <- function(mvrmObj,labels){
    all.labels <- c("alpha","calpha","cbeta","delta","beta","gamma","sigma2")
    mtch<-match(labels,all.labels)
	R<-NULL
    if (any(mtch==1) && mvrmObj$LD > 0) R<-matrix(unlist(read.table(paste(mvrmObj$DIR,"BNSP.alpha.txt",sep=""))),ncol=mvrmObj$LD,dimnames=list(c(),colnames(mvrmObj$Z)[-1]))
    if (any(mtch==2)) R<-matrix(unlist(read.table(paste(mvrmObj$DIR,"BNSP.calpha.txt",sep=""))),ncol=1,dimnames=list(c(),c("c_alpha")))
    if (any(mtch==3)) R<-matrix(unlist(read.table(paste(mvrmObj$DIR,"BNSP.cbeta.txt",sep=""))),ncol=1,dimnames=list(c(),c("c_beta")))
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
	if (!is.function(centre)) stop("centre must be a function, usually mean or median")
	if (length(quantiles)==1) quantiles <- c(quantiles,1-quantiles)
    if (length(quantiles) > 2) stop("up to two quantiles")
    if (!is.null(quantiles)) quantiles <- unique(sort(quantiles))
    grid<-round(grid)
	MEAN<-0
	if (model=="mean" || model=="both") MEAN<-1
    STDEV<-0
    if (model=="stdev" || model=="both") STDEV<-1    
    if (missing(term)) term<-1    
    countX<-0
    if (MEAN){
        if (is.character(term)){    
            int.labelX<-grep(term,mvrmObj$labelsX,fixed=TRUE)
            k<-1
            while (length(int.labelX)==0 && (nchar(term)>k)){
                int.labelX<-grep(substr(term,1,nchar(term)-k),mvrmObj$labelsX,fixed=TRUE)
                k<-k+1 
	        }
	    }
	    if (is.numeric(term)) int.labelX<-term    
        if (length(int.labelX)==0) stop(cat(cat("`term` in mean model should be an integer between 1 and", length(mvrmObj$labelsX),"or one of: "), cat(mvrmObj$labelsX,sep=", ",fill=TRUE))) 
        if (int.labelX > length(mvrmObj$labelsX)) stop(cat(cat("`term` in mean model should be an integer between 1 and", length(mvrmObj$labelsX),"Or one of:"), cat(mvrmObj$labelsX,sep=", ",fill=TRUE)))         
        T<-0        
        if (length(mvrmObj$togetherX)>0) 
            for (i in 1:length(mvrmObj$togetherX)) 
                if (int.labelX %in% mvrmObj$togetherX[[i]])
                    T<-i         
        if (T > 0) int.labelX<-mvrmObj$togetherX[[i]]                
        countX<-mvrmObj$countX[int.labelX[length(int.labelX)]]
    }    
    countZ<-0
    if (STDEV){
	    if (is.character(term)){    
            int.labelZ<-grep(term,mvrmObj$labelsZ,fixed=TRUE)
            k<-1
            while (length(int.labelZ)==0 && (nchar(term)>k)){
                int.labelZ<-grep(substr(term,1,nchar(term)-k),mvrmObj$labelsZ,fixed=TRUE)
                k<-k+1 
	        }
	    }
	    if (is.numeric(term)) int.labelZ<-term    
	    if (length(int.labelZ)==0) stop(cat(cat("`term` in stdev model should be an integer between 1 and", length(mvrmObj$labelsZ),"or one of: "), cat(mvrmObj$labelsZ,sep=", ",fill=TRUE))) 
        if (int.labelZ > length(mvrmObj$labelsZ)) stop(cat(cat("`term` in stdev model should be an integer between 1 and", length(mvrmObj$labelsZ),"or one of: "), cat(mvrmObj$labelsZ,sep=", ",fill=TRUE))) 
        T<-0        
        if (length(mvrmObj$togetherZ)>0) 
            for (i in 1:length(mvrmObj$togetherZ)) 
                if (int.labelZ %in% mvrmObj$togetherZ[[i]])
                    T<-i         
        if (T > 0) int.labelZ<-mvrmObj$togetherZ[[i]]                
        countZ<-mvrmObj$countZ[int.labelZ[length(int.labelZ)]]        
	}	
    if (model=="both" && (!countZ==countX && min(countZ,countX)>0)) 
        stop("don't know how to arrange ggplots and 3d-plots in a single grid; plot one of `model` at a time")    
    if (MEAN){
        int.label<-int.labelX
        count<-countX
        vars<-mvrmObj$varsX[[int.label[length(int.label)]]]
        label<-mvrmObj$labelsX[int.label[length(int.label)]]
        is.D<-mvrmObj$is.Dx[[int.label[length(int.label)]]]
        formula.term<-mvrmObj$formula.termsX[int.label]
		if (!grepl("knots",formula.term) && mvrmObj$which.SpecX[[int.label[length(int.label)]]]>0){
    	    z<-length(formula.term) #use this to fix both sm and s smooths 
     		formula.term[z]<-substr(formula.term[z],1,nchar(formula.term[z])-1)
     		formula.term[z]<-paste(formula.term[z],", knots=knots)")
		}
        ML<-which(mvrmObj$assignX %in% int.label)
        if ((intercept || sum(is.D)) & MEAN) ML<-c(1,ML)
        if (count==1){  
			if (!is.D){			
	            min1<-min(with(mvrmObj$data,eval(as.name(vars[1]))))
			    max1<-max(with(mvrmObj$data,eval(as.name(vars[1]))))
			    newR1<-seq(min1,max1,length.out=grid)
                newData<-data.frame(newR1)
                colnames(newData)<-vars                      
                whichKnots <- mvrmObj$which.SpecX[[int.label]] 
                if (length(mvrmObj$Xknots)>0) {Dstar<-data.frame(knots=mvrmObj$Xknots[[whichKnots]])}else{Dstar<-NULL}                
                DsM<-DM(formula=reformulate(formula.term),data=newData,n=NROW(newData),knots=Dstar)$X
                if (! 1%in%ML) DsM<-DsM[,-1]               
		    }
	        if (is.D){ 
	            newData<-data.frame(unique(with(mvrmObj$data,eval(as.name(vars[1])))))
                colnames(newData)<-vars                                      
                DsM<-DM(formula=reformulate(formula.term),data=newData,n=NROW(newData),knots=NULL)$X
	        }      
		}
	    if (count==2){    
	        if (sum(is.D)){	            
	            dv<-which(is.D==1)
	            cv<-which(is.D==0)	        	        
	            min1<-min(with(mvrmObj$data,eval(as.name(vars[cv]))))
				max1<-max(with(mvrmObj$data,eval(as.name(vars[cv]))))
                newR1<-seq(min1,max1,length.out=grid)
                newR2<-unique(with(mvrmObj$data,eval(as.name(vars[dv]))))
                newData<-expand.grid(newR1,newR2)
                colnames(newData)<-vars[c(cv,dv)] 
                whichKnots <- mvrmObj$which.SpecX[[int.label[length(int.label)]]] 
                Dstar<-mvrmObj$Xknots[[whichKnots]]
                DsM<-DM(formula=reformulate(formula.term),data=newData,n=NROW(newData),knots=Dstar)$X                                              
			}
            if (sum(is.D)==0){
				min1<-min(with(mvrmObj$data,eval(as.name(vars[1]))))
				min2<-min(with(mvrmObj$data,eval(as.name(vars[2]))))
				max1<-max(with(mvrmObj$data,eval(as.name(vars[1]))))
				max2<-max(with(mvrmObj$data,eval(as.name(vars[2]))))
                newR1<-seq(min1,max1,length.out=grid)
                newR2<-seq(min2,max2,length.out=grid)
                newData<-expand.grid(newR1,newR2)
                colnames(newData)<-vars                      
                whichKnots <- mvrmObj$which.SpecX[[int.label]] 
                Dstar<-mvrmObj$Xknots[[whichKnots]]
                DsM<-DM(formula=reformulate(formula.term),data=newData,n=NROW(newData),knots=Dstar)$X
		    }                    
        }
        fitM<-matrix(0,nrow=mvrmObj$nSamples,ncol=NROW(DsM))
		etaFN <- file.path(paste(mvrmObj$DIR,"BNSP.beta.txt",sep="")) 
        eFile<-file(etaFN,open="r")
		for (i in 1:mvrmObj$nSamples){
            eta<-scan(eFile,what=numeric(),n=mvrmObj$LG+1,quiet=TRUE)
            fitM[i,]<-as.matrix(DsM)%*%matrix(c(eta[ML]))
            if (centreEffects) fitM[i,]<-fitM[i,]-mean(fitM[i,])
		}
		close(eFile)
		if (count==1 && !is.D){
		    centreM<-apply(fitM,2,centre)
		    QM<-NULL
		    if (!is.null(quantiles)) QM<-drop(t(apply(fitM,2,quantile,probs=quantiles,na.rm=TRUE)))
		    newX<-newData
		    dataM<-data.frame(cbind(newX,centreM,QM))		
		    nms <- c(vars,"centreM")
		    if (!is.null(quantiles)) nms<-c(nms,"QLM","QUM")
		    colnames(dataM) <- nms
		    
		    dataRug<-data.frame(b=with(mvrmObj$data,eval(as.name(vars))))
		    plotElM<-c(geom_line(aes_string(x=as.name(vars),y=centreM),col=4,alpha=0.5),
                       geom_rug(mapping=aes(x=dataRug$b),data=dataRug,alpha=0.3),                       
                       list(ylab(label)))
            
            if (!is.null(quantiles))
                plotElM<-c(geom_ribbon(data=dataM,aes_string(x=vars, ymin="QLM", ymax="QUM"),alpha=0.2),plotElM)		    
		    ggM<-ggplot(data=dataM)
	        plotM<-ggM + plotElM + ggtitle("mean") + plotOptions
	    }
	    if (count==1 && is.D){  
			lvs<-levels(with(mvrmObj$data,eval(as.name(vars[1]))))  
			if (is.null(lvs)) lvs<-unique(with(mvrmObj$data,eval(as.name(vars[1]))))  
	        df<-data.frame(x=rep(lvs,each=mvrmObj$nSamples),y=c(fitM))	        
            plotElM<-c(geom_boxplot(),list(xlab(label),ylab("")))
            ggM<-ggplot(aes(x=factor(x),y=df$y),data=df)
	        plotM<-ggM + plotElM + ggtitle("mean") + plotOptions
	    }	   
	    if (count==2 && sum(is.D)==1){		    
		    centreM<-apply(fitM,2,centre)
		    QM<-NULL
		    if (!is.null(quantiles)) QM<-drop(t(apply(fitM,2,quantile,probs=quantiles,na.rm=TRUE)))		      
		    disc.var<-vars[which(is.D==1)]
		    cont.var<-vars[which(is.D==0)]			    		    	                        
		    lvs<-levels(with(mvrmObj$data,eval(as.name(disc.var))))		    		    		   
		    nms<-paste(disc.var,lvs[-1],sep="")	    		    
		    dataM<-data.frame(newData,centreM,QM)   		    
		    nms <- c(colnames(newData),"centreM")		    
		    if (!is.null(quantiles)) nms<-c(nms,"QLM","QUM")
		    colnames(dataM) <- nms		    
		    DG<-data.frame(b=with(mvrmObj$data,eval(as.name(cont.var))),c=with(mvrmObj$data,eval(as.name(disc.var))))		    		    
		    plotElM<-c(geom_line(aes_string(x=cont.var,y=centreM,
		               group=disc.var,colour=disc.var,linetype=disc.var),alpha=0.80),
                       geom_rug(mapping=aes(x=DG$b,group=c,colour=c,linetype=c),data=DG,alpha=0.3),
                       list(ylab(label)))
            if (!is.null(quantiles))
                plotElM<-c(geom_ribbon(data=dataM,aes_string(x=cont.var,ymin="QLM",ymax="QUM",
                                       group=disc.var,fill=as.name(disc.var)),alpha=0.2),plotElM)
		    ggM<-ggplot(data=dataM)
	        plotM<-ggM + plotElM + ggtitle("mean") + plotOptions
	    }
	    if (count==2 && sum(is.D)==0){
		    centreM<-apply(fitM,2,centre)		    
		    if (static){		        
		        defaultList<-list(x=as.numeric(newR1),y=as.numeric(newR2),z=matrix(centreM,length(newR1),length(newR2)),colvar=matrix(centreM,length(newR1),length(newR2)))
                along="xy"; 
                space=0.6; 
                optionalList<-list(xlab=vars[1],ylab=vars[2],zlab=label,along=along,space=space,add=FALSE,bty="g",main="mean")
                allOptions<-c(defaultList,plotOptions,optionalList)
                if (STDEV==1) par(mfrow=c(1,2))
                do.call("ribbon3D",allOptions[!duplicated(names(allOptions))])
		    }
            if (!static){             
                a<-as.matrix(cbind(newData,centreM))
                if (is.null(plotOptions$col)) col=rainbow(16,2/3)
                else col<-plotOptions$col
                plotCentreM <- centreM
                if (min(centreM)<=0) plotCentreM <- centreM + abs(min(centreM)) + 1
                ra<-ceiling(length(col)*plotCentreM/max(plotCentreM))
                defaultList<-list(x=a,col=col[ra])
                optionalList<-list(size=0.4,bg=1,axisLabels=c(vars[1],label,vars[2]),main="mean")
                allOptions<-c(defaultList,plotOptions,optionalList)
                plotM<-do.call("scatterplot3js",allOptions[!duplicated(names(allOptions))])
			}
		}
	}
	if (STDEV){
        int.label<-int.labelZ
        count<-countZ
        vars<-mvrmObj$varsZ[[int.label[length(int.label)]]]
        label<-mvrmObj$labelsZ[int.label[length(int.label)]]
        is.D<-mvrmObj$is.Dz[[int.label[length(int.label)]]]
        formula.term<-mvrmObj$formula.termsZ[int.label]
        if (!grepl("knots",formula.term) && mvrmObj$which.SpecZ[[int.label[length(int.label)]]]>0){
    	    z<-length(formula.term) #use this to fix both sm and s smooths 
     		formula.term[z]<-substr(formula.term[z],1,nchar(formula.term[z])-1)
     		formula.term[z]<-paste(formula.term[z],", knots=knots)")
		}
        VL<-which(mvrmObj$assignZ %in% int.label)
        if (count==1){              
	        if (!is.D){
	            min1<-min(with(mvrmObj$data,eval(as.name(vars[1]))))
			    max1<-max(with(mvrmObj$data,eval(as.name(vars[1]))))
			    newR1<-seq(min1,max1,length.out=grid)
                newData<-data.frame(newR1)
                colnames(newData)<-vars                      
                whichKnots <- mvrmObj$which.SpecZ[[int.label]] 
                if (length(mvrmObj$Zknots)>0) {Dstar<-data.frame(mvrmObj$Zknots[[whichKnots]])}else{Dstar<-NULL}
                DsV<-DM(formula=reformulate(formula.term),data=newData,n=NROW(newData),knots=Dstar)$X[,-1]
		    }	        
	        if (is.D){ 
	            newData<-data.frame(unique(with(mvrmObj$data,eval(as.name(vars[1])))))
                colnames(newData)<-vars                                      
                DsV<-DM(formula=reformulate(formula.term),data=newData,n=NROW(newData),knots=NULL)$X[,-1]
			}      
		}
		if (count==2){
	        if (sum(is.D)){
	            dv<-which(is.D==1)
	            cv<-which(is.D==0)	        	        
	            min1<-min(with(mvrmObj$data,eval(as.name(vars[cv]))))
				max1<-max(with(mvrmObj$data,eval(as.name(vars[cv]))))
                newR1<-seq(min1,max1,length.out=grid)
                newR2<-unique(with(mvrmObj$data,eval(as.name(vars[dv]))))
                newData<-expand.grid(newR1,newR2)
                colnames(newData)<-vars[c(cv,dv)]                      
                whichKnots <- mvrmObj$which.SpecZ[[int.label[length(int.label)]]] 
                Dstar<-mvrmObj$Zknots[[whichKnots]]
                DsV<-DM(formula=reformulate(formula.term),data=newData,n=NROW(newData),knots=Dstar)$X
                DsV<-DsV[,-1]
			}		        
            if (sum(is.D)==0){
                min1<-min(with(mvrmObj$data,eval(as.name(vars[1]))))
				min2<-min(with(mvrmObj$data,eval(as.name(vars[2]))))
				max1<-max(with(mvrmObj$data,eval(as.name(vars[1]))))
				max2<-max(with(mvrmObj$data,eval(as.name(vars[2]))))
                newR1<-seq(min1,max1,length.out=grid)
                newR2<-seq(min2,max2,length.out=grid)
                newData<-expand.grid(newR1,newR2)
                colnames(newData)<-vars
                whichKnots <- mvrmObj$which.SpecZ[[int.label]] 
                Dstar<-mvrmObj$Zknots[[whichKnots]]
                DsV<-DM(formula=reformulate(formula.term),data=newData,n=NROW(newData),knots=Dstar)$X    
                DsV<-DsV[,-1]
		    }                                
        }
        fitV<-matrix(0,nrow=mvrmObj$nSamples,ncol=NROW(DsV))
		alphaFN <- file.path(paste(mvrmObj$DIR,"BNSP.alpha.txt",sep="")) 
        aFile<-file(alphaFN,open="r")
		sigma2FN <- file.path(paste(mvrmObj$DIR,"BNSP.sigma2.txt",sep="")) 
        s2File<-file(sigma2FN,open="r")
        for (i in 1:mvrmObj$nSamples){
        	alpha<-scan(aFile,what=numeric(),n=mvrmObj$LD,quiet=TRUE)
		    if (intercept) s2<-scan(s2File,what=numeric(),n=1,quiet=TRUE)
		    else s2=1
            fitV[i,]<-sqrt(s2*exp(as.matrix(DsV)%*%matrix(c(alpha[VL-1]))))
            if (centreEffects) fitV[i,]<-fitV[i,]/mean(fitV[i,])
		}
		close(aFile)
		close(s2File)		
		if (count==1 && !is.D){
		    centreV<-apply(fitV,2,centre)
		    QV<-NULL
		    if (!is.null(quantiles)) QV<-drop(t(apply(fitV,2,quantile,probs=quantiles,na.rm=TRUE)))
		    newV<-newData#DsV[,vars]+mean(with(mvrmObj$data,eval(as.name(vars))))
		    dataV<-data.frame(newV,centreV,QV)
            nms<-c(vars,"centreV")
		    if (!is.null(quantiles)) nms<-c(nms,"QLV","QUV")
		    colnames(dataV) <- nms        
		    dataRug<-data.frame(b=with(mvrmObj$data,eval(as.name(vars))))
		    plotElV<-c(geom_line(aes_string(x=as.name(vars),y=centreV),col=4,alpha=0.5),                     
                       geom_rug(mapping=aes(x=dataRug$b),data=dataRug,alpha=0.3),
                       list(ylab(label)))
            if (!is.null(quantiles))
                plotElV<-c(geom_ribbon(data=dataV,aes_string(x=vars, ymin="QLV", ymax="QUV"),alpha =0.2),plotElV)   
		    ggV<-ggplot(data=dataV)
		    plotV<-ggV + plotElV + ggtitle("st dev") + plotOptions
		}
		if (count==1 && is.D){	        	        
	        lvs<-levels(with(mvrmObj$data,eval(as.name(vars[1]))))  
			if (is.null(lvs)) lvs<-unique(with(mvrmObj$data,eval(as.name(vars[1]))))  
	        df<-data.frame(x=rep(lvs,each=mvrmObj$nSamples),y=c(fitV))
            plotElV<-c(geom_boxplot(),list(xlab(label),ylab("")))
            ggV<-ggplot(data=df,aes(x=factor(x),y=df$y))
	        plotV<-ggV + plotElV + ggtitle("st dev") + plotOptions
	    }
	    if (count==2 && sum(is.D)==1){		    						    		    
		    centreV<-apply(fitV,2,centre)
		    QV<-NULL
		    if (!is.null(quantiles)) QV<-drop(t(apply(fitV,2,quantile,probs=quantiles,na.rm=TRUE)))		
		    disc.var<-vars[which(is.D==1)]
		    cont.var<-vars[which(is.D==0)]
		    lvs<-levels(with(mvrmObj$data,eval(as.name(disc.var))))
		    nms<-paste(disc.var,lvs[-1],sep="")		    		    
		    dataV<-data.frame(newData,centreV,QV)
		    nms <- c(colnames(newData),"centreV")
		    if (!is.null(quantiles)) nms<-c(nms,"QLV","QUV")
		    colnames(dataV) <- nms
   		    DG<-data.frame(b=with(mvrmObj$data,eval(as.name(cont.var))),c=with(mvrmObj$data,eval(as.name(disc.var))))
		    plotElV<-c(geom_line(aes_string(x=cont.var,y=centreV,
		               group=disc.var,linetype=disc.var,colour=disc.var),alpha=0.80),
                       geom_rug(mapping=aes(x=DG$b,group=DG$c,colour=DG$c,linetype=DG$c),data=DG,alpha=0.3),
                       list(ylab(label)))
            if (!is.null(quantiles))
                plotElV<-c(geom_ribbon(data=dataV,aes_string(x=cont.var,ymin="QLV",ymax="QUV",
                                       group=disc.var,fill=as.name(disc.var)),alpha=0.2),plotElV)
		    ggV<-ggplot(data=dataV)
	        plotV<-ggV + plotElV + ggtitle("st dev") + plotOptions
	    }
		if (count==2 && sum(is.D)==0){
		    centreV<-apply(fitV,2,centre)
		    if (static){		        				
				defaultList<-list(x=as.numeric(newR1),y=as.numeric(newR2),z=matrix(centreV,length(newR1),length(newR2)),colvar=matrix(centreV,length(newR1),length(newR2)))
                along="xy"; 
                space=0.6; 
                optionalList<-list(xlab=vars[1],ylab=vars[2],zlab=label,along=along,space=space,add=FALSE,bty="g",main="st dev")
                allOptions<-c(defaultList,plotOptions,optionalList)
                do.call("ribbon3D",allOptions[!duplicated(names(allOptions))])
		    }
            if (!static){                
                a<-as.matrix(cbind(newData,centreV))
                if (is.null(plotOptions$col)) col=rainbow(16,2/3)
                else col<-plotOptions$col
                ra<-ceiling(length(col)*centreV/max(centreV))
                defaultList<-list(x=a,col=col[ra])
                optionalList<-list(size=0.4,bg=1,axisLabels=c(vars[1],label,vars[2]),main="st dev")
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
    if (missing(newdata) || is.null(newdata)) newdata<-object$data
    if (length(newdata)==0) stop("no data found")
    newdata <- as.data.frame(newdata)
    terms.reform<-NULL
    k<-0
    for (i in 1:length(object$formula.termsX)){
        term<-object$formula.termsX[i]
        if (!i %in%  which(unlist(object$which.SpecX)==-99)){ 
            k<-k+1
            if (!grepl("knots",term)){
                term<-substr(term,1,nchar(term)-1)
     		    term<-paste(term,",knots= knots[[",k,"]])")
			}
	    }
        terms.reform<-c(terms.reform,term)
    }
    formula2<-reformulate(terms.reform)
    if (length(object$data)>0){
        nd<-object$data[0,match(colnames(newdata),colnames(object$data)),drop=FALSE]
        nd[1:NROW(newdata),] <- newdata
    }else{nd<-newdata}
    X<-DM(formula=formula2,data=nd,n=NROW(newdata),knots=object$Xknots,predInd=TRUE,meanVector=object$storeMeanVectorX,indicator=object$storeIndicatorX,mvrmObj=object)$X
	fitM<-matrix(0,nrow=object$nSamples,ncol=NROW(newdata))
    etaFN <- file.path(paste(object$DIR,"BNSP.beta.txt",sep="")) 
    eFile<-file(etaFN,open="r")
	for (i in 1:object$nSamples){
        eta<-scan(eFile,what=numeric(),n=object$LG+1,quiet=TRUE)
        fitM[i,]<-c(as.matrix(X)%*%matrix(c(eta)))
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
	    terms.reform<-NULL
        k<-0
        if (length(object$formula.termsZ)>0){
            for (i in 1:length(object$formula.termsZ)){
                term<-object$formula.termsZ[i]
                if (!i %in%  which(unlist(object$which.SpecZ)==-99)){ 
                    k<-k+1
                    if (!grepl("knots",term)){
                       term<-substr(term,1,nchar(term)-1)
     		           term<-paste(term,",knots= knots[[",k,"]])")
			        }
	            }
                terms.reform<-c(terms.reform,term)
            }
	    }
        if (is.null(terms.reform)) {formula2<-~1}else{formula2<-reformulate(terms.reform)}                
        Z<-DM(formula=formula2,data=nd,n=NROW(newdata),knots=object$Zknots,predInd=TRUE,meanVector=object$storeMeanVectorZ,indicator=object$storeIndicatorZ,mvrmObj=object)$X
        Z<-Z[,-1]
        fitV<-matrix(0,nrow=object$nSamples,ncol=NROW(newdata))     
		alphaFN <- file.path(paste(object$DIR,"BNSP.alpha.txt",sep="")) 
        aFile<-file(alphaFN,open="r")
		sigma2FN <- file.path(paste(object$DIR,"BNSP.sigma2.txt",sep="")) 
        s2File<-file(sigma2FN,open="r")
        for (i in 1:object$nSamples){
        	alpha<-scan(aFile,what=numeric(),n=object$LD,quiet=TRUE)
		    s2<-scan(s2File,what=numeric(),n=1,quiet=TRUE)
            fitV[i,]<-sqrt(s2*exp(as.matrix(Z)%*%matrix(c(alpha))))            
		}
		close(aFile)
		close(s2File)	
		fitSD<-apply(fitV,2,mean)
		predictions<-cbind(fit,fit-qnorm(1-(1-level)/2)*fitSD,fit+qnorm(1-(1-level)/2)*fitSD,fitSD)
		colnames(predictions) <- c("fit", "lwr", "upr", "fit.sd")
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
        cat("\nMean model - marginal inclusion probabilities\n")
        print(apply(G,2,mean), digits=5)
    }
    if (x$LD > 0){
        cat("\nVariance model - marginal inclusion probabilities\n")
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
    if (nModels>0){
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
      TrueDim<-min(nModels,dim(g)[1])
      print(g[1:TrueDim,])
      cat("Displaying", TrueDim, "models of the",NROW(g),"visited\n")
      cat(TrueDim,"models account for",sub(" ","",paste(g$cumulative[TrueDim],"%")),"of the posterior mass\n")
    }
}
