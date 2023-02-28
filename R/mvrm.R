# s, te, ti, sm, DM, match.call.defaults, quiet  

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
    if (d > 2) stop("Up to bivariate covariates supported; other arguments are k, knots and bs")
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
    type<-1
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
            knots<-unique(round(quantile(x2,knots,type=type),5))           
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
        if (!is.null(knots)) nknots<-c(dim(knots))
        if (is.null(knots) && length(nknots)==1) nknots<-rep(nknots,2)
        if (is.null(knots)){
			knots1<-seq(from = 0, to = 1, length = nknots[1] + 2)[-c(1, nknots[1] + 2)]
            knots1<-seq(from = 0, to = 1, length = nknots[1])
            knots2<-seq(from = 0, to = 1, length = nknots[2] + 2)[-c(1, nknots[2] + 2)]
            knots2<-seq(from = 0, to = 1, length = nknots[2])
            if (dim(x2)[2] == 2){
                knots1<-unique(round(quantile(x2[,1],knots1,type=type),5))
                knots2<-unique(round(quantile(x2[,2],knots2,type=type),5))
                knots<-as.matrix(expand.grid(knots1,knots2))
			}else if (is.D[1] == 1){
			    knots1<-unique(x2[,-dim(x2)[2]])
			    knots2<-unique(round(quantile(x2[,dim(x2)[2],drop=FALSE],knots2,type=type),5))
			    knots <- cbind(knots1[rep(1:nrow(knots1), length(knots2)), ], rep(knots2, each = nrow(knots1)))
			}else if (is.D[2] == 1){
				knots1<-unique(round(quantile(x2[,1],knots1,type=type),5))
			    knots2<-unique(x2[,-1])
			    knots <- cbind(rep(knots1, nrow(knots2)), knots2[rep(1:nrow(knots2), length(knots1)), ])
			}
			colnames(knots)<-colnames(x2)
		}
        if (bs=="rd"){
            for (i in 1:NROW(knots)){
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

sinusoid<-function(..., harmonics = 1, amplitude = 1, period = 1, knots = NULL, breaks = NULL){    
    pf <- parent.frame()
    vars<-as.list(substitute(list(...)))[-1]
    d<-length(vars)
    if (d > 1) stop("univariate effects supported; other arguments are harmonics, amplitude, period")
    if (amplitude > 1 && !is.null(breaks)) stop("cannot do varying amplitude with breaks")
    term<-NULL
    for (i in 1:d){
	    term[i]<-deparse(vars[[i]],backtick=TRUE,width.cutoff=500)
	    term[i] <- attr(terms(reformulate(term[i])),"term.labels")
    }
    if (length(unique(term))!=d) stop("repeated variables are not permitted")
    if (harmonics <= 0) stop("harmonics should be a positive integer")    
    if (harmonics!=round(harmonics)) {warning("number of harmonics should be a positive integer; it has been rounded")
		                              harmonics <- round(harmonics)}    
    label<-paste("sinusoid(",term[1],sep="")
    if (d>1) for (i in 2:d) label<-paste(label,",",term[i],sep="")
    label<-paste(label,")",sep="")
    is.D<-rep(0,d)
    X<-NULL
    Xnames<-NULL
    for (k in 1:harmonics){ 
        X<-cbind(X,sin(2 * k * pi * eval(vars[[1]],pf)/period))            
        X<-cbind(X,cos(2 * k * pi * eval(vars[[1]],pf)/period))
        Xnames<-c(Xnames,paste(paste("sin(",2*k,sep=""),"pi", vars, "/ p)",sep=" "),paste(paste("cos(",2*k,sep=""),"pi", vars, "/ p)",sep=" "))
	}                    
    X<-data.frame(X)
    colnames(X)[1:NCOL(X)]<-Xnames
    Dynamic <- NULL
    if (amplitude > 1){         
        XK<-sm(eval(parse(text=vars),pf),k=amplitude,knots=knots)        
        Design.Xamp<-XK$X 
        knots<-XK$knots
        Dynamic<-sweep(Design.Xamp, MARGIN=1, X[,1], `*`)
        interMat2<-Dynamic/sqrt(2)
        for(j in 2:NCOL(X)){ 
            interMat<-sweep(Design.Xamp, MARGIN=1, X[,j], `*`)
            interMat2<-interMat2+interMat/sqrt(2)
            Dynamic<-cbind(Dynamic,interMat)
		}
        X<-interMat2
        colnames(X)[1:NCOL(X)] <- paste("sm(amplitude)",seq(1:(amplitude+1)),sep=".")
	}
    XK<-list(X=X,knots=knots,count=d,vars=unlist(term),is.D=is.D,label=label,
             amplitude=amplitude,harmonics=harmonics,Dynamic=Dynamic,breaks=breaks,
             period=period)
    return(XK)
}

DM<-function(formula,data,n,knots=NULL,predInd=FALSE,meanVector,indicator,mvrmObj,centre=centre){
    is.M <- NULL
    count <- NULL
    k <- 0
    Rknots <- list()
    vars <- list()
    is.D <- list()
    which.Spec<-list()
    specials <- c('sm','s','te','ti','sinusoid')
    trms<-terms.formula(formula,specials=specials)
    attr(trms,"intercept") <- 1
    nFactors <- dim(attr(trms,"factors"))[2]
    if (attr(trms,"response")){
        y <- with(data,eval(attr(trms,"variables")[[2]]))
        n <- length(y)
	}else{
	    y <- NULL
	    n <- n
	}
	labels <- colnames(attr(trms,"factors"))
	formula.terms <- attr(trms,"term.labels")
	Design.X <- matrix(1,ncol=1,nrow=n)
    colnames(Design.X) <- "(Intercept)"
    assign <- 0
    assign2 <- 0 # for 1-hot encoding    
    assign3 <- 0 # for breaking up couples of sinusoidal terms
    repsX <- NULL # repeats for fixing differences between NG abd NG2
    isSin <- 0    
    amplitude <- 0
    harmonics <- 0
    startSin <- 0
    Dynamic <- NULL
    breaks <- NULL
    period <- 0  
    if (!is.null(nFactors)) isSin<-rep(0,nFactors)
    if (!is.null(nFactors)){
        for (i in 1:nFactors){
            trms2<-drop.terms(trms, dropx = -i)
            is.spec<-sum(unlist(lapply(attr(trms2,"specials"),sum))) > 0
            is.sm<-!is.null(attr(trms2,"specials")$sm)
            is.sin<-!is.null(attr(trms2,"specials")$sinusoid)
            is.s<-!is.null(attr(trms2,"specials")$s) || !is.null(attr(trms2,"specials")$te) ||!is.null(attr(trms2,"specials")$ti)
            is.M.i<-0
            if (!is.spec && (length(with(data,eval(attr(trms2,"variables"))[[1]]))/n > 1)) is.M.i <-1
            is.M <- c(is.M, is.M.i)
            repsX <- c(repsX,i)
            if (!is.spec){
                Design.Xt<-model.matrix(trms2,data=data)
                remove.cols<-which(c(apply(Design.Xt,2,sd)==0))
                if (length(remove.cols)>0 && !predInd) Design.Xt<-Design.Xt[,-remove.cols,drop=FALSE]
                count<-c(count,1)
                vars[[i]]<-as.character(attr(trms2,"variables")[[2]])
                LF<-which(vars[[i]]=="as.factor")
                if (length(LF)>0) {vars[[i]]<-vars[[i]][-which(vars[[i]]=="as.factor")];labels[i]<-vars[[i]]}
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
                count<-c(count,XK$count)
                vars[[i]]<-XK$vars
                is.D[[i]]<-XK$is.D
                labels[i]<-XK$label
                Rknots[[k<-k+1]]<- XK$knots
                which.Spec[[i]]<-k                
			}
			if (is.sin){
			    XK<-with(data,eval(trms2[1][[2]]))
                Design.Xt<-XK$X                
                count<-c(count,XK$count)
                vars[[i]]<-XK$vars                
                is.D[[i]]<-XK$is.D
                labels[i]<-XK$label                
                amplitude<-XK$amplitude 
                harmonics<-XK$harmonics
                breaks<-XK$breaks
                period<-XK$period
                repsX<-c(repsX,rep(i, harmonics - 1))
                isSin[i]<-1 #is it a sinusoidal term?                
                startSin<-NCOL(Design.X)
                which.Spec[[i]]<--99
                Dynamic<-XK$Dynamic 
                if (amplitude > 1){ 
                    Rknots[[k<-k+1]]<- XK$knots
                    which.Spec[[i]]<-k                    
				}
                if (amplitude == 1) is.M[i] <- 1 #convention to avoid centering sinusoidals with fixed amplitude
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
			    #remove.cols<-which(c(apply(Design.Xt,2,sd)==0))
			    colnames(Design.Xt)<-paste(labels[i],1:NCOL(Design.Xt),sep=".")
			}			
			Design.X<-cbind(Design.X,Design.Xt)
			assign<-c(assign,rep(i,NCOL(Design.Xt)))
			ifelse(is.spec,assign2<-c(assign2,rep(max(assign2)+1,NCOL(Design.Xt))),
			               assign2<-c(assign2,max(assign2)+c(1:NCOL(Design.Xt))))
			               
			ifelse(is.sin && amplitude==1,
			       assign3<-c(assign3,max(assign3)+rep(seq(from=1,to=NCOL(Design.Xt)/2),each=2)),
                   assign3<-c(assign3,rep(max(assign3)+1,NCOL(Design.Xt))))
                          
        }
	}
	remove.cols <- which(colnames(Design.X) == "(Intercept)")
    if (length(remove.cols) > 1) Design.X<-Design.X[,-remove.cols[-1],drop=FALSE]
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
	if (centre) Design.X[,!indicator]<-Design.X[,!indicator]-matrix(1,nrow=n)%*%matrix(meanVector[!indicator],nrow=1)
	if (!centre) meanVector <- 0*meanVector
	rDyn <- NULL
	if (!is.null(Dynamic)) rDyn <- as.matrix(Dynamic)
    return(list(y=y,X=Design.X,assign=assign,Rknots=Rknots,meanVector=meanVector,
                indicator=indicator,labels=labels,count=count,vars=vars,is.D=is.D,
                which.Spec=which.Spec,formula.terms=formula.terms,assign2=assign2,
                assign3=assign3,repsX=repsX,isSin=isSin,Dynamic=rDyn,
                DSP=c(amplitude,startSin,harmonics),breaks=breaks,period=period))
}

match.call.defaults <- function(...) {
    call <- evalq(match.call(expand.dots = FALSE), parent.frame(1))
    formals <- evalq(formals(), parent.frame(1))
    for(k in setdiff(names(formals), names(call)))
        call[k] <- list( formals[[k]] )
    match.call(sys.function(sys.parent()), call)
}

quiet <- function(x) {
    sink(tempfile())
    on.exit(sink())
    invisible(force(x))
}
