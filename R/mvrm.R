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
        #print("x2");print("x2")
        #print("mm");print(mm)
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
