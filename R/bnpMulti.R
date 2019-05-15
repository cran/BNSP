mmvrm <- function(formula,data=list(),sweeps,burn=0,thin=1,seed,StorageDir,
                  c.betaPrior="IG(0.5,0.5*n*p)", pi.muPrior="Beta(1,1)", c.alphaPrior="IG(1.1,1.1)",
                  sigma2Prior="HN(2)", pi.sigmaPrior="Beta(1,1)", mu.RPrior="N(0,1)",
                  sigma.RPrior="HN(1)",corr.Model=c("common",nClust=1),DP.concPrior="Gamma(5,2)",
                  tuneAlpha,tuneSigma2,tuneCb,tuneCa,tuneR,tuneSigma2R,tau,FT=1,...){
    #Samples etc
    if (thin <= 0) thin <- 1
    thin <- as.integer(thin)
    sweeps <- as.integer(sweeps)
    burn <- as.integer(burn)
    if (missing(sweeps)) stop("provide sweeps argument")
    nSamples<-0
    LASTsw<-0
    if (sweeps > 0 && (sweeps-burn) > 0){
		nSamples <- length(seq(1,(sweeps-burn),by=thin))
		LASTsw<-seq(burn,by=thin,length.out=nSamples)[nSamples]
	}
	LASTWB<-1
	#print(LASTWB)
    # Match call
    call <- match.call(expand.dots = FALSE)
    call2 <- match.call.defaults()
    #Data
    if (!is.list(data) && !is.data.frame(data))
        data <- as.data.frame(data)
    #Formula & data dimensions
    p <- length(as.Formula(formula))[1]
    if (length(as.Formula(formula))[2] > 2) stop("more than two regression models provided")
    if (length(as.Formula(formula))[2] == 1) formula <- as.Formula(formula, ~1)
    formula.save <- formula(formula)
    formula.v<-formula(as.Formula(formula),lhs=0,rhs=2)
    formula<-formula(as.Formula(formula),lhs=1,rhs=1)
    # Responses, design matrices, indicators
    Y<-NULL
    varsY<-list()
    for (i in 1:p){
        trms<-terms.formula(formula(as.Formula(formula.save,~1),lhs=i,rhs=3))
        Y<-cbind(Y,with(data,eval(attr(trms,"variables")[[2]])))
        varsY[[i]] <- attr(trms,"variables")[[2]]
    }
    # Null Deviance
    nullDeviance <- 0
    for (i in 1:p) nullDeviance <- nullDeviance - 2 * logLik(lm(Y[,i] ~ 1))
    # Response etc
    Y<-c(t(Y))
    n<-length(Y)/p
    #
    XYK<-DM(formula=formula,data=data,n=n)
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
    #print(X)
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
    if (LD==0) MVLD <- 1
    if (LD > 0) MVLD<-max(vecLD)
    assignZ<-ZK$assign
    labelsZ<-ZK$labels
    countZ<-ZK$count
    varsZ<-ZK$vars
    is.Dz<-ZK$is.D
    which.SpecZ<-ZK$which.Spec
    formula.termsZ<-ZK$formula.terms
    togetherZ<-ZK$together
    #print(Z)
    ##Dim fix
    MVLD<-1
    if (LD > 0) MVLD <-max(vecLD)
    #Prior for pi.mu
    if (!length(pi.muPrior)==1 && !length(pi.muPrior)==NG && !length(pi.muPrior)==(p*NG))
        stop("pi.muPrior has incorrect dimension")
    pimu<-NULL
    for (k in 1:length(pi.muPrior)){
        sp<-strsplit(pi.muPrior[k],"Beta\\(")
        sp<-strsplit(sp[[1]][2],"\\)")
        sp<-strsplit(sp[[1]][1],",")
        pimu<-c(pimu,as.numeric(sp[[1]]))
    }
    if (length(pi.muPrior)==1) pimu<-rep(pimu,each=p*NG)
    if (length(pi.muPrior)==NG) pimu<-rep(pimu,p)
    #Prior for c.beta
    sp<-strsplit(c.betaPrior,"IG\\(")
    sp<-strsplit(sp[[1]][2],"\\)")
    sp<-strsplit(sp[[1]][1],",")
    cetaParams<-c(as.numeric(sp[[1]][1]),eval(parse(text=sp[[1]][2])))
    #Prior for pi.sigma
    if (!length(pi.sigmaPrior)==1 && !length(pi.sigmaPrior)==ND && !length(pi.sigmaPrior)==(p*ND))
        stop("pi.sigmaPrior has incorrect dimension")
    pisigma<-NULL
    for (k in 1:length(pi.sigmaPrior)){
        sp<-strsplit(pi.sigmaPrior[k],"Beta\\(")
        sp<-strsplit(sp[[1]][2],"\\)")
        sp<-strsplit(sp[[1]][1],",")
        pisigma<-c(pisigma,as.numeric(sp[[1]]))
    }
    if (length(pi.sigmaPrior)==1) pisigma<-rep(pisigma,each=p*ND)
    if (length(pi.sigmaPrior)==ND) pisigma<-rep(pisigma,p)
    #Prior for c.alpha
    if (!length(c.alphaPrior)==1 && !length(c.alphaPrior)==p)
        stop("c.alphaPrior has incorrect dimension")
    specials<-c("HN","IG")
    calphaParams<-NULL
    HNca<-vector()
    for (k in 1:length(c.alphaPrior)){
        sp<-strsplit(c.alphaPrior[k],"\\(")
        if (sp[[1]][1] %in% specials){
            if (match(sp[[1]][1],specials)==1) HNca[k]<-1
            if (match(sp[[1]][1],specials)==2) HNca[k]<-0
        } else stop("unrecognised prior for c.alpha_0k")
        sp<-strsplit(sp[[1]][2],"\\)")
        sp<-strsplit(sp[[1]][1],",")
        calphaParams<-c(calphaParams,as.numeric(sp[[1]]))
    }
    if (length(c.alphaPrior)==1){
        calphaParams<-rep(calphaParams,p)
        HNca<-rep(HNca,p)
	}
    #Prior for sigma2_{zero k}
    if (!length(sigma2Prior)==1 && !length(sigma2Prior)==p)
        stop("sigma2Prior has incorrect dimension")
    specials<-c("HN","IG")
    szkParams<-NULL
    HNszk<-vector()
    for (k in 1:length(sigma2Prior)){
        sp<-strsplit(sigma2Prior[k],"\\(")
        if (sp[[1]][1] %in% specials){
            if (match(sp[[1]][1],specials)==1) HNszk[k]<-1
            if (match(sp[[1]][1],specials)==2) HNszk[k]<-0
        } else stop("unrecognised prior for sigma^2_k")
        sp<-strsplit(sp[[1]][2],"\\)")
        sp<-strsplit(sp[[1]][1],",")
        szkParams<-c(szkParams,as.numeric(sp[[1]]))
    }
    if (length(sigma2Prior)==1){
        szkParams<-rep(szkParams,p)
        HNszk<-rep(HNszk,p)
	}
	#Prior for muR
	sp<-strsplit(mu.RPrior,"N\\(")
    sp<-strsplit(sp[[1]][2],"\\)")
    sp<-strsplit(sp[[1]][1],",")
    Rparams<-c(as.numeric(unlist(sp)))
    #Prior for sigmaR
    sp<-strsplit(sigma.RPrior,"HN\\(")
    sp<-strsplit(sp[[1]][2],"\\)")
    sp<-strsplit(sp[[1]][1],",")
    Rparams<-c(Rparams,as.numeric(unlist(sp)))
    #Cor model
    corModels<-c("common","groupC","groupV")
    mcm<-match(corr.Model[1],corModels)
    if (is.na(mcm)) stop("unrecognised correlation model")
    H <- G <- 1
    if (p==2) mcm=1
    if (mcm==2){
        H<-as.numeric(corr.Model[2])
        if (H==1) {mcm<-1; warning("Common correlations model specified with nClust = 1")}
        if (is.na(H) || (!H%%1==0) || H==0) {H <- p*(p-1)/2; warning(cat("mispecified number of clusters. nClust set to ",H,"\n"))}
    }
    if (mcm==3){
        G<-as.numeric(corr.Model[2])
        if (G==1) {mcm<-1; warning("Common correlations model specified with nClust = 1")}
        if (is.na(G) || (!G%%1==0) || G==0) {G <- p; warning(cat("mispecified number of clusters. nClust set to ",G,"\n"))}
        H<-G*(G-1)/2+G #min(d,G*(G-1)/2+G) #min(G,abs(p-G))
	}
    #Prior for alpha DP
    if (mcm > 1){
        sp<-strsplit(DP.concPrior,"Gamma\\(")
        sp<-strsplit(sp[[1]][2],"\\)")
        sp<-strsplit(sp[[1]][1],",")
        DPparams<-as.numeric(sp[[1]])
        DPparams <- c(DPparams,0.01)
    }
    #Seed
    if (missing(seed)) seed<-as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31)
    # Storage directory & files
    WF <- 1
    if (!missing(StorageDir)){
        StorageDir <- path.expand(StorageDir)
        ncharwd <- nchar(StorageDir)}
    if (!missing(StorageDir)) if (!(substr(StorageDir,ncharwd,ncharwd)=="/")) StorageDir <- paste(StorageDir,"/",sep="")
    if (!missing(StorageDir)) if (!dir.exists(StorageDir)) dir.create(StorageDir, recursive = TRUE)
    if (missing(StorageDir)) stop("provide a storage directory via argument StorageDir")
    FL <- c("gamma", "cbeta", "delta", "alpha", "R", "theta", "muR", "sigma2R", "calpha", "sigma2", "beta", "test",
            "compAlloc", "nmembers", "deviance", "DPconc",
            "compAllocV","nmembersV",
            "DE")
    for (i in 1:length(FL)){
        oneFile <- paste(StorageDir, paste("BNSP",FL[i], "txt",sep="."),sep="/")
        if (file.exists(oneFile)) file.remove(oneFile)
	}
    #Tuning Parameters
    if (missing(tuneR)) tuneR<-100*(p+2)
    tuneR[which(tuneR<p+2)]<-p+2
    if (missing(tuneSigma2R)) tuneSigma2R<-1
    if (missing(tuneCa)) tuneCa<-rep(1,p)
    if (!length(tuneCa)==p) tuneCa<-rep(mean(tuneCa),p)
    if (missing(tuneSigma2)) tuneSigma2<-rep(1,p)
    if (!length(tuneSigma2)==p) tuneSigma2<-rep(mean(tuneSigma2),p)
    if (missing(tuneCb)) tuneCb<-10
    if (missing(tuneAlpha)) tuneAlpha<-rep(5,ND*p)
    if (!length(tuneAlpha)==(ND*p)) tuneAlpha<-rep(mean(tuneAlpha),ND*p)
	if (missing(tau)) tau = 0.01
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
    #Deviance
    deviance <- c(0,0)
    #Call C
    if (mcm==1) out<-.C("mult",
            as.integer(seed),as.character(StorageDir),as.integer(WF),
            as.integer(sweeps),as.integer(burn),as.integer(thin),
            as.integer(n),as.integer(p),as.double(Y),as.double(t(as.matrix(X))),as.double(as.matrix(Z)),
            as.integer(LG),as.integer(LD),
            as.double(blockSizeProbG),as.integer(maxBSG),as.double(blockSizeProbD),as.integer(maxBSD),
            as.integer(NG),as.integer(ND),as.integer(vecLG),as.integer(vecLD),
            as.integer(cusumVecLG),as.integer(cusumVecLD),as.integer(MVLD),
            as.double(tuneSigma2R),as.double(tuneCa),as.double(tuneSigma2),as.double(tuneCb),as.double(tuneAlpha),as.double(tuneR),
            as.double(pimu),as.double(cetaParams),as.double(pisigma),
            as.integer(HNca),as.double(calphaParams),as.double(Rparams),
            as.integer(HNszk),as.double(szkParams),as.double(tau),as.integer(FT),as.double(deviance),
            as.integer(c(0)),as.integer(c(0)),as.integer(c(0)),as.double(c(0)),as.double(c(0)),as.double(c(0)),
            as.double(c(0)),as.double(c(0)),as.double(c(0)),as.double(c(0)),
            as.integer(c(LASTsw)),as.double(c(0)),as.integer(LASTWB))
    if (mcm==2) out<-.C("multg",
            as.integer(seed),as.character(StorageDir),as.integer(WF),
            as.integer(sweeps),as.integer(burn),as.integer(thin),
            as.integer(n),as.integer(p),as.double(Y),as.double(t(as.matrix(X))),as.double(as.matrix(Z)),
            as.integer(LG),as.integer(LD),
            as.double(blockSizeProbG),as.integer(maxBSG),as.double(blockSizeProbD),as.integer(maxBSD),
            as.integer(NG),as.integer(ND),as.integer(vecLG),as.integer(vecLD),
            as.integer(cusumVecLG),as.integer(cusumVecLD),as.integer(MVLD),
            as.double(tuneSigma2R),as.double(tuneCa),as.double(tuneSigma2),as.double(tuneCb),as.double(tuneAlpha),as.double(tuneR),
            as.double(pimu),as.double(cetaParams),as.double(pisigma),
            as.integer(HNca),as.double(calphaParams),as.double(Rparams),
            as.integer(HNszk),as.double(szkParams),as.double(tau),as.integer(FT),as.double(deviance),as.integer(H),
            as.double(DPparams), as.integer(c(0)),as.integer(c(0)),as.integer(c(0)),as.double(c(0)),as.double(c(0)),
            as.double(c(0)),as.double(c(0)),as.double(c(0)),as.double(c(0)),
            as.double(c(0)),as.integer(c(0)),as.double(c(0)),
            as.integer(c(LASTsw)),as.double(c(0)),as.integer(LASTWB))
    if (mcm==3) out<-.C("multgv",
            as.integer(seed),as.character(StorageDir),as.integer(WF),
            as.integer(sweeps),as.integer(burn),as.integer(thin),
            as.integer(n),as.integer(p),as.double(Y),as.double(t(as.matrix(X))),as.double(as.matrix(Z)),
            as.integer(LG),as.integer(LD),
            as.double(blockSizeProbG),as.integer(maxBSG),as.double(blockSizeProbD),as.integer(maxBSD),
            as.integer(NG),as.integer(ND),as.integer(vecLG),as.integer(vecLD),
            as.integer(cusumVecLG),as.integer(cusumVecLD),as.integer(MVLD),
            as.double(tuneSigma2R),as.double(tuneCa),as.double(tuneSigma2),as.double(tuneCb),as.double(tuneAlpha),as.double(tuneR),
            as.double(pimu),as.double(cetaParams),as.double(pisigma),
            as.integer(HNca),as.double(calphaParams),as.double(Rparams),
            as.integer(HNszk),as.double(szkParams),as.double(tau),as.integer(FT),as.double(deviance),as.integer(G),
            as.double(DPparams),as.integer(c(0)),as.integer(c(0)),as.integer(c(0)),as.double(c(0)),as.double(c(0)),
            as.double(c(0)),as.double(c(0)),as.double(c(0)),as.double(c(0)),
            as.double(c(0)),as.integer(c(0)),as.double(c(0)),
            as.integer(c(LASTsw)),as.double(c(0)),as.integer(LASTWB))
    #Output
    loc<-25
    fit <- list(call=call,call2=call2,formula=formula.save,seed=seed,p=p,d=p*(p-1)/2,
                data=data,Y=Y,X=X,Xknots=Xknots,Z=Z,Zknots=Zknots,LG=LG,LD=LD,ND=ND,NG=NG,
                mcpar=c(as.integer(burn+1),as.integer(seq(from=burn+1,by=thin,length.out=nSamples)[nSamples]),as.integer(thin)),
                nSamples=nSamples,storeMeanVectorX=storeMeanVectorX,storeMeanVectorZ=storeMeanVectorZ,
                tuneSigma2R=c(tuneSigma2R,out[[loc+0]][1]),tuneCa=c(tuneCa,out[[loc+1]][1:p]),tuneSigma2=c(tuneSigma2,out[[loc+2]][1:p]),
                tuneCb=c(tuneCb,out[[loc+3]][1]),tuneAlpha=c(tuneAlpha,out[[loc+4]][1:(p*ND)]),tuneR=c(tuneR,out[[loc+5]][1]),
                DIR=StorageDir,deviance=c(out[[41]][1:2]),nullDeviance=nullDeviance,
                storeIndicatorX=storeIndicatorX,storeIndicatorZ=storeIndicatorZ,assignX=assignX,
                assignZ=assignZ,labelsX=labelsX,labelsZ=labelsZ,countX=countX,countZ=countZ,
                varsX=varsX,varsZ=varsZ,varsY=varsY,is.Dx=is.Dx,is.Dz=is.Dz,which.SpecX=which.SpecX,
                which.SpecZ=which.SpecZ,formula.termsX=formula.termsX,formula.termsZ=formula.termsZ,
                togetherX=togetherX,togetherZ=togetherZ,HNca=HNca,H=H,G=G,mcm=mcm,out=out,totalSweeps=sweeps)
    class(fit) <- 'mmvrm'
    return(fit)
}

continue <- function(object,sweeps,discard=FALSE,...){
	#Sweeps
	if (missing(sweeps)) stop("provide sweeps argument")
    sweeps <- as.integer(sweeps)
    nSamples <-0
    if (sweeps > 0){
        nSamples <- length(seq(1,sweeps,by=object$mcpar[3]))
        LASTsw<-seq(0,by=object$mcpar[3],length.out=nSamples)[nSamples]
	}
    if (nSamples<=0) stop("problem with sweeps argument")
    totalSweeps<-object$totalSweeps
    LASTWB <- floor(totalSweeps/50)+1
    #print(LASTWB)
    totalSweeps<-object$totalSweeps + sweeps
    #Files
    FL <- c("gamma", "cbeta", "delta", "alpha", "R", "theta", "muR", "sigma2R", "calpha", "sigma2", "beta", "test",
            "compAlloc", "nmembers", "deviance", "DPconc",
            "compAllocV","nmembersV",
            "DE")
    gamma <- paste(object$DIR, paste("BNSP",FL[1], "txt",sep="."),sep="/")
    gamma <- scan(gamma,what=numeric(),n=object$p*object$LG,quiet=TRUE,skip=object$nSamples-1)
    cbeta <- paste(object$DIR, paste("BNSP",FL[2], "txt",sep="."),sep="/")
    cbeta <- scan(cbeta,what=numeric(),n=1,quiet=TRUE,skip=object$nSamples-1)
    delta <- paste(object$DIR, paste("BNSP",FL[3], "txt",sep="."),sep="/")
    delta <- scan(delta,what=numeric(),n=object$p*object$LD,quiet=TRUE,skip=object$nSamples-1)
    alpha <- paste(object$DIR, paste("BNSP",FL[4], "txt",sep="."),sep="/")
    alpha <- scan(alpha,what=numeric(),n=object$p*object$LD,quiet=TRUE,skip=object$nSamples-1)
    R <- paste(object$DIR, paste("BNSP",FL[5], "txt",sep="."),sep="/")
    R <- scan(R,what=numeric(),n=object$p^2,quiet=TRUE,skip=object$nSamples-1)
    muR <- paste(object$DIR, paste("BNSP",FL[7], "txt",sep="."),sep="/")
    muR <- scan(muR,what=numeric(),n=1,quiet=TRUE,skip=object$nSamples-1)
    sigma2R <- paste(object$DIR, paste("BNSP",FL[8], "txt",sep="."),sep="/")
    sigma2R <- scan(sigma2R,what=numeric(),n=1,quiet=TRUE,skip=object$nSamples-1)
    calpha <- paste(object$DIR, paste("BNSP",FL[9], "txt",sep="."),sep="/")
    calpha <- scan(calpha,what=numeric(),n=object$p,quiet=TRUE,skip=object$nSamples-1)
    sigma2 <- paste(object$DIR, paste("BNSP",FL[10], "txt",sep="."),sep="/")
    sigma2 <- scan(sigma2,what=numeric(),n=object$p,quiet=TRUE,skip=object$nSamples-1)
    DE <- paste(object$DIR, paste("BNSP",FL[19], "txt",sep="."),sep="/")
    DE <- scan(DE,what=numeric(),n=2*object$p^2,quiet=TRUE,skip=0)
    compAlloc <- paste(object$DIR, paste("BNSP",FL[13], "txt",sep="."),sep="/")
    if (file.exists(compAlloc)) compAlloc <- scan(compAlloc,what=numeric(),n=object$d,quiet=TRUE,skip=object$nSamples-1)
    DPconc <- paste(object$DIR, paste("BNSP",FL[16], "txt",sep="."),sep="/")
    if (file.exists(DPconc)) DPconc <- scan(DPconc,what=numeric(),n=1,quiet=TRUE,skip=object$nSamples-1)
    compAllocV <- paste(object$DIR, paste("BNSP",FL[17], "txt",sep="."),sep="/")
    if (file.exists(compAllocV)) compAllocV <- scan(compAllocV,what=numeric(),n=object$p,quiet=TRUE,skip=object$nSamples-1)
    if (discard==TRUE){
        for (i in 1:length(FL)){
            oneFile <- paste(object$DIR, paste("BNSP",FL[i], "txt",sep="."),sep="/")
            if (file.exists(oneFile)) file.remove(oneFile)
	    }
	}
	oneFile <- paste(object$DIR, paste("BNSP",FL[length(FL)], "txt",sep="."),sep="/")
    if (file.exists(oneFile)) file.remove(oneFile)
    #Deviance & nSamples
    deviance<-c(0,0)
    if (discard==FALSE) deviance<-as.double(object$out[[41]])
    if (discard==TRUE) {object$nSamples<-0; object$mcpar[1]<-1}
    #Tuning
    loc<-25
    tuneSigma2R<-object$out[[loc+0]][1]
    tuneCa<-object$out[[loc+1]][1:object$p]
    tuneSigma2<-object$out[[loc+2]][1:object$p]
    tuneCb<-object$out[[loc+3]][1]
    tuneAlpha<-object$out[[loc+4]][1:(object$p*object$ND)]
    tuneR<-object$out[[loc+5]][1]
	#Call C
    if (object$mcm==1) out<-.C("mult",
            as.integer(object$out[[1]]),as.character(object$out[[2]]),as.integer(object$out[[3]]),
            as.integer(sweeps),as.integer(0),as.integer(object$out[[6]]),
            as.integer(object$out[[7]]),as.integer(object$out[[8]]),as.double(object$out[[9]]),as.double(object$out[[10]]),as.double(object$out[[11]]),
            as.integer(object$out[[12]]),as.integer(object$out[[13]]),
            as.double(object$out[[14]]),as.integer(object$out[[15]]),as.double(object$out[[16]]),as.integer(object$out[[17]]),
            as.integer(object$out[[18]]),as.integer(object$out[[19]]),as.integer(object$out[[20]]),as.integer(object$out[[21]]),
            as.integer(object$out[[22]]),as.integer(object$out[[23]]),as.integer(object$out[[24]]),
            as.double(object$out[[25]]),as.double(object$out[[26]]),as.double(object$out[[27]]),as.double(object$out[[28]]),as.double(object$out[[29]]),
            as.double(object$out[[30]]), as.double(object$out[[31]]),as.double(object$out[[32]]),as.double(object$out[[33]]),
            as.integer(object$out[[34]]),as.double(object$out[[35]]),as.double(object$out[[36]]),
            as.integer(object$out[[37]]),as.double(object$out[[38]]),as.double(object$out[[39]]),as.integer(object$out[[40]]),as.double(deviance),
            as.integer(c(1)),as.integer(gamma),as.integer(delta),as.double(alpha),as.double(sigma2),as.double(R),
            as.double(muR),as.double(sigma2R),as.double(cbeta),as.double(calpha),
            as.integer(c(LASTsw)),as.double(c(DE)),as.integer(LASTWB))
    if (object$mcm==2) out<-.C("multg",
            as.integer(object$out[[1]]),as.character(object$out[[2]]),as.integer(object$out[[3]]),
            as.integer(sweeps),as.integer(0),as.integer(object$out[[6]]),
            as.integer(object$out[[7]]),as.integer(object$out[[8]]),as.double(object$out[[9]]),as.double(object$out[[10]]),as.double(object$out[[11]]),
            as.integer(object$out[[12]]),as.integer(object$out[[13]]),
            as.double(object$out[[14]]),as.integer(object$out[[15]]),as.double(object$out[[16]]),as.integer(object$out[[17]]),
            as.integer(object$out[[18]]),as.integer(object$out[[19]]),as.integer(object$out[[20]]),as.integer(object$out[[21]]),
            as.integer(object$out[[22]]),as.integer(object$out[[23]]),as.integer(object$out[[24]]),
            as.double(object$out[[25]]),as.double(object$out[[26]]),as.double(object$out[[27]]),as.double(object$out[[28]]),as.double(object$out[[29]]),as.double(object$out[[30]]),
            as.double(object$out[[31]]),as.double(object$out[[32]]),as.double(object$out[[33]]),
            as.integer(object$out[[34]]),as.double(object$out[[35]]),as.double(object$out[[36]]),
            as.integer(object$out[[37]]),as.double(object$out[[38]]),as.double(object$out[[39]]),as.integer(object$out[[40]]),as.double(deviance),as.integer(object$out[[42]]),as.double(object$out[[43]]),
            as.integer(c(1)),as.integer(gamma),as.integer(delta),as.double(alpha),as.double(sigma2),as.double(R),
            as.double(muR),as.double(sigma2R),as.double(cbeta),as.double(calpha),
            as.integer(compAlloc),as.double(DPconc),as.integer(c(LASTsw)),as.double(c(DE)),as.integer(LASTWB))
    if (object$mcm==3) out<-.C("multgv",
            as.integer(object$out[[1]]),as.character(object$out[[2]]),as.integer(object$out[[3]]),
            as.integer(sweeps),as.integer(0),as.integer(object$out[[6]]),
            as.integer(object$out[[7]]),as.integer(object$out[[8]]),as.double(object$out[[9]]),as.double(object$out[[10]]),as.double(as.matrix(object$out[[11]])),
            as.integer(object$out[[12]]),as.integer(object$out[[13]]),
            as.double(object$out[[14]]),as.integer(object$out[[15]]),as.double(object$out[[16]]),as.integer(object$out[[17]]),
            as.integer(object$out[[18]]),as.integer(object$out[[19]]),as.integer(object$out[[20]]),as.integer(object$out[[21]]),
            as.integer(object$out[[22]]),as.integer(object$out[[23]]),as.integer(object$out[[24]]),
            as.double(object$out[[25]]),as.double(object$out[[26]]),as.double(object$out[[27]]),as.double(object$out[[28]]),as.double(object$out[[29]]),
            as.double(object$out[[30]]),as.double(object$out[[31]]),as.double(object$out[[32]]),as.double(object$out[[33]]),
            as.integer(object$out[[34]]),as.double(object$out[[35]]),as.double(object$out[[36]]),
            as.integer(object$out[[37]]),as.double(object$out[[38]]),as.double(object$out[[39]]),as.integer(object$out[[40]]),as.double(deviance),as.integer(object$out[[42]]),as.double(object$out[[43]]),
            as.integer(c(1)),as.integer(gamma),as.integer(delta),as.double(alpha),as.double(sigma2),as.double(R),
            as.double(muR),as.double(sigma2R),as.double(cbeta),as.double(calpha),
            as.integer(compAlloc),as.double(DPconc),as.integer(c(LASTsw)),as.double(c(DE)),as.integer(LASTWB))
    #Output
    fit <- list(call=object$call,call2=object$call2,formula=object$formula,seed=object$seed,p=object$p,d=object$d,
                data=object$data,Y=object$Y,X=object$X,Xknots=object$Xknots,Z=object$Z,Zknots=object$Zknots,LG=object$LG,LD=object$LD,ND=object$ND,NG=object$NG,
                mcpar=c(as.integer(object$mcpar[1]),as.integer(seq(from=object$mcpar[1],by=object$mcpar[3],length.out=(nSamples+object$nSamples))[nSamples+object$nSamples]),as.integer(object$mcpar[3])),
                nSamples=nSamples+object$nSamples,
                storeMeanVectorX=object$storeMeanVectorX,storeMeanVectorZ=object$storeMeanVectorZ,
                tuneSigma2R=c(tuneSigma2R,out[[loc+0]][1]),tuneCa=c(tuneCa,out[[loc+1]][1:object$p]),tuneSigma2=c(tuneSigma2,out[[loc+2]][1:object$p]),
                tuneCb=c(tuneCb,out[[loc+3]][1]),tuneAlpha=c(tuneAlpha,out[[loc+4]][1:(object$p*object$ND)]),tuneR=c(tuneR,out[[loc+5]][1]),
                DIR=object$DIR,deviance=c(out[[41]][1:2]),nullDeviance=object$nullDeviance,
                storeIndicatorX=object$storeIndicatorX,storeIndicatorZ=object$storeIndicatorZ,assignX=object$assignX,
                assignZ=object$assignZ,labelsX=object$labelsX,labelsZ=object$labelsZ,countX=object$countX,countZ=object$countZ,
                varsX=object$varsX,varsZ=object$varsZ,varsY=object$varsY,is.Dx=object$is.Dx,is.Dz=object$is.Dz,which.SpecX=object$which.SpecX,
                which.SpecZ=object$which.SpecZ,formula.termsX=object$formula.termsX,formula.termsZ=object$formula.termsZ,
                togetherX=object$togetherX,togetherZ=object$togetherZ,HNca=object$HNca,H=object$H,G=object$G,mcm=object$mcm,out=out,
                totalSweeps=totalSweeps)
    class(fit) <- 'mmvrm'
    return(fit)
}

mmvrm2mcmc <- function(mmvrmObj,labels){
    mvrmObj <- mmvrmObj
    all.labels <- c("alpha","calpha","cbeta","delta","beta","gamma","sigma2")
    more.labels1 <- c("muR","sigma2R","R","theta")
    more.labels2 <- c("compAlloc","nmembers","DPconc")
    more.labels3 <- c("compAllocV","nmembersV","deviance")
    all.labels<-c(all.labels,more.labels1,more.labels2,more.labels3)
    if (missing(labels)) labels <- all.labels
    mtch<-match(labels,all.labels)
	p<-mvrmObj$p
	R<-NULL
    if (any(mtch==1) && mvrmObj$LD > 0){
        file <- paste(mmvrmObj$DIR,"BNSP.alpha.txt",sep="")
        if (file.exists(file)) R<-cbind(R,matrix(unlist(read.table(file)),ncol=mvrmObj$LD*p,
           dimnames=list(c(),paste("alpha",rep(colnames(mvrmObj$Z)[-1],p),rep(mvrmObj$varsY,each=mvrmObj$LD),sep="."))))
    }
    if (any(mtch==2)){
        file <- paste(mmvrmObj$DIR,"BNSP.calpha.txt",sep="")
        if (file.exists(file)) R<-cbind(R,matrix(unlist(read.table(file)),ncol=p,
                                  dimnames=list(c(),paste(rep("c_alpha",p),mvrmObj$varsY,sep="."))))
    }
    if (any(mtch==3)){
       file <- paste(mmvrmObj$DIR,"BNSP.cbeta.txt",sep="")
       if (file.exists(file)) R<-cbind(R,matrix(unlist(read.table(file)),ncol=1,dimnames=list(c(),c("c_beta"))))
    }
    if (any(mtch==4) && mvrmObj$LD > 0){
        file <- paste(mmvrmObj$DIR,"BNSP.delta.txt",sep="")
        if (file.exists(file)) R<-cbind(R,matrix(unlist(read.table(file)),ncol=mvrmObj$LD*p,
           dimnames=list(c(),paste(rep(paste("delta",seq(1,mvrmObj$LD),sep="_"),p),rep(mvrmObj$varsY,each=mvrmObj$LD),sep="."))))
    }
    if (any(mtch==5)){
        file <- paste(mmvrmObj$DIR,"BNSP.beta.txt",sep="")
        if (file.exists(file)) R<-cbind(R,matrix(unlist(read.table(file)),ncol=p*(mvrmObj$LG+1),
           dimnames=list(c(),paste("beta",rep(colnames(mvrmObj$X),p),rep(mvrmObj$varsY,each=mvrmObj$LG+1),sep="."))))
    }
    if (any(mtch==6) && mvrmObj$LG > 0){
        file <- paste(mmvrmObj$DIR,"BNSP.gamma.txt",sep="")
        if (file.exists(file)) R<-cbind(R,matrix(unlist(read.table(file)),ncol=p*mvrmObj$LG,
           dimnames=list(c(),paste(rep(paste("gamma",seq(1,mvrmObj$LG),sep="_"),p),rep(mvrmObj$varsY,each=mvrmObj$LG),sep="."))))
    }
    if (any(mtch==7)){
        file <- paste(mmvrmObj$DIR,"BNSP.sigma2.txt",sep="")
        if (file.exists(file)) R<-cbind(R,matrix(unlist(read.table(file)),ncol=p,
                          dimnames=list(c(),paste(rep("sigma^2",p),mvrmObj$varsY,sep="."))))
	}
	if (any(mtch==8)){
		names <- paste("muR of cluster",seq(1,mvrmObj$H),sep=" ")
		if (mvrmObj$H==1) names <- "muR"
	    file <- paste(mmvrmObj$DIR,"BNSP.muR.txt",sep="")
        if (file.exists(file)) R<-cbind(R,matrix(unlist(read.table(file)),nrow=mvrmObj$nSamples,
                dimnames=list(c(),names)))
    }
    if (any(mtch==9)){
        file <- paste(mmvrmObj$DIR,"BNSP.sigma2R.txt",sep="")
        if (file.exists(file)) R<-cbind(R,matrix(unlist(read.table(file)),ncol=1,dimnames=list(c(),c("sigma2R"))))
	}
	subset <- rep(seq(1,p),each=p) < rep(seq(1,p),p)
	cor.names <- paste(rep(seq(1,p),each=p),rep(seq(1,p),p),sep="")
    if (any(mtch==10)){
        file <- paste(mmvrmObj$DIR,"BNSP.R.txt",sep="")
        if (file.exists(file)) R<-cbind(R,matrix(unlist(read.table(file)),ncol=p*p,
                 dimnames=list(c(),paste("cor",cor.names,sep=".")))[,subset,drop=FALSE])
    }
    if (any(mtch==11)){
        file <- paste(mmvrmObj$DIR,"BNSP.theta.txt",sep="")
        if (file.exists(file)) R<-cbind(R,matrix(unlist(read.table(file)),ncol=mvrmObj$d,
                          dimnames=list(c(),paste("theta",cor.names[subset],sep="."))))
	}
	if (any(mtch==12)){
	    file <- paste(mmvrmObj$DIR,"BNSP.compAlloc.txt",sep="")
	    if (file.exists(file)) R<-cbind(R,matrix(unlist(read.table(file)),ncol=mmvrmObj$d,
	                           dimnames=list(c(),paste("clustering of cor",cor.names[subset],sep="."))))
	}
    if (any(mtch==13)){
        file <- paste(mmvrmObj$DIR,"BNSP.nmembers.txt",sep="")
        if (file.exists(file))R<-cbind(R,matrix(unlist(read.table(file)),ncol=mmvrmObj$H,dimnames=list(c(),paste("elements in cor cluster ",seq(1,mmvrmObj$H)))))
    }
    #if (any(mtch==14)) R<-cbind(R,matrix(unlist(read.table(paste(mmvrmObj$DIR,"BNSP.probs.txt",sep=""))),ncol=mmvrmObj$H,dimnames=list(c(),seq(1,mmvrmObj$H))))
    if (any(mtch==14)){
        file <- paste(mmvrmObj$DIR,"BNSP.DPconc.txt",sep="")
        if (file.exists(file)) R<-cbind(R,matrix(unlist(read.table(file)),ncol=1,dimnames=list(c(),c("DPconc"))))
    }
    if (any(mtch==15)){
        file <- paste(mmvrmObj$DIR,"BNSP.compAllocV.txt",sep="")
        if (file.exists(file)) R<-cbind(R,matrix(unlist(read.table(file)),ncol=mmvrmObj$p,dimnames=list(c(),paste("clustering of ",mmvrmObj$varsY))))
	}
    if (any(mtch==16)){
        file <- paste(mmvrmObj$DIR,"BNSP.nmembersV.txt",sep="")
        if (file.exists(file)) R<-cbind(R,matrix(unlist(read.table(file)),ncol=mmvrmObj$G,dimnames=list(c(),paste("elements in var cluster ",seq(1,mmvrmObj$G)))))
	}
	if (any(mtch==17)){
       file <- paste(mmvrmObj$DIR,"BNSP.deviance.txt",sep="")
       if (file.exists(file)) R<-cbind(R,matrix(unlist(read.table(file)),ncol=2,dimnames=list(c(),c("marginal deviance","full deviance"))))
    }
	if (!is.null(R)){
	    attr(R, "mcpar") <- mvrmObj$mcpar
        attr(R, "class") <- "mcmc"
	}
	return(R)
}

plot2 <- function(x, model="mean", term, response=1, intercept=TRUE, grid=30,
                  centre=mean, quantiles=c(0.1, 0.9), static=TRUE,
                  centreEffects=FALSE, plotOptions=list(), ...){
    mvrmObj<-x
    p<-mvrmObj$p
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
                if (length(mvrmObj$Xknots)>0 && whichKnots>0) {Dstar<-data.frame(knots=mvrmObj$Xknots[[whichKnots]])}else{Dstar<-NULL}
                DsM<-DM(formula=reformulate(formula.term),data=newData,n=NROW(newData),knots=Dstar,
                        meanVector=mvrmObj$storeMeanVectorX[ML],indicator=mvrmObj$storeIndicatorX[ML])$X
                if (! 1%in%ML) DsM<-DsM[,-1]
		    }
	        if (is.D){
	            newData<-data.frame(unique(with(mvrmObj$data,eval(as.name(vars[1])))))
                colnames(newData)<-vars
                DsM<-DM(formula=reformulate(formula.term),data=newData,n=NROW(newData),knots=NULL,
                        meanVector=mvrmObj$storeMeanVectorX[ML],indicator=mvrmObj$storeIndicatorX[ML])$X
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
                DsM<-DM(formula=reformulate(formula.term),data=newData,n=NROW(newData),knots=Dstar,
                        meanVector=mvrmObj$storeMeanVectorX[ML],indicator=mvrmObj$storeIndicatorX[ML])$X
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
                DsM<-DM(formula=reformulate(formula.term),data=newData,n=NROW(newData),knots=Dstar,
                        meanVector=mvrmObj$storeMeanVectorX[ML],indicator=mvrmObj$storeIndicatorX[ML])$X
		    }
        }
        fitM<-matrix(0,nrow=mvrmObj$nSamples,ncol=NROW(DsM))
		etaFN <- file.path(paste(mvrmObj$DIR,"BNSP.beta.txt",sep=""))
        eFile<-file(etaFN,open="r")
		for (i in 1:mvrmObj$nSamples){
            eta<-scan(eFile,what=numeric(),n=p*(mvrmObj$LG+1),quiet=TRUE)
            fitM[i,]<-as.matrix(DsM)%*%matrix(c(eta[ML+(response-1)*(mvrmObj$LG+1)]))
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
        VL<-which(mvrmObj$assignZ %in% int.label)-1
        if (count==1){
	        if (!is.D){
	            min1<-min(with(mvrmObj$data,eval(as.name(vars[1]))))
			    max1<-max(with(mvrmObj$data,eval(as.name(vars[1]))))
			    newR1<-seq(min1,max1,length.out=grid)
                newData<-data.frame(newR1)
                colnames(newData)<-vars
                whichKnots <- mvrmObj$which.SpecZ[[int.label]]
                if (length(mvrmObj$Zknots)>0 && whichKnots>0) {Dstar<-data.frame(mvrmObj$Zknots[[whichKnots]])}else{Dstar<-NULL}
                DsV<-DM(formula=reformulate(formula.term),data=newData,n=NROW(newData),knots=Dstar,
                        meanVector=mvrmObj$storeMeanVectorZ[c(1,VL+1)],indicator=mvrmObj$storeIndicatorZ[c(1,VL+1)])$X[,-1]
		    }
	        if (is.D){
	            newData<-data.frame(unique(with(mvrmObj$data,eval(as.name(vars[1])))))
                colnames(newData)<-vars
                DsV<-DM(formula=reformulate(formula.term),data=newData,n=NROW(newData),knots=NULL,
                        meanVector=mvrmObj$storeMeanVectorZ[c(1,VL+1)],indicator=mvrmObj$storeIndicatorZ[c(1,VL+1)])$X[,-1]
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
                DsV<-DM(formula=reformulate(formula.term),data=newData,n=NROW(newData),knots=Dstar,
                        meanVector=mvrmObj$storeMeanVectorZ[c(1,VL+1)],indicator=mvrmObj$storeIndicatorZ[c(1,VL+1)])$X
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
                DsV<-DM(formula=reformulate(formula.term),data=newData,n=NROW(newData),knots=Dstar,
                        meanVector=mvrmObj$storeMeanVectorZ[c(1,VL+1)],indicator=mvrmObj$storeIndicatorZ[c(1,VL+1)])$X
                DsV<-DsV[,-1]
		    }
        }
        fitV<-matrix(0,nrow=mvrmObj$nSamples,ncol=NROW(DsV))
		alphaFN <- file.path(paste(mvrmObj$DIR,"BNSP.alpha.txt",sep=""))
        aFile<-file(alphaFN,open="r")
		sigma2FN <- file.path(paste(mvrmObj$DIR,"BNSP.sigma2.txt",sep=""))
        s2File<-file(sigma2FN,open="r")
        for (i in 1:mvrmObj$nSamples){
        	alpha<-scan(aFile,what=numeric(),n=p*mvrmObj$LD,quiet=TRUE)
		    if (intercept) s2<-scan(s2File,what=numeric(),n=p,quiet=TRUE)[response]
		    else s2=1
            fitV[i,]<-sqrt(s2*exp(as.matrix(DsV)%*%matrix(c(alpha[VL+(response-1)*mvrmObj$LD]))))
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

quiet <- function(x) {
    sink(tempfile())
    on.exit(sink())
    invisible(force(x))
}

plot.mmvrm<-function(x, model, term, response, intercept=TRUE, grid=30,
                    centre=mean, quantiles=c(0.1, 0.9), static=TRUE,
                    centreEffects=FALSE,plotOptions=list(),nrow,ask=FALSE,...)
{
    oldpar <- NULL
    on.exit(par(oldpar))

    if (missing(response)) response <- c(1:x$p)
    if (max(response) > x$p) stop("argument response exceeds the number of responses");
    if (missing(model)) model<-"both"
    MEAN <- 0; STDEV <- 0
    if ((model=="both" || model=="mean") && (x$NG > 0)) MEAN <- 1
    if ((model=="both" || model=="stdev") && (x$ND > 0)) STDEV <- 1
    if (MEAN==0 && STDEV==0) stop("no terms to plot; only intercepts in the model");
    if (missing(term)) {termM<-1:x$NG; termSD<-1:x$ND}
    if (!missing(term)) {termM <- termSD <- term[[1]]; if (length(term)==2) termSD<-term[[2]]}
    if (missing(nrow)) nrow <- length(response)
    my_plots <- list()
    count <- 1
    for (r in response) {
		if (MEAN)
        for (i in termM){
            plotOptions <- list(ggtitle(paste("mean of", x$varsY[[r]])),plotOptions)

            my_plots[[count]] <- plot2(x, model="mean", term=i, response=r, intercept=intercept, grid=grid,
                  centre=centre, quantiles=quantiles, static=static,
                  centreEffects=centreEffects, plotOptions=plotOptions)

            if (ask==TRUE) print(my_plots[[count]])

            if (count==1)
                oldpar <- c(oldpar, par(ask=ask))

            count <- count + 1
        }
        if (STDEV)
        for (i in termSD){
            plotOptions <- list(ggtitle(paste("stdev of", x$varsY[[r]])),plotOptions)

            my_plots[[count]] <- plot2(x, model = "stdev", term=i, response=r, intercept=intercept, grid=grid,
                  centre=centre, quantiles=quantiles, static=static,
                  centreEffects=centreEffects, plotOptions=plotOptions)

            if (ask==TRUE) print(my_plots[[count]])

            if (count==1)
                oldpar <- c(oldpar, par(ask=ask))

            count <- count + 1
        }
	}
	if (ask==FALSE) quiet(print(grid.arrange(grobs = my_plots,nrow=nrow)))
}

plotCorr <- function(x, term="R", centre=mean, quantiles=c(0.1, 0.9), ...){
    mmvrmObj<-x
    if (mmvrmObj$p == 1) stop("doesn't apply for univariate response")
    if (length(quantiles)==1) quantiles <- c(quantiles,1-quantiles)
    if (length(quantiles) > 2) stop("up to two quantiles")
    if (!is.null(quantiles)) quantiles <- unique(sort(quantiles))
    if (!is.null(quantiles)) par(mfrow=c(1,2))
    mt<-match(term,c("R","muR"))
    if (is.na(mt)) stop("unrecognised term");
    if (mt==1) R<-mmvrm2mcmc(mmvrmObj,"R")
    if (mt==2){
        muR<-mmvrm2mcmc(mmvrmObj,"muR")
        if (mmvrmObj$mcm == 1) return(plot(tanh(muR)))
        compAlloc <- matrix(0,nrow=mmvrmObj$nSamples,ncol=mmvrmObj$d)
        compAllocFile <- paste(mmvrmObj$DIR, paste("BNSP","compAlloc","txt",sep="."),sep="/")
        if (file.exists(compAllocFile)) compAlloc<-mmvrm2mcmc(mmvrmObj,"compAlloc")
        R <- sapply(1:mmvrmObj$d,
                    function(i) muR[cbind(1:mmvrmObj$nSamples,compAlloc[,i]+1)])
        if (x$out[[40]]==1) R <- tanh(R)
	}
    ec <- diag(rep(1,x$p))
    ec[lower.tri(ec)] <- apply(R,2,centre)
    ec[upper.tri(ec)]<-t(ec)[upper.tri(ec)]
    colnames(ec) <- x$varsY
    corrplot.mixed(ec,lower.col="black")
    if (!is.null(quantiles)){
        q<-apply(R,2,quantile,probs=quantiles)
        ci<-diag(1,mmvrmObj$p)
        ci[lower.tri(ci)] <- q[2,]
        ci[upper.tri(ci)]<-t(ci)[upper.tri(ci)]
        ci[lower.tri(ci)] <- q[1,]
        colnames(ci) <- x$varsY
        rownames(ci) <- x$varsY
        corrplot(ci,col="black",method="number")
    }
}

histCorr <- function(x, term="R", plotOptions=list(), ...){
    mmvrmObj<-x
    if (mmvrmObj$p == 1) stop("doesn't apply for univariate response")
    mt<-match(term,c("R","muR"))
    if (is.na(mt)) stop("unrecognised term");
    if (mt==1) R<-mmvrm2mcmc(mmvrmObj,"R")
    if (mt==2){
        muR<-mmvrm2mcmc(mmvrmObj,"muR")
        if (mmvrmObj$mcm == 1) return(plot(tanh(muR)))
        compAlloc <- matrix(0,nrow=mmvrmObj$nSamples,ncol=mmvrmObj$d)
        compAllocFile <- paste(mmvrmObj$DIR, paste("BNSP","compAlloc","txt",sep="."),sep="/")
        if (file.exists(compAllocFile)) compAlloc<-mmvrm2mcmc(mmvrmObj,"compAlloc")
        R <- sapply(1:mmvrmObj$d,
                    function(i) muR[cbind(1:mmvrmObj$nSamples,compAlloc[,i]+1)])
        if (x$out[[40]]==1) R <- tanh(R)
	}
    r<-rep(rep(seq(1,mmvrmObj$p-1),times=seq(mmvrmObj$p-1,1)),each=mmvrmObj$nSample)
    c<-rep(unlist(mapply(seq,seq(2,mmvrmObj$p),mmvrmObj$p)),each=mmvrmObj$nSample)
    df<-data.frame(cor=c(R),r=r,c=c)
    pp<-ggplot(df) + geom_histogram(aes(x=cor),binwidth=0.01) + facet_wrap(r~c) + plotOptions #facet_grid(r~c)
    return(pp)
}

predict.mmvrm <- function(object,newdata,interval=c("none","credible","prediction"),level=0.95,nSamples=100, ...){
    if (missing(newdata) || is.null(newdata)) newdata<-object$data
    if (length(newdata)==0) stop("no data found")
    #if (length(response) > object$p) stop(paste("`response' is a vector of up to", object$p,"integers"))
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
    npred<-NROW(newdata)
    X<-DM(formula=formula2,data=nd,n=npred,knots=object$Xknots,predInd=TRUE,meanVector=object$storeMeanVectorX,indicator=object$storeIndicatorX,mvrmObj=object)$X
	fitM<-matrix(0,nrow=object$nSamples,ncol=npred*object$p)
    etaFN <- file.path(paste(object$DIR,"BNSP.beta.txt",sep=""))
    eFile<-file(etaFN,open="r")
	for (i in 1:object$nSamples){
        eta<-scan(eFile,what=numeric(),n=object$p*(object$LG+1),quiet=TRUE)
        for (k in 1:object$p)
            fitM[i,(1+(k-1)*npred):(k*npred)]<-c(as.matrix(X)%*%matrix(c(eta[(1+(k-1)*(object$LG+1)):(k*(object$LG+1))])))
	}
	close(eFile)
	predictions<-fit<-apply(fitM,2,mean)
	interval <- match.arg(interval)
	if (interval=="credible"){
		QMb<-apply(fitM,2,quantile,probs=c((1-level)/2,1-(1-level)/2),na.rm=TRUE)
		predictions<-cbind(fit,t(QMb))
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
        Z<-DM(formula=formula2,data=nd,n=npred,knots=object$Zknots,predInd=TRUE,meanVector=object$storeMeanVectorZ,indicator=object$storeIndicatorZ,mvrmObj=object)$X
        Z<-Z[,-1]
        fitV<-matrix(0,nrow=object$nSamples,ncol=npred*object$p)
		alphaFN <- file.path(paste(object$DIR,"BNSP.alpha.txt",sep=""))
        aFile<-file(alphaFN,open="r")
		sigma2FN <- file.path(paste(object$DIR,"BNSP.sigma2.txt",sep=""))
        s2File<-file(sigma2FN,open="r")
        for (i in 1:object$nSamples){
        	alpha<-scan(aFile,what=numeric(),n=object$p*object$LD,quiet=TRUE)
		    s2<-scan(s2File,what=numeric(),n=object$p,quiet=TRUE)
		    for (k in 1:object$p)
                fitV[i,(1+(k-1)*npred):(k*npred)]<-sqrt(s2[k]*exp(as.matrix(Z)%*%matrix(c(alpha[(1+(k-1)*object$LD):(k*object$LD)]))))
		}
		close(aFile)
		close(s2File)
		fitSD<-apply(fitV,2,mean)
		predictions<-cbind(fit,fit-qnorm(1-(1-level)/2)*fitSD,fit+qnorm(1-(1-level)/2)*fitSD,fitSD)
		colnames(predictions) <- c("fit", "lwr", "upr", "fit.sd")
	}
	return(data.frame(predictions))
}

print.mmvrm <- function(x,  digits = 5, ...) {
    cat("\nCall:\n")
    print(x$call)
    cat("\n")
    cat(x$nSamples,"posterior samples\n")
    G<-as.data.frame(mmvrm2mcmc(x,"gamma"))
    D<-as.data.frame(mmvrm2mcmc(x,"delta"))
    colnames(G)<- rep(colnames(x$X)[-1],x$p)
    colnames(D)<- rep(colnames(x$Z)[-1],x$p)
    if (x$LG > 0){
        cat("\nMean model - marginal inclusion probabilities\n")
        cat("\n")
        for (k in 1:x$p){
            cat(x$varsY[[k]])
            cat("\n")
            print(apply(G[(1+(k-1)*x$LG):(k*x$LG)],2,mean), digits=5)
            cat("\n")
		}
    }
    if (x$LD > 0){
        cat("\nVariance model - marginal inclusion probabilities\n")
        cat("\n")
        for (k in 1:x$p){
            cat(x$varsY[[k]])
            cat("\n")
            print(apply(D[(1+(k-1)*x$LD):(k*x$LD)],2,mean), digits=5)
            cat("\n")
		}
    }
}

summary.mmvrm <- function(object, nModels = 5, digits = 5, printTuning = FALSE, ...) {
    cat("\nSpecified model for the mean and variance:\n")
    print(object$formula, showEnv = FALSE)
    cat("\nSpecified priors:\n")
    upTo<-13
    if (object$p > 1 && object$mcm==1) upTo<-16
    if (object$p > 1 && object$mcm>1) upTo<-17
    prior<-NULL
    for (k in 9:upTo) prior<-c(prior,noquote(paste(c(names(as.list(object$call2))[k],object$call2[[k]]),collapse="=")))
    print(noquote(sub("Prior","",prior)))
    cat("\nTotal posterior samples:",object$nSamples,"; burn-in:",object$mcpar[1]-1,"; thinning:",object$mcpar[3],"\n")
    cat("\nFiles stored in",object$DIR,"\n")
    deviance <- matrix(c(object$nullDeviance,object$deviance[1:2]/object$nSamples),ncol=1)
    rownames(deviance) <- c("Null deviance:", "Mean posterior deviance (marginal):", "Mean posterior deviance:")
    colnames(deviance) <- c("")
    print(deviance, digits = digits)
    if ((nModels > 0) && ((object$LG+object$LD) > 0)){
      cat("\nJoint mean/variance model posterior probabilities\n\n")
      G<-as.data.frame(mmvrm2mcmc(object,"gamma"))
      D<-as.data.frame(mmvrm2mcmc(object,"delta"))
      if (object$LG > 0) colnames(G)<-paste(rep(sub(" ",".",paste("mean",colnames(object$X)[-1])),object$p))#,rep(seq(1,p),each=object$LG),sep=".")
      if (object$LD > 0) colnames(D)<- paste(rep(sub(" ",".",paste("var",colnames(object$Z)[-1])),object$p))#,rep(seq(1,p),each=object$LD),sep=".")
      for (k in 1:object$p){
          if (object$LG > 0 && object$LD > 0) Joint<-cbind(G[(1+(k-1)*object$LG):(k*object$LG)],D[(1+(k-1)*object$LD):(k*object$LD)])
          if (object$LG > 0 && object$LD == 0) Joint<-G
          if (object$LG == 0 && object$LD > 0) Joint<-D
          g<-count(Joint)
          g<-g[order(g$freq,decreasing=TRUE),]
          rownames(g)<-seq(1,NROW(g))
          g$prob<-100*g$freq/sum(g$freq)
          g$cumulative<-cumsum(g$prob)
          TrueDim<-min(nModels,dim(g)[1])
          cat(object$varsY[[k]]) #cat(paste(paste("Response model",k),": ",sep=""))
          cat("\n")
          print(g[1:TrueDim,])
          cat("Displaying", TrueDim, "models of the",NROW(g),"visited\n")
          cat(TrueDim,"models account for",sub(" ","",paste(g$cumulative[TrueDim],"%")),"of the posterior mass\n\n")
      }
    }
    if (printTuning){
		cat("\nTuning parameters: start and end values\n\n")
	    dm<-c("start","end")
	    names(object$tuneSigma2R)<-dm
	    names(object$tuneCb)<-dm
	    names(object$tuneR)<-dm
	    ind.f3<-lapply(seq(1:object$p),function(x)seq(x,2*object$p,by=object$p))
	    sigma2<-lapply(ind.f3,function(x){object$tuneSigma2[x]})
	    for (k in 1:object$p) names(sigma2[[k]])<-dm
	    c.alpha<-lapply(ind.f3,function(x){object$tuneCa[x]})
	    for (k in 1:object$p) names(c.alpha[[k]])<-dm
	    if (object$ND > 0){
	        tot<-object$p*object$ND
	        ind.h<-lapply(seq(1:tot),function(x)seq(x,2*tot,by=tot))
	        tuneAlpha<-lapply(ind.h,function(x){object$tuneAlpha[x]})
	        for (k in 1:tot) names(tuneAlpha[[k]])<-dm
	    }
	    pT1<-list(c.beta = object$tuneCb, sigma2 = sigma2)
	    if (object$ND > 0) {pT2<-list(Alpha=tuneAlpha, c.alpha = c.alpha); pT1<-c(pT1,pT2)}
	    if (object$p > 1) {pT2<-list(sigma2R = object$tuneSigma2R, R = object$tuneR); pT1<-c(pT1,pT2)}
	    print(pT1)
	}
}

clustering <- function(object, ...){
    R <- list()
    if (object$mcm==1) stop("Common correlations model")
    if (object$mcm > 1){
        simMatC<-matrix(0,object$d,object$d)
        compAllocFP <- file.path(paste(object$DIR,"BNSP.compAlloc.txt",sep=""))
        compAlloc<-file(compAllocFP,open="r")
        for (j in 1:object$nSamples){
            ca<-scan(compAlloc,what=numeric(),n=object$d,quiet=TRUE)
            for (i in 1:object$d) {
                temp<-which(ca==ca[i])
                simMatC[i,temp]<-simMatC[i,temp]+1
            }
        }
	    close(compAlloc)
	    subset <- rep(seq(1,object$p),each=object$p) < rep(seq(1,object$p),object$p)
	    cor.names <- paste("cor",paste(rep(seq(1,object$p),each=object$p),rep(seq(1,object$p),object$p),sep=""),sep=".")[subset]
	    colnames(simMatC) <- cor.names
	    rownames(simMatC) <- cor.names
	    R[[1]]<-simMatC/object$nSamples
    }
    if (object$mcm==3){
        simMatV<-matrix(0,object$p,object$p)
        compAllocVFP <- file.path(paste(object$DIR,"BNSP.compAllocV.txt",sep=""))
        compAllocV<-file(compAllocVFP,open="r")
        for (j in 1:object$nSamples){
            ca<-scan(compAllocV,what=numeric(),n=object$p,quiet=TRUE)
            for (i in 1:object$p) {
                temp<-which(ca==ca[i])
                simMatV[i,temp]<-simMatV[i,temp]+1
            }
        }
	    close(compAllocV)
	    colnames(simMatV) <- object$varsY
	    rownames(simMatV) <- object$varsY
	    R[[2]]<-simMatV/object$nSamples
    }
    return(R)
}
