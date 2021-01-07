lmrm <- function(formula,data=list(),id,time,
                 sweeps,burn=0,thin=1,seed,StorageDir,
                 c.betaPrior="IG(0.5,0.5*n*p)", pi.muPrior="Beta(1,1)", c.alphaPrior="IG(1.1,1.1)",
                 pi.phiPrior="Beta(1,1)", c.psiPrior="HN(2)",
                 sigmaPrior="HN(2)", pi.sigmaPrior="Beta(1,1)", 
                 corr.Model=c("common",nClust=1), DP.concPrior="Gamma(5,2)",
                 c.etaPrior="IG(0.5,0.5*samp)", pi.nuPrior="Beta(1,1)", pi.fiPrior="Beta(1,1)",
                 c.omegaPrior="IG(1.1,1.1)",sigmaCorPrior="HN(2)",
                 tuneCa,tuneSigma2,tuneCb,tuneAlpha,tuneSigma2R,tuneR,tuneCpsi,tuneCbCor,tuneOmega,tuneComega,
                 tau,FT=1,...){
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
    # Match call
    call <- match.call(expand.dots = FALSE)
    call2 <- match.call.defaults()
    #Data
    if (!is.list(data) && !is.data.frame(data))
        data <- as.data.frame(data)
    #Formula & data dimensions
    p <- length(as.Formula(formula))[1]    
    if (length(as.Formula(formula))[2] != 5) stop("ambiguous definition of regression model")    
    formula.m<-formula(as.Formula(formula),lhs=0,rhs=1)
    formula.v<-formula(as.Formula(formula),lhs=0,rhs=2)
    formula.d<-formula(as.Formula(formula),lhs=0,rhs=3)
    formula.corm<-formula(as.Formula(formula),lhs=0,rhs=4)
    formula.corv<-formula(as.Formula(formula),lhs=0,rhs=5)
    # Responses, design matrices, indicators
    Y<-NULL
    varsY<-list()
    for (i in 1:p){
        trms<-terms.formula(formula(as.Formula(formula,~1),lhs=i,rhs=0))
        Y<-cbind(Y,with(data,eval(attr(trms,"variables")[[2]])))
        varsY[[i]] <- as.character(attr(trms,"variables")[[2]])
    }
    # Null Deviance
    nullDeviance <- 0
    for (i in 1:p) nullDeviance <- nullDeviance - 2 * logLik(lm(Y[,i] ~ 1))
    # Response etc
    if (missing(id)) stop("id needed")    
    id<-eval(substitute(id), data)
    if (is.null(id)) stop("unrecognised id")
    n<-length(unique(id))
    #print("n");print(n)
    niVec<-table(id)[1:n] 
    niMax<-max(niVec)
    N<-sum(niVec)
    cusumniVec<-c(0,cumsum(niVec))
    cusumC<-c(0,cumsum(niVec*(niVec-1)/2))     
    # Times
    varTime<-as.character(substitute(list(time))[-1])
    time2<-eval(substitute(time), data)
    LUT<-length(unique(time2)) 
    SUT<-sort(unique(time2))
    FUT<-table(time2)
    intime<-time2
    for (i in 1:length(intime)) intime[i] <- which(intime[i]==SUT)-1   
    intime2 <- rep(intime,each=p)
    #Desing matrix C (dependence)
    dataNew<-as.data.frame(cbind(data,lag=rnorm(dim(data)[1])))
    bb<-DM(formula=formula.d,data=dataNew,n=dim(data)[1])
    vars <- unlist(bb$vars)
    if ("lag" %in% vars)
        vars <- vars[-which(vars=="lag")]   
    C<-NULL
    for (i in 1:n) 
        if (niVec[i] > 1) C<-rbind(C,cbind(cusumniVec[i]+rep(seq(2,niVec[i],1),seq(1,niVec[i]-1,1)), 
                                           cusumniVec[i]+unlist(sapply(1:(niVec[i]-1), function(i) seq(1,i,1)))                                     
                                     ))
    lag<-time2[C[,1]]-time2[C[,2]]
    lag<-data.frame(lag=lag)    
    if (length(vars)>0){
        for (k in 1:length(vars)){
            V<-data[,vars[k]]
            lag<-data.frame(lag,V[C[,1]])
            colnames(lag)[k+1] <- vars[[k]]
	    }
	}
	#print(lag)
	if (niMax>1){    
        XYK<-DM(formula=formula.d,data=lag,n=dim(lag)[1])
        C<-as.matrix(XYK$X)
    }else{
    	XYK<-NULL
	    C<-NULL
	}
    Cknots<-XYK$Rknots
    storeMeanVectorC<-XYK$meanVector
    storeIndicatorC<-XYK$indicator
    LK<-NCOL(C)
    if (!is.null(C)) C<-t(C)
    vecLK<-table(XYK$assign)
    NK<-length(vecLK)
    cusumVecLK<-c(0,cumsum(vecLK))
    assignC<-XYK$assign
    labelsC<-XYK$labels
    countC<-XYK$count
    varsC<-XYK$vars
    is.Dc<-XYK$is.D
    which.SpecC<-XYK$which.Spec
    formula.termsC<-XYK$formula.terms
    #Desing matrix X (mean)
    XYK<-DM(formula=formula.m,data=data,n=N)
    X<-as.matrix(XYK$X)
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
    #print("2")
    #Desing matrix Z (variances)
    ZK<-DM(formula=formula.v,data=data,n=N)
    Z<-as.matrix(ZK$X)
    Zknots<-ZK$Rknots
    storeMeanVectorZ<-ZK$meanVector
    storeIndicatorZ<-ZK$indicator
    LD<-NCOL(Z)-1
    vecLD<-table(ZK$assign)[-1]
    ND<-length(vecLD)
    cusumVecLD<-c(0,cumsum(vecLD))
    MVLD <- 1
    if (LD > 0 || LK > 0) MVLD<-max(vecLD,vecLK)
    assignZ<-ZK$assign
    labelsZ<-ZK$labels
    countZ<-ZK$count
    varsZ<-ZK$vars
    is.Dz<-ZK$is.D
    which.SpecZ<-ZK$which.Spec
    formula.termsZ<-ZK$formula.terms
    #Design matrix Xc (mean of correlations)
    d<-p*(p-1)/2
    t <- rep(SUT,each=d)
    dataCor<-data.frame(t)
    if (p > 1){
		lab<-DM(formula=formula.corm,data=data,n=1)        
        if (length(lab$vars) > 0) colnames(dataCor)<-lab$vars[[1]]   
		XYK<-DM(formula=formula.corm,data=dataCor,n=(LUT*d))
		Xc<-as.matrix(XYK$X)
	}else{
	    XYK<-NULL
	    Xc<-NULL
	}    
    Xcknots<-XYK$Rknots
    storeMeanVectorXc<-XYK$meanVector
    storeIndicatorXc<-XYK$indicator
    LGc<-NCOL(Xc)-1
    vecLGc<-table(XYK$assign)[-1]
    NGc<-length(vecLGc)
    cusumVecLGc<-c(0,cumsum(vecLGc))
    assignXc<-XYK$assign
    labelsXc<-XYK$labels
    countXc<-XYK$count
    varsXc<-XYK$vars
    is.Dxc<-XYK$is.D
    which.SpecXc<-XYK$which.Spec
    formula.termsXc<-XYK$formula.terms
    #Design matrix Zc (variance of correlations)
    if (p > 1){
		XYK<-DM(formula=formula.corv,data=dataCor,n=(LUT*d))    
        Zc<-as.matrix(XYK$X)
    }else{
		XYK<-NULL
		Zc<-NULL
	}
    Zcknots<-XYK$Rknots
    storeMeanVectorZc<-XYK$meanVector
    storeIndicatorZc<-XYK$indicator
    LDc<-NCOL(Zc)-1
    vecLDc<-table(XYK$assign)[-1]
    MVLD<-max(vecLDc,MVLD)
    NDc<-length(vecLDc)
    cusumVecLDc<-c(0,cumsum(vecLDc))
    assignZc<-XYK$assign
    labelsZc<-XYK$labels
    countZc<-XYK$count
    varsZc<-XYK$vars
    is.Dzc<-XYK$is.D
    which.SpecZc<-XYK$which.Spec
    formula.termsZc<-XYK$formula.terms   
    #Initialize covariance & correlation matrix
    #print("3")
    LASTR<-LASTD<-LASTE<-rep(1,LUT)    
    if (p>1){        
        LASTR<-LASTD<-LASTE<-NULL
        for (t in SUT){
            Res<-NULL
            for (i in 1:p){            
				#print(head(X))
				#print(grep("(",colnames(X),fixed=TRUE)) 
				#RM<-grep("(",colnames(X),fixed=TRUE)
				#print(RM)    			
                RM<-grep("(",colnames(X),fixed=TRUE)[-1]
				#print(RM)
				Xinit<-X[time2==t,]
                if (length(RM) > 0) Xinit<-Xinit[,-RM]
                #print(head(Xinit))
                lm1<-lm(Y[time2==t,i] ~ Xinit)
                Res<-cbind(Res,residuals(lm1))
			}
            CR<-0.9*cov(Res)+0.1*diag(1,p)            
            LASTR<-c(LASTR,c(cov2cor(CR)))            
            D<-matrix(0,p,p)
            diag(D)<-sqrt(diag(CR))
            LASTD<-c(LASTD,c(D))
            LASTE<-c(LASTE,c(CR))
            #LASTmuR<-mean(LASTR[upper.tri(LASTR)])
            #LASTsigma2R<-1
            #if (p > 2) LASTsigma2R<-var(LASTR[upper.tri(LASTR)])
            #print(CR)
            #print(D)
            #print(cov2cor(CR))            
		}
    }  
    LASTAll<-c(LASTR,LASTD,LASTE)    
    #print("4")
    #print(LASTAll)
    #print("2")     
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
    if (length(pi.muPrior)==1) pimu<-rep(pimu,p*NG)
    if (length(pi.muPrior)==NG) pimu<-rep(pimu,p)
    #Prior for pi.phi
    if (!length(pi.phiPrior)==1 && !length(pi.phiPrior)==NK && !length(pi.phiPrior)==(p*p*NK))
        stop("pi.phiPrior has incorrect dimension")
    piphi<-NULL
    for (k in 1:length(pi.phiPrior)){
        sp<-strsplit(pi.phiPrior[k],"Beta\\(")
        sp<-strsplit(sp[[1]][2],"\\)")
        sp<-strsplit(sp[[1]][1],",")
        piphi<-c(piphi,as.numeric(sp[[1]]))
    }
    if (length(pi.phiPrior)==1) piphi<-rep(piphi,p*p*NK)
    if (length(pi.phiPrior)==NK) piphi<-rep(piphi,p*p)
    #Prior for c.beta
    if (!length(c.betaPrior)==1)
        stop("c.betaPrior has incorrect dimension")
    #c.betaPrior<-sub("samp","p*n",c.betaPrior)
    sp<-strsplit(c.betaPrior,"IG\\(")
    sp<-strsplit(sp[[1]][2],"\\)")
    sp<-strsplit(sp[[1]][1],",")
    cetaParams<-c(as.numeric(sp[[1]][1]),eval(parse(text=sp[[1]][2])))
    #print(cetaParams)
    #Prior for c.eta
    if (!length(c.etaPrior)==1)
        stop("c.etaPrior has incorrect dimension")
    c.etaPrior<-sub("samp","LUT*d",c.etaPrior)
    sp<-strsplit(c.etaPrior,"IG\\(")
    sp<-strsplit(sp[[1]][2],"\\)")
    sp<-strsplit(sp[[1]][1],",")
    cetaCorParams<-c(as.numeric(sp[[1]][1]),eval(parse(text=sp[[1]][2])))
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
    if (length(pi.sigmaPrior)==1) pisigma<-rep(pisigma,p*ND)
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
        } else stop("unrecognised prior for c.alpha")
        sp<-strsplit(sp[[1]][2],"\\)")
        sp<-strsplit(sp[[1]][1],",")
        calphaParams<-c(calphaParams,as.numeric(sp[[1]]))
    }
    if (length(c.alphaPrior)==1){
        calphaParams<-rep(calphaParams,p)
        HNca<-rep(HNca,p)
	}
	#Prior for cPsi
    if (!length(c.psiPrior)==1 && !length(c.psiPrior)==(p*p))
        stop("c.psiPrior has incorrect dimension")
    specials<-c("HN","IG")
    cpsiParams<-NULL
    HNcpsi<-vector()
    for (k in 1:length(c.psiPrior)){
        sp<-strsplit(c.psiPrior[k],"\\(")
        if (sp[[1]][1] %in% specials){
            if (match(sp[[1]][1],specials)==1) HNcpsi[k]<-1
            if (match(sp[[1]][1],specials)==2) HNcpsi[k]<-0
        } else stop("unrecognised prior for c.psi")
        sp<-strsplit(sp[[1]][2],"\\)")
        sp<-strsplit(sp[[1]][1],",")
        cpsiParams<-c(cpsiParams,as.numeric(sp[[1]]))
	}
	if (length(c.psiPrior)==1){
        cpsiParams<-rep(cpsiParams,(p*p))
        HNcpsi<-rep(HNcpsi,(p*p))
	}    
    #Prior for sigma2_{zero k}
    if (!length(sigmaPrior)==1 && !length(sigmaPrior)==p)
        stop("sigmaPrior has incorrect dimension")
    specials<-c("HN","IG")
    sigmaParams<-NULL
    HNsg<-vector()
    for (k in 1:length(sigmaPrior)){
        sp<-strsplit(sigmaPrior[k],"\\(")
        if (sp[[1]][1] %in% specials){
            if (match(sp[[1]][1],specials)==1) HNsg[k]<-1
            if (match(sp[[1]][1],specials)==2) HNsg[k]<-0
        } else stop("unrecognised prior for sigma^2")
        sp<-strsplit(sp[[1]][2],"\\)")
        sp<-strsplit(sp[[1]][1],",")
        sigmaParams<-c(sigmaParams,as.numeric(sp[[1]]))
    }
    if (length(sigmaPrior)==1){
        sigmaParams<-rep(sigmaParams,p)
        HNsg<-rep(HNsg,p)
	}
    #Prior for pi.nu
    if (!length(pi.nuPrior)==NGc) if (!length(pi.nuPrior)==1) 
        stop("pi.nuPrior has incorrect dimension")
    pinu<-NULL
    for (k in 1:length(pi.nuPrior)){
        sp<-strsplit(pi.nuPrior[k],"Beta\\(")
        sp<-strsplit(sp[[1]][2],"\\)")
        sp<-strsplit(sp[[1]][1],",")
        pinu<-c(pinu,as.numeric(sp[[1]]))
    }
    if (length(pi.nuPrior)==1) pinu<-rep(pinu,NGc)
    #Prior for pi.fi
    if (!length(pi.fiPrior)==NDc) if (!length(pi.fiPrior)==1) 
        stop("pi.fiPrior has incorrect dimension")
    pifi<-NULL
    for (k in 1:length(pi.fiPrior)){
        sp<-strsplit(pi.fiPrior[k],"Beta\\(")
        sp<-strsplit(sp[[1]][2],"\\)")
        sp<-strsplit(sp[[1]][1],",")
        pifi<-c(pifi,as.numeric(sp[[1]]))
    }
    if (length(pi.fiPrior)==1) pifi<-rep(pifi,ND)    
    #Prior for c.omega
    if (!length(c.omegaPrior)==1)
        stop("c.omegaPrior has incorrect dimension")
    specials<-c("HN","IG")
    sp<-strsplit(c.omegaPrior,"\\(")
    if (sp[[1]][1] %in% specials){ 
        if (match(sp[[1]][1],specials)==1) HNco<-1
        if (match(sp[[1]][1],specials)==2) HNco<-0
    } else stop("unrecognised prior for c.omega")
    sp<-strsplit(sp[[1]][2],"\\)")
    sp<-strsplit(sp[[1]][1],",")
    comegaParams<-as.numeric(sp[[1]])
	#Prior for sigma2Cor
	if (!length(sigmaCorPrior)==1)
        stop("sigmaCorPrior has incorrect dimension")
    HNscor<-0    
    specials<-c("HN","IG")
    sp<-strsplit(sigmaCorPrior,"\\(")
    if (sp[[1]][1] %in% specials){ 
        if (match(sp[[1]][1],specials)==1 && p > 1) HNscor<-1
        if (match(sp[[1]][1],specials)==2) HNscor<-0
    } else stop("unrecognised prior for sigma2Cor")
    sp<-strsplit(sp[[1]][2],"\\)")
    sp<-strsplit(sp[[1]][1],",")
    sigmaCorParams<-as.numeric(sp[[1]])
    #Cor model
    corModels<-c("common","groupC","groupV")
    mcm<-4+match(corr.Model[1],corModels)
    if (p==1) mcm=5
    if (p==2) mcm=5
    if (is.na(mcm)) stop("unrecognised correlation model")
    H <- G <- 1
    if (mcm==6){
        H<-as.numeric(corr.Model[2])
        if (H==1) {mcm<-5; warning("Common correlations model specified with nClust = 1")}
        if (is.na(H) || (!H%%1==0) || H==0) {H <- p*(p-1)/2; warning(cat("mispecified number of clusters. nClust set to ",H,"\n"))}
    }
    if (mcm==7){
        G<-as.numeric(corr.Model[2])
        if (G==1) {mcm<-5; warning("Common correlations model specified with nClust = 1")}
        if (is.na(G) || (!G%%1==0) || G==0) {G <- p; warning(cat("mispecified number of clusters. nClust set to ",G,"\n"))}
        H<-G*(G-1)/2+G #min(d,G*(G-1)/2+G) #min(G,abs(p-G))
	}
    #Prior for alpha DP
    if (mcm > 5){
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
    FL <- c("gamma", "cbeta", "delta", "alpha", "sigma2", "calpha", "beta", "psi", "ksi", "cpsi", "R", "nu", "fi", "omega",
            "sigma2R","ceta","comega","eta", "deviance", 
            "compAlloc", "nmembers", "DPconc",
            "compAllocV","nmembersV",
            "DE", "test","nu.ls","eta.ls","nmembers.ls","clusters","probs","tune")
    for (i in 1:length(FL)){
        oneFile <- paste(StorageDir, paste("BNSP",FL[i], "txt",sep="."),sep="/")
        if (file.exists(oneFile)) file.remove(oneFile)
	}
    #Tuning Parameters            
    if (missing(tuneCa)) tuneCa<-rep(1,p)
    if (!length(tuneCa)==p) tuneCa<-rep(mean(tuneCa),p)        
    if (missing(tuneSigma2)) tuneSigma2<-rep(1,p)
    if (!length(tuneSigma2)==p) tuneSigma2<-rep(mean(tuneSigma2),p)                 
    if (missing(tuneCb)) tuneCb<-100
    if (missing(tuneAlpha)) tuneAlpha<-rep(5,ND*p)
    if (!length(tuneAlpha)==(ND*p)) tuneAlpha<-rep(mean(tuneAlpha),ND*p)
    if (missing(tuneSigma2R)) tuneSigma2R<-1
    if (missing(tuneR)) tuneR<-rep(40*(p+2)^3,LUT)
    tuneR[which(tuneR<p+2)]<-p+2
    if (!length(tuneR)==LUT) tuneR<-rep(mean(tuneR),LUT)
    if (missing(tuneCpsi)) tuneCpsi<-rep(5,p*p)
    if (!length(tuneCpsi)==(p*p)) tuneCpsi<-rep(mean(tuneCpsi),p*p)
	if (missing(tuneCbCor)) tuneCbCor<-10
	if (missing(tuneOmega)) tuneOmega<-rep(5,NDc)
	if (!length(tuneOmega)==NDc) tuneOmega<-rep(mean(tuneOmega),NDc)
    if (missing(tuneComega)) tuneComega<-1
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
    #if (missing(blockSizeProbK)){
	blockSizeProbK <- rep(0,LK)
    blockSizeProbK[1:5]<-c(10,25,30,25,10)
    #}
    #if (missing(blockSizeProbGc)){ 
    blockSizeProbGc <- rep(0,LGc) 
    blockSizeProbGc[1:5]<-c(10,25,30,25,10)
	#}
    #if (missing(blockSizeProbD)){
	blockSizeProbDc <- rep(0,LDc)
    blockSizeProbDc[1:5]<-c(10,25,30,25,10)
    #}
    maxBSG <- max(which(blockSizeProbG>0))
    maxBSD <- max(which(blockSizeProbD>0))
    maxBSK <- max(which(blockSizeProbK>0))
    maxBSGc <- max(which(blockSizeProbGc>0))
    maxBSDc <- max(which(blockSizeProbDc>0))
    
    #print(as.double(c(blockSizeProbG, blockSizeProbD, blockSizeProbK, blockSizeProbGc, blockSizeProbDc)))
    #print(as.integer(c(maxBSG,maxBSD,maxBSK,maxBSGc,maxBSDc)))
    
    #Deviance
    deviance <- c(0,0)  
    #tol
    tol <- sqrt(.Machine$double.eps)  
    #Call C
    if (mcm==5){
        out<-.C("longmult",
        as.integer(seed), as.character(StorageDir), as.integer(WF),
        as.integer(c(sweeps,burn,thin,n,p,N,MVLD)),
        as.integer(niVec), as.integer(cusumniVec), as.integer(intime2), as.integer(intime), 
        as.integer(c(niMax,LUT)), as.integer(FUT),
        as.integer(cusumC), as.double(C), as.double(c(t(Y))), as.double(t(X)), as.double(Z), 
        as.double(Xc),as.double(Zc),
        as.integer(c(LG,LD,LK,LGc,LDc)), as.integer(c(NG,ND,NK,NGc,NDc)), 
        as.integer(vecLG), as.integer(vecLD), as.integer(vecLK), as.integer(vecLGc),as.integer(vecLDc), 
        as.integer(cusumVecLG), as.integer(cusumVecLD), as.integer(cusumVecLK), as.integer(cusumVecLGc),as.integer(cusumVecLDc),        
        as.double(c(blockSizeProbG, blockSizeProbD, blockSizeProbK, blockSizeProbGc, blockSizeProbDc)), 
        as.integer(c(maxBSG,maxBSD,maxBSK,maxBSGc,maxBSDc)),                        
        as.double(tuneCa), as.double(tuneSigma2), as.double(tuneCb), as.double(tuneAlpha), as.double(tuneSigma2R), 
        as.double(tuneR), as.double(tuneCpsi), as.double(tuneCbCor), as.double(tuneOmega), as.double(tuneComega),                                    
        as.double(pimu), as.double(cetaParams), as.double(pisigma), as.integer(HNca), as.double(calphaParams),            
        as.integer(HNsg), as.double(sigmaParams), as.double(piphi), as.integer(HNcpsi), as.double(cpsiParams),
        as.double(cetaCorParams),as.integer(HNco),as.double(comegaParams),as.double(pinu),as.double(pifi),
        as.integer(HNscor),as.double(sigmaCorParams),
        as.double(c(tau,tol)), as.integer(FT), as.double(deviance),            
        as.integer(c(0,LASTsw,LASTWB)),as.double(LASTAll))}
    if (mcm==6){
        out<-.C("longmultg",
        as.integer(seed), as.character(StorageDir), as.integer(WF),
        as.integer(c(sweeps,burn,thin,n,p,N,MVLD)),
        as.integer(niVec), as.integer(cusumniVec), as.integer(intime2), as.integer(intime), 
        as.integer(c(niMax,LUT)), as.integer(FUT),
        as.integer(cusumC), as.double(C), as.double(c(t(Y))), as.double(t(X)), as.double(Z), 
        as.double(Xc),as.double(Zc),
        as.integer(c(LG,LD,LK,LGc,LDc)), as.integer(c(NG,ND,NK,NGc,NDc)), 
        as.integer(vecLG), as.integer(vecLD), as.integer(vecLK), as.integer(vecLGc),as.integer(vecLDc), 
        as.integer(cusumVecLG), as.integer(cusumVecLD), as.integer(cusumVecLK), as.integer(cusumVecLGc),as.integer(cusumVecLDc),
        as.double(c(blockSizeProbG, blockSizeProbD, blockSizeProbK, blockSizeProbGc, blockSizeProbDc)), 
        as.integer(c(maxBSG,maxBSD,maxBSK,maxBSGc,maxBSDc)),                        
        as.double(tuneCa), as.double(tuneSigma2), as.double(tuneCb), as.double(tuneAlpha), as.double(tuneSigma2R), 
        as.double(tuneR), as.double(tuneCpsi), as.double(tuneCbCor), as.double(tuneOmega), as.double(tuneComega),  
        as.double(pimu), as.double(cetaParams), as.double(pisigma), as.integer(HNca), as.double(calphaParams),            
        as.integer(HNsg), as.double(sigmaParams), as.double(piphi), as.integer(HNcpsi), as.double(cpsiParams),
        as.double(cetaCorParams),as.integer(HNco),as.double(comegaParams),as.double(pinu),as.double(pifi),
        as.integer(HNscor),as.double(sigmaCorParams),
        as.double(tau), as.integer(FT), as.double(deviance),            
        as.integer(c(0,LASTsw,LASTWB)),as.double(LASTAll),
        as.integer(H), as.double(DPparams))}   
    if (mcm==7){
        out<-.C("longmultgv",
        as.integer(seed), as.character(StorageDir), as.integer(WF),
        as.integer(c(sweeps,burn,thin,n,p,N,MVLD)),
        as.integer(niVec), as.integer(cusumniVec), as.integer(intime2), as.integer(intime), 
        as.integer(c(niMax,LUT)), as.integer(FUT),
        as.integer(cusumC), as.double(C), as.double(c(t(Y))), as.double(t(X)), as.double(Z), 
        as.double(Xc),as.double(Zc),
        as.integer(c(LG,LD,LK,LGc,LDc)), as.integer(c(NG,ND,NK,NGc,NDc)), 
        as.integer(vecLG), as.integer(vecLD), as.integer(vecLK), as.integer(vecLGc),as.integer(vecLDc), 
        as.integer(cusumVecLG), as.integer(cusumVecLD), as.integer(cusumVecLK), as.integer(cusumVecLGc),as.integer(cusumVecLDc),          
        as.double(c(blockSizeProbG, blockSizeProbD, blockSizeProbK, blockSizeProbGc, blockSizeProbDc)), 
        as.integer(c(maxBSG,maxBSD,maxBSK,maxBSGc,maxBSDc)),                        
        as.double(tuneCa), as.double(tuneSigma2), as.double(tuneCb), as.double(tuneAlpha), as.double(tuneSigma2R), 
        as.double(tuneR), as.double(tuneCpsi), as.double(tuneCbCor), as.double(tuneOmega), as.double(tuneComega),  
        as.double(pimu), as.double(cetaParams), as.double(pisigma), as.integer(HNca), as.double(calphaParams),            
        as.integer(HNsg), as.double(sigmaParams), as.double(piphi), as.integer(HNcpsi), as.double(cpsiParams),
        as.double(cetaCorParams),as.integer(HNco),as.double(comegaParams),as.double(pinu),as.double(pifi),
        as.integer(HNscor),as.double(sigmaCorParams),
        as.double(tau), as.integer(FT), as.double(deviance),            
        as.integer(c(0,LASTsw,LASTWB)), as.double(LASTAll),
        as.integer(G), as.double(DPparams))}   
    #Output  
    if (!is.null(C)) C<-t(C)
    loc1<-32 
    loc2<-41 
    tuneSigma2Ra<-out[[loc1+4]][1]
    tuneRa<-out[[loc1+5]][1:LUT]
    fit <- list(call=call,call2=call2,formula=formula,seed=seed,p=p,d=p*(p-1)/2,
                data=data,lag=lag,Y=Y,
                X=X,Xknots=Xknots,LG=LG,NG=NG,
                Z=Z,Zknots=Zknots,LD=LD,ND=ND,
                storeMeanVectorX=storeMeanVectorX,
                storeMeanVectorZ=storeMeanVectorZ,
                storeIndicatorX=storeIndicatorX,
                storeIndicatorZ=storeIndicatorZ,
                assignX=assignX,
                assignZ=assignZ,
                labelsX=labelsX,
                labelsZ=labelsZ,
                countX=countX,
                countZ=countZ,
                varsY=varsY,                
                varsX=varsX,
                varsZ=varsZ,
                is.Dx=is.Dx,
                is.Dz=is.Dz,
                which.SpecX=which.SpecX,
                which.SpecZ=which.SpecZ,
                formula.termsX=formula.termsX,
                formula.termsZ=formula.termsZ,                
                nSamples=nSamples,
                totalSweeps=sweeps,
                mcpar=c(as.integer(burn+1),as.integer(seq(from=burn+1,by=thin,length.out=nSamples)[nSamples]),as.integer(thin)),                                
                mcm=mcm,H=H,G=G,                
                tuneCa=c(tuneCa,out[[loc1+0]][1:p]),        
                tuneSigma2=c(tuneSigma2,out[[loc1+1]][1:p]),
                tuneCb=c(tuneCb,out[[loc1+2]][1]),                                
                tuneAlpha=c(tuneAlpha,out[[loc1+3]][1:(p*ND)]),                
                tuneSigma2R=c(tuneSigma2R,tuneSigma2Ra),
                tuneR=c(tuneR,tuneRa),
                deviance=c(out[[loc2]][1:2]),
                nullDeviance=nullDeviance,
                DIR=StorageDir,  
                out=out,
                LUT=LUT,
                SUT=SUT,
                C=C,Cknots=Cknots,LK=LK,NK=NK,
                Xc=Xc,Xcknots=Xcknots,LGc=LGc,NGc=NGc,
                Zc=Zc,Zcknots=Zcknots,LDc=LDc,NDc=NDc,                
                storeMeanVectorC=storeMeanVectorC,
                storeMeanVectorXc=storeMeanVectorXc,
                storeMeanVectorZc=storeMeanVectorZc,                               
                storeIndicatorC=storeIndicatorC,
                storeIndicatorXc=storeIndicatorXc,
                storeIndicatorZc=storeIndicatorZc, 
                assignC=assignC,          
                assignXc=assignXc,
                assignZc=assignZc,
                labelsC=labelsC,                
                labelsXc=labelsXc,
                labelsZc=labelsZc,                                               
                countC=countC,                            
                countXc=countXc,
                countZc=countZc,
                varsC=varsC,                
                varsXc=varsXc,
                varsZc=varsZc,
                varTime=varTime,  
                is.Dc=is.Dc,
                is.Dxc=is.Dxc,
                is.Dzc=is.Dzc,  
                which.SpecC=which.SpecC,                
                which.SpecXc=which.SpecXc,
                which.SpecZc=which.SpecZc, 
                formula.termsC=formula.termsC,
                formula.termsXc=formula.termsXc,
                formula.termsZc=formula.termsZc,                 
                tuneCpsi=c(tuneCpsi,out[[38]][1:(p*p)]),                
	            tuneCbCor=c(tuneCbCor,out[[39]][1]),	            
	            tuneOmega=c(tuneOmega,out[[40]][1:NDc]),	            
                tuneComega=c(tuneComega,out[[41]][1]),
                HNca=HNca,
                HNcpsi=HNcpsi,
                HNscor=HNscor,
                HNco=HNco,
                niVec=niVec,
                intime=intime)                                                    
    class(fit) <- 'mvrm'
    return(fit)
} 

#psiAll <- mvrm2mcmc(x,"psi")
#alphaAll <- mvrm2mcmc(x,"alpha") 
#sigma2All <- mvrm2mcmc(x,"sigma2")
#RAll <- mvrm2mcmc(x,"R")    

postSigma<-function(x, time, psi, alpha, sigma2, R, samples, ...){
    p<-x$p
    ifelse(missing(time), t<-x$SUT, t<-time)
    T<-length(t)
    postS<-matrix(0,nrow=T*p,ncol=T*p)    

    # Vars
    varsd<-as.list(substitute(list(...)))[-1]
    if (x$varTime %in% names(varsd)) 
        varsd <-varsd[-which(names(varsd)==x$varTime)]
    if ("lag" %in% names(varsd)) 
        varsd <-varsd[-which(names(varsd)=="lag")]

    comp.des.mats<-1
    if (comp.des.mats==1){
		comp.des.mats<-0
        # Zmat
        vars<-x$varsZ
        wv<-NULL 
        newdata<-t
        for (v in 1:length(vars)){
	    	if (vars[[v]] %in% names(varsd)){
                min1<-eval(varsd[[which(vars[[v]]==names(varsd))]])
                newdata<-cbind(newdata,rep(min1,T))
                wv<-c(wv,v)
		    }
        }
        if (! (length(wv)+1) == length(vars)) stop("insufficient input on covariates")
        newdata <- as.data.frame(newdata) 
        colnames(newdata)<-c(x$varTime,vars[wv])
        terms.reform<-NULL
        k<-0
        for (i in 1:length(x$formula.termsZ)){
            term<-x$formula.termsZ[i]
            if (!i %in%  which(unlist(x$which.SpecZ)==-99)){
                k<-k+1
                if (!grepl("knots",term)){
                    term<-substr(term,1,nchar(term)-1)
     		        term<-paste(term,",knots= knots[[",k,"]])")
			    }
	        }
            terms.reform<-c(terms.reform,term)
        }
        formula2<-reformulate(terms.reform)
        if (length(x$data)>0){
            nd<-x$data[0,match(colnames(newdata),colnames(x$data)),drop=FALSE]
            for (j in 1:dim(nd)[2])
                nd[,j]<-drop(nd[,j])  
            nd[1:NROW(newdata),] <- newdata
        }else{nd<-newdata}
        Zmat<-DM(formula=formula2,data=nd,n=NROW(nd),knots=x$Zknots,meanVector=x$storeMeanVectorZ,indicator=x$storeIndicatorZ)$X
        Zmat<-as.matrix(Zmat)[,-1]
    
        # Cmat  
        vars<-x$varsC
        wv<-NULL 
        C<-cbind(rep(seq(2,T,1),seq(1,T-1,1)), unlist(sapply(1:(T-1), function(i) seq(1,i,1))))
        lag<-t[C[,1]]-t[C[,2]]
        newdata<-lag
        for (v in 1:length(vars)){
	    	if (vars[[v]] %in% names(varsd)){
                min1<-eval(varsd[[which(vars[[v]]==names(varsd))]])
                newdata<-cbind(newdata,rep(min1,T))
                wv<-c(wv,v)
		    }
        }
        if (! (length(wv)+1) == length(vars)) stop("insufficient input on covariates")
        newdata<-as.data.frame(newdata)
        colnames(newdata)<-c("lag",vars[wv])
        newdata <- as.data.frame(newdata) 
        terms.reform<-NULL
        k<-0
        for (i in 1:length(x$formula.termsC)){
            term<-x$formula.termsC[i]
            if (!i %in%  which(unlist(x$which.SpecC)==-99)){
                k<-k+1
                if (!grepl("knots",term)){
                    term<-substr(term,1,nchar(term)-1)
     		        term<-paste(term,",knots= knots[[",k,"]])")
			    }
	        }
            terms.reform<-c(terms.reform,term)
        }
        formula2<-reformulate(terms.reform)
        if (length(x$data)>0){
	    	dataNew<-as.data.frame(cbind(x$data,lag=rnorm(dim(x$data)[1])))
            nd<-dataNew[0,match(colnames(newdata),colnames(dataNew)),drop=FALSE]
            for (j in 1:dim(nd)[2])
                nd[,j]<-drop(nd[,j])
            nd[1:NROW(newdata),] <- newdata
        }else{nd<-newdata}    
        Cmat<-DM(formula=formula2,data=nd,n=NROW(newdata),knots=x$Cknots,meanVector=x$storeMeanVectorC,indicator=x$storeIndicatorC)$X   
        Cmat<-as.matrix(Cmat)
    }
    
    # #
    Sij<-matrix(0,p,p)
    RA<-diag(rep(1,p))
    Li<-diag(rep(1,p*T))
    Di<-matrix(0,p*T,p*T)
    
    if (missing(samples)) samples <- 1:x$nSamples
    
    for (sw in samples){		        
	    buildPhi<-NULL
		for (r1 in 1:p){
	        for (r2 in 1:p){
		        Pair<-(r1-1)*p + r2
		        psiA<-psi[sw,(1+((Pair-1)*x$LK)):(Pair*x$LK)]		        
		        buildPhi<-cbind(buildPhi,Cmat%*%matrix(psiA))
            }
		}				        
		move<-1
        if (T > 1){
            for (j in 2:T){
                for (k in 1:(j-1)){
                    Li[(j*p-p+1):(j*p),(k*p-p+1):(k*p)] <- -matrix(buildPhi[move,],p,p,byrow=TRUE) 		
                    move <- move + 1
			    }
		    }
		}
		Sr<-NULL
		for (r in 1:p){
			s2<-sigma2[sw,r]
            alphaA<-alpha[sw,(1+((r-1)*x$LD)):(r*x$LD)]                                                            
            Sr<-cbind(Sr,sqrt(s2*exp(Zmat%*%matrix(alphaA))))
        }
        for (j in 1:T){
            diag(Sij) <- Sr[j,]  
            pick.time<-which(t[j]==x$SUT)#j ; x$intime[cusumniVec[subject]+j]+1
            RA[lower.tri(RA)] <- R[sw,(1+(pick.time-1)*p*(p-1)/2):(pick.time*p*(p-1)/2)]
            RA[upper.tri(RA)] <- t(RA)[upper.tri(RA)] 
            Di[((j-1)*p+1):(j*p),((j-1)*p+1):(j*p)] <- Sij %*% RA %*% Sij
        }       
        SigmaSW <- solve(Li)%*%Di%*%solve(t(Li))
        postS <- postS + SigmaSW
	}
	return(postS/length(samples))
}

postSigma2<-function(x, subject, psi, alpha, sigma2, R, samples, ...){
    niVec<-x$niVec
    p<-x$p
    T<-niVec[subject]
    postS<-matrix(0,nrow=T*p,ncol=T*p)    
    cusumniVec<-c(0,cumsum(niVec))
    cusumC<-c(0,cumsum(niVec*(niVec-1)/2))                
    Cmat<-x$C[(cusumC[subject]+1):cusumC[subject+1],]      
    Zmat<-x$Z[(cusumniVec[subject]+1):cusumniVec[subject+1],-1]        
    Sij<-matrix(0,p,p)
    RA<-diag(rep(1,p))
    Li<-diag(rep(1,p*T))
    Di<-matrix(0,p*T,p*T)
    if (missing(samples)) samples <- 1:x$nSamples
    for (sw in samples){		        
	    buildPhi<-NULL
		for (r1 in 1:p){
	        for (r2 in 1:p){
		        Pair<-(r1-1)*p + r2
		        psiA<-psi[sw,(1+((Pair-1)*x$LK)):(Pair*x$LK)]		        
		        buildPhi<-cbind(buildPhi,Cmat%*%matrix(psiA))
            }
		}				        
		move<-1
        for (j in 2:T){
            for (k in 1:(j-1)){
                Li[(j*p-p+1):(j*p),(k*p-p+1):(k*p)] <- -matrix(buildPhi[move,],p,p,byrow=TRUE) 		
                move <- move + 1
			}
		}
		Sr<-NULL
		for (r in 1:p){
			s2<-sigma2[sw,r]
            alphaA<-alpha[sw,(1+((r-1)*x$LD)):(r*x$LD)]                                                            
            Sr<-cbind(Sr,sqrt(s2*exp(Zmat%*%matrix(alphaA))))
        }
        for (j in 1:T){
            diag(Sij) <- Sr[j,]  
            pick.time<-x$intime[cusumniVec[subject]+j]+1              
            RA[lower.tri(RA)] <- R[sw,(1+(pick.time-1)*p*(p-1)/2):(pick.time*p*(p-1)/2)]
            RA[upper.tri(RA)] <- t(RA)[upper.tri(RA)]                                    
            Di[((j-1)*p+1):(j*p),((j-1)*p+1):(j*p)] <- Sij %*% RA %*% Sij
        }       
        SigmaSW <- solve(Li)%*%Di%*%solve(t(Li))
        postS <- postS + SigmaSW
	}
	return(postS/length(samples))
}

plot3.bcmg<-function(x, model="mean", centre=mean, quantiles=c(0.10,0.90), plotEmptyCluster=FALSE, plotOptions=list(), ...){
    t1<-unique(x$Xc[,2]) # this line is not needed
    t1<-x$SUT
    if (model=="mean"){
        meanReg<-array(0,dim=c(x$nSamples,x$H,x$LUT)) 
        uXc<-unique(x$Xc)
        file1 <- paste(x$DIR,"BNSP.eta.ls.txt",sep="")
        file2 <- paste(x$DIR,"BNSP.clusters.txt",sep="")
        if (!file.exists(file1) && !file.exists(file2)) quiet(ls.mvrm(x))
        eta<-array(unlist(read.table("BNSP.eta.ls.txt")),dim=c(x$nSamples,x$H,(x$LGc+1)))
        tabs<-table(array(unlist(read.table("BNSP.clusters.txt"))))
        labs<-array(0,x$H)
        labs[as.numeric(names(tabs))]<-tabs
        for (i in 1:x$nSamples)
            meanReg[i,,] <- eta[i,,]%*%t(uXc)
        if (x$out[[60]]==1) meanReg <- tanh(meanReg)
        centreM<- apply(meanReg,c(2,3),centre)
        QM<-NULL; SQM<-matrix(NA,ncol=2,nrow=NROW(centreM))
        if (!is.null(quantiles)){ 
            QM<-apply(meanReg,c(2,3),quantile,probs=quantiles,na.rm=TRUE)
            SQM<-cbind(QLM=c(t(QM[1,,])),QUM=c(t(QM[2,,])))
	    }
        dataM<-data.frame(group=factor(rep(1:x$H,each=x$LUT)),size=rep(labs,each=x$LUT),t=rep(t1,x$H),centre=c(t(centreM)),QLM=SQM[,1],QUM=SQM[,2])		
        dataM<-dataM[with(dataM, order(-size, t)),]
        #dataM$group<-factor(rep(seq(1,x$H,1),each=x$LUT))# not needed
        labs2<-unique(dataM$size)
        if (!plotEmptyCluster) dataM<-dataM[dataM$size > 0,]
        plotElM<-geom_line(aes_string(x="t",y="centre", group="group",col="group"),alpha=1.0,cex=1.3)
        if (!is.null(quantiles))  
            plotElM<-c(geom_ribbon(aes_string(x="t",ymin="QLM",ymax="QUM",group="group",fill="group",col="group"),alpha =0.2),plotElM)
        ggM<-ggplot(data=dataM)
        plotM <- ggM + plotElM + ylab(expression(mu[t])) + guides(fill=FALSE) + 
                 scale_color_discrete(name = "group [size]", labels = paste("[",labs2,"]",sep="")) + plotOptions  
    }
    if (model=="stdev"){
	    fitV<-matrix(0,nrow=x$nSamples,ncol=length(t1))
	    uZ<-unique(x$Z)
		alphaFN <- file.path(paste(x$DIR,"BNSP.alpha.txt",sep="")) 
        aFile<-file(alphaFN,open="r")
		sigma2FN <- file.path(paste(x$DIR,"BNSP.sigma2.txt",sep="")) 
        s2File<-file(sigma2FN,open="r")
        for (i in 1:x$nSamples){
        	alpha<-scan(aFile,what=numeric(),n=x$LD,quiet=T)
		    s2<-scan(s2File,what=numeric(),n=1,quiet=T)
            fitV[i,]<-sqrt(s2*exp(uZ[,-1]%*%matrix(alpha)))
		}
		close(aFile)
		close(s2File)	
		centreM<-apply(fitV,2,centre)
        QM<-NULL; SQM<-matrix(NA,ncol=2,nrow=NROW(centreM))
        if (!is.null(quantiles)){ 
            QM<-apply(fitV,2,quantile,probs=quantiles,na.rm=TRUE)
            SQM<-cbind(QLM=QM[1,],QUM=QM[2,])
	    }
        dataM<-data.frame(t=t1,centre=centreM,QLM=SQM[,1],QUM=SQM[,2])		
        plotElM<-geom_line(aes_string(x="t",y="centre"),alpha=1.0,cex=1.3)
        if (!is.null(quantiles))  
            plotElM<-c(geom_ribbon(aes_string(x="t",ymin="QLM",ymax="QUM"),alpha =0.2),plotElM)
        ggM<-ggplot(data=dataM)
        plotM <- ggM + plotElM + ylab(expression(sigma[t])) + plotOptions  
	}
	return(plotM)
}

