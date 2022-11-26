# mvrm, continue, mvrm2mcmc, plotCorr, histCorr, predict.mvrm, print.mvrm, summary.mvrm, clustering
mvrm <- function(formula,distribution="normal",
                 data=list(),centre=TRUE,
                 sweeps,burn=0,thin=1,seed,StorageDir,
                 c.betaPrior="IG(0.5, 0.5 * n * p)",                  
                 pi.muPrior="Beta(1, 1)",                  
                 c.alphaPrior="IG(1.1, 1.1)",                                  
                 sigmaPrior="HN(2)",                  
                 pi.sigmaPrior="Beta(1, 1)",
                 c.psiPrior = "HN(1)",
                 phiPrior="HN(2)", 
                 pi.omegaPrior="Beta(1, 1)",                  
                 mu.RPrior="N(0, 1)",
                 sigma.RPrior="HN(1)",
                 corr.Model=c("common",nClust=1),
                 DP.concPrior="Gamma(5, 2)",
                 tuneCbeta,
                 tuneCalpha,
                 tuneAlpha,               
                 tuneSigma2,
                 tuneCpsi,
                 tunePsi,
                 tunePhi,
                 tuneR,
                 tuneSigma2R,                 
                 tuneHar,
                 tau,FT=1,
                 compDeviance=FALSE,...){
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
    if (nSamples<0) stop("problem with sweeps & burn arguments")
	LASTWB<-1	
	#Distribution
    distributions <- c('normal','t')
    fam<-match(distribution,distributions) 
    if (is.na(fam)) stop("unrecognised distribution; should be one of: normal or t")
    # Match call
    call <- match.call(expand.dots = FALSE)
    call2 <- match.call.defaults()
    #Data
    if (!is.list(data) && !is.data.frame(data))
        data <- as.data.frame(data)
    #Formula & data dimensions
    p <- length(as.Formula(formula))[1]
    if (fam==1 && length(as.Formula(formula))[2] > 2) stop("more than two regression models provided")
    if (fam==1 && length(as.Formula(formula))[2] == 1) formula <- as.Formula(formula, ~1)    
    if (fam==2 && length(as.Formula(formula))[2] > 3) stop("more than three regression models provided")    
    if (fam==2 && length(as.Formula(formula))[2] == 1) formula <- as.Formula(formula, ~1, ~1)
    if (fam==2 && length(as.Formula(formula))[2] == 2) stop("ambiguous definition of regression model")        
    formula.m<-formula(as.Formula(formula),lhs=1,rhs=1)
    formula.v<-formula(as.Formula(formula),lhs=0,rhs=2)
    if (fam==2) formula.s<-formula(as.Formula(formula),lhs=0,rhs=3)    
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
    for (i in 1:p) 
        nullDeviance <- nullDeviance - 2 * logLik(lm(Y[,i] ~ 1))
    # Sample size
    n<-length(c(Y))/p
    #X
    XYK<-DM(formula=formula.m,data=data,n=n,centre=centre)
    X<-as.matrix(XYK$X)
    Xknots<-XYK$Rknots
    storeMeanVectorX<-XYK$meanVector
    storeIndicatorX<-XYK$indicator
    LG<-NCOL(X)-1
    labelsX<-XYK$labels
    countX<-XYK$count
    varsX<-XYK$vars
    which.SpecX<-XYK$which.Spec
    formula.termsX<-XYK$formula.terms
    is.Dx<-XYK$is.D
    assignX<-XYK$assign
    NG<-max(assignX)
    assignX3<-XYK$assign3
    vecLG<-table(assignX3)[-1]    
    cusumVecLG<-c(0,cumsum(vecLG))    
    NG2<-max(assignX3)
    repsX <- XYK$repsX
    DynamicSinPar <-XYK$DSP
    amplitude <- DynamicSinPar[1]
    nHar <- DynamicSinPar[3]
    Dynamic <- XYK$Dynamic
    isSin <- XYK$isSin
    varSin <- varsX[[min(which(isSin==1))]]
    if (amplitude==1) isSin <- isSin[repsX]
    #Z
    ZK<-DM(formula=formula.v,data=data,n=n,centre=centre)
    Z<-as.matrix(ZK$X)
    Zknots<-ZK$Rknots
    storeMeanVectorZ<-ZK$meanVector
    storeIndicatorZ<-ZK$indicator
    LD<-NCOL(Z)-1
    labelsZ<-ZK$labels
    countZ<-ZK$count
    varsZ<-ZK$vars
    which.SpecZ<-ZK$which.Spec
    formula.termsZ<-ZK$formula.terms
    is.Dz<-ZK$is.D    
    assignZ<-ZK$assign
    #vecLD<-table(assignZ)[-1]
    #cusumVecLD<-c(0,cumsum(vecLD))
    ND<-max(assignZ)#length(vecLD)
    repeats<-rep(1,ND)
    repeats[which(is.Dz==1)]<-table(assignZ)[which(is.Dz==1)+1]
    oneHE<-1; if (ND > 0) oneHE<-rep(c(1:ND),repeats)
    assignZ2<-ZK$assign2        
    vecLD<-table(assignZ2)[-1]
    cusumVecLD<-c(0,cumsum(vecLD))
    ND2<-max(assignZ2)#length(vecLD) 
    isDz2<-as.numeric(is.Dz)[oneHE]
    if (is.na(isDz2[1])) isDz2<-1
    #W
    W<-Wknots<-storeMeanVectorW<-storeIndicatorW<-assignW<-labelsW<-countW<-varsW<-is.Dw<-which.SpecW<-formula.termsW<-NULL
    LC<-vecLC<-NC<-cusumVecLC<-0
    if (fam==2){
        WK<-DM(formula=formula.s,data=data,n=n,centre=centre)
        W<-as.matrix(WK$X)
        Wknots<-WK$Rknots
        storeMeanVectorW<-WK$meanVector
        storeIndicatorW<-WK$indicator
        LC<-NCOL(W)-1
        vecLC<-table(WK$assign)[-1]
        NC<-length(vecLC)
        cusumVecLC<-c(0,cumsum(vecLC))
        assignW<-WK$assign
        labelsW<-WK$labels
        countW<-WK$count
        varsW<-WK$vars
        is.Dw<-WK$is.D
        which.SpecW<-WK$which.Spec
        formula.termsW<-WK$formula.terms    
    }
    MVLD <- 1
    if (LD > 0 || LC > 0) MVLD<-max(vecLD,vecLC)    
    #Initialize covariance & correlation matrix
    LASTDE <- LASTR <- LASTmuR <- LASTsigma2R <- 1       
    Res<-NULL
    find.sm <- grep("sm(",colnames(X),fixed=TRUE)
    Xmain <- X
    if (length(find.sm) > 0) Xmain<-X[,-find.sm]
    for (i in 1:p)
        Res<-cbind(Res,residuals(lm(Y[,i] ~ Xmain - 1)))
    CR<-0.9*cov(Res)+0.1*diag(mean(eigen(cov(Res))$values),p)
    D<-matrix(0,p,p)
    diag(D)<-sqrt(diag(CR))
    LASTsigma2zk<-diag(CR)
    LASTDE<-c(c(D),c(CR))
    LASTR<-cov2cor(CR)
    if (p > 1) LASTmuR<-mean(LASTR[upper.tri(LASTR)])
    if (p > 2) LASTsigma2R<-var(LASTR[upper.tri(LASTR)])
    #Prior for pi.mu
    if (!length(pi.muPrior)==1 && !length(pi.muPrior)==NG && !length(pi.muPrior)==(p*NG))
        stop("pi.muPrior has incorrect dimension")
    if (length(pi.muPrior)==NG && NG2!=NG) pi.muPrior <- pi.muPrior[repsX]    
    if (length(pi.muPrior)==(p*NG) && NG2!=NG){ 
	    for (i in 1:(p-1)) repsX<-c(repsX,repsX+max(repsX))	    
	    pi.muPrior <- pi.muPrior[repsX]
	}    
    pimu<-NULL
    for (k in 1:length(pi.muPrior)){
        sp<-strsplit(pi.muPrior[k],"Beta\\(")
        sp<-strsplit(sp[[1]][2],"\\)")
        sp<-strsplit(sp[[1]][1],",")
        pimu<-c(pimu,as.numeric(sp[[1]]))
    }
    if (length(pi.muPrior)==1) pimu<-rep(pimu,p*NG2)
    if (length(pi.muPrior)==NG2) pimu<-rep(pimu,p)
    piHar<-1 #pi of the harminocs for forming the dynamic matrix
    if (sum(isSin) && amplitude > 1){ 
        loc <- 1+2*(which(isSin==1)-1)
        loc <-c(loc,loc+1)
        expand.loc <- rep(loc,p) + rep(seq(from=0,length.out=p,by=2*NG),each=2)
        piHar <- pimu[loc]  
	}    
    #Prior for c.beta
    sp<-strsplit(c.betaPrior,"IG\\(")
    sp<-strsplit(sp[[1]][2],"\\)")
    sp<-strsplit(sp[[1]][1],",")
    cetaParams<-c(as.numeric(sp[[1]][1]),eval(parse(text=sp[[1]][2])))
    #Prior for pi.sigma
    if (ND > 0){
        if (!length(pi.sigmaPrior)==1 && !length(pi.sigmaPrior)==ND && !length(pi.sigmaPrior)==(p*ND))
            stop("pi.sigmaPrior has incorrect dimension")
        if (length(pi.sigmaPrior)==(p*ND)) 
            pi.sigmaPrior<-pi.sigmaPrior[oneHE+rep(seq(from=0,to=p*ND-1,by=ND),each=length(oneHE))]
        if (length(pi.sigmaPrior)==ND) pi.sigmaPrior<-pi.sigmaPrior[oneHE]    
        pisigma<-NULL
        for (k in 1:length(pi.sigmaPrior)){
            sp<-strsplit(pi.sigmaPrior[k],"Beta\\(")
            sp<-strsplit(sp[[1]][2],"\\)")
            sp<-strsplit(sp[[1]][1],",")
            pisigma<-c(pisigma,as.numeric(sp[[1]]))
        }
        if (length(pi.sigmaPrior)==1) pisigma<-rep(pisigma,p*ND2)
        if (length(pi.sigmaPrior)==ND2) pisigma<-rep(pisigma,p)
	}else{pisigma<-1}
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
	#Prior for pi.omega
    if (!length(pi.omegaPrior)==1 && !length(pi.omegaPrior)==NC && !length(pi.omegaPrior)==(p*NC))
        stop("pi.omegaPrior has incorrect dimension")
    piomega<-NULL
    for (k in 1:length(pi.omegaPrior)){
        sp<-strsplit(pi.omegaPrior[k],"Beta\\(")
        sp<-strsplit(sp[[1]][2],"\\)")
        sp<-strsplit(sp[[1]][1],",")
        piomega<-c(piomega,as.numeric(sp[[1]]))
    }
    if (length(pi.omegaPrior)==1) piomega<-rep(piomega,p*NC)
    if (length(pi.omegaPrior)==NC) piomega<-rep(piomega,p)
    #Prior for c.psi
    if (!length(c.psiPrior)==1 && !length(c.psiPrior)==p)
        stop("c.psiPrior has incorrect dimension")
    specials<-c("HN","IG")
    cpsiParams<-NULL
    HNcp<-vector()
    for (k in 1:length(c.psiPrior)){
        sp<-strsplit(c.psiPrior[k],"\\(")
        if (sp[[1]][1] %in% specials){
            if (match(sp[[1]][1],specials)==1) HNcp[k]<-1
            if (match(sp[[1]][1],specials)==2) HNcp[k]<-0
        } else stop("unrecognised prior for c.psi")
        sp<-strsplit(sp[[1]][2],"\\)")
        sp<-strsplit(sp[[1]][1],",")
        cpsiParams<-c(cpsiParams,as.numeric(sp[[1]]))
    }
    if (length(c.psiPrior)==1){
        cpsiParams<-rep(cpsiParams,p)
        HNcp<-rep(HNcp,p)
	}	
    #Prior for phi2_k 
    if (!length(phiPrior)==1 && !length(phiPrior)==p)
        stop("phiPrior has incorrect dimension")
    specials<-c("HN","IG")
    phiParams<-NULL
    HNphi<-vector()
    for (k in 1:length(phiPrior)){
        sp<-strsplit(phiPrior[k],"\\(")
        if (sp[[1]][1] %in% specials){
            if (match(sp[[1]][1],specials)==1) HNphi[k]<-1
            if (match(sp[[1]][1],specials)==2) HNphi[k]<-0
        } else stop("unrecognised prior for phi")
        sp<-strsplit(sp[[1]][2],"\\)")
        sp<-strsplit(sp[[1]][1],",")
        phiParams<-c(phiParams,as.numeric(sp[[1]]))
    }
    if (length(phiPrior)==1){
        phiParams<-rep(phiParams,p)
        HNphi<-rep(HNphi,p)
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
    corModels<-c("common","groupC","groupV","uni")
    mcm<-match(corr.Model[1],corModels)
    if (p==1 && fam==1) mcm=4
    if (p==1 && fam==2) mcm=1
    if (p==2) mcm=1
    if (is.na(mcm)) stop("unrecognised correlation model")
    H <- G <- 1
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
	if (fam==2) mcm <- mcm + 7
    #Prior for alpha DP
    sp<-strsplit(DP.concPrior,"Gamma\\(")
    sp<-strsplit(sp[[1]][2],"\\)")
    sp<-strsplit(sp[[1]][1],",")
    DPparams<-as.numeric(sp[[1]])
    DPparams <- c(DPparams,0.01)
    #Seed
    if (missing(seed)) seed<-as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31)
    # Storage directory & files
    if (!missing(StorageDir)){
        StorageDir <- path.expand(StorageDir)
        ncharwd <- nchar(StorageDir)}
    if (!missing(StorageDir)) if (!(substr(StorageDir,ncharwd,ncharwd)=="/")) StorageDir <- paste(StorageDir,"/",sep="")
    if (!missing(StorageDir)) if (!dir.exists(StorageDir)) dir.create(StorageDir, recursive = TRUE)
    if (missing(StorageDir)) stop("provide a storage directory via argument StorageDir")
    FL <- c("gamma", "cbeta", "delta", "alpha", "R", "muR", "sigma2R", "calpha", "sigma2", "beta", 
            "ksi", "psi", "cpsi", "phi2",
            "compAlloc", "nmembers", "deviance", "DPconc",
            "compAllocV", "nmembersV",
            "DE",
            "nu", "fi", "omega", "ceta","comega","eta", "test", "nu.ls", "eta.ls", "nmembers.ls", "clusters",  
            "probs", "tune")
    for (i in 1:length(FL)){
        oneFile <- paste(StorageDir, paste("BNSP",FL[i], "txt",sep="."),sep="/")
        if (file.exists(oneFile)) file.remove(oneFile)
	}
    #Tuning Parameters
	if (missing(tuneCbeta)) tuneCbeta<-20
    if (missing(tuneCalpha)) tuneCalpha<-rep(1,p)
    if (!length(tuneCalpha)==p) tuneCalpha<-rep(mean(tuneCalpha),p)    
    if (ND > 0){
        if (missing(tuneAlpha) || !length(tuneAlpha)==ND || !length(tuneAlpha)==(ND*p)){
		    tuneAlpha<-rep(5,ND)
		    tuneAlpha[which(is.Dz==1)]<-1
	    }    
	    if (length(tuneAlpha)==(ND*p)) 
            tuneAlpha<-tuneAlpha[oneHE+rep(seq(from=0,to=p*ND-1,by=ND),each=length(oneHE))]                
        if (length(tuneAlpha)==ND) tuneAlpha<-tuneAlpha[oneHE]
        if (length(tuneAlpha)==ND2) tuneAlpha<-rep(tuneAlpha,p)
    }else{tuneAlpha<-1}
    if (missing(tuneSigma2)) tuneSigma2<-rep(0.5,p)
    if (!length(tuneSigma2)==p) tuneSigma2<-rep(mean(tuneSigma2),p)
    if (missing(tuneCpsi)) tuneCpsi<-rep(1,p)
    if (!length(tuneCpsi)==p) tuneCpsi<-rep(mean(tuneCpsi),p)
	if (missing(tunePsi)) tunePsi<-rep(5,NC*p)
    if (!length(tunePsi)==(NC*p)) tunePsi<-rep(mean(tunePsi),NC*p)
    if (missing(tunePhi)) tunePhi<-rep(0.5,p)
    if (!length(tunePhi)==p) tunePhi<-rep(mean(tunePhi),p)
    if (missing(tuneR)) tuneR<-40*(p+2)^3
    tuneR[which(tuneR<p+2)]<-p+2
    if (missing(tuneSigma2R)) tuneSigma2R<-0.25
    if (nHar > 0){
        ifelse(missing(tuneHar), tuneHar<-rep(100,nHar), tuneHar<-rep(mean(tuneHar),nHar))
	}else{tuneHar<-1}
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
    #if (missing(blockSizeProbC)){
	blockSizeProbC <- rep(0,LC)
    blockSizeProbC[1:5]<-c(10,25,30,25,10)
    #}
    maxBSG <- max(which(blockSizeProbG>0))
    maxBSD <- max(which(blockSizeProbD>0))
    maxBSC <- max(which(blockSizeProbC>0))
    maxBSGDC <- c(maxBSG,maxBSD,maxBSC)
    #Deviance
    deviance <- c(0,0)
    #Cont
    cont <- 0    
    #Call C
    if (mcm==4) {out<-.C("mvrmC",
        as.integer(seed),as.character(StorageDir),
        as.integer(sweeps),as.integer(burn),as.integer(thin),
        as.double(c(t(Y))),as.double(X),as.double(Z),as.integer(n),as.integer(LG),as.integer(LD),
        as.double(blockSizeProbG),as.integer(maxBSG),as.double(blockSizeProbD),as.integer(maxBSD),
        as.double(tuneCalpha),as.double(tuneSigma2),as.double(tuneCbeta),as.double(tuneAlpha),
        as.integer(NG2),as.integer(ND2),as.integer(vecLG),as.integer(vecLD),
        as.integer(cusumVecLG),as.integer(cusumVecLD),as.integer(MVLD),
        as.double(cetaParams),as.integer(HNca),as.double(calphaParams),as.double(pimu),as.double(pisigma),
        as.integer(HNsg),as.double(sigmaParams),as.double(deviance),as.integer(isDz2),
        as.integer(c(cont)),as.integer(c(0)),as.integer(c(0)),as.double(c(0)),as.double(LASTsigma2zk),
        as.double(c(0)),as.double(c(0)),as.integer(LASTWB),
        as.integer(isSin),as.double(Dynamic),as.integer(DynamicSinPar),as.double(tuneHar),as.double(piHar))}        
    if (mcm==1) {out<-.C("mult",
        as.integer(seed),as.character(StorageDir),
        as.integer(sweeps),as.integer(burn),as.integer(thin),
        as.integer(n),as.integer(c(p,mcm)),
        as.double(c(t(Y))),as.double(t(X)),as.double(Z),
        as.integer(LG),as.integer(LD),
        as.double(blockSizeProbG),as.integer(maxBSG),as.double(blockSizeProbD),as.integer(maxBSD),
        as.integer(NG2),as.integer(ND2),as.integer(vecLG),as.integer(vecLD),
        as.integer(cusumVecLG),as.integer(cusumVecLD),as.integer(MVLD),
        as.double(tuneCalpha),as.double(tuneSigma2),as.double(tuneCbeta),as.double(tuneAlpha),
        as.double(tuneSigma2R),as.double(tuneR),
        as.double(pimu),as.double(cetaParams),as.double(pisigma),
        as.integer(HNca),as.double(calphaParams),as.double(Rparams),
        as.integer(HNsg),as.double(sigmaParams),as.double(tau),as.integer(FT),as.double(deviance),
        as.integer(isDz2),
        as.integer(c(cont)),as.integer(c(0)),as.integer(c(0)),as.double(c(0)),as.double(LASTsigma2zk),
        as.double(c(LASTR)),
        as.double(c(LASTmuR)),as.double(c(LASTsigma2R)),as.double(c(0)),as.double(c(0)),
        as.integer(c(LASTsw)),as.double(c(LASTDE)),as.integer(LASTWB),as.integer(isSin))}
    if (mcm==2) {out<-.C("multg",
        as.integer(seed),as.character(StorageDir),
        as.integer(sweeps),as.integer(burn),as.integer(thin),
        as.integer(n),as.integer(c(p,mcm)),
        as.double(c(t(Y))),as.double(t(X)),as.double(Z),
        as.integer(LG),as.integer(LD),
        as.double(blockSizeProbG),as.integer(maxBSG),as.double(blockSizeProbD),as.integer(maxBSD),
        as.integer(NG2),as.integer(ND2),as.integer(vecLG),as.integer(vecLD),
        as.integer(cusumVecLG),as.integer(cusumVecLD),as.integer(MVLD),
        as.double(tuneCalpha),as.double(tuneSigma2),as.double(tuneCbeta),as.double(tuneAlpha),
        as.double(tuneSigma2R),as.double(tuneR), as.double(pimu),as.double(cetaParams),as.double(pisigma),
        as.integer(HNca),as.double(calphaParams),as.double(Rparams),
        as.integer(HNsg),as.double(sigmaParams),as.double(tau),as.integer(FT),as.double(deviance),as.integer(H),
        as.double(DPparams), as.integer(isDz2),
        as.integer(c(cont)),as.integer(c(0)),as.integer(c(0)),as.double(c(0)),as.double(LASTsigma2zk),
        as.double(LASTR),as.double((rnorm(n=H,mean=LASTmuR,sd=0.01))),as.double(LASTsigma2R/H),
        as.double(c(0)), as.double(c(0)),as.integer(c(0)),as.double(c(0)),
        as.integer(c(LASTsw)),as.double(c(LASTDE)),as.integer(LASTWB),as.integer(isSin))}
    if (mcm==3) {out<-.C("multgv",
        as.integer(seed),as.character(StorageDir),
        as.integer(sweeps),as.integer(burn),as.integer(thin),
        as.integer(n),as.integer(c(p,mcm)),
        as.double(c(t(Y))),as.double(t(X)),as.double(Z),
        as.integer(LG),as.integer(LD),
        as.double(blockSizeProbG),as.integer(maxBSG),as.double(blockSizeProbD),as.integer(maxBSD),
        as.integer(NG2),as.integer(ND2),as.integer(vecLG),as.integer(vecLD),
        as.integer(cusumVecLG),as.integer(cusumVecLD),as.integer(MVLD),
        as.double(tuneCalpha),as.double(tuneSigma2),as.double(tuneCbeta),as.double(tuneAlpha),
        as.double(tuneSigma2R),as.double(tuneR), as.double(pimu),as.double(cetaParams),as.double(pisigma),
        as.integer(HNca),as.double(calphaParams),as.double(Rparams),
        as.integer(HNsg),as.double(sigmaParams),as.double(tau),as.integer(FT),as.double(deviance),as.integer(G),
        as.double(DPparams), as.integer(isDz2), 
        as.integer(c(cont)),as.integer(c(0)),as.integer(c(0)),as.double(c(0)),as.double(LASTsigma2zk),
        as.double(LASTR),as.double((rnorm(n=H,mean=LASTmuR,sd=0.01))),as.double(LASTsigma2R/H),
        as.double(c(0)),
        as.double(c(0)),as.integer(c(0)),as.double(c(0)),
        as.integer(c(LASTsw)),as.double(c(LASTDE)),as.integer(LASTWB),as.integer(isSin))}
     if (mcm==8) {out<-.C("multT",
        as.integer(seed),as.character(StorageDir),
        as.integer(c(sweeps,burn,thin,n,p,mcm,MVLD,compDeviance)),
        as.double(c(t(Y))),as.double(t(X)),as.double(Z),as.double(W),
        as.integer(c(LG,LD,LC)),
        as.double(blockSizeProbG),as.double(blockSizeProbD),as.double(blockSizeProbC),
        as.integer(maxBSGDC),as.integer(c(NG2,ND2,NC)),as.integer(vecLG),as.integer(vecLD),as.integer(vecLC), 
        as.integer(cusumVecLG),as.integer(cusumVecLD),as.integer(cusumVecLC),   
        as.double(tuneCalpha),as.double(tuneSigma2),as.double(tuneCbeta),as.double(tuneAlpha),
        as.double(tuneSigma2R),as.double(tuneR),as.double(tunePsi),as.double(tunePhi),as.double(tuneCpsi),        
        as.double(pimu),as.double(cetaParams),as.double(pisigma),
        as.integer(HNca),as.double(calphaParams),as.integer(HNcp),as.double(cpsiParams),
        as.integer(HNphi),as.double(phiParams),as.double(Rparams),
        as.integer(HNsg),as.double(sigmaParams),as.double(piomega),
        as.double(tau),as.integer(FT),as.double(deviance),as.double(Res),
        as.integer(isDz2),        
        as.integer(c(cont)),as.integer(c(0)),as.integer(c(0)),as.double(c(0)),
        as.double(c(LASTsigma2zk)),as.double(c(LASTR)),as.double(c(LASTDE)), 
        as.double(c(LASTmuR)),as.double(c(LASTsigma2R)),
        as.integer(c(LASTsw)),as.integer(LASTWB),
        as.double(c(0)),as.double(c(0)),
        as.integer(c(0)),as.double(c(0,0,0)),as.integer(isSin))}    
    if (mcm==9) {out<-.C("multgT",
        as.integer(seed),as.character(StorageDir),
        as.integer(c(sweeps,burn,thin,n,p,mcm,MVLD,compDeviance)),
        as.double(c(t(Y))),as.double(t(X)),as.double(Z),as.double(W),
        as.integer(c(LG,LD,LC)),
        as.double(blockSizeProbG),as.double(blockSizeProbD),as.double(blockSizeProbC),
        as.integer(maxBSGDC),as.integer(c(NG,ND2,NC)),as.integer(vecLG),as.integer(vecLD),as.integer(vecLC),                
        as.integer(cusumVecLG),as.integer(cusumVecLD),as.integer(cusumVecLC),
        as.double(tuneCalpha),as.double(tuneSigma2),as.double(tuneCbeta),as.double(tuneAlpha),
        as.double(tuneSigma2R),as.double(tuneR),as.double(tunePsi),as.double(tunePhi), as.double(tuneCpsi),        
        as.double(pimu),as.double(cetaParams),as.double(pisigma),
        as.integer(HNca),as.double(calphaParams),as.integer(HNcp),as.double(cpsiParams),
        as.integer(HNphi),as.double(phiParams),as.double(Rparams),
        as.integer(HNsg),as.double(sigmaParams),as.double(piomega),        
        as.double(tau),as.integer(FT),as.double(deviance),as.double(Res),                
        as.integer(H),as.double(DPparams),        
        as.integer(isDz2),
        as.integer(c(cont)),as.integer(c(0)),as.integer(c(0)),as.double(c(0)),        
        as.double(c(LASTsigma2zk)),as.double(c(LASTR)),as.double(c(LASTDE)),                        
        as.double(rnorm(n=H,mean=LASTmuR,sd=0.01),LASTsigma2R/H),        
        as.integer(c(LASTsw)),as.integer(LASTWB),                
        as.double(c(0)),as.double(c(0)),
        as.integer(c(0)),as.double(c(0)),
        as.integer(c(0)),as.double(c(0,0,0)),as.integer(isSin))}    
    if (mcm==10) {out<-.C("multgvT",
        as.integer(seed),as.character(StorageDir),
        as.integer(c(sweeps,burn,thin,n,p,mcm,MVLD,compDeviance)),
        as.double(c(t(Y))),as.double(t(X)),as.double(Z),as.double(W),
        as.integer(c(LG,LD,LC)),        
        as.double(blockSizeProbG),as.double(blockSizeProbD),as.double(blockSizeProbC),
        as.integer(maxBSGDC),as.integer(c(NG,ND2,NC)),as.integer(vecLG),as.integer(vecLD),as.integer(vecLC),          
        as.integer(cusumVecLG),as.integer(cusumVecLD),as.integer(cusumVecLC),
        as.double(tuneCalpha),as.double(tuneSigma2),as.double(tuneCbeta),as.double(tuneAlpha),
        as.double(tuneSigma2R),as.double(tuneR),as.double(tunePsi),as.double(tunePhi), as.double(tuneCpsi),        
        as.double(pimu),as.double(cetaParams),as.double(pisigma),as.integer(HNca),as.double(calphaParams),        
        as.integer(HNcp),as.double(cpsiParams),as.integer(HNphi),as.double(phiParams),        
        as.double(Rparams),as.integer(HNsg),as.double(sigmaParams),as.double(piomega),        
        as.double(tau),as.integer(FT),as.double(deviance),as.double(Res),        
        as.integer(G), as.double(DPparams),        
        as.integer(isDz2),        
        as.integer(c(cont)),as.integer(c(0)),as.integer(c(0)),as.double(c(0)),
        as.double(c(LASTsigma2zk)),as.double(c(LASTR)),as.double(c(LASTDE)),        
        as.double(rnorm(n=H,mean=LASTmuR,sd=0.01),LASTsigma2R/H),
        as.integer(c(LASTsw)),as.integer(LASTWB),
        as.double(c(0)),as.double(c(0)),        
        as.integer(c(0)),as.double(c(0)),
        as.integer(c(0)),as.double(c(0,0,0)),as.integer(isSin))}    
    #Output
    if (mcm<4){
        loc1<-24
        loc2<-40 
        tuneSigma2Ra<-out[[loc1+4]][1]
        tuneRa<-out[[loc1+5]][1] 
        tunePsia<-tunePsi
        tunePhia<-tunePhi
        tuneCpsia<-tuneCpsi 
        DevCalcs<-2  
    }
    if (mcm==4){
		loc1<-16 
		loc2<-34
		tuneSigma2Ra<-tuneSigma2R 
		tuneRa<-tuneR
		tunePsia<-tunePsi
        tunePhia<-tunePhi
        tuneCpsia<-tuneCpsi
		DevCalcs<-2
    }  
    if (mcm>7){
		loc1<-20
		loc2<-44
		tuneSigma2Ra<-out[[loc1+4]][1]
        tuneRa<-out[[loc1+5]][1]
        tunePsia<-out[[loc1+6]][1:(p*NC)]
        tunePhia<-out[[loc1+7]][1:p]
        tuneCpsia<-out[[loc1+8]][1:p]
		DevCalcs<-1
	}      
    fit <- list(call=call,call2=call2,formula=formula,seed=seed,p=p,d=p*(p-1)/2,
                data=data,Y=Y,
                X=X,Xknots=Xknots,LG=LG,NG=NG,NG2=NG2,
                Z=Z,Zknots=Zknots,LD=LD,ND=ND,ND2=ND2,
                W=W,Wknots=Wknots,LC=LC,NC=NC,
                storeMeanVectorX=storeMeanVectorX,
                storeMeanVectorZ=storeMeanVectorZ,
                storeMeanVectorW=storeMeanVectorW,
                storeIndicatorX=storeIndicatorX,
                storeIndicatorZ=storeIndicatorZ,
                storeIndicatorW=storeIndicatorW,
                assignX=assignX,
                assignZ=assignZ,
                assignW=assignW,
                labelsX=labelsX,
                labelsZ=labelsZ,
                labelsW=labelsW,
                countX=countX,
                countZ=countZ,
                countW=countW,
                varsY=varsY,
                varsX=varsX,
                varsZ=varsZ,
                varsW=varsW,
                is.Dx=is.Dx,
                is.Dz=is.Dz,
                is.Dw=is.Dw,
                which.SpecX=which.SpecX,
                which.SpecZ=which.SpecZ,
                which.SpecW=which.SpecW,
                formula.termsX=formula.termsX,
                formula.termsZ=formula.termsZ,
                formula.termsW=formula.termsW,
                nSamples=nSamples,
                totalSweeps=sweeps,
                mcpar=c(as.integer(burn+1),as.integer(seq(from=burn+1,by=thin,length.out=nSamples)[nSamples]),as.integer(thin)),
                mcm=mcm,H=H,G=G,
                tuneCalpha=c(tuneCalpha,out[[loc1+0]][1:p]),            
                tuneSigma2=c(tuneSigma2,out[[loc1+1]][1:p]),
                tuneCbeta=c(tuneCbeta,out[[loc1+2]][1]),
                tuneAlpha=c(tuneAlpha,out[[loc1+3]][1:(p*ND2)]),
                tuneSigma2R=c(tuneSigma2R,tuneSigma2Ra),
                tuneR=c(tuneR,tuneRa),
                tunePsi=c(tunePsi,tunePsia),                
                tunePhi=c(tunePhi,tunePhia),
                tuneCpsi=c(tuneCpsi,tuneCpsia),
                tuneHar=c(tuneHar,out[[47]][1:nHar]),
                deviance=c(out[[loc2]][1:DevCalcs]),
                nullDeviance=nullDeviance,                            
                DIR=StorageDir,
                out=out,
                LUT=1, SUT=1, LGc=0, LDc=0, NGc=0, NDc=0, NK=0, FT=FT,
                qCont=0,
                HNca=HNca,HNsg=HNsg,HNcp=HNcp,HNphi=HNphi,nHar=nHar,varSin=varSin)
    class(fit) <- 'mvrm'
    return(fit)
}

continue <- function(object,sweeps,burn=0,thin,discard=FALSE,...){
	if (object$qCont==1) stop("current object has been continued; there is a more recent one")
	eval.parent(substitute(object$qCont<-1))
	cont<-1
	#Sweeps
	if (missing(sweeps)) stop("provide sweeps argument")
    sweeps <- as.integer(sweeps)
    burn <- as.integer(burn)
    if (missing(thin) || !discard) thin <- object$mcpar[3]
    if (thin <= 0) thin <- 1 
    thin <- as.integer(thin)    
    nSamples <- 0
    LASTsw<-0
    if (sweeps > 0 && (sweeps-burn) > 0){
		nSamples <- length(seq(1,(sweeps-burn),by=thin))
		LASTsw<-seq(burn,by=thin,length.out=nSamples)[nSamples]
	}
    if (nSamples<=0) stop("problem with sweeps & burn arguments")
    LASTWB <- floor(object$totalSweeps/50)+1
    #Files
    FL <- c("gamma", "cbeta", "delta", "alpha", "R", "muR", "sigma2R", "calpha", "sigma2", "beta", #10  
            "compAlloc", "nmembers", "deviance", "DPconc", #4
            "compAllocV","nmembersV", 
            "DE", #mcm 1,2,3,4
            "psi", "ksi", "cpsi", "nu", "fi", "omega", "ceta", "comega", "eta","tune","probs", #mcm 5,6,7            
            "phi2", #mcm 8,9,10     
            "test")            
    gamma <- paste(object$DIR, paste("BNSP",FL[1], "txt",sep="."),sep="/")
    gamma <- scan(gamma,what=numeric(),n=object$p*object$LG,quiet=TRUE,skip=object$nSamples-1)
    cbeta <- paste(object$DIR, paste("BNSP",FL[2], "txt",sep="."),sep="/")
    cbeta <- scan(cbeta,what=numeric(),n=1,quiet=TRUE,skip=object$nSamples-1)
    delta <- paste(object$DIR, paste("BNSP",FL[3], "txt",sep="."),sep="/")
    delta <- scan(delta,what=numeric(),n=object$p*object$LD,quiet=TRUE,skip=object$nSamples-1)
    alpha <- paste(object$DIR, paste("BNSP",FL[4], "txt",sep="."),sep="/")
    alpha <- scan(alpha,what=numeric(),n=object$p*object$LD,quiet=TRUE,skip=object$nSamples-1)
    R <- paste(object$DIR, paste("BNSP",FL[5], "txt",sep="."),sep="/")
    if (file.exists(R)) R <- scan(R,what=numeric(),n=object$p^2*object$LUT,quiet=TRUE,skip=object$nSamples-1)
    muR <- paste(object$DIR, paste("BNSP",FL[6], "txt",sep="."),sep="/")
    if (file.exists(muR)) muR <- scan(muR,what=numeric(),n=object$H,quiet=TRUE,skip=object$nSamples-1)
    sigma2R <- paste(object$DIR, paste("BNSP",FL[7], "txt",sep="."),sep="/")
    if (file.exists(sigma2R)) sigma2R <- scan(sigma2R,what=numeric(),n=1,quiet=TRUE,skip=object$nSamples-1)
    calpha <- paste(object$DIR, paste("BNSP",FL[8], "txt",sep="."),sep="/")
    calpha <- scan(calpha,what=numeric(),n=object$p,quiet=TRUE,skip=object$nSamples-1)
    sigma2 <- paste(object$DIR, paste("BNSP",FL[9], "txt",sep="."),sep="/")
    sigma2 <- scan(sigma2,what=numeric(),n=object$p,quiet=TRUE,skip=object$nSamples-1)    
    compAlloc <- paste(object$DIR, paste("BNSP",FL[11], "txt",sep="."),sep="/")
    if (file.exists(compAlloc)) compAlloc <- scan(compAlloc,what=numeric(),n=object$d,quiet=TRUE,skip=object$nSamples-1)
    DPconc <- paste(object$DIR, paste("BNSP",FL[14], "txt",sep="."),sep="/")
    if (file.exists(DPconc)) DPconc <- scan(DPconc,what=numeric(),n=1,quiet=TRUE,skip=object$nSamples-1)
    compAllocV <- paste(object$DIR, paste("BNSP",FL[15], "txt",sep="."),sep="/")
    if (file.exists(compAllocV)) compAllocV <- scan(compAllocV,what=numeric(),n=object$p,quiet=TRUE,skip=object$nSamples-1)
    DE <- paste(object$DIR, paste("BNSP",FL[17], "txt",sep="."),sep="/")
    if (file.exists(DE)) DE <- scan(DE,what=numeric(),n=2*object$p^2*object$LUT,quiet=TRUE,skip=0)
    psi <- paste(object$DIR, paste("BNSP",FL[18], "txt",sep="."),sep="/")
    if (file.exists(psi)) psi <- scan(psi,what=numeric(),n=object$p*object$p*object$LK,quiet=TRUE,skip=object$nSamples-1)                       
    ksi <- paste(object$DIR, paste("BNSP",FL[19], "txt",sep="."),sep="/")
    if (file.exists(ksi)) ksi <- scan(ksi,what=numeric(),n=object$p*object$p*object$LK,quiet=TRUE,skip=object$nSamples-1)    
    cpsi <- paste(object$DIR, paste("BNSP",FL[20], "txt",sep="."),sep="/")
    if (file.exists(cpsi)) cpsi <- scan(cpsi,what=numeric(),n=object$p^2,quiet=TRUE,skip=object$nSamples-1)        
    gammaCor <- paste(object$DIR, paste("BNSP",FL[21], "txt",sep="."),sep="/")
    gammaCor <- if (file.exists(gammaCor)) scan(gammaCor,what=numeric(),n=object$LGc*object$H,quiet=TRUE,skip=object$nSamples-1)  
    deltaCor <- paste(object$DIR, paste("BNSP",FL[22], "txt",sep="."),sep="/")
    deltaCor <- if (file.exists(deltaCor)) scan(deltaCor,what=numeric(),n=object$LDc,quiet=TRUE,skip=object$nSamples-1)  
    omega <- paste(object$DIR, paste("BNSP",FL[23], "txt",sep="."),sep="/")
    omega <- if (file.exists(omega)) scan(omega,what=numeric(),n=object$LDc,quiet=TRUE,skip=object$nSamples-1)
    ceta <- paste(object$DIR, paste("BNSP",FL[24], "txt",sep="."),sep="/")
    ceta <- if (file.exists(ceta)) scan(ceta,what=numeric(),n=1,quiet=TRUE,skip=object$nSamples-1)    
    comega <- paste(object$DIR, paste("BNSP",FL[25], "txt",sep="."),sep="/")
    comega <- if (file.exists(comega)) scan(comega,what=numeric(),n=1,quiet=TRUE,skip=object$nSamples-1)
    LASTAll<-c(R,DE,gamma,delta,alpha,sigma2,ksi,psi,cbeta,calpha,cpsi,gammaCor)
    if (object$mcm==6) LASTAll<-c(LASTAll,compAlloc,DPconc)
    if (object$mcm==7) LASTAll<-c(LASTAll,compAllocV,DPconc)
    LASTAll<-c(LASTAll,deltaCor,comega,sigma2R,omega,ceta)
    #LASTAll<-c(LASTAll,deltaCor)
    #LASTAll<-c(LASTAll,comega)
    #LASTAll<-c(LASTAll,sigma2R)
    #LASTAll<-c(LASTAll,omega)
    #LASTAll<-c(LASTAll,ceta)
    #Discard
    if (discard==TRUE){
        for (i in 1:length(FL)){
            oneFile <- paste(object$DIR, paste("BNSP",FL[i], "txt",sep="."),sep="/")
            if (file.exists(oneFile)) file.remove(oneFile)
	    }
	}
    #Deviance 
    deviance<-c(0,0)
    if (discard==FALSE) deviance<-as.double(object$deviance)    
    #Tuning
    tuneCalpha<-object$tuneCalpha[(object$p+1):(2*object$p)]
    tuneSigma2<-object$tuneSigma2[(object$p+1):(2*object$p)]    
    tuneCbeta<-object$tuneCbeta[2]
    tuneAlpha<-object$tuneAlpha[(object$p*object$ND2+1):(2*object$p*object$ND2)]
    tuneSigma2R<-object$tuneSigma2R[2]
    tuneR<-object$tuneR[2]    
    tunePsi<-NULL
    if (object$NC > 0){
        tunePsi<-object$tunePsi[(object$p*object$NC+1):(2*object$p*object$NC)]            
    }   
    tunePhi<-object$tunePhi[(object$p+1):(2*object$p)]     
    tuneCpsi<-object$tuneCpsi[(object$p+1):(2*object$p)]
	#Call C
	if (object$mcm==4) out<-.C("mvrmC",
        as.integer(object$out[[1]]),as.character(object$out[[2]]),
        as.integer(sweeps),as.integer(burn),as.integer(thin),            
        as.double(object$out[[6]]),as.double(object$out[[7]]),as.double(object$out[[8]]),
        as.integer(object$out[[9]]),as.integer(object$out[[10]]),as.integer(object$out[[11]]),
        as.double(object$out[[12]]),as.integer(object$out[[13]]),as.double(object$out[[14]]),
        as.integer(object$out[[15]]),as.double(object$out[[16]]),as.double(object$out[[17]]),
        as.double(object$out[[18]]),as.double(object$out[[19]]),as.integer(object$out[[20]]),
        as.integer(object$out[[21]]),as.integer(object$out[[22]]),as.integer(object$out[[23]]),
        as.integer(object$out[[24]]),as.integer(object$out[[25]]),as.integer(object$out[[26]]),
        as.double(object$out[[27]]),as.integer(object$out[[28]]),as.double(object$out[[29]]),
        as.double(object$out[[30]]),as.double(object$out[[31]]),as.integer(object$out[[32]]),
        as.double(object$out[[33]]),as.double(deviance),as.integer(object$out[[35]]),
        as.integer(cont),as.integer(gamma),as.integer(delta),as.double(alpha),as.double(sigma2),
        as.double(cbeta),as.double(calpha),as.integer(LASTWB),
        as.integer(object$out[[44]]),as.double(object$out[[45]]),as.integer(object$out[[46]]),
        as.double(object$out[[47]]),as.double(object$out[[48]]))
    if (object$mcm==1) out<-.C("mult",
        as.integer(object$out[[1]]),as.character(object$out[[2]]),
        as.integer(sweeps),as.integer(burn),as.integer(thin),
        as.integer(object$out[[6]]),as.integer(object$out[[7]]),as.double(object$out[[8]]),
        as.double(object$out[[9]]),as.double(object$out[[10]]),
        as.integer(object$out[[11]]),as.integer(object$out[[12]]),
        as.double(object$out[[13]]),as.integer(object$out[[14]]),as.double(object$out[[15]]),
        as.integer(object$out[[16]]),
        as.integer(object$out[[17]]),as.integer(object$out[[18]]),as.integer(object$out[[19]]),
        as.integer(object$out[[20]]),
        as.integer(object$out[[21]]),as.integer(object$out[[22]]),as.integer(object$out[[23]]),            
        as.double(object$out[[24]]),as.double(object$out[[25]]),as.double(object$out[[26]]),
        as.double(object$out[[27]]),
        as.double(object$out[[28]]),as.double(object$out[[29]]),
        as.double(object$out[[30]]),as.double(object$out[[31]]),as.double(object$out[[32]]),
        as.integer(object$out[[33]]),as.double(object$out[[34]]),as.double(object$out[[35]]),
        as.integer(object$out[[36]]),as.double(object$out[[37]]),as.double(object$out[[38]]),
        as.integer(object$out[[39]]),as.double(deviance),as.integer(object$out[[41]]),
        as.integer(cont),as.integer(gamma),as.integer(delta),as.double(alpha),as.double(sigma2),
        as.double(R),
        as.double(muR),as.double(sigma2R),as.double(cbeta),as.double(calpha),
        as.integer(c(LASTsw)),as.double(c(DE)),as.integer(LASTWB),as.integer(object$out[[55]]))
    if (object$mcm==2) out<-.C("multg",
        as.integer(object$out[[1]]),as.character(object$out[[2]]),
        as.integer(sweeps),as.integer(burn),as.integer(thin),
        as.integer(object$out[[6]]),as.integer(object$out[[7]]),as.double(object$out[[8]]),
        as.double(object$out[[9]]),as.double(object$out[[10]]),
        as.integer(object$out[[11]]),as.integer(object$out[[12]]),
        as.double(object$out[[13]]),as.integer(object$out[[14]]),as.double(object$out[[15]]),
        as.integer(object$out[[16]]),
        as.integer(object$out[[17]]),as.integer(object$out[[18]]),as.integer(object$out[[19]]),
        as.integer(object$out[[20]]),
        as.integer(object$out[[21]]),as.integer(object$out[[22]]),as.integer(object$out[[23]]),
        as.double(object$out[[24]]),as.double(object$out[[25]]),as.double(object$out[[26]]),
        as.double(object$out[[27]]),
        as.double(object$out[[28]]),as.double(object$out[[29]]),            
        as.double(object$out[[30]]),as.double(object$out[[31]]),as.double(object$out[[32]]),
        as.integer(object$out[[33]]),as.double(object$out[[34]]),as.double(object$out[[35]]),
        as.integer(object$out[[36]]),as.double(object$out[[37]]),as.double(object$out[[38]]),
        as.integer(object$out[[39]]),as.double(deviance),as.integer(object$out[[41]]),
        as.double(object$out[[42]]),as.integer(object$out[[43]]),
        as.integer(cont),as.integer(gamma),as.integer(delta),as.double(alpha),as.double(sigma2),
        as.double(R),
        as.double(muR),as.double(sigma2R),as.double(cbeta),as.double(calpha),
        as.integer(compAlloc),as.double(DPconc),as.integer(c(LASTsw)),as.double(c(DE)),
        as.integer(LASTWB),as.integer(object$out[[59]]))
    if (object$mcm==3) out<-.C("multgv",
        as.integer(object$out[[1]]), as.character(object$out[[2]]), as.integer(sweeps),
        as.integer(burn), as.integer(thin), as.integer(object$out[[6]]), as.integer(object$out[[7]]), 
        as.double(object$out[[8]]), as.double(object$out[[9]]), as.double(as.matrix(object$out[[10]])),
        as.integer(object$out[[11]]), as.integer(object$out[[12]]), as.double(object$out[[13]]), 
        as.integer(object$out[[14]]),
        as.double(object$out[[15]]), as.integer(object$out[[16]]), as.integer(object$out[[17]]), 
        as.integer(object$out[[18]]),
        as.integer(object$out[[19]]), as.integer(object$out[[20]]), as.integer(object$out[[21]]), 
        as.integer(object$out[[22]]),
        as.integer(object$out[[23]]), as.double(object$out[[24]]), as.double(object$out[[25]]), 
        as.double(object$out[[26]]),
        as.double(object$out[[27]]), as.double(object$out[[28]]), as.double(object$out[[29]]), 
        as.double(object$out[[30]]),
        as.double(object$out[[31]]), as.double(object$out[[32]]), as.integer(object$out[[33]]), 
        as.double(object$out[[34]]),
        as.double(object$out[[35]]), as.integer(object$out[[36]]), as.double(object$out[[37]]), 
        as.double(object$out[[38]]),
        as.integer(object$out[[39]]), as.double(deviance), as.integer(object$out[[41]]), 
        as.double(object$out[[42]]),as.integer(object$out[[43]]),
        as.integer(cont), as.integer(gamma), as.integer(delta), as.double(alpha), as.double(sigma2), 
        as.double(R),
        as.double(muR), as.double(sigma2R), as.double(cbeta), as.double(calpha), as.integer(compAllocV), 
        as.double(DPconc), 
        as.integer(c(LASTsw)), as.double(c(DE)), as.integer(LASTWB),as.integer(object$out[[59]]))
    if (object$mcm==5)
        out<-.C("longmult",
        as.integer(object$out[[1]]), as.character(object$out[[2]]), 
        as.integer(c(sweeps,burn,thin,object$out[[3]][4:8])), as.integer(object$out[[4]]), as.integer(object$out[[5]]), 
        as.integer(object$out[[6]]), as.integer(object$out[[7]]), as.integer(object$out[[8]]), as.integer(object$out[[9]]),
        as.integer(object$out[[10]]), as.double(object$out[[11]]), as.double(object$out[[12]]), as.double(object$out[[13]]), 
        as.double(object$out[[14]]), as.double(object$out[[15]]),as.double(object$out[[16]]), as.integer(object$out[[17]]), 
        as.integer(object$out[[18]]), as.integer(object$out[[19]]), as.integer(object$out[[20]]), as.integer(object$out[[21]]), 
        as.integer(object$out[[22]]), as.integer(object$out[[23]]), as.integer(object$out[[24]]), as.integer(object$out[[25]]), 
        as.integer(object$out[[26]]), as.integer(object$out[[27]]), as.integer(object$out[[28]]), as.double(object$out[[29]]), 
        as.integer(object$out[[30]]), as.double(object$out[[31]]), as.double(object$out[[32]]), as.double(object$out[[33]]), 
        as.double(object$out[[34]]), as.double(object$out[[35]]), as.double(object$out[[36]]), as.double(object$out[[37]]), 
        as.double(object$out[[38]]), as.double(object$out[[39]]), as.double(object$out[[40]]), as.double(object$out[[41]]), 
        as.double(object$out[[42]]), as.double(object$out[[43]]), as.integer(object$out[[44]]), as.double(object$out[[45]]),            
        as.integer(object$out[[46]]), as.double(object$out[[47]]), as.double(object$out[[48]]), as.integer(object$out[[49]]), 
        as.double(object$out[[50]]), as.double(object$out[[51]]), as.integer(object$out[[52]]), as.double(object$out[[53]]), 
        as.double(object$out[[54]]), as.integer(object$out[[55]]), as.double(object$out[[56]]),
        as.double(object$out[[57]]), as.integer(object$out[[58]]), as.double(deviance), as.integer(object$out[[60]]),
        as.integer(c(cont,LASTsw,LASTWB)), as.double(LASTAll))
    if (object$mcm==6)
        out<-.C("longmultg",
        as.integer(object$out[[1]]), as.character(object$out[[2]]), 
        as.integer(c(sweeps,burn,thin,object$out[[3]][4:8])), as.integer(object$out[[4]]), as.integer(object$out[[5]]), 
        as.integer(object$out[[6]]), as.integer(object$out[[7]]), as.integer(object$out[[8]]), as.integer(object$out[[9]]),
        as.integer(object$out[[10]]), as.double(object$out[[11]]), as.double(object$out[[12]]), as.double(object$out[[13]]), 
        as.double(object$out[[14]]), as.double(object$out[[15]]),as.double(as.matrix(16)), as.integer(object$out[[17]]), 
        as.integer(object$out[[18]]), as.integer(object$out[[19]]), as.integer(object$out[[20]]), as.integer(object$out[[21]]), 
        as.integer(object$out[[22]]), as.integer(object$out[[23]]), as.integer(object$out[[24]]), as.integer(object$out[[25]]), 
        as.integer(object$out[[26]]), as.integer(object$out[[27]]), as.integer(object$out[[28]]), as.double(object$out[[29]]), 
        as.integer(object$out[[30]]), as.double(object$out[[31]]), as.double(object$out[[32]]), as.double(object$out[[33]]), 
        as.double(object$out[[34]]), as.double(object$out[[35]]), as.double(object$out[[36]]), as.double(object$out[[37]]), 
        as.double(object$out[[38]]), as.double(object$out[[39]]), as.double(object$out[[40]]), as.double(object$out[[41]]), 
        as.double(object$out[[42]]), as.double(object$out[[43]]), as.integer(object$out[[44]]), as.double(object$out[[45]]),
        as.integer(object$out[[46]]), as.double(object$out[[47]]), as.double(object$out[[48]]), as.integer(object$out[[49]]), 
        as.double(object$out[[50]]), as.double(object$out[[51]]), as.integer(object$out[[52]]), as.double(object$out[[53]]),
        as.double(object$out[[54]]), as.integer(object$out[[55]]), as.double(object$out[[56]]),
        as.double(object$out[[57]]), as.integer(object$out[[58]]), as.double(deviance),as.integer(object$out[[60]]), 
        as.integer(c(cont,LASTsw,LASTWB)),
        as.double(LASTAll), as.integer(object$out[[63]]), as.double(object$out[[64]]))
    if (object$mcm==7)
        out<-.C("longmultgv",
        as.integer(object$out[[1]]), as.character(object$out[[2]]), 
        as.integer(c(sweeps,burn,thin,object$out[[3]][4:8])), as.integer(object$out[[4]]), as.integer(object$out[[5]]), 
        as.integer(object$out[[6]]), as.integer(object$out[[7]]), as.integer(object$out[[8]]), as.integer(object$out[[9]]),
        as.integer(object$out[[10]]), as.double(object$out[[11]]), as.double(object$out[[12]]), as.double(object$out[[13]]), 
        as.double(object$out[[14]]), as.double(object$out[[15]]),as.double(as.matrix(16)), as.integer(object$out[[17]]), 
        as.integer(object$out[[18]]), as.integer(object$out[[19]]), as.integer(object$out[[20]]), as.integer(object$out[[21]]), 
        as.integer(object$out[[22]]), as.integer(object$out[[23]]), as.integer(object$out[[24]]), as.integer(object$out[[25]]), 
        as.integer(object$out[[26]]), as.integer(object$out[[27]]), as.integer(object$out[[28]]), as.double(object$out[[29]]), 
        as.integer(object$out[[30]]), as.double(object$out[[31]]), as.double(object$out[[32]]), as.double(object$out[[33]]), 
        as.double(object$out[[34]]), as.double(object$out[[35]]), as.double(object$out[[36]]), as.double(object$out[[37]]), 
        as.double(object$out[[38]]), as.double(object$out[[39]]), as.double(object$out[[40]]), as.double(object$out[[41]]), 
        as.double(object$out[[42]]), as.double(object$out[[43]]), as.integer(object$out[[44]]), as.double(object$out[[45]]),
        as.integer(object$out[[46]]), as.double(object$out[[47]]), as.double(object$out[[48]]), as.integer(object$out[[49]]), 
        as.double(object$out[[50]]), as.double(object$out[[51]]), as.integer(object$out[[52]]), as.double(object$out[[53]]),
        as.double(object$out[[54]]), as.integer(object$out[[55]]), as.double(object$out[[56]]),
        as.double(object$out[[57]]), as.integer(object$out[[58]]), as.double(deviance),as.integer(object$out[[60]]), 
        as.integer(c(cont,LASTsw,LASTWB)),
        as.double(LASTAll), as.integer(object$out[[63]]), as.double(object$out[[64]]))
    #Output
    if (object$mcm<4){
        loc1<-25
        loc2<-41 
        tuneSigma2Ra<-out[[loc1+4]][1]
        tuneRa<-out[[loc1+5]][1] 
        tunePsia<-tunePsi
        tunePhia<-tunePhi
        tuneCpsia<-tuneCpsi   
    }
    if (object$mcm==4){
		loc1<-16 
		loc2<-34
		tuneSigma2Ra<-tuneSigma2R 
		tuneRa<-tuneR
		tunePsia<-tunePsi
        tunePhia<-tunePhi
        tuneCpsia<-tuneCpsi
		DevCalcs<-2
    }
    if (object$mcm %in% c(5,6,7)){
	    loc1<-31
	    loc2<-60
	    tuneSigma2Ra<-out[[loc1+4]][1]
        tuneR<-object$tuneR[(1+object$LUT):(2*object$LUT)]    
        tuneRa<-out[[loc1+5]][1:object$LUT]
        tuneCpsi<-object$tuneCpsi[(object$p*object$p+1):(2*object$p*object$p)]
        tuneCpsia<-out[[loc1+6]][1:object$p*object$p]
        DevCalcs<-2
	}
	#nSamples & mcpar
    totalSweeps <- object$totalSweeps + sweeps
	if (discard==TRUE) object$nSamples <- 0
	nSamples <- nSamples + object$nSamples	
	totalBurn <- object$mcpar[1] + burn 
	if (discard==TRUE) totalBurn <- object$totalSweeps + burn 
	#Call
    call <- object$call
    call2 <- object$call2
    sweeps.pos <- which.max(pmatch(names(call), "sweeps"))
    call[[sweeps.pos]] <- as.integer(totalSweeps)
    burn.pos <- which.max(pmatch(names(call), "burn"))
    call[[burn.pos]] <- totalBurn
    thin.pos <- which.max(pmatch(names(call), "thin"))
    call[[thin.pos]] <- thin
    sweeps.pos <- which.max(pmatch(names(call2), "sweeps"))
    call2[[sweeps.pos]] <- totalSweeps
    burn.pos <- which.max(pmatch(names(call2), "burn"))
    call2[[burn.pos]] <- totalBurn
    thin.pos <- which.max(pmatch(names(call2), "thin"))
    call2[[thin.pos]] <- thin
    #fit
    fit <- list(call=call,call2=call2,formula=object$formula,seed=object$seed,p=object$p,d=object$d,
                data=object$data,Y=object$Y,
                X=object$X,Xknots=object$Xknots,LG=object$LG,NG=object$NG,NG2=object$NG2,
                Z=object$Z,Zknots=object$Zknots,LD=object$LD,ND=object$ND,ND2=object$ND2,
                W=object$W,Wknots=object$Wknots,LC=object$LC,NC=object$NC,
                storeMeanVectorX=object$storeMeanVectorX,
                storeMeanVectorZ=object$storeMeanVectorZ,
                storeMeanVectorW=object$storeMeanVectorW,
                storeIndicatorX=object$storeIndicatorX,                
                storeIndicatorZ=object$storeIndicatorZ,
                storeIndicatorW=object$storeIndicatorW,
                assignX=object$assignX,
                assignZ=object$assignZ,
                assignW=object$assignW,
                labelsX=object$labelsX,
                labelsZ=object$labelsZ,
                labelsW=object$labelsW,
                countX=object$countX,
                countZ=object$countZ,
                countW=object$countW,
                varsY=object$varsY,
                varsX=object$varsX,
                varsZ=object$varsZ,
                varsW=object$varsW,
                is.Dx=object$is.Dx,
                is.Dz=object$is.Dz,
                is.Dw=object$is.Dw,
                which.SpecX=object$which.SpecX,
                which.SpecZ=object$which.SpecZ,
                which.SpecW=object$which.SpecW,
                formula.termsX=object$formula.termsX,
                formula.termsZ=object$formula.termsZ,
                formula.termsW=object$formula.termsW,
                nSamples=nSamples,
                totalSweeps=totalSweeps,
                mcpar=c(totalBurn,as.integer(seq(from=totalBurn,by=thin,length.out=nSamples)[nSamples]),thin),
                mcm=object$mcm,H=object$H,G=object$G,
                tuneCalpha=c(tuneCalpha,out[[loc1+0]][1:object$p]),    
                tuneSigma2=c(tuneSigma2,out[[loc1+1]][1:object$p]),
                tuneCbeta=c(tuneCbeta,out[[loc1+2]][1]),
                tuneAlpha=c(tuneAlpha,out[[loc1+3]][1:(object$p*object$ND2)]),
                tuneSigma2R=c(tuneSigma2R,tuneSigma2Ra),
                tuneR=c(tuneR,tuneRa),
                tunePsi=c(tunePsi,tunePsia),
                tunePhi=c(tunePhi,tunePhia),
                tuneCpsi=c(tuneCpsi,tuneCpsia),
                deviance=c(out[[loc2]][1:DevCalcs]),
                nullDeviance=object$nullDeviance,                            
                DIR=object$DIR,
                out=out,
                LUT=object$LUT,
                SUT=object$SUT, 
                LGc=object$LGc,
                LDc=object$LDc,                
                NGc=object$NGc, 
                NDc=object$NDc, 
                NK=object$NK,
                FT=object$FT,
                qCont=0,
                HNca=object$HNca,HNsg=object$HNsg,HNcp=object$HNcp,HNphi=object$HNphi)
    if (object$mcm %in% c(5,6,7)){
        tuneCbCor<-object$tuneCbCor[2]
	    tuneOmega<-object$tuneOmega[(object$NDc+1):(object$NDc*2)]	            
        tuneComega<-object$tuneComega[2]
		fit2 <- list(
		        C=object$C,Cknots=object$Cknots,LK=object$LK,
                Xc=object$Xc,Xcknots=object$Xcknots,
                Zc=object$Zc,Zcknots=object$Zcknots,                
                storeMeanVectorC=object$storeMeanVectorC,
                storeMeanVectorXc=object$storeMeanVectorXc,
                storeMeanVectorZc=object$storeMeanVectorZc,                               
                storeIndicatorC=object$storeIndicatorC,
                storeIndicatorXc=object$storeIndicatorXc,
                storeIndicatorZc=object$storeIndicatorZc, 
                assignC=object$assignC,          
                assignXc=object$assignXc,
                assignZc=object$assignZc,
                labelsC=object$labelsC,                
                labelsXc=object$labelsXc,
                labelsZc=object$labelsZc,                                               
                countC=object$countC,                            
                countXc=object$countXc,
                countZc=object$countZc,
                varsC=object$varsC,                
                varsXc=object$varsXc,
                varsZc=object$varsZc,  
                is.Dc=object$is.Dc,
                is.Dxc=object$is.Dxc,
                is.Dzc=object$is.Dzc,  
                which.SpecC=object$which.SpecC,                
                which.SpecXc=object$which.SpecXc,
                which.SpecZc=object$which.SpecZc, 
                formula.termsC=object$formula.termsC,
                formula.termsXc=object$formula.termsXc,
                formula.termsZc=object$formula.termsZc, 
	            tuneCbCor=c(tuneCbCor,out[[loc1+7]][1]),	            
	            tuneOmega=c(tuneOmega,out[[loc1+8]][1:object$NDc]),	            
                tuneComega=c(tuneComega,out[[loc1+9]][1]),
                HNcpsi=object$HNcpsi,HNscor=object$HNscor,HNco=object$HNco,
                lag=object$lag,varTime=object$varTime,niVec=object$niVec,intime=object$intime)
        fit<-c(fit,fit2)        
	}
    class(fit) <- 'mvrm'
    return(fit)
}

mvrm2mcmc <- function(mvrmObj,labels){
    labels1 <- c("alpha","calpha","cbeta","delta","beta","gamma","sigma2") #7
    labels2 <- c("muR","sigma2R","R") #3
    labels3 <- c("compAlloc","nmembers","DPconc") #3
    labels4 <- c("compAllocV","nmembersV","deviance")#3
    labels5 <- c("psi","ksi","cpsi","nu","fi","omega","comega","eta","ceta")#9
    labels6 <- c("probs","tune")#2
    labels7 <- c("phi2")#1
    labels8 <- c("Hbeta")#1
    all.labels<-c(labels1,labels2,labels3,labels4,labels5,labels6,labels7,labels8)
    if (missing(labels)) labels <- all.labels
    mtch<-match(labels,all.labels)
	p<-mvrmObj$p
	R<-NULL
    if (any(mtch==1) && mvrmObj$LD > 0){
        file <- paste(mvrmObj$DIR,"BNSP.alpha.txt",sep="")
        if (file.exists(file)){ 
           if (p > 1) names1<-paste("alpha",rep(colnames(mvrmObj$Z)[-1],p),rep(mvrmObj$varsY,each=mvrmObj$LD),sep=".")
           if (p == 1) names1<-paste("alpha",colnames(mvrmObj$Z)[-1],sep=".")
           R<-cbind(R,matrix(unlist(read.table(file)),ncol=mvrmObj$LD*p,dimnames=list(c(),names1)))
	   }
    }
    if (any(mtch==2)){
        file <- paste(mvrmObj$DIR,"BNSP.calpha.txt",sep="")
        if (file.exists(file)){ 
            if (p > 1) names1<-paste(rep("c_alpha",p),mvrmObj$varsY,sep=".")
            if (p == 1) names1<-"c_alpha"
            R<-cbind(R,matrix(unlist(read.table(file)),ncol=p,dimnames=list(c(),names1)))
		}
    }
    if (any(mtch==3)){
       file <- paste(mvrmObj$DIR,"BNSP.cbeta.txt",sep="")
       if (file.exists(file)) 
           R<-cbind(R,matrix(unlist(read.table(file)),ncol=1,dimnames=list(c(),c("c_beta"))))
    }
    if (any(mtch==4) && mvrmObj$LD > 0){
        file <- paste(mvrmObj$DIR,"BNSP.delta.txt",sep="")
        if (file.exists(file)){ 
            if (p > 1) names1<-paste("delta",rep(colnames(mvrmObj$Z)[-1],p),rep(mvrmObj$varsY,each=mvrmObj$LD),sep=".")
            if (p == 1) names1<-paste("delta",colnames(mvrmObj$Z)[-1],sep=".")
            R<-cbind(R,matrix(unlist(read.table(file)),ncol=mvrmObj$LD*p,dimnames=list(c(),names1)))
	    }
    }
    if (any(mtch==5)){
        file <- paste(mvrmObj$DIR,"BNSP.beta.txt",sep="")
        if (file.exists(file)){ 
            if (p > 1) names1<-paste("beta",rep(colnames(mvrmObj$X),p),rep(mvrmObj$varsY,each=mvrmObj$LG+1),sep=".")
            if (p == 1) names1<-paste("beta",colnames(mvrmObj$X),sep=".")
            R<-cbind(R,matrix(unlist(read.table(file)),ncol=p*(mvrmObj$LG+1),dimnames=list(c(),names1)))
		}
    }
    if (any(mtch==6) && mvrmObj$LG > 0){
        file <- paste(mvrmObj$DIR,"BNSP.gamma.txt",sep="")
        if (file.exists(file)){ 
            if (p > 1) names1<-paste("gamma",rep(colnames(mvrmObj$X)[-1],p),rep(mvrmObj$varsY,each=mvrmObj$LG),sep=".")
            if (p == 1) names1<-paste("gamma",colnames(mvrmObj$X)[-1],sep=".")
            R<-cbind(R,matrix(unlist(read.table(file)),ncol=p*mvrmObj$LG,dimnames=list(c(),names1)))
		}
    }
    if (any(mtch==7)){
        file <- paste(mvrmObj$DIR,"BNSP.sigma2.txt",sep="")
        if (file.exists(file)) 
            if (p > 1) names1<-paste(rep("sigma2",p),mvrmObj$varsY,sep=".")
            if (p == 1) names1<-"sigma2"
            R<-cbind(R,matrix(unlist(read.table(file)),ncol=p,dimnames=list(c(),names1)))
	}
	if (any(mtch==8) && p > 1 && mvrmObj$LUT==1){
		names <- paste("muR_clust",seq(1,mvrmObj$H),sep="_")
		if (mvrmObj$H==1) names <- "muR"
	    file <- paste(mvrmObj$DIR,"BNSP.muR.txt",sep="")
        if (file.exists(file)) R<-cbind(R,matrix(unlist(read.table(file)),nrow=mvrmObj$nSamples,
                dimnames=list(c(),names)))
    }
    if (any(mtch==9) && p > 1){
        file <- paste(mvrmObj$DIR,"BNSP.sigma2R.txt",sep="")
        if (file.exists(file)) R<-cbind(R,matrix(unlist(read.table(file)),ncol=1,dimnames=list(c(),c("sigma2R"))))
	}
	subset <- rep(seq(1,p),each=p) < rep(seq(1,p),p)
	cor.names <- paste(rep(seq(1,p),each=p),rep(seq(1,p),p),sep="")
    if (any(mtch==10) && p > 1){
        file <- paste(mvrmObj$DIR,"BNSP.R.txt",sep="")
        if (file.exists(file)) R<-cbind(R,matrix(unlist(read.table(file)),ncol=mvrmObj$LUT*p*p,
                 dimnames=list(c(),paste("cor",rep(cor.names,mvrmObj$LUT),sep=".")))[,subset,drop=FALSE])
    }
	if (any(mtch==11)){
	    file <- paste(mvrmObj$DIR,"BNSP.compAlloc.txt",sep="")
	    if (file.exists(file)) R<-cbind(R,matrix(unlist(read.table(file)),ncol=mvrmObj$d,
	                           dimnames=list(c(),paste("clustering of cor",cor.names[subset],sep="."))))
	}
    if (any(mtch==12)){
        file <- paste(mvrmObj$DIR,"BNSP.nmembers.txt",sep="")
        if (file.exists(file))R<-cbind(R,matrix(unlist(read.table(file)),ncol=mvrmObj$H,dimnames=list(c(),paste("cor cluster ",seq(1,mvrmObj$H)))))
    }    
    if (any(mtch==13)){
        file <- paste(mvrmObj$DIR,"BNSP.DPconc.txt",sep="")
        if (file.exists(file)) R<-cbind(R,matrix(unlist(read.table(file)),ncol=1,dimnames=list(c(),c("DPconc"))))
    }
    if (any(mtch==14)){
        file <- paste(mvrmObj$DIR,"BNSP.compAllocV.txt",sep="")
        if (file.exists(file)) R<-cbind(R,matrix(unlist(read.table(file)),ncol=p,dimnames=list(c(),paste("clustering of ",mvrmObj$varsY))))
	}	
    if (any(mtch==15)){
        file <- paste(mvrmObj$DIR,"BNSP.nmembersV.txt",sep="")
        if (file.exists(file)) R<-cbind(R,matrix(unlist(read.table(file)),ncol=mvrmObj$G,dimnames=list(c(),paste("var cluster ",seq(1,mvrmObj$G)))))
	}
	if (any(mtch==16)){
       file <- paste(mvrmObj$DIR,"BNSP.deviance.txt",sep="")
       if (file.exists(file)) R<-cbind(R,matrix(unlist(read.table(file)),ncol=2,dimnames=list(c(),c("marginal deviance","full deviance"))))
    }
    if (any(mtch==17) && mvrmObj$mcm < 8){
        file <- paste(mvrmObj$DIR,"BNSP.psi.txt",sep="")
        if (file.exists(file)){ 
            if (p > 1){ 
                e <- cbind(rep(unlist(mvrmObj$varsY),each=p),rep(unlist(mvrmObj$varsY),p))
                e <- paste(e[,1],e[,2],sep=".")
                names1<-paste("psi",rep(colnames(mvrmObj$C),p*p),rep(e,each=mvrmObj$LK),sep=".")
			}
            if (p == 1) names1<-paste("psi",colnames(mvrmObj$C),sep=".")
            R<-cbind(R,matrix(unlist(read.table(file)),ncol=p*p*mvrmObj$LK,dimnames=list(c(),names1)))
		}
    }    
    if (any(mtch==17) && mvrmObj$LC > 0 && mvrmObj$mcm > 7){
        file <- paste(mvrmObj$DIR,"BNSP.psi.txt",sep="")
        if (file.exists(file)){ 
           if (p > 1) names1<-paste("psi",rep(colnames(mvrmObj$W)[-1],p),rep(mvrmObj$varsY,each=mvrmObj$LC),sep=".")
           if (p == 1) names1<-paste("psi",colnames(mvrmObj$W)[-1],sep=".")
           R<-cbind(R,matrix(unlist(read.table(file)),ncol=mvrmObj$LC*p,dimnames=list(c(),names1)))
	   }
    }
    if (any(mtch==18) && mvrmObj$mcm < 8){
        file <- paste(mvrmObj$DIR,"BNSP.ksi.txt",sep="")
        if (file.exists(file)){ 
            if (p > 1){ 
                e <- cbind(rep(unlist(mvrmObj$varsY),each=p),rep(unlist(mvrmObj$varsY),p))
                e <- paste(e[,1],e[,2],sep=".")
                names1<-paste("ksi",rep(colnames(mvrmObj$C),p*p),rep(e,each=mvrmObj$LK),sep=".")
			}
            if (p == 1) names1<-paste("ksi",colnames(mvrmObj$C),sep=".")
            R<-cbind(R,matrix(unlist(read.table(file)),ncol=p*p*mvrmObj$LK,dimnames=list(c(),names1)))
		}
    }
    if (any(mtch==18) && mvrmObj$LC > 0 && mvrmObj$mcm > 7){
        file <- paste(mvrmObj$DIR,"BNSP.ksi.txt",sep="")
        if (file.exists(file)){ 
            if (p > 1) names1<-paste("ksi",rep(colnames(mvrmObj$W)[-1],p),rep(mvrmObj$varsY,each=mvrmObj$LC),sep=".")
            if (p == 1) names1<-paste("ksi",colnames(mvrmObj$W)[-1],sep=".")
            R<-cbind(R,matrix(unlist(read.table(file)),ncol=mvrmObj$LC*p,dimnames=list(c(),names1)))
	    }
    }
    if (any(mtch==19) && mvrmObj$mcm < 8){
        file <- paste(mvrmObj$DIR,"BNSP.cpsi.txt",sep="")
        if (file.exists(file)){ 
            if (p > 1){ 
                e <- cbind(rep(unlist(mvrmObj$varsY),each=p),rep(unlist(mvrmObj$varsY),p))
                e <- paste(e[,1],e[,2],sep=".")
                names1<-paste("cpsi",e,sep=".")
			}
            if (p == 1) names1<-paste("cpsi")
            R<-cbind(R,matrix(unlist(read.table(file)),ncol=p*p,dimnames=list(c(),names1)))
		}
    }
    if (any(mtch==19) && mvrmObj$mcm > 7){
        file <- paste(mvrmObj$DIR,"BNSP.cpsi.txt",sep="")
        if (file.exists(file)){ 
            if (p > 1) names1<-paste(rep("c_psi",p),mvrmObj$varsY,sep=".")
            if (p == 1) names1<-"c_psi"
            R<-cbind(R,matrix(unlist(read.table(file)),ncol=p,dimnames=list(c(),names1)))
		}
    }    
    if (any(mtch==20) && mvrmObj$LGc > 0){
        file <- paste(mvrmObj$DIR,"BNSP.nu.txt",sep="")
        if (file.exists(file)){ 
            if (mvrmObj$H==1) names1<-paste("nu",colnames(mvrmObj$Xc)[-1],sep=".")
            if (mvrmObj$H>1) names1<-paste(paste("nu",rep(1:mvrmObj$H,each=mvrmObj$LGc),sep="."),rep(seq(1:mvrmObj$LGc),mvrmObj$H),sep=".")
            R<-cbind(R,matrix(unlist(read.table(file)),ncol=mvrmObj$LGc*mvrmObj$H,dimnames=list(c(),names1)))
		}
    }    
    if (any(mtch==21) && mvrmObj$LDc > 0){
        file <- paste(mvrmObj$DIR,"BNSP.fi.txt",sep="")
        if (file.exists(file)){ 
            names1<-paste("fi",colnames(mvrmObj$Zc)[-1],sep=".")
            R<-cbind(R,matrix(unlist(read.table(file)),ncol=mvrmObj$LDc,dimnames=list(c(),names1)))
		}
    }
    if (any(mtch==22) && mvrmObj$LDc > 0){
        file <- paste(mvrmObj$DIR,"BNSP.omega.txt",sep="")
        if (file.exists(file)){ 
            names1<-paste("omega",colnames(mvrmObj$Zc)[-1],sep=".")
            R<-cbind(R,matrix(unlist(read.table(file)),ncol=mvrmObj$LDc,dimnames=list(c(),names1)))
		}
    }
    if (any(mtch==23) && p > 1){
        file <- paste(mvrmObj$DIR,"BNSP.comega.txt",sep="")
        if (file.exists(file)){ 
            names1<-"c_omega"
            R<-cbind(R,matrix(unlist(read.table(file)),ncol=1,dimnames=list(c(),names1)))
		}
    }
    if (any(mtch==24) && p > 1){
        file <- paste(mvrmObj$DIR,"BNSP.eta.txt",sep="")
        if (file.exists(file)){ 
            if (mvrmObj$H==1) names1<-paste("eta",colnames(mvrmObj$Xc),sep=".")
            if (mvrmObj$H>1) names1<-paste(paste("eta",rep(1:mvrmObj$H,each=(mvrmObj$LGc+1)),sep="."),rep(seq(1:(mvrmObj$LGc+1)),mvrmObj$H),sep=".")
            R<-cbind(R,matrix(unlist(read.table(file)),ncol=((mvrmObj$LGc+1)*mvrmObj$H),dimnames=list(c(),names1)))
		}
    }
    if (any(mtch==25) && p > 1){
       file <- paste(mvrmObj$DIR,"BNSP.ceta.txt",sep="")
       if (file.exists(file)) 
           R<-cbind(R,matrix(unlist(read.table(file)),ncol=1,dimnames=list(c(),c("c_eta"))))
    }
    if (any(mtch==26)){ 
        file <- paste(mvrmObj$DIR,"BNSP.probs.txt",sep="")
        if (file.exists(file))
            R<-cbind(R,matrix(unlist(read.table(file)),ncol=mvrmObj$H,dimnames=list(c(),seq(1,mvrmObj$H))))
    }
    if (any(mtch==27)){ 
        file <- paste(mvrmObj$DIR,"BNSP.tune.txt",sep="")
        if (file.exists(file)){
			
			if (p > 1) names1<-paste("alpha",rep(mvrmObj$varsY,each=mvrmObj$ND),rep(mvrmObj$varsZ,p),sep=".") 
            if (p == 1) names1<-paste("alpha",mvrmObj$varsZ,sep=".")
            
            if (p > 1) names2<-paste(rep("sigma2",p),mvrmObj$varsY,sep=".")
            if (p == 1) names2<-"sigma2"
            			
			if (p > 1) names3<-paste(rep("c_alpha",p),mvrmObj$varsY,sep=".")
            if (p == 1) names3<-"c_alpha"
            
            if (p > 1){ 
                e <- cbind(rep(unlist(mvrmObj$varsY),each=p),rep(unlist(mvrmObj$varsY),p))
                e <- paste(e[,1],e[,2],sep=".")
                names4<-paste("cpsi",e,sep=".")
			}
            if (p == 1) names4<-paste("cpsi")
            
            names5<-NULL
            if (p > 1) names5<-paste("R",seq(1,mvrmObj$LUT),sep="_")
            
            names6<-NULL
            if (p > 1 && mvrmObj$NDc > 0) {names6<-"omega"}
            names7<-NULL
            if (p > 1) {names7<-"c_eta"}
            
            tune.names<-c("c_beta",names1,names2,names3[mvrmObj$HNca==1],names4[mvrmObj$HNcpsi==1],names5,names6,
                          "sigms2R"[mvrmObj$HNscor==1],names7,"c_omega"[mvrmObj$HNco==1])
            
            R<-cbind(R,matrix(unlist(read.table(file)),ncol=length(tune.names),dimnames=list(c(),tune.names)))
		}
    }
    if (any(mtch==28)){
        file <- paste(mvrmObj$DIR,"BNSP.phi2.txt",sep="")
        if (file.exists(file)){ 
            if (p > 1) names1<-paste(rep("phi2",p),mvrmObj$varsY,sep=".")
            if (p == 1) names1<-"phi2"
            R<-cbind(R,matrix(unlist(read.table(file)),ncol=p,dimnames=list(c(),names1)))
		}
	}
	if (any(mtch==29)){
        file <- paste(mvrmObj$DIR,"BNSP.Hbeta.txt",sep="")
        if (file.exists(file)){ 
            
            names1 <- NULL
            for (k in 1:mvrmObj$nHar)
                names1<-c(names1,paste(paste("sin(",2*k,sep=""),"pi", mvrmObj$varSin, "/ p)",sep=" "),paste(paste("cos(",2*k,sep=""),"pi", mvrmObj$varSin, "/ p)",sep=" "))
            
            if (p > 1) names1<-paste(rep(names1,p),rep(mvrmObj$varsY,2*mvrmObj$nHar),sep=".")
            
            R<-cbind(R,matrix(unlist(read.table(file)),ncol=2*mvrmObj$nHar*p,dimnames=list(c(),names1)))
		}
	}
	
	if (!is.null(R) && !any(mtch==27)){
	    attr(R, "mcpar") <- mvrmObj$mcpar
        attr(R, "class") <- "mcmc"
	}
	if (!is.null(R) && any(mtch==27)){
	    attr(R, "mcpar") <- c(1,mvrmObj$mcpar[2],50)
        attr(R, "class") <- "mcmc"
	}
	return(R)
}

plotCorr <- function(x, term="R", centre=mean, quantiles=c(0.1, 0.9), ...){
    mvrmObj<-x
    if (mvrmObj$p == 1) stop("doesn't apply for univariate response")
    if (length(quantiles)==1) quantiles <- c(quantiles,1-quantiles)
    if (length(quantiles) > 2) stop("up to two quantiles")
    if (!is.null(quantiles)) quantiles <- unique(sort(quantiles))
    if (!is.null(quantiles)) par(mfrow=c(1,2))
    mt<-match(term,c("R","muR"))
    if (is.na(mt)) stop("unrecognised term");
    if (mt==1) R<-mvrm2mcmc(mvrmObj,"R")
    if (mt==2){
        muR<-mvrm2mcmc(mvrmObj,"muR")
        if (mvrmObj$mcm == 1) return(plot(tanh(muR)))
        compAlloc <- matrix(0,nrow=mvrmObj$nSamples,ncol=mvrmObj$d)
        compAllocFile <- paste(mvrmObj$DIR, paste("BNSP","compAlloc","txt",sep="."),sep="/")
        if (file.exists(compAllocFile)) compAlloc<-mvrm2mcmc(mvrmObj,"compAlloc")
        R <- sapply(1:mvrmObj$d,
                    function(i) muR[cbind(1:mvrmObj$nSamples,compAlloc[,i]+1)])
        if (x$FT==1) R <- tanh(R)
	}
    ec <- diag(rep(1,x$p))
    ec[lower.tri(ec)] <- apply(R,2,centre)
    ec[upper.tri(ec)]<-t(ec)[upper.tri(ec)]
    colnames(ec) <- x$varsY
    corrplot.mixed(ec,lower.col="black")
    if (!is.null(quantiles)){
        q<-apply(R,2,quantile,probs=quantiles)
        ci<-diag(1,mvrmObj$p)
        ci[lower.tri(ci)] <- q[2,]
        ci[upper.tri(ci)]<-t(ci)[upper.tri(ci)]
        ci[lower.tri(ci)] <- q[1,]
        colnames(ci) <- x$varsY
        rownames(ci) <- x$varsY
        corrplot(ci,col="black",method="number")
    }
}

histCorr <- function(x, term="R", plotOptions=list(), ...){
    mvrmObj<-x
    if (mvrmObj$p == 1) stop("doesn't apply for univariate response")
    mt<-match(term,c("R","muR"))
    if (is.na(mt)) stop("unrecognised term");
    if (mt==1) R<-mvrm2mcmc(mvrmObj,"R")
    if (mt==2){
        muR<-mvrm2mcmc(mvrmObj,"muR")
        if (mvrmObj$mcm == 1) return(plot(tanh(muR)))
        compAlloc <- matrix(0,nrow=mvrmObj$nSamples,ncol=mvrmObj$d)
        compAllocFile <- paste(mvrmObj$DIR, paste("BNSP","compAlloc","txt",sep="."),sep="/")
        if (file.exists(compAllocFile)) compAlloc<-mvrm2mcmc(mvrmObj,"compAlloc")
        R <- sapply(1:mvrmObj$d,
                    function(i) muR[cbind(1:mvrmObj$nSamples,compAlloc[,i]+1)])
        if (x$FT==1) R <- tanh(R)
	}
    r<-rep(rep(seq(1,mvrmObj$p-1),times=seq(mvrmObj$p-1,1)),each=mvrmObj$nSample)
    c<-rep(unlist(mapply(seq,seq(2,mvrmObj$p),mvrmObj$p)),each=mvrmObj$nSample)
    df<-data.frame(cor=c(R),r=r,c=c)
    pp<-ggplot(df) + geom_histogram(aes(x=cor),binwidth=0.01) + facet_wrap(r~c) + plotOptions #facet_grid(r~c)
    return(pp)
}

predict.mvrm <- function(object,newdata,interval=c("none","credible","prediction"),level=0.95,
                         nSamples=100, ind.preds=FALSE,...){
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
    
    #no.match<-all(is.na(match(colnames(newdata),colnames(object$data))))
    #if (length(object$data)>0 && !no.match){
    #    nd<-object$data[0,match(colnames(newdata),colnames(object$data)),drop=FALSE]
    #    for (j in 1:dim(nd)[2])
    #        nd[,j]<-drop(nd[,j])
    #    nd[1:NROW(newdata),] <- newdata
    #}else{nd<-newdata}
    
    nd<-newdata
    npred<-NROW(newdata)
    X<-DM(formula=formula2,data=nd,n=npred,knots=object$Xknots,predInd=TRUE,meanVector=object$storeMeanVectorX,indicator=object$storeIndicatorX,mvrmObj=object,centre=TRUE)$X
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
        Z<-DM(formula=formula2,data=nd,n=npred,knots=object$Zknots,predInd=TRUE,meanVector=object$storeMeanVectorZ,indicator=object$storeIndicatorZ,mvrmObj=object,centre=TRUE)$X
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
	returns <- data.frame(predictions)
    if (ind.preds) returns <- list(returns,fitM) 		
	return(returns)
}

print.mvrm <- function(x,  digits = 5, ...) {
    cat("\nCall:\n")
    print(x$call)
    cat("\n")
    cat(x$nSamples,"posterior samples\n")
    G<-as.data.frame(mvrm2mcmc(x,"gamma"))
    D<-as.data.frame(mvrm2mcmc(x,"delta"))
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

summary.mvrm <- function(object, nModels = 5, digits = 5, printTuning = FALSE, ...) {
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
      G<-as.data.frame(mvrm2mcmc(object,"gamma"))
      D<-as.data.frame(mvrm2mcmc(object,"delta"))
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
	    names(object$tuneCbeta)<-dm
	    sigma2<-c(t(matrix(object$tuneSigma2,ncol=2)))
	    names(sigma2)<-rep(dm,object$p)	    
	    names(object$tuneSigma2R)<-dm	    
	    R<-c(t(matrix(object$tuneR,ncol=2)))
	    names(R)<-rep(dm,object$LUT)	    
	    #rearrange<-lapply(seq(1:object$p),function(x)seq(x,2*object$p,by=object$p))
	    #sigma2<-lapply(rearrange,function(x){object$tuneSigma2[x]})
	    #for (k in 1:object$p) names(sigma2[[k]])<-dm	    	    	    
	    #c.alpha<-lapply(rearrange,function(x){object$tuneCa[x]})
	    #for (k in 1:object$p) names(c.alpha[[k]])<-dm
	    c.alpha<-c(t(matrix(object$tuneCalpha,ncol=2)))
	    names(c.alpha)<-rep(dm,object$p)	    
	    if (object$ND > 0){
	        tot<-object$p*object$ND
	        #rearrange<-lapply(seq(1:tot),function(x)seq(x,2*tot,by=tot))
	        #tuneAlpha<-lapply(rearrange,function(x){object$tuneAlpha[x]})
	        #for (k in 1:tot) names(tuneAlpha[[k]])<-dm
	        tuneAlpha<-c(t(matrix(object$tuneAlpha,ncol=2)))
	        names(tuneAlpha)<-rep(dm,tot)
	    }
	    pT1<-list(c.beta=object$tuneCbeta,sigma2=sigma2)
	    if (object$ND > 0) {pT2<-list(Alpha=tuneAlpha, c.alpha = c.alpha); pT1<-c(pT1,pT2)}
	    if (object$p > 1) {pT2<-list(sigma2R = object$tuneSigma2R, R =  R); pT1<-c(pT1,pT2)}
	    print(pT1)
	}
}

clustering <- function(object, ...){
    R <- list()
    if (object$mcm==1 || object$mcm==5 || object$mcm==8) stop("common correlations model has no clustering")
    if (object$mcm %in% c(2,3,6,7,9,10)){
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
    if (object$mcm %in%c(3,7,10)){
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
