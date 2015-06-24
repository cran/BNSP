bnpglmf <- function(formula,data,sampler="slice",StorageDir,ncomp,sweeps,burn,thin=1,seed,
                    V,Vdf,Mu.mu,Sigma.mu,Alpha.alpha,Beta.alpha,Turnc.alpha,Xpred,...){
    # Match call
    call <- match.call()
    # Data environment & design matrix
    if (missing(data))
        data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    X<-model.matrix(formula,data=data)[,-1]
    Y <- model.response(mf, "any")
    # Dimensions
    n<-NROW(X)
    p<-NCOL(X)
    # Prior parameters
    if (missing(seed)) seed<-as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31)
    if (missing(V)) V<-diag(1,p+1)/(p+1)
    if (missing(Vdf)) Vdf<-p+1
    if (missing(Mu.mu) & p > 1) Mu.mu<-c(mean(Y),apply(X,2,mean))
    if (missing(Mu.mu) & p == 1) Mu.mu<-c(mean(Y),mean(X))
    if (missing(Sigma.mu) & p > 1) Sigma.mu<-diag(c((max(Y)-min(Y))^2,(apply(X,2,max)-apply(X,2,min))^2))
    if (missing(Sigma.mu) & p == 1) Sigma.mu<-diag(c((max(Y)-min(Y))^2,(max(X)-min(X))^2))
    if (missing(Alpha.alpha)) Alpha.alpha<-2
    if (missing(Beta.alpha)) Beta.alpha<-4
    if (missing(Turnc.alpha)) Turnc.alpha<-0.25
    if (p > 1){
        xbar<-apply(X,2,mean)
        xsd<-apply(X,2,sd)
    }else if (p == 1){
        xbar<-mean(X)
        xsd<-sd(X)
    }
    # Predictions
    if (missing(Xpred)){
        npred <- 0
        Xpred <- 1
    }else{
        npred <- length(Xpred)/p
    }
    meanReg <- array(0,npred)
    medianReg <- array(0,npred)
    q1Reg <- array(0,npred)
    q3Reg <- array(0,npred)
    modeReg <- array(0,npred)
    # Range of predictions
    maxy <- 2000
    #Sampler
    sampler.indicator <- match(sampler,c("slice","truncated"))
    if (is.na(sampler.indicator)){
        stop(c(sampler," 'sampler' not recognized"))
    }
    # Storage directory & files
    WF <- 1
    if (!missing(StorageDir)){
        StorageDir <- path.expand(StorageDir)
        ncharwd <- nchar(StorageDir)}
    if (!missing(StorageDir)) if (!(substr(StorageDir,ncharwd,ncharwd)=="/")) StorageDir <- paste(StorageDir,"/",sep="")
    if (!missing(StorageDir)) if (!file.exists(StorageDir)) dir.create(StorageDir)
    if (missing(StorageDir)) {WF <- 0; StorageDir <- paste(getwd(),"/",sep="")}
    on.exit(if (WF==0) file.remove(paste(StorageDir,"BNSP.Th.txt",sep=""), paste(StorageDir,"BNSP.Sigmah.txt",sep=""),
    paste(StorageDir,"BNSP.muh.txt",sep=""), paste(StorageDir,"BNSP.alpha.txt",sep=""),
    paste(StorageDir,"BNSP.compAlloc.txt",sep=""), paste(StorageDir,"BNSP.nmembers.txt",sep=""),
    paste(StorageDir,"BNSP.Updated.txt",sep="")))
    #Call C
    out<-.C("OneResLtntFx", as.integer(seed), as.double(unlist(c(X))), as.integer(Y),
            as.integer(sweeps), as.integer(burn), as.integer(thin), as.integer(ncomp),
            as.integer(n), as.integer(p),
            as.double(V), as.double(Vdf),
            as.double(Mu.mu), as.double(Sigma.mu),
            as.double(Alpha.alpha),as.double(Beta.alpha),as.double(Turnc.alpha),
            as.double(xbar), as.double(xsd), as.double(mean(Y)), as.double(sd(Y)),
            as.integer(sampler.indicator),
            as.integer(npred),as.double(Xpred),as.integer(maxy),
            as.double(meanReg),as.double(medianReg),as.double(q1Reg),as.double(q3Reg),as.double(modeReg),
            as.character(StorageDir),as.integer(WF))
    #Output
    location<-25
    meanReg <- out[[location+0]][1:npred]
    medianReg <- out[[location+1]][1:npred]
    q1Reg <- out[[location+2]][1:npred]
    q3Reg <- out[[location+3]][1:npred]
    modeReg <- out[[location+4]][1:npred]
    fit <- list(call=call,seed=seed,meanReg=meanReg,medianReg=medianReg,q1Reg=q1Reg,q3Reg=q3Reg,modeReg=modeReg)
    class(fit) <- 'bnp'
    return(fit)
}
