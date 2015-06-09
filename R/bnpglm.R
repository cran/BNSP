bnpglm <- function(formula,family,data,offset,sampler="slice",StorageDir,ncomp,sweeps,burn,thin=1,seed,prec,
                   V,Vdf,Mu.nu,Sigma.nu,Mu.mu,Sigma.mu,Alpha.xi,Beta.xi,Alpha.alpha,Beta.alpha,Turnc.alpha,
                   Xpred,offsetPred,...){
    # Match call
    call <- match.call()
    # Family
    family.indicator <- match(c(family),c("poisson","binomial","negative binomial","beta binomial","generalized poisson"))
    if (is.na(family.indicator)){
        stop('family must be character, and one of "poisson", "binomial", "negative binomial", "beta binomial" and "generalized poisson"')
    }
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
    # Dimensions
    n<-NROW(X)
    p<-NCOL(X)
    # Prior parameters
    if (missing(seed)) seed<-as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31)
    if (missing(V)) V<-diag(1,p)/p
    if (missing(Vdf)) Vdf<-p
    if (missing(Mu.nu)) Mu.nu<-rep(0,p)
    if (missing(Sigma.nu)) Sigma.nu<-diag(1,p)*100
    if (missing(Mu.mu) & p > 1) Mu.mu<-apply(X,2,mean)
    if (missing(Mu.mu) & p == 1) Mu.mu<-mean(X)
    if (missing(Sigma.mu) & p > 1) Sigma.mu<-diag((apply(X,2,max)-apply(X,2,min))^2)
    if (missing(Sigma.mu) & p == 1) Sigma.mu<-(max(X)-min(X))^2
    if (missing(Alpha.alpha)) Alpha.alpha<-2
    if (missing(Beta.alpha)) Beta.alpha<-4
    if (missing(Turnc.alpha)) Turnc.alpha<-0.25
    if (p > 1){
        xbar<-apply(X,2,mean)
        xsd<-apply(X,2,sd)
    }
    else if (p == 1){
        xbar<-mean(X)
        xsd<-sd(X)
    }
    # Family specific responses and offset terms
    if (family.indicator==1 | family.indicator==3 | family.indicator==5){
        Y <- model.response(mf, "any")
        offset <- as.vector(model.offset(mf))
        if (!is.null(offset)) {
           if (length(offset) != NROW(Y))
               stop(gettextf("number of offsets is %d should equal %d (number of observations)",
                   length(offset), NROW(Y)), domain = NA)
        }
        if (missing(offset)) offset<-rep(1,n)
    } else if (family.indicator==2 | family.indicator==4){
        Y1 <- model.response(mf, "any")
        if (NCOL(Y1)==1){
	        if (any(Y1 < 0 | Y1 > 1)) stop("y values must be 0 <= y <= 1")
	        offset <- array(1,n)
	        Y<-Y1
	    } else if (NCOL(Y1) == 2){
            offset <- Y1[, 1] + Y1[, 2]
	        Y<-Y1[,1]
	      } else
	         stop(paste("For the binomial family, y must be",
			             "a vector of 0 and 1's or a 2 column",
			             "matrix where col 1 is no. successes",
			             "and col 2 is no. failures"))
    }
    # Family specific prior parameters
    if (family.indicator==1){
        if (!missing(Alpha.xi)) if (!length(Alpha.xi)==1) stop(paste("For Poisson mixtures, argument Alpha.xi must be of length 1"))
        if (!missing(Alpha.xi)) if (Alpha.xi < 0) stop(paste("For Poisson mixtures, argument Alpha.xi must be positive"))
        if (missing(Alpha.xi)) Alpha.xi<-1.0
        if (!missing(Beta.xi)) if (!length(Beta.xi)==1) stop(paste("For Poisson mixtures, argument Beta.xi must be of length 1"))
        if (!missing(Beta.xi)) if (Beta.xi < 0) stop(paste("For Poisson mixtures, argument Beta.xi must be positive"))
        if (missing(Beta.xi)) Beta.xi<-0.1
    } else if (family.indicator==2){
        if (!missing(Alpha.xi)) if (!length(Alpha.xi)==1) stop(paste("For Binomial mixtures, argument Alpha.xi must be of length 1"))
        if (!missing(Alpha.xi)) if (Alpha.xi < 0) stop(paste("For Binomial mixtures, argument Alpha.xi must be positive"))
        if (missing(Alpha.xi)) Alpha.xi<-1.0
        if (!missing(Beta.xi)) if (!length(Beta.xi)==1) stop(paste("For Binomial mixtures, argument Beta.xi must be of length 1"))
        if (!missing(Beta.xi)) if (Beta.xi < 0) stop(paste("For Binomial mixtures, argument Beta.xi must be positive"))
        if (missing(Beta.xi)) Beta.xi<-1.0
    } else if (family.indicator==3){
        if (!missing(Alpha.xi)) if (!length(Alpha.xi)==2) stop(paste("For Negative Binimial mixtures, argument Alpha.xi must be of length 2"))
        if (!missing(Alpha.xi)) if (sum(Alpha.xi < 0) >0) stop(paste("For Negative Binimial mixtures, vector Alpha.xi must have positive elements"))
        if (missing(Alpha.xi)) Alpha.xi<-c(1.0,1.0)
        if (!missing(Beta.xi)) if (!length(Beta.xi)==2) stop(paste("For Negative Binimial mixtures, argument Beta.xi must be of length 2"))
        if (!missing(Beta.xi)) if (sum(Beta.xi < 0) >0) stop(paste("For Negative Binimial mixtures, vector Beta.xi must have positive elements"))
        if (missing(Beta.xi)) Beta.xi<-c(0.1,0.1)
    } else if (family.indicator==4){
        if (!missing(Alpha.xi)) if (!length(Alpha.xi)==2) stop(paste("For Beta Binimial mixtures, argument Alpha.xi must be of length 2"))
        if (!missing(Alpha.xi)) if (sum(Alpha.xi < 0) >0) stop(paste("For Beta Binimial mixtures, vector Alpha.xi must have positive elements"))
        if (missing(Alpha.xi)) Alpha.xi<-c(1.0,1.0)
        if (!missing(Beta.xi)) if (!length(Beta.xi)==2) stop(paste("For Beta Binimial mixtures, argument Beta.xi must be of length 2"))
        if (!missing(Beta.xi)) if (sum(Beta.xi < 0) >0) stop(paste("For Beta Binimial mixtures, vector Beta.xi must have positive elements"))
        if (missing(Beta.xi)) Beta.xi<-c(0.1,0.1)
    } else if (family.indicator==5){
        if (!missing(Alpha.xi)) if (!length(Alpha.xi)==2) stop(paste("For Generalized Poisson mixtures, argument Alpha.xi must be of length 2"))
        if (!missing(Alpha.xi)) if (sum(Alpha.xi < 0) >0) stop(paste("For Generalized Poisson mixtures, vector Alpha.xi must have positive elements"))
        if (missing(Alpha.xi)) Alpha.xi<-c(1.0,1.0)
        if (!missing(Beta.xi)) if (!length(Beta.xi)==2) stop(paste("For Generalized Poisson mixtures, argument Beta.xi must be of length 2"))
        if (!missing(Beta.xi)) if (sum(Beta.xi < 0) >0) stop(paste("For Generalized Poisson mixtures, vector Beta.xi must have positive elements"))
        if (missing(Beta.xi)) Beta.xi<-c(0.1,0.1)
    }
    #Precision
    if (!missing(prec)) if (length(prec)==1) prec<-c(prec,prec)
    if (missing(prec)) prec <- c(10,10)
    # Predictions
    if (missing(Xpred)){
        npred <- 0
        Xpred <- 1
    }
    else{
        npred <- length(Xpred)/p
    }
    if ((!missing(offsetPred)) & (npred > 0)){
           if (length(offsetPred) != npred)
               stop(gettextf("number of prediction offsets is %d, but should equal %d: the number of prediction scenarios",
                   length(offsetPred), npred), domain = NA)
    }
    if (missing(offsetPred) & (npred > 0)) offsetPred <- rep(1,npred)
    if (missing(offsetPred) & (npred == 0)) offsetPred <- 1.0
    meanReg <- array(0,npred)
    medianReg <- array(0,npred)
    q1Reg <- array(0,npred)
    q3Reg <- array(0,npred)
    modeReg <- array(0,npred)
    # Family specific predictions
    if (family.indicator==1 | family.indicator==3 | family.indicator==5){
        maxy <- 2000
    } else if (family.indicator==2 | family.indicator==4){
        maxy <- max(offsetPred)+1
    }
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
    paste(StorageDir,"BNSP.SigmahI.txt",sep=""), paste(StorageDir,"BNSP.nuh.txt",sep=""), paste(StorageDir,"BNSP.muh.txt",sep=""),
    paste(StorageDir,"BNSP.xih.txt",sep=""), paste(StorageDir,"BNSP.alpha.txt",sep=""),
    paste(StorageDir,"BNSP.compAlloc.txt",sep=""), paste(StorageDir,"BNSP.nmembers.txt",sep=""),
    paste(StorageDir,"BNSP.Updated.txt",sep="")))
    #Call C
    out<-.C("OneResLtnt", as.integer(seed), as.double(unlist(c(X))), as.integer(Y),
            as.double(offset), as.integer(sweeps), as.integer(burn), as.integer(thin), as.integer(ncomp),
            as.integer(n), as.integer(p),
            as.double(V), as.double(Vdf),
            as.double(Mu.nu), as.double(Sigma.nu),
            as.double(Mu.mu), as.double(Sigma.mu),
            as.double(Alpha.xi),as.double(Beta.xi),
            as.double(Alpha.alpha),as.double(Beta.alpha),as.double(Turnc.alpha),
            as.double(xbar), as.double(xsd), as.double(sum(Y)/sum(offset)), as.double(prec),
            as.integer(family.indicator), as.integer(sampler.indicator),
            as.integer(npred),as.double(Xpred),as.double(offsetPred),as.integer(maxy),
            as.double(meanReg),as.double(medianReg),as.double(q1Reg),as.double(q3Reg),as.double(modeReg),
            as.character(StorageDir),as.integer(WF))
    #Output
    location<-32
    meanReg <- out[[location+0]][1:npred]
    medianReg <- out[[location+1]][1:npred]
    q1Reg <- out[[location+2]][1:npred]
    q3Reg <- out[[location+3]][1:npred]
    modeReg <- out[[location+4]][1:npred]
    fit <- list(call=call,seed=seed,meanReg=meanReg,medianReg=medianReg,q1Reg=q1Reg,q3Reg=q3Reg,modeReg=modeReg)
    class(fit) <- 'bnp'
    return(fit)
}
