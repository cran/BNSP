bnpglm <- function(formula,family=poisson,data,offset,sampler="slice",WorkingDir,ncomp,sweeps,burn,seed,
                   V,Vdf,Mu.nu,Sigma.nu,Mu.mu,Sigma.mu,Alpha.gamma,Beta.gamma,Alpha.alpha,Beta.alpha,Turnc.alpha,
                   Xpred,offsetPred,...){
    # Match call and family
    call <- match.call()
    if (is.character(family))
        family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family))
        family <- family()
    if (is.null(family$family)) {
        print(family)
        stop("'family' not recognized")
    }
    family.indicator <- match(family$family,c("poisson","binomial"))
    if (is.na(family.indicator)){
        stop("only the poisson and binomial families are supported for now")
    }
    # Data environment, design matrix
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
    # Family specific
    if (family.indicator==1){
        Y <- model.response(mf, "any")
        offset <- as.vector(model.offset(mf))
        if (!is.null(offset)) {
           if (length(offset) != NROW(Y))
               stop(gettextf("number of offsets is %d should equal %d (number of observations)",
                   length(offset), NROW(Y)), domain = NA)
        }
        if (missing(offset)) offset<-rep(1,n)
        if (missing(Alpha.gamma)) Alpha.gamma<-1.0
        if (missing(Beta.gamma)) Beta.gamma<-0.1
    } else if (family.indicator==2){
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
        if (missing(Alpha.gamma)) Alpha.gamma<-1.0
        if (missing(Beta.gamma)) Beta.gamma<-1.0
    }
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
    if (family.indicator==1){
        m.pred.mean <- max(offsetPred) * max(Y/offset)
        maxy <- round(m.pred.mean)
        while (dpois(maxy,m.pred.mean)>0.00000001) maxy<-maxy+1
    } else if (family.indicator==2){
        maxy <- max(offsetPred)+1
    }
    sampler.indicator <- match(sampler,c("slice","truncated"))
    if (is.na(sampler.indicator)){
        stop(c(sampler," 'sampler' not recognized"))
    }
    # Working directory
    WF <- 1
    if (!missing(WorkingDir)){
        WorkingDir <- path.expand(WorkingDir)
        ncharwd <- nchar(WorkingDir)}
    if (!missing(WorkingDir)) if (!(substr(WorkingDir,ncharwd,ncharwd)=="/")) WorkingDir <- paste(WorkingDir,"/",sep="")
    if (!missing(WorkingDir)) if (!file.exists(WorkingDir)) dir.create(WorkingDir)
    if (missing(WorkingDir)) {WF <- 0; WorkingDir <- paste(getwd(),"/",sep="")}
    on.exit(if (WF==0) file.remove(paste(WorkingDir,"BNSP.Th.txt",sep=""), paste(WorkingDir,"BNSP.Sigmah.txt",sep=""),
    paste(WorkingDir,"BNSP.SigmahI.txt",sep=""), paste(WorkingDir,"BNSP.nuh.txt",sep=""), paste(WorkingDir,"BNSP.muh.txt",sep=""),
    paste(WorkingDir,"BNSP.gammah.txt",sep=""), paste(WorkingDir,"BNSP.alpha.txt",sep=""),
    paste(WorkingDir,"BNSP.compAlloc.txt",sep=""), paste(WorkingDir,"BNSP.nmembers.txt",sep=""),
    paste(WorkingDir,"BNSP.Updated.txt",sep=""),paste(WorkingDir,"BNSP.PD.txt",sep="")))
    out<-.C("OneResLtnt", as.integer(seed), as.double(unlist(c(X))), as.integer(Y),
            as.double(offset), as.integer(sweeps), as.integer(burn), as.integer(ncomp),
            as.integer(n), as.integer(p),
            as.double(V), as.double(Vdf),
            as.double(Mu.nu), as.double(Sigma.nu),
            as.double(Mu.mu), as.double(Sigma.mu),
            as.double(Alpha.gamma),as.double(Beta.gamma),
            as.double(Alpha.alpha),as.double(Beta.alpha),as.double(Turnc.alpha),
            as.double(xbar), as.double(xsd), as.double(sum(Y)/sum(offset)),
            as.integer(family.indicator), as.integer(sampler.indicator),
            as.integer(npred),as.double(Xpred),as.double(offsetPred),as.integer(maxy),
            as.double(meanReg),as.double(medianReg),as.double(q1Reg),as.double(q3Reg),as.double(modeReg),
            as.character(WorkingDir),as.integer(WF))
    location<-30
    meanReg <- out[[location+0]][1:npred]
    medianReg <- out[[location+1]][1:npred]
    q1Reg <- out[[location+2]][1:npred]
    q3Reg <- out[[location+3]][1:npred]
    modeReg <- out[[location+4]][1:npred]
    fit <- list(call=call,seed=seed,meanReg=meanReg,medianReg=medianReg,q1Reg=q1Reg,q3Reg=q3Reg,modeReg=modeReg)
    class(fit) <- 'bnp'
    return(fit)
}
