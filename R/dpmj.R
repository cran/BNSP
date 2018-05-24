dpmj <- function(formula,Fcdf,data,offset,sampler="truncated",Xpred,offsetPred,
                 StorageDir,ncomp,sweeps,burn,thin=1,seed,
                 H,Hdf,d,D,Alpha.xi,Beta.xi,Alpha.alpha,Beta.alpha,Trunc.alpha,...){
    # Match call
    call <- match.call()
    # Seed
    if (missing(seed)) seed<-as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31)
    # Family
    Fcdf.indicator <- match(c(Fcdf),c("poisson","binomial","negative binomial","beta binomial","generalized poisson"))
    if (is.na(Fcdf.indicator)){
        stop('Fcdf must be character, and one of "poisson", "binomial", "negative binomial", "beta binomial", "generalized poisson"')
    }
    # Design matrix
    if (missing(data))
        data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    X<-as.matrix(model.matrix(formula,data=data)[,-1])
    # Dimensions
    n<-NROW(X)
    p<-NCOL(X)
    totran<-p+1
    # Binary covariates
    binary <- numeric(p)
    for (i in 1:p) binary[i] <- length(unique(X[,i]))
    binary <- as.numeric(binary == 2)
    NBC <- sum(binary) #number of binary covariates                
    NDV <- NBC + 1 #number of discrete variables
    NCC <- p - NBC #number of continuous covariates
    # Rearrange design matrix
    X <- as.matrix(X[,c(which(binary==1),which(binary==0))],ncol=p)
    # Prior parameters
      # Concentration parameter
        if (missing(Alpha.alpha)) Alpha.alpha<-2
        if (missing(Beta.alpha)) Beta.alpha<-5
        if (missing(Trunc.alpha)) Trunc.alpha<-0.25
      # mu_h ~ N(d,D)   
        if (missing(d)){ 
            d<-apply(as.matrix(X),2,mean)
            if (NBC > 0) d[1:NBC] <- -qnorm(1-d[1:NBC])
		}
        if (missing(D)){ 
            elements<-((apply(as.matrix(X),2,max)-apply(as.matrix(X),2,min))^2)/8
            if (NBC > 0) elements[1:NBC]<-5     
            D<-diag(elements,p)
        }
      # Wishart(h,H)
        if (missing(Hdf)) Hdf<-totran+2 
        if (missing(H)){ 
            elements<-((apply(as.matrix(X),2,max)-apply(as.matrix(X),2,min))^2)/8
            if (NBC > 0) elements[1:NBC]<-1     
            H<-diag(c(1,elements),totran)
        }
      # Family specific prior parameters
		if (Fcdf.indicator==1){
            if (!missing(Alpha.xi)) if (!length(Alpha.xi)==1) stop(paste("For Poisson mixtures, argument Alpha.xi must be of length 1"))
            if (!missing(Alpha.xi)) if (Alpha.xi < 0) stop(paste("For Poisson mixtures, argument Alpha.xi must be positive"))
            if (missing(Alpha.xi)) Alpha.xi<-1.0
            if (!missing(Beta.xi)) if (!length(Beta.xi)==1) stop(paste("For Poisson mixtures, argument Beta.xi must be of length 1"))
            if (!missing(Beta.xi)) if (Beta.xi < 0) stop(paste("For Poisson mixtures, argument Beta.xi must be positive"))
            if (missing(Beta.xi)) Beta.xi<-0.1
        } else if (Fcdf.indicator==2){
            if (!missing(Alpha.xi)) if (!length(Alpha.xi)==1) stop(paste("For Binomial mixtures, argument Alpha.xi must be of length 1"))
            if (!missing(Alpha.xi)) if (Alpha.xi < 0) stop(paste("For Binomial mixtures, argument Alpha.xi must be positive"))
            if (missing(Alpha.xi)) Alpha.xi<-1.0
            if (!missing(Beta.xi)) if (!length(Beta.xi)==1) stop(paste("For Binomial mixtures, argument Beta.xi must be of length 1"))
            if (!missing(Beta.xi)) if (Beta.xi < 0) stop(paste("For Binomial mixtures, argument Beta.xi must be positive"))
            if (missing(Beta.xi)) Beta.xi<-1.0
        } else if (Fcdf.indicator==3){
            if (!missing(Alpha.xi)) if (!length(Alpha.xi)==2) stop(paste("For Negative Binimial mixtures, argument Alpha.xi must be of length 2"))
            if (!missing(Alpha.xi)) if (sum(Alpha.xi < 0) >0) stop(paste("For Negative Binimial mixtures, vector Alpha.xi must have positive elements"))
            if (missing(Alpha.xi)) Alpha.xi<-c(1.0,1.0)
            if (!missing(Beta.xi)) if (!length(Beta.xi)==2) stop(paste("For Negative Binimial mixtures, argument Beta.xi must be of length 2"))
            if (!missing(Beta.xi)) if (sum(Beta.xi < 0) >0) stop(paste("For Negative Binimial mixtures, vector Beta.xi must have positive elements"))
            if (missing(Beta.xi)) Beta.xi<-c(0.1,0.1)
        } else if (Fcdf.indicator==4){
            if (!missing(Alpha.xi)) if (!length(Alpha.xi)==2) stop(paste("For Beta Binimial mixtures, argument Alpha.xi must be of length 2"))
            if (!missing(Alpha.xi)) if (sum(Alpha.xi < 0) >0) stop(paste("For Beta Binimial mixtures, vector Alpha.xi must have positive elements"))
            if (missing(Alpha.xi)) Alpha.xi<-c(1.0,1.0)
            if (!missing(Beta.xi)) if (!length(Beta.xi)==2) stop(paste("For Beta Binimial mixtures, argument Beta.xi must be of length 2"))
            if (!missing(Beta.xi)) if (sum(Beta.xi < 0) >0) stop(paste("For Beta Binimial mixtures, vector Beta.xi must have positive elements"))
            if (missing(Beta.xi)) Beta.xi<-c(0.1,0.1)
        } else if (Fcdf.indicator==5){
            if (!missing(Alpha.xi)) if (!length(Alpha.xi)==2) stop(paste("For Generalized Poisson mixtures, argument Alpha.xi must be of length 2"))
            if (!missing(Alpha.xi)) if (sum(Alpha.xi < 0) >0) stop(paste("For Generalized Poisson mixtures, vector Alpha.xi must have positive elements"))
            if (missing(Alpha.xi)) Alpha.xi<-c(1.0,1.0)
            if (!missing(Beta.xi)) if (!length(Beta.xi)==2) stop(paste("For Generalized Poisson mixtures, argument Beta.xi must be of length 2"))
            if (!missing(Beta.xi)) if (sum(Beta.xi < 0) >0) stop(paste("For Generalized Poisson mixtures, vector Beta.xi must have positive elements"))
            if (missing(Beta.xi)) Beta.xi<-c(0.1,1.0)
        } 
    # Sample mean and sd
    xbar<-apply(as.matrix(X),2,mean)
    xsd<-apply(as.matrix(X),2,sd)
    # Family specific responses and offset terms
    if (Fcdf.indicator==1 | Fcdf.indicator==3 | Fcdf.indicator==5){
        Y <- model.response(mf, "any")
        offset <- as.vector(model.offset(mf))
        if (!is.null(offset)){
           if (length(offset) != NROW(Y))
               stop(gettextf("number of offsets is %d should equal %d (number of observations)",
                   length(offset), NROW(Y)), domain = NA)
        }
        if (is.null(offset)) offset<-rep(1,n)
    } else if (Fcdf.indicator==2 | Fcdf.indicator==4){
        Y1 <- model.response(mf, "any")
        if (NCOL(Y1)==1){
	        if (any(Y1 < 0 | Y1 > 1)) stop("y values must be 0 <= y <= 1")
	        offset <- array(1,n)
	        Y<-Y1
	    } else if (NCOL(Y1) == 2){
            offset <- Y1[, 1] + Y1[, 2]
	        Y<-Y1[,1]
	      } else
	         stop(paste("For the binomial Fcdf, y must be",
			             "a vector of 0 and 1's or a 2 column",
			             "matrix where col 1 is no. successes",
			             "and col 2 is no. failures"))
    }
    # Predictions
    if (missing(Xpred)){
        npred <- 0
        Xpred <- 1
    } else{
        npred <- NROW(Xpred)
    }
    if (npred > 0){
           if (NCOL(Xpred) != p)
               stop(gettextf("Xpred matrix includes %d covariates, but it should include %d: the number covariates in the model",
                   NCOL(Xpred), p), domain = NA)
    }
    if ((!missing(offsetPred)) & (npred > 0)){
           if (length(offsetPred) != 1)
               stop(gettextf("the length of offsetPred %d, but it should be 1",
                   length(offsetPred)), domain = NA)
    }
    if (missing(offsetPred) & (npred > 0)) offsetPred <- round(mean(offset))
    if (missing(offsetPred) & (npred == 0)) offsetPred <- 1.0
    c1 <- cbind(rep(1,npred),!is.na(Xpred))
    Xpred[is.na(Xpred)]<-0
    meanReg <- array(0,npred)
    medianReg <- array(0,npred)
    q1Reg <- array(0,npred)
    q3Reg <- array(0,npred)
    modeReg <- array(0,npred)
    # Family specific predictions
    if (Fcdf.indicator==1 | Fcdf.indicator==3 | Fcdf.indicator==5){
        maxy <- max(Y)+1000 # maxy <- max(y)+max(10,floor(0.1*max(y)))
    } else if (Fcdf.indicator==2 | Fcdf.indicator==4){
        maxy <- max(offsetPred)+1
    }
    denReg <- array(0,npred*maxy)
    denVar <- array(0,npred*maxy)
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
    if (!missing(StorageDir)) if (!file.exists(StorageDir)) DIRC<-dir.create(StorageDir,recursive = TRUE)
    #if (!DIRC) stop("selected directory does not exist")
    if (missing(StorageDir)) {WF <- 0; StorageDir <- paste(getwd(),"/",sep="")}
    on.exit(if (WF==0) file.remove(paste(StorageDir,"BNSP.Sigmah.txt",sep=""),
    paste(StorageDir,"BNSP.muh.txt",sep=""), paste(StorageDir,"BNSP.xih.txt",sep=""),
    paste(StorageDir,"BNSP.alpha.txt",sep=""), paste(StorageDir,"BNSP.compAlloc.txt",sep=""),
    paste(StorageDir,"BNSP.nmembers.txt",sep=""), paste(StorageDir,"BNSP.Updated.txt",sep=""),
    paste(StorageDir,"BNSP.MeanReg.txt",sep=""), paste(StorageDir,"BNSP.MedianReg.txt",sep=""),
    paste(StorageDir,"BNSP.Q05Reg.txt",sep=""), paste(StorageDir,"BNSP.Q10Reg.txt",sep=""),
    paste(StorageDir,"BNSP.Q15Reg.txt",sep=""), paste(StorageDir,"BNSP.Q20Reg.txt",sep=""),
    paste(StorageDir,"BNSP.Q25Reg.txt",sep=""), paste(StorageDir,"BNSP.Q75Reg.txt",sep=""),
    paste(StorageDir,"BNSP.Q80Reg.txt",sep=""), paste(StorageDir,"BNSP.Q85Reg.txt",sep=""),
    paste(StorageDir,"BNSP.Q90Reg.txt",sep=""), paste(StorageDir,"BNSP.Q95Reg.txt",sep="")))
    #Call C
    out<-.C("OneResLtnt",as.integer(seed),as.double(unlist(c(X))),as.integer(Y),as.double(offset), 
            as.integer(sweeps), as.integer(burn), as.integer(thin), as.integer(ncomp),
            as.integer(n), as.integer(p), as.integer(NBC),
            as.double(H), as.double(Hdf),
            as.double(d), as.double(D),
            as.double(Alpha.xi),as.double(Beta.xi),
            as.double(Alpha.alpha),as.double(Beta.alpha),as.double(Trunc.alpha),
            as.double(xbar), as.double(xsd), as.double(sum(Y)/sum(offset)),
            as.integer(Fcdf.indicator), as.integer(sampler.indicator),
            as.integer(npred),as.double(Xpred),as.double(offsetPred),as.integer(maxy), as.integer(c1),
            as.double(meanReg),as.double(medianReg),as.double(q1Reg),as.double(q3Reg),as.double(modeReg),
            as.double(denReg),as.double(denVar),
            as.character(StorageDir),as.integer(WF))
    #Output
    location<-31
    meanReg <- out[[location+0]][1:npred]
    medianReg <- out[[location+1]][1:npred]
    q1Reg <- out[[location+2]][1:npred]
    q3Reg <- out[[location+3]][1:npred]
    modeReg <- out[[location+4]][1:npred]
    denReg <- matrix(out[[location+5]][1:(maxy*npred)],nrow=npred,ncol=maxy,byrow=TRUE)
    denVar <- matrix(out[[location+6]][1:(maxy*npred)],nrow=npred,ncol=maxy,byrow=TRUE)
    denVar <- denVar - denReg^2
    fit <- list(call=call,seed=seed,meanReg=meanReg,medianReg=medianReg,q1Reg=q1Reg,q3Reg=q3Reg,modeReg=modeReg,
                denReg=denReg,denVar=denVar)
    class(fit) <- 'bnp'
    return(fit)
}

print.bnp<-function(x,digits=max(3,getOption('digits')-3), ...){
   #Print formula
   cat('\nCall: ',deparse(x$call),'\n\n')
}
