chol <- function(x, mod = TRUE, p = 1, ...){
    S<-base::chol(x)
    list <- list(S)
    if (mod){
		list<-list()
		nblocks <- dim(x)[2] / p
		if (!nblocks == trunc(nblocks)) stop("dim of matrix and p non-conformable")  
        S<-t(S)
        T<-solve(S %*% diag(1/diag(S)))
        D<-T %*% x %*% t(T)
		if (p > 1){		    		    
            for (l in (0:(p-2))){
                E<-diag(1,dim(x)[2])
                indeces <- unlist(lapply(seq(0,(nblocks-1)),function(k) ((2+l) : p) + k * (dim(x)[2] + 1) * p + l * dim(x)[2]))
                E[indeces] <- -T[indeces]
		        T <- E %*% T 
                #D <- E %*% D %*% t(E)  
                D <- T %*% x %*% t(T)   		        		         
		    }
		}
		list<-list(L=T,D=D)
    }
    return(list) 
}
