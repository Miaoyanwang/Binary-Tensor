library(rTensor)
library(MASS)
library(Matrix)
library(rootSolve)
library(gtools)
######################### main function for low-rank tensor estimation from binary observations #####
### Currently the algorithm only handels order-3 binary tensors with logit link. 
### Input: Y -- input binary tensor.
### r: the prespecified rank.
### alpha: upper bound for the infinity norm.
### sigma: scale parameter in the logit distribution. By default set to 1, which corresponds to standard logit link. 
### option = 1: alternating glm algorithm; option = 2: MM algorithm. In general, MM algorithm is faster for large-scale binary tensors. 
### random.ini: random initilization (TRUE) or spectral initilization (FALSE)
### const: whether or not impose infinity norm bound on the parameter estimations. if const=TRUE, then the argument alpha will be ignored
### nrand: number of random initilization. If random.ini=FALSE, then this argment will be ignored.

### output: 
### A,B,C -- factor matrices along each of the three modes. Each factor matrix is of dimension d_k-by-r, where d_k is the dimension at mode k and r is the pre-specified rank
### para -- estimated canonical (real-valued) tensor from multilinear glm model. 
#########################
binary_fit=function(Y,r,alpha=10,sigma=1,option=2,random.ini=FALSE,const=TRUE,nrand=5){
    Y0=Y
    missing=which(is.na(Y)==T)
    nonmissing=which(is.na(Y)==F)
    Y[missing]=mean(Y,na.rm=T)
    
    Z=2*Y-1
    W=alpha/sigma*Z
    
    if(random.ini==TRUE){
        
        A=normalize(matrix(rnorm(dim(Y)[1]*r,0,1),nrow=dim(Y)[1]))
        B=normalize(matrix(rnorm(dim(Y)[2]*r,0,1),nrow=dim(Y)[2]))
          
        Predictor=KhatriRao(A,B)
        C=unfold(Z,3)%*%ginv(t(as.matrix(Predictor)))
        M=tensorize(A,B,C)
        
        if((max(abs(M))>alpha/sigma)&(const==TRUE)){
            C=C/max(abs(M))*(alpha/sigma)
        }
        
        start=list("A"=A,"B"=B,"C"=C)
        
        fit0=binary_fit_fix(Y0,r,alpha/sigma,start,option,const)
        cost=fit0$cost
        
        est=tensorize(fit0$A,fit0$B,fit0$C)
        return(list("A"=fit0$A,"B"=fit0$B,"C"=fit0$C*sigma,"cost"=fit0$cost,"para"=est))
        
    }
    
    
    if(r<=min(dim(Y))){ ## if rank is less than the minimum dimension
        
        ## spectral initilization from unfoldings
        A=as.matrix(svd(unfold(W,1))$u[,1:r])
        B=as.matrix(svd(unfold(W,2))$u[,1:r])
        Predictor=KhatriRao(A,B)
        C=unfold(Z,3)%*%ginv(t(as.matrix(Predictor)))
        M=tensorize(A,B,C)
        
      
        if((max(abs(M))>alpha/sigma)&(const==TRUE)){
            C=C/max(abs(M))*(alpha/sigma)
        }
        
        start=list("A"=A,"B"=B,"C"=C)
        ## binary tensor decomposition
        fit0=binary_fit_fix(Y0,r,alpha/sigma,start,option,const)
        cost0=fit0$cost
        
        ## random initilization
        for(n in 1:nrand){
            start=CP_decomp(W,r)
            start=list("A"=start$A,"B"=start$B,"C"=start$C)
            M=tensorize(start$A,start$B,start$C)
            if((max(abs(M))>alpha/sigma)&(const==TRUE)){
                start$C=start$C/max(abs(M))*(alpha/sigma)
            }
            ## binary tensor decomposition
            fit=binary_fit_fix(Y0,r,alpha/sigma,start,option,const)
            cost=fit$cost
            if(max(cost)>max(cost0)){
                fit0=fit
                cost0=cost
            }
        }
        
    }
    else{ ## rank exceeds the minimum dimension
        
        ## spectral initilization from real-valued tensor decomposition
        start=CP_decomp(W,r)
        start=list("A"=start$A,"B"=start$B,"C"=start$C)
        M=tensorize(start$A,start$B,start$C)
        if((max(abs(M))>alpha/sigma)&(const==TRUE)){
            start$C=start$C/max(abs(M))*(alpha/sigma)
        }
        ## binary tensor decomposition
        fit0=binary_fit_fix(Y0,r,alpha/sigma,start,option,const)
        cost0=fit0$cost
        
        ## random initilization
        for(n in 1:nrand){
            start=CP_decomp(W,r)
            start=list("A"=start$A,"B"=start$B,"C"=start$C)
            ## binary tensor decomposition
            fit=binary_fit_fix(Y0,r,alpha/sigma,start,const)
            cost=fit$cost
            if((max(abs(M))>alpha/sigma)&(const==TRUE)){
                fit0=fit
                cost0=cost
            }
        }
    }
    
    est=tensorize(fit0$A,fit0$B,fit0$C)
    ## output
    if(sigma==0){
        return(list("A"=fit0$A,"B"=fit0$B,"C"=fit0$C,"cost"=fit0$cost,"para"=est))
    }
    else
    return(list("A"=fit0$A,"B"=fit0$B,"C"=fit0$C*sigma,"cost"=fit0$cost,"para"=est))
}


######### subroutine for binary tensor decomposition with initilization start. 
## Option=1 alternating glm; option=2 miminization majorization
binary_fit_fix=function(Y,r,alpha,start,option=2,const=TRUE){
    missing=which(is.na(Y)==T)
    nonmissing=which(is.na(Y)==F)
    Y[missing]=mean(Y,na.rm=T)
    
    Z=2*Y-1
    W=alpha*Z
    bound=NULL
    
    A=start$A
    B=start$B
    C=start$C
    
    M=tensorize(A,B,C)
    cost=objective(Z[nonmissing],M[nonmissing])
    cost_new=cost*2
    cost_traj=cost
    

    d=dim(Y)
    nloop=0
    
    while((abs((cost_new-cost)/cost)>=0.001)&(nloop<50)){ ## modified for simulation
    ##while(nloop<14){
        nloop=nloop+1
        A0=A
        B0=B
        C0=C
        cost=cost_new
        
        if(option==1){ ## alternating glm
        Predictor=KhatriRao(C,B)
        for(n in 1:d[1]){
            res=glm_modify(unfold(Y,1)[n,],as.matrix(Predictor),start=A[n,])
            A[n,]=res
            #res=glm(unfold(Y,1)[n,]~-1+as.matrix(Predictor),family=binomial(link='logit'))
            #A[n,]=coef(res)
        }
        
        Predictor=KhatriRao(A,C)
        for(n in 1:d[2]){
            res=glm_modify(unfold(Y,2)[n,],as.matrix(Predictor),start=B[n,])
            B[n,]=res
            #res=glm(unfold(Y,2)[n,]~-1+as.matrix(Predictor),family=binomial(link='logit'))
            #B[n,]=coef(res)
        }
        
        Predictor=KhatriRao(A,B)
        for(n in 1:d[3]){
            res=glm_modify(unfold(Y,3)[n,],as.matrix(Predictor),start=C[n,])
            C[n,]=res
            #res=glm(unfold(Y,3)[n,]~-1+as.matrix(Predictor),family=binomial(link='logit'))
            #C[n,]=coef(res)
        }
        
        res=rescale(A,B,C)
        A=res$A
        B=res$B
        C=res$C
        }
        else if (option==2){ ## minimization majorization
            H=h(Z*M) ## pointwise
            W=M+4*Z*H ## pointwise logistic 
            W[missing]=M[missing]
            
            Predictor=KhatriRao(C,B)
            A=unfold(W,1)%*%ginv(t(as.matrix(Predictor)))
            A=normalize(A)
            
            Predictor=KhatriRao(A,C)
            B=unfold(W,2)%*%ginv(t(as.matrix(Predictor)))
            B=normalize(B)
            
            Predictor=KhatriRao(A,B)
            C=unfold(W,3)%*%ginv(t(as.matrix(Predictor)))
        }
        
        M=tensorize(A,B,C)
        cost_new=objective(Z[nonmissing],M[nonmissing])
        
        if(is.na(max(abs(M)))==T){
            if(is.na(cost_new)==T)
            cost_new=-Inf
            return(list("M"=M,"A"=A,"B"=B,"C"=C,"cost"=cost_traj))
        }
        
     
        if((max(abs(M))>alpha)&(const==TRUE)){ ## project into the feasible space
            print("project back to the feasible space")
            fun_upper=function(x)constraint_max(x,A,B,C,A0,B0,C0,alpha)
            fun_lower=function(x)constraint_min(x,A,B,C,A0,B0,C0,alpha)
            find_root=uniroot.all(fun_upper,c(-10^(-4),1))
            upper=find_root[1]
            find_root=uniroot.all(fun_lower,c(-10^(-4),1))
            lower=find_root[1]
            
            if(is.na(upper)==T) upper=0
            if(is.na(lower)==T) lower=0
            
            bound=c(bound,min(upper,lower))
            if(min(upper,lower)>0){
                temp=optimize(function(x)search_new(x,A,B,C,Z,A0,B0,C0,nonmissing),lower=-10^(-4),upper=min(upper,lower),maximum=TRUE)
                if(temp$maximum<0) gamma=0
                else gamma=temp$maximum
            }
            else gamma=0
            
            A=gamma*A+(1-gamma)*A0
            B=gamma*B+(1-gamma)*B0
            C=gamma*C+(1-gamma)*C0
            res=rescale(A,B,C)
            A=res$A
            B=res$B
            C=res$C
            
            M=tensorize(A,B,C)
            
            cost_new=objective(Z[nonmissing],M[nonmissing])
        }
        
        cost_traj=c(cost_traj,cost_new) 
    }
    
    scale=apply(C,2,function(x)sqrt(sum(x^2)))
    ind=sort(scale,index=T,decreasing=T)$ix
    A=A[,ind]
    B=B[,ind]
    C=C[,ind]
    
    return(list("A"=A,"B"=B,"C"=C,"cost"=cost_traj,"bound"=bound))
}


glm_modify=function(y,x,start){
    
    ## initial coefficent
    ini_loglik=sum(log(inv.logit((2*y-1)*(x%*%start))))
    
    ## Option 1: glm fittig with default initilization
    fit1 = glm(y~-1+x,family=binomial(link="logit"),control = list(maxit = 50))
    ##return(coef(fit1))
    ## Option 2: glm with user specified initilization
    fit2= glm(y~-1+x,family=binomial(link="logit"),control = list(maxit = 50),start=start)
    
    ## report the result whichever gives the highest likelihood
    value1=logLik(fit1)
    value2=logLik(fit2)
    if(is.na(value1)) value1=-Inf
    if(is.na(value2)) value2=-Inf
    
    if(max(value1,value2)<ini_loglik) return (start)
    else if(value1>value2) return(coef(fit1))
    else return(coef(fit2))
}


##################### construct CP tensor using factor matrices X, Y, Z ###################################
tensorize=function(X,Y,Z){
    r=dim(X)[2]
    tensor=0
    if(is.matrix(X)==0){
        tensor=X%o%Y%o%Z
        return(tensor)
    }
    
    for(i in 1:r){
        tensor=tensor+X[,i]%o%Y[,i]%o%Z[,i]
    }
    return(tensor)
}


#####################  unfold a tensor along a given mode ###################################

unfold=function(tensor,mode=1){
    d1=dim(tensor)[1]
    d2=dim(tensor)[2]
    d3=dim(tensor)[3]
    
    
    if(mode==1){
        unfold=array(tensor,dim=c(d1,d2*d3))
    }
    
    if(mode==2){
        unfold=NULL
        for(i in 1:d2){
            unfold=cbind(unfold,c(t(tensor[,i,])))
        }
        unfold=t(unfold)
    }
    
    if (mode==3){
        unfold=NULL
        for(i in 1:d3){
            unfold=cbind(unfold,c(t(tensor[,,i])))
        }
        unfold=t(unfold)
        
    }
    
    return(unfold)
}

tensordensity=function(p1,p2,p3){
    r=dim(p1)[2]
    prop=array(1,dim=c(dim(p1)[1],dim(p2)[1],dim(p3)[1]))
    for(i in 1:r){
        prob0=1-tensorize(p1[,i],p2[,i],p3[,i])
        prop=prob0*prop
    }
    return(1-prop)
}
randomtensor=function(E){
    Y=E
    d=dim(E)
    for(i in 1:d[1]){
        for(j in 1:d[2]){
            for(k in 1:d[3]){
                Y[i,j,k]=rbinom(1,1,E[i,j,k])   
            }
        }
    }
    return(Y)
}


mismatch=function(tensor1,tensor2){
    tensor1[tensor1>=0.5]=1
    tensor1[tensor1<0.5]=0
    
    tensor2[tensor2>=0.5]=1
    tensor2[tensor2<0.5]=0
    return(mean(abs(tensor1-tensor2)))    
}

BIC=function(Y,rankrange,alpha,sigma,option,random.ini=TRUE,const=FALSE,nrand=0){
    BIC_list=logL=pe=NULL
    d=dim(Y)
    for (rank in rankrange){
        print(sprintf("fitting binary tensor with rank %d",rank))
        fit=binary_fit(Y,rank,alpha,sigma,option,random.ini,const,nrand)
        logL=c(logL,-2*objective(2*Y-1,fit$para))
        pe=c(pe,rank*(sum(d)-length(d)+1))
        BIC_list=c(BIC_list,-2*objective(2*Y-1,fit$para)+(rank*(sum(d)-length(d)+1)*log10(prod(d))))
    }
    return(list("BIC"=BIC_list,"log"=logL,"df"=pe*log10(prod(d)),"rank"=rankrange[which.min(BIC_list)],"AIC"=rankrange[which.min(logL+pe*2)]))
}

## classical CP decomposition
CP_decomp=function(Z,r){
    start=cp(as.tensor(Z),r)
    A=normalize(start$U[[1]])
    B=normalize(start$U[[2]])
    AtimesB=KhatriRao(A,B)
    C=unfold(Z,3)%*%ginv(t(as.matrix(AtimesB)))
    
    return(list("A"=A,"B"=B,"C"=C))
}

## objective function (proportional to log-likelihood) 
objective=function(Z,M){
    obj=sum(Fun(Z*M))
    if(obj<(-10^(10))) return(-10^(10))
    else if(obj>(10^(10))) return(10^(10)) 
    else return(obj)
}


## propotional to log-likelihood 
Fun=function(x,option="logit",sigma=1){
    if(option=="logit"){ ## logit link
        x[x>500]=500 ## threshold the entries
        return(log10(exp(x)/(1+exp(x))))
    }
    else if(option=="probit"){ ## probit link
        return(log10(pnorm(x,0,sigma)))
    }
}


################################### normalize each column of X to be unit-one. ###################################
normalize=function(X){
    d=dim(X)[2]
    for(i in 1:d){
        X[,i]=X[,i]/sqrt(sum(X[,i]^2))
    }
    return(X)
}

## rescale each column of the factor matrices to be unit-one. 
rescale=function(A,B,C){
    d=dim(A)[2]
    for(i in 1:d){
        C[,i]=C[,i]*sqrt(sum(A[,i]^2))*sqrt(sum(B[,i]^2))
        A[,i]=A[,i]/sqrt(sum(A[,i]^2))
        B[,i]=B[,i]/sqrt(sum(B[,i]^2))
        
    }
    return(list("A"=A,"B"=B,"C"=C))     
}


## find the maximum among the positive tensor entries
constraint_max=function(gamma,A,B,C,A0,B0,C0,alpha){
    nvector=length(gamma)
    res=NULL
    for(k in 1:nvector){
        M=tensorize(((1-gamma[k])*A0+gamma[k]*A),((1-gamma[k])*B0+gamma[k]*B),((1-gamma[k])*C0+gamma[k]*C))-alpha
        res=c(res,max(M))
    }
    return(res)
}

## find the minimum among the negative tensor entries
constraint_min=function(gamma,A,B,C,A0,B0,C0,alpha){
    nvector=length(gamma)
    res=NULL
    for(k in 1:nvector){
        M=-tensorize(((1-gamma[k])*A0+gamma[k]*A),((1-gamma[k])*B0+gamma[k]*B),((1-gamma[k])*C0+gamma[k]*C))-alpha
        res=c(res,max(M))
    }
    return(res)
}


## specify the cost function
search_new=function(alpha,A,B,C,Z,A0,B0,C0,nonmissing){
    M=tensorize(((1-alpha)*A0+alpha*A),((1-alpha)*B0+alpha*B),((1-alpha)*C0+alpha*C))
    return(objective(M[nonmissing],Z[nonmissing]))
}

##################### auxillary function for binary-valued tensor decomposition ##################### 
h=function(x){
    return(1/(1+exp(x)))
}


### simulate low-rank binary tensor
## input: d 
low_rank=function(d,rank,alpha,sigma=1){
    factor=list()
    for(i in 1:length(d)){
        factor[[i]]=matrix(runif(d[i]*rank,-1,1),nrow=d[i],ncol=rank)
    }
    parameter=tensorize(factor[[1]],factor[[2]],factor[[3]])
    scale=max(abs(parameter))/alpha
    
    parameter=parameter/scale
    A=factor[[1]]
    B=factor[[2]]
    C=factor[[3]]/scale
    
    if(rank>=2){
        scale_A=apply(A,2,function(x) sqrt(sum(x^2)))
        scale_B=apply(B,2,function(x) sqrt(sum(x^2)))
        scale_C=apply(C,2,function(x) sqrt(sum(x^2)))
        scale=scale_A*scale_B*scale_C
        
        
        ind=sort(scale,decreasing=T,index=T)$ix
        A=normalize(A)[,ind]
        B=normalize(B)[,ind]
        C=normalize(C)[,ind]%*%diag(scale[ind])
        parameter=tensorize(A,B,C)
    }
    
    dtotal=d[1]*d[2]*d[3]
    data=parameter+array(rnorm(dtotal,0,sigma),dim=d)
    data=1*(data>0)
    
    return(list("Y"=data,"A"=A,"B"=B,"C"=C,"parameter"=parameter))
}
######################################## 

error_bound=function(tensor1,tensor2){
    error=sqrt(sum((tensor1-tensor2)^2))
    return(error)
}

block=function(d,nc,option){
    mu=array(0,dim=rep(nc,3))
    if(option=="add"){
        x=runif(nc,-1,1)
        y=runif(nc,-1,1)
        z=runif(nc,-1,1)
        mu=x%o%rep(1,nc)%o%rep(1,nc)+rep(1,nc)%o%y%o%rep(1,nc)+rep(1,nc)%o%rep(1,nc)%o%z
        return(pnorm(mu,0,1))
    }
    if(option=="mul"){
        x=runif(nc,-1,1)
        y=runif(nc,-1,1)
        z=runif(nc,-1,1)
        mu=x%o%y%o%z
        return(pnorm(mu,0,1))
    }
    if(option=="com"){
        mu=array(runif(nc^3,-1,1),dim=rep(nc,3))
        return(pnorm(mu,0,1))
    }
    
}
