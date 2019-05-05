source("binary_tensor.R")
library(lattice)
d=20
rank=3
alpha=10
sigma=1

## simulate a binary tensor with a rank-3 signal tensor
simulate=low_rank(c(d,d,d),rank,alpha,sigma) ## simulate a rank-3 signal tensor with infinity norm bound alpha. 
ttensor=simulate$parameter ## extract the real-valued canonical tensor

Y=array(0,dim=rep(d,3))
for(i in 1:d){
    for(j in 1:d){
        for(k in 1:d){
            Y[i,j,k]=rbinom(1,1,inv.logit(ttensor[i,j,k]))
        }
    }
}
### generate Bernoulli tensor from logistic model

## fit the model with input data Y
option=1
BIC_select=BIC(Y,1:5,alpha,sigma,option,random.ini=FALSE,const=TRUE,nrand=3) ## select rank from BIC 
fit=binary_fit(Y,BIC_select$rank,alpha,sigma,option,random.ini=FALSE,const=TRUE,nrand=3) ## estimate the low-rank tensor from the binary observations. 

######################### main function for low-rank tensor estimation from binary observations #####
### Currently the algorithm only handels order-3 binary tensors with logit link. 
### Input: 
### Y -- input binary tensor.
### r -- the prespecified rank.
### alpha -- upper bound for the infinity norm.
### sigma: scale parameter in the logit distribution. By default set to 1, which corresponds to standard logit link. 
### option = 1: alternating glm algorithm; option = 2: MM algorithm. In general, MM algorithm is faster for large-scale binary tensors. 
### random.ini: random initilization (TRUE) or spectral initilization (FALSE)
### const: whether or not impose infinity norm bound on the parameter estimations. if const=TRUE, then the argument alpha will be ignored
### nrand: number of random initilizations in additional to the spectral initilization


### output: 
### A,B,C -- factor matrix along each of the three modes. Each factor matrix is of dimension d_k-by-r, where d_k is the dimension at mode k and r is the pre-specified rank
### para -- estimated canonical (real-valued) tensor from low-rank binary tensor model. 
#########################
plot(fit$para,ttensor) ## assess the estimation accuracy for the continous-valued tensor parameter
plot(inv.logit(fit$para),inv.logit(ttensor)) #assess the estimation accuracy for the Bernoulli probability 
error=sqrt(sum(ttensor-fit$para)^2/(d^3))
