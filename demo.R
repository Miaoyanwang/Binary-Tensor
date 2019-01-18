source("binary_tensor.R")
library(lattice)
d=20
rank=3
alpha=10
sigma=1

## simulate a binary tensor with a rank-3 signal tensor
simulate=low_rank(c(d,d,d),rank,alpha,sigma) ## simulate a rank-3 signal tensor with infinity norm bound alpha. 
ttensor=simulate$parameter ## extract the real-valued canonical tensor
Y=ttensor+array(rnorm(d^3,0,sigma),dim=c(d,d,d)) ## add Gaussian noise and generate the corresponding binary tensor 
Y=1*(Y>0)

option=1
fit=binary_fit(Y,rank,alpha,sigma,option) ## estimate the low-rank tensor from the binary observations. 

######################### main function for low-rank tensor estimation from binary observations #####
### Currently the algorithm only handels order-3 binary tensors with logit link. 
### Input: 
### Y -- input binary tensor.
### r -- the prespecified rank.
### alpha -- upper bound for the infinity norm.
### sigma -- scale parameter in the logit distribution. By default set to 1, which corresponds to standard logit link. 
### option = 1: alternating glm algorithm; option = 2: MM (minimization-majorization) algorithm. In general, MM algorithm is faster but less accurate...Preferable for large-scale binary tensors. 

### output: 
### A,B,C -- factor matrix along each of the three modes. Each factor matrix is of dimension d_k-by-r, where d_k is the dimension at mode k and r is the pre-specified rank
### para -- estimated canonical (real-valued) tensor from low-rank binary tensor model. 
#########################

error=sqrt(sum(ttensor-fit$para)^2/(d^3))
table(Y-(fit$para>0)) ## observed binary tensor Y and estimated binary tensor from low-rank model
