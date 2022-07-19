#### Gibbs Sampler with normalization 

data <- read.csv("Lakechaddata.csv", na.strings = NA, header = TRUE)
data1 <-as.matrix(data)

V = as.matrix(data1[1:16,7]) 
V = scale(V)
X = as.matrix(data1[1:16,8:13])
X = X[,-4]
T = as.matrix(data1[1:16,15])
T = scale(T)
X = cbind(X,T)
X = scale(X)
b = cbind(rep(1,16))
X = cbind(b,X)

thn =100
n=10000
q=2 
r=0.1
k =length(V) 

alpha_initial = c(0,0,0,0,0,0,0)

I = diag(length(alpha_initial)) 
alpha = matrix(0,nrow=length(alpha_initial),ncol=n) 
sigma_2 = matrix(0,nrow=1,ncol=n) 

alpha[,1] =t(rmvnorm(1,mean=alpha_initial,sigma=I)) 

sigma_2[1,1] =rinvgamma(1,shape = q,rate = r) 

for (i in 2:n) {
  
  omega_alpha =solve(1/sigma_2[,i-1]*t(X)%*%X+diag(7))
  delta_alpha =1/sigma_2[,i-1]*t(V)%*%X+t(alpha[,(i-1)])%*%diag(7)
  temp =t(rmvnorm(1,mean=delta_alpha%*%omega_alpha,sigma=omega_alpha))
  alpha[,i] =temp[1:7,1]
  
  q = q + 0.5*k
  r = r + 0.5*crossprod(V - X%*%alpha[,i])
  sigma_2[i] = rinvgamma(1,shape = q,rate = r)
  
}


#### Gibbs Sampler without normalization 

data1 <-as.matrix(data)

V = as.matrix(data1[1:16,7]) 
X = as.matrix(data1[1:16,8:13])
X = X[,-4]
T = as.matrix(data1[1:16,15])
X = cbind(X,T)
b = cbind(rep(1,16))
X = cbind(b,X)

thn =100
n=10000
q=2 
r=0.1
k =length(V) 

alpha_initial = c(0,0,0,0,0,0,0)

I = diag(length(alpha_initial)) 
alpha = matrix(0,nrow=length(alpha_initial),ncol=n) 
sigma_2 = matrix(0,nrow=1,ncol=n) 

alpha[,1] =t(rmvnorm(1,mean=alpha_initial,sigma=I)) 

sigma_2[1,1] =rinvgamma(1,shape = q,rate = r) 

for (i in 2:n) {
  
  omega_alpha =solve(1/sigma_2[,i-1]*t(X)%*%X+diag(7))
  delta_alpha =1/sigma_2[,i-1]*t(V)%*%X+t(alpha[,(i-1)])%*%diag(7)
  temp =t(rmvnorm(1,mean=delta_alpha%*%omega_alpha,sigma=omega_alpha))
  alpha[,i] =temp[1:7,1]
  
  q = q + 0.5*k
  r = r + 0.5*crossprod(V - X%*%alpha[,i])
  sigma_2[i] = rinvgamma(1,shape = q,rate = r)
 
}

