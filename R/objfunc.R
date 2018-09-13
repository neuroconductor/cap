objfunc <-
function(Y,X,gamma,beta)
{
  n<-length(Y)
  p<-ncol(Y[[1]])
  Tvec<-rep(NA,n)
  
  q<-ncol(X)
  
  # Estimate covariance matrix for each subject
  Sigma<-array(NA,c(p,p,n))
  for(i in 1:n)
  {
    Tvec[i]<-nrow(Y[[i]])
    
    Sigma[,,i]<-t(scale(Y[[i]],center=TRUE,scale=FALSE))%*%(scale(Y[[i]],center=TRUE,scale=FALSE))/nrow(Y[[i]])
  }
  
  Q1<-sum((X%*%beta)*Tvec)/2
  S<-matrix(0,p,p)
  for(i in 1:n)
  {
    S<-S+Sigma[,,i]*(Tvec[i]/(exp(t(X[i,])%*%beta)[1,1]))
  }
  Q2<-(t(gamma)%*%S%*%gamma/2)[1,1]
  
  re<-Q1+Q2
  
  return(re)
}
