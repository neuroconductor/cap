MatReg_QC_beta <-
function(Y,X,gamma,max.itr=1000,tol=1e-4,trace=FALSE,score.return=TRUE)
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
  
  if(score.return)
  {
    score<-rep(NA,n)
    for(i in 1:n)
    {
      score[i]<-t(gamma)%*%Sigma[,,i]%*%gamma
    }
  }
  
  beta0<-rep(0,q)
  
  if(trace)
  {
    beta.trace<-beta0
    
    obj<-objfunc(Y=Y,X=X,gamma=gamma,beta=beta0)
  }
  
  s<-0
  diff<-100
  while(s<=max.itr&diff>tol)
  {
    s<-s+1
    
    # update beta
    # Q1<-t(X)%*%diag(Tvec)%*%X
    Q1<-matrix(0,q,q)
    Q2<-rep(0,q)
    for(i in 1:n)
    {
      # Q1<-Q1+Tvec[i]*X[i,]%*%t(X[i,])
      
      Q1<-Q1+(Tvec[i]*(t(gamma)%*%Sigma[,,i]%*%gamma)[1,1]*exp(-t(X[i,])%*%beta0)[1,1])*(X[i,]%*%t(X[i,]))
      
      Q2<-Q2+Tvec[i]*(1-(t(gamma)%*%Sigma[,,i]%*%gamma)[1,1]*(exp(-t(X[i,])%*%beta0)[1,1]))*X[i,]
    }
    # beta.new<-beta0-solve(Q1)%*%Q2
    beta.new<-beta0-ginv(Q1)%*%Q2
    
    if(trace)
    {
      beta.trace<-cbind(beta.trace,beta.new)
      
      obj<-c(obj,objfunc(Y=Y,X=X,gamma=gamma,beta=beta.new))
    }
    
    diff<-max(abs(beta.new-beta0))
    
    beta0<-beta.new
    
    # print(diff)
  }
  
  if(trace)
  {
    if(score.return)
    {
      re<-list(gamma=c(gamma),beta=c(beta.new),convergence=(s<max.itr),score=score,beta.trace=beta.trace,obj=obj)  
    }else
    {
      re<-list(gamma=c(gamma),beta=c(beta.new),convergence=(s<max.itr),beta.trace=beta.trace,obj=obj)
    }
  }else
  {
    if(score.return)
    {
      re<-list(gamma=c(gamma),beta=c(beta.new),convergence=(s<max.itr),score=score) 
    }else
    {
      re<-list(gamma=c(gamma),beta=c(beta.new),convergence=(s<max.itr))
    }
  }
  
  return(re)
}
