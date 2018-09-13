MatReg_QC_RE <-
function(Y,X,Phi0,max.itr=1000,tol=1e-4,trace=FALSE,score.return=TRUE)
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
  
  beta0<-rep(0,q)
  
  if(trace)
  {
    gamma.trace<-NULL
    beta.trace<-beta0
  }
  
  s<-0
  diff<-100
  while(s<=max.itr&diff>tol)
  {
    s<-s+1
    
    # update gamma
    A<-matrix(0,p,p)
    for(i in 1:n)
    {
      A<-A+Tvec[i]*Sigma[,,i]/(exp(X[i,]%*%beta0)[1,1])
    }
    B<-apply(Sigma,c(1,2),mean)
    B.inv<-solve(B)
    P<-Phi0%*%ginv(t(Phi0)%*%B.inv%*%Phi0)%*%t(Phi0)%*%B.inv
    B.inv.eigen<-eigen(B.inv)
    B.inv.rt<-B.inv.eigen$vectors%*%diag(sqrt(B.inv.eigen$values))%*%solve(B.inv.eigen$vectors)
    re.eigen<-eigen(B.inv.rt%*%(diag(rep(1,p))-P)%*%A%*%B.inv.rt)
    re.eigen.vec<-Re(re.eigen$vectors)
    
    x.tmp<-B.inv.rt%*%re.eigen.vec
    # diag(t(x.tmp)%*%B%*%x.tmp)
    
    # which.min(diag(t(x.tmp)%*%A%*%x.tmp))
    # gamma.new<-x.tmp[,which.min(diag(t(x.tmp)%*%A%*%x.tmp))]
    
    gamma.new.mat<-x.tmp
    
    # update beta
    # Q1<-t(X)%*%diag(Tvec)%*%X
    beta.new.mat<-NULL
    objfunc.tmp<-NULL
    for(j in 1:ncol(gamma.new.mat))
    {
      gamma.new<-gamma.new.mat[,j]
      
      Q1<-matrix(0,q,q)
      Q2<-rep(0,q)
      for(i in 1:n)
      {
        # Q1<-Q1+Tvec[i]*X[i,]%*%t(X[i,])
        
        Q1<-Q1+(Tvec[i]*(t(gamma.new)%*%Sigma[,,i]%*%gamma.new)[1,1]*exp(-t(X[i,])%*%beta0)[1,1])*(X[i,]%*%t(X[i,]))
        
        Q2<-Q2+Tvec[i]*(1-(t(gamma.new)%*%Sigma[,,i]%*%gamma.new)[1,1]*(exp(-t(X[i,])%*%beta0)[1,1]))*X[i,]
      }
      beta.new.mat<-cbind(beta.new.mat,beta0-ginv(Q1)%*%Q2)
      
      objfunc.tmp<-c(objfunc.tmp,objfunc(Y,X,gamma.new,beta.new.mat[,j]))
    }
    beta.new<-beta.new.mat[,which.min(objfunc.tmp)]
    gamma.new<-gamma.new.mat[,which.min(objfunc.tmp)]
    
    if(trace)
    {
      gamma.trace<-cbind(gamma.trace,gamma.new)
      beta.trace<-cbind(beta.trace,beta.new)
    }
    
    diff<-max(abs(beta.new-beta0))
    
    beta0<-beta.new
    
    # print(diff)
  }
  
  # scale gamma
  gamma.new<-gamma.new/sqrt(sum(gamma.new^2))
  if(gamma.new[1]<0)
  {
    gamma.new<--gamma.new
  }
  beta.new<-MatReg_QC_beta(Y,X,gamma=gamma.new,max.itr=max.itr,tol=tol,trace=trace,score.return=score.return)$beta
  
  if(score.return)
  {
    score<-rep(NA,n)
    for(i in 1:n)
    {
      score[i]<-t(gamma.new)%*%Sigma[,,i]%*%gamma.new
    }
  }
  
  if(trace)
  {
    colnames(v.trace)<-NULL
    colnames(beta.trace)<-NULL
    
    if(score.return)
    {
      re<-list(gamma=c(gamma.new),beta=c(beta.new),convergence=(s<max.itr),score=score,gamma.trace=gamma.trace,beta.trace=beta.trace)
    }else
    {
      re<-list(gamma=c(gamma.new),beta=c(beta.new),convergence=(s<max.itr),gamma.trace=gamma.trace,beta.trace=beta.trace)
    }
    
  }else
  {
    if(score.return)
    {
      re<-list(gamma=c(gamma.new),beta=c(beta.new),convergence=(s<max.itr),score=score)
    }else
    {
      re<-list(gamma=c(gamma.new),beta=c(beta.new),convergence=(s<max.itr))
    }
  }
  
  return(re)
}
