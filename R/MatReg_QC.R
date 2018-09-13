MatReg_QC <-
function(Y,X,method=c("CAP","CAP-C","CAP-C1"),max.itr=1000,tol=1e-4,trace=FALSE,gamma0=NULL,score.return=TRUE)
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
  
  # common PCA based method
  if(method[1]=="CAP-C")
  {
    # find common eigenvectors and subject-specific eigenvalues
    Ymat<-NULL
    Group<-NULL
    for(i in 1:n)
    {
      Tvec[i]<-nrow(Y[[i]])
      
      Ymat<-rbind(Ymat,Y[[i]])
      
      Group<-c(Group,rep(i,Tvec[i]))
    }
    re.FCPCA<-FCPCA(Data=Ymat,Group=Group)
    
    # eigenvalues
    lambda<-re.FCPCA$lambda
    # common eigenvectors
    phi<-re.FCPCA$loadings.common
  }
  
  # SVD based method
  if(method[1]=="CAP")
  {
    Sigma.bar<-apply(Sigma,c(1,2),mean,na.rm=TRUE)
    Sbar.svd<-eigen(Sigma.bar)
    Ds<-diag(sqrt(Sbar.svd$values))
    Ds.inv<-diag(1/sqrt(Sbar.svd$values))
    Sbar.u<-Sbar.svd$vectors
    
    theta0<-gamma0
    if(is.null(theta0))
    {
      theta0<-rep(1/sqrt(p),p)
    }
    v0<-c(Sbar.u%*%Ds.inv%*%theta0)
  }else
  {
    v0<-gamma0
    if(is.null(v0))
    {
      if(method[1]=="CAP-C1")
      {
        v0<-rep(1/sqrt(p),p) 
      }else
        if(method[1]=="CAP-C")
        {
          v0<-phi[,p]
        }
    }
  }
  beta0<-rep(0,q)
  
  if(trace)
  {
    v.trace<-v0
    beta.trace<-beta0
    
    obj<-objfunc(Y=Y,X=X,gamma=v0,beta=beta0)
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
      
      Q1<-Q1+(Tvec[i]*(t(v0)%*%Sigma[,,i]%*%v0)[1,1]*exp(-t(X[i,])%*%beta0)[1,1])*(X[i,]%*%t(X[i,]))
      
      Q2<-Q2+Tvec[i]*(1-(t(v0)%*%Sigma[,,i]%*%v0)[1,1]*(exp(-t(X[i,])%*%beta0)[1,1]))*X[i,]
    }
    # beta.new<-beta0-solve(Q1)%*%Q2
    beta.new<-beta0-ginv(Q1)%*%Q2
    
    # update gamma
    if(method[1]=="CAP-C1")
    {
      Q3<-exp(X%*%beta.new)
      Q4<-matrix(0,p,p)
      for(i in 1:n)
      {
        Q4<-Q4+Sigma[,,i]*(Tvec[i]/Q3[i])
      }
      v.new<-svd(Q4)$u[,p]
    }else
      if(method[1]=="CAP-C")
      {
        nu<-apply(apply(apply(lambda,2,function(x){x/exp(X%*%beta.new)}),2,function(x){return(x*Tvec)}),2,sum)
        de<-apply(lambda,2,sum)
        cpc.idx<-which.min(nu/de)
        v.new<-phi[,cpc.idx]
      }else
        if(method[1]=="CAP")
        {
          V1<-matrix(0,p,p)
          for(i in 1:n)
          {
            V1<-V1+Sigma[,,i]*((Tvec*exp(-X%*%beta.new))[i])
          }
          svd.new<-eigen(Ds.inv%*%t(Sbar.u)%*%V1%*%Sbar.u%*%Ds.inv)
          theta.new<-svd.new$vectors[,p]
          v.new<-c(Sbar.u%*%Ds.inv%*%theta.new)
        }
        
    if(trace)
    {
      v.trace<-cbind(v.trace,v.new)
      beta.trace<-cbind(beta.trace,beta.new)
      
      obj<-c(obj,objfunc(Y=Y,X=X,gamma=v.new,beta=beta.new))
    }
    
    v.diff<-max(abs(v.new-v0))
    beta.diff<-max(abs(beta.new-beta0))
    
    diff<-max(c(v.diff,beta.diff))
    
    beta0<-beta.new
    v0<-v.new
  }
  
  if(v.new[1]<0)
  {
    v.new<--v.new
  }
  
  if(score.return)
  {
    score<-rep(NA,n)
    for(i in 1:n)
    {
      score[i]<-t(v.new)%*%Sigma[,,i]%*%v.new
    }
  }
  
  if(trace)
  {
    colnames(v.trace)<-NULL
    colnames(beta.trace)<-NULL
    
    if(score.return)
    {
      re<-list(gamma=c(v.new),beta=c(beta.new),convergence=(s<max.itr),score=score,gamma.trace=v.trace,beta.trace=beta.trace,obj=obj)
    }else
    {
      re<-list(gamma=c(v.new),beta=c(beta.new),convergence=(s<max.itr),gamma.trace=v.trace,beta.trace=beta.trace,obj=obj)
    }
    
  }else
  {
    if(score.return)
    {
      re<-list(gamma=c(v.new),beta=c(beta.new),convergence=(s<max.itr),score=score)
    }else
    {
      re<-list(gamma=c(v.new),beta=c(beta.new),convergence=(s<max.itr))
    }
  }
  
  return(re)
}
