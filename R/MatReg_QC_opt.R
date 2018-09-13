MatReg_QC_opt <-
function(Y,X,method=c("CAP","CAP-C","CAP-C1"),max.itr=1000,tol=1e-4,trace=FALSE,score.return=TRUE,gamma0.mat=NULL,ninitial=NULL)
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
  
  # set initial values
  if(is.null(gamma0.mat))
  {
    gamma0.mat<-matrix(NA,p,p+1+5)
    for(j in 1:p)
    {
      gamma0.mat[,j]<-rep(0,p)
      gamma0.mat[j,j]<-1
    }
    gamma0.mat[,p+1]<-rep(1,p)/sqrt(sum(rep(1,p)^2))
    
    set.seed(500)
    gamma.tmp<-matrix(rnorm(5*p,mean=0,sd=1),nrow=p)
    gamma0.mat[,(p+2):(p+1+5)]<-apply(gamma.tmp,2,function(x){return(x/sqrt(sum(x^2)))})
  }
  if(is.null(ninitial))
  {
    ninitial<-min(ncol(gamma0.mat),10)
  }else
  {
    if(ninitial>ncol(gamma0.mat))
    {
      ninitial<-ncol(gamma0.mat)
    }
  }
  set.seed(500)
  gamma0.mat<-matrix(gamma0.mat[,sort(sample(1:ncol(gamma0.mat),ninitial,replace=FALSE))],ncol=ninitial)
  
  if(method[1]=="CAP-C1")
  {
    gamma0.mat<-apply(gamma0.mat,2,function(x){return(x/sqrt(sum(x^2)))})
    
    re<-vector("list",ncol(gamma0.mat))
    obj.func<-rep(NA,ncol(gamma0.mat))
    for(kk in 1:ncol(gamma0.mat))
    {
      try(re[[kk]]<-MatReg_QC(Y,X,method=method[1],max.itr=max.itr,tol=tol,trace=trace,gamma0=gamma0.mat[,kk],score.return=score.return))
      
      try(obj.func[kk]<-objfunc(Y,X,re[[kk]]$gamma,re[[kk]]$beta))
    }
    
    opt.idx<-which.min(obj.func)
    re.opt<-re[[opt.idx]]
    
    if(method[1]=="CAP-C")
    {
      dis<-rep(NA,1,p)
      for(j in 1:p)
      {
        if(phi[1,j]<0)
        {
          phi.tmp<--phi[,j]
        }else
        {
          phi.tmp<-phi[,j]
        }
        dis[j]<-sqrt(sum(re.opt$gamma-phi.tmp)^2)
      }
      re.opt$PC.idx<-which.min(dis)
    }
    
    return(re.opt)
  }
  if(method[1]=="CAP")
  {
    theta0.mat<-gamma0.mat
    
    re<-vector("list",ncol(gamma0.mat))
    obj.func<-rep(NA,ncol(gamma0.mat))
    
    re.scale<-vector("list",ncol(gamma0.mat))
    obj.func.scale<-rep(NA,ncol(gamma0.mat))
    
    for(kk in 1:ncol(gamma0.mat))
    {
      try(re[[kk]]<-MatReg_QC(Y,X,method="CAP",max.itr=max.itr,tol=tol,trace=trace,gamma0=theta0.mat[,kk],score.return=score.return))
      
      try(obj.func[kk]<-objfunc(Y,X,re[[kk]]$gamma,re[[kk]]$beta))
      
      try(re.scale[[kk]]<-MatReg_QC_beta(Y,X,gamma=re[[kk]]$gamma/sqrt(sum((re[[kk]]$gamma)^2)),max.itr=max.itr,tol=tol,trace=trace,score.return=score.return))
      try(obj.func.scale[kk]<-objfunc(Y,X,re.scale[[kk]]$gamma,re.scale[[kk]]$beta))
    }
    
    opt.idx<-which.min(obj.func)
    opt.idx.scale<-which.min(obj.func.scale)
    
    # re.opt<-list(unscale=re[[opt.idx]],scale=re.scale[[opt.idx]])
    re.opt<-re.scale[[opt.idx]]
    
    return(re.opt)
  }
  if(method[1]=="CAP-C")
  {
    optmat<-matrix(NA,p,p)
    colnames(optmat)<-paste0("Dim",1:p)
    rownames(optmat)<-paste0("BetaDim",1:p)
    beta.est<-matrix(NA,q,p)
    colnames(beta.est)<-paste0("Dim",1:p)
    rownames(beta.est)<-colnames(X)
    for(j in 1:p)
    {
      beta.tmp<-MatReg_QC_beta(Y,X,gamma=phi[,j])$beta
      optmat[j,]<-(apply(lambda,2,function(x){return(sum(x*Tvec*exp(-X%*%beta.tmp)))})/apply(lambda,2,sum))*n/2+sum(Tvec*X%*%beta.tmp)/2
      
      beta.est[,j]<-beta.tmp
    }
    
    min.idx<-apply(optmat,1,which.min)
    sidx<-which(apply(cbind(min.idx,1:p),1,function(x){x[1]==x[2]})==TRUE)
    if(length(sidx)>0)
    {
      svar<-rep(NA,length(sidx))
      for(j in 1:length(sidx))
      {
        svar[j]<-sum(exp(X%*%beta.est[,sidx[j]]))
      }
      opt.idx<-sidx[which.max(svar)]
      beta.opt<-beta.est[,opt.idx]
      gamma.opt<-phi[,opt.idx]
      if(gamma.opt[1]<0)
      {
        gamma.opt<--gamma.opt
      }
      
      if(score.return)
      {
        re.opt<-list(gamma=gamma.opt,beta=beta.opt,PC.index=opt.idx,score=c(lambda[,opt.idx]),minmix=TRUE)
      }else
      {
        re.opt<-list(gamma=gamma.opt,beta=beta.opt,PC.index=opt.idx,minmix=TRUE)
      }
    }else
    {
      optvec<-rep(NA,p)
      for(j in 1:p)
      {
        optvec[j]<-optmat[j,min.idx[j]]
      }
      opt.idx<-which.min(optvec)
      beta.opt<-beta.est[,opt.idx]
      gamma.opt<-phi[,opt.idx]
      if(gamma.opt[1]<0)
      {
        gamma.opt<--gamma.opt
      }
      
      if(score.return)
      {
        re.opt<-list(gamma=gamma.opt,beta=beta.opt,PC.index=opt.idx,score=c(lambda[,opt.idx]),minmix=FALSE)
      }else
      {
        re.opt<-list(gamma=gamma.opt,beta=beta.opt,PC.index=opt.idx,minmix=FALSE)
      }
    }
    
    return(re.opt)
  }
}
