MatReg_QC_opt2 <-
function(Y,X,Phi0=NULL,method=c("CAP","CAP-C","CAP-C1"),CAP.OC=FALSE,max.itr=1000,tol=1e-4,trace=FALSE,score.return=TRUE,gamma0.mat=NULL,ninitial=NULL) 
{
  if(is.null(Phi0))
  {
    return(MatReg_QC_opt(Y,X,method=method,max.itr=max.itr,tol=tol,trace=trace,score.return=score.return,gamma0.mat=gamma0.mat,ninitial=ninitial))
  }else
  {
    n<-length(Y)
    q<-ncol(X)
    p<-ncol(Y[[1]])
    
    p0<-ncol(Phi0)
    # estimate beta0
    beta0<-rep(NA,p0)
    for(j in 1:p0)
    {
      beta0[j]<-MatReg_QC_beta(Y,X,gamma=Phi0[,j],max.itr=max.itr,tol=tol,trace=FALSE)$beta[1]
    }
    Ytmp<-vector("list",length=n)
    for(i in 1:n)
    {
      Y2tmp<-Y[[i]]-Y[[i]]%*%(Phi0%*%t(Phi0))
      
      Y2tmp.svd<-svd(Y2tmp)
      
      Ytmp[[i]]<-Y2tmp.svd$u%*%diag(c(Y2tmp.svd$d[1:(p-p0)],exp(beta0)))%*%t(Y2tmp.svd$v)
      
      # Y2tmp0<-Y2tmp.svd$u%*%diag(c(rep(0,p-p0),exp(beta0)))%*%t(Y2tmp.svd$v)
      # Ytmp[[i]]<-Y2tmp+Y2tmp0
    }
    
    if(method=="CAP")
    {
      if(CAP.OC==FALSE)
      {
        re.tmp<-MatReg_QC_opt(Ytmp,X,method=method,max.itr=max.itr,tol=tol,trace=FALSE,score.return=score.return,gamma0.mat=gamma0.mat,ninitial=ninitial)
      }else
      {
        re.tmp<-MatReg_QC_RE(Ytmp,X,Phi0=Phi0,max.itr=max.itr,tol=tol,trace=FALSE,score.return=score.return)
      }
    }else
    {
      re.tmp<-MatReg_QC_opt(Ytmp,X,method=method,max.itr=max.itr,tol=tol,trace=FALSE,score.return=score.return,gamma0.mat=gamma0.mat,ninitial=ninitial) 
    }
    
    re<-re.tmp
    re$orthogonal<-c(t(re.tmp$gamma)%*%Phi0)
    
    return(re)
  }
}
