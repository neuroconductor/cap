cap_beta <-
function(Y,X,gamma=NULL,beta=NULL,method=c("asmp","LLR"),boot=FALSE,sims=1000,boot.ci.type=c("bca","perc"),conf.level=0.95,verbose=TRUE)
{
  n<-length(Y)
  p<-ncol(Y[[1]])
  Tvec<-rep(NA,n)
  
  q<-ncol(X)
  
  if(is.null(colnames(X)))
  {
    colnames(X)<-c("Intercept",paste0("X",1:(q-1)))
  }
  
  for(i in 1:n)
  {
    Tvec[i]<-nrow(Y[[i]])
  }
  
  if(boot)
  {
    if(is.null(gamma))
    {
      stop("Error! Need gamma value.")
    }else
    {
      beta.boot<-matrix(NA,q,sims)
      
      for(b in 1:sims)
      {
        idx.tmp<-sample(1:n,n,replace=TRUE)
        
        Ytmp<-Y[idx.tmp]
        Xtmp<-matrix(X[idx.tmp,],ncol=q)
        
        beta.boot[,b]<-MatReg_QC_beta(Ytmp,Xtmp,gamma=gamma)$beta
        
        if(verbose)
        {
          print(paste0("Bootstrap sample ",b))
        }
      }
      
      beta.est<-apply(beta.boot,1,mean,na.rm=TRUE)
      beta.se<-apply(beta.boot,1,sd,na.rm=TRUE)
      beta.stat<-beta.est/beta.se
      pv<-(1-pnorm(abs(beta.stat)))*2
      
      if(boot.ci.type[1]=="bca")
      {
        beta.ci<-t(apply(beta.boot,1,BC.CI,sims=sims,conf.level=conf.level))
      }
      if(boot.ci.type[1]=="perc")
      {
        beta.ci<-t(apply(beta.boot,1,quantile,probs=c((1-conf.level)/2,1-(1-conf.level)/2)))
      }
      
      re<-data.frame(Estiamte=beta.est,SE=beta.se,statistic=beta.stat,pvalue=pv,LB=beta.ci[,1],UB=beta.ci[,2])
      rownames(re)<-colnames(X)
      
      return(list(Inference=re,beta.boot=beta.boot))
    }
  }else
  {
    if(is.null(beta)&is.null(gamma)==FALSE)
    {
      beta<-MatReg_QC_beta(Y,X,gamma=gamma)$beta
    }else
      if(is.null(gamma))
      {
        stop("Error! Need gamma value.")
      }
    
    if(method[1]=="asmp")
    {
      beta.var<-2*solve(t(X)%*%X)/min(Tvec)
      
      beta.se<-sqrt(diag(beta.var))
      beta.stat<-beta/beta.se
      pv<-(1-pnorm(abs(beta.stat)))*2
      
      LB<-beta-beta.se*qnorm((1-conf.level)/2,lower.tail=FALSE)
      UB<-beta+beta.se*qnorm((1-conf.level)/2,lower.tail=FALSE)
      
      re<-data.frame(Estimate=beta,SE=beta.se,statistic=beta.stat,pvalue=pv,LB=LB,UB=UB)
      rownames(re)<-colnames(X)
    }else
      if(method[1]=="LLR")
      {
        stat=pv<-rep(NA,q)
        for(j in 1:q)
        {
          Xtmp<-matrix(X[,-j],nrow=n)
          beta0<-MatReg_QC_beta(Y,Xtmp,gamma=gamma)$beta
          stat[j]<-2*((-objfunc(Y,X,gamma,beta))-(-objfunc(Y,Xtmp,gamma,beta0)))
          pv[j]<-1-pchisq(stat[j],df=1)
        }
        re<-data.frame(Estimate=beta,statistic=stat,pvalue=pv)
        rownames(re)<-colnames(X)
      }
    
    return(re)
  }
}
