diag.level <-
function(Y,Phi)
{
  if(is.null(ncol(Phi))|ncol(Phi)==1)
  {
    stop("dimension of Phi is less than 2")
  }else
  {
    n<-length(Y)
    p<-ncol(Y[[1]])
    
    Tvec<-rep(NA,n)
    
    ps<-ncol(Phi)
    
    dl.sub<-matrix(NA,n,ps)
    colnames(dl.sub)<-paste0("Dim",1:ps)
    dl.sub[,1]<-1
    for(i in 1:n)
    {
      cov.tmp<-cov(Y[[i]])
      Tvec[i]<-nrow(Y[[i]])
      
      for(j in 2:ps)
      {
        phi.tmp<-Phi[,1:j]
        mat.tmp<-t(phi.tmp)%*%cov.tmp%*%phi.tmp
        dl.sub[i,j]<-det(diag(diag(mat.tmp)))/det(mat.tmp)
      }
    }
    
    pmean<-apply(dl.sub,2,function(y){return(prod(apply(cbind(y,Tvec),1,function(x){return(x[1]^(x[2]/sum(Tvec)))})))})
    
    re<-list(avg.level=pmean,sub.level=dl.sub)
    return(re)
  }
}
