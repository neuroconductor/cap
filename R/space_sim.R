space_sim <-
function(M1,M2,thred=1e-4)
{
  M1.svd<-svd(M1)
  M2.svd<-svd(M2)
  
  k1<-length(which(abs(M1.svd$d)>thred))
  k2<-length(which(abs(M2.svd$d)>thred))
  
  M1.new<-matrix(matrix(M1.svd$u[,1:k1],ncol=k1)%*%diag(M1.svd$d[1:k1])%*%t(matrix(M1.svd$v[,1:k1],ncol=k1)),ncol=k1)
  M2.new<-matrix(matrix(M2.svd$u[,1:k1],ncol=k2)%*%diag(M2.svd$d[1:k2])%*%t(matrix(M2.svd$v[,1:k2],ncol=k2)),ncol=k2)
  
  S<-t(M1.new)%*%M2.new%*%t(M2.new)%*%M1.new
  
  re<-list(similarity=sum(diag(S))/min(k1,k2),k1=k1,k2=k2)
  
  return(re)
}
