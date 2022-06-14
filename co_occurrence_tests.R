require(CooccurrenceAffinity)

comb <- function(n, k){
  if (k > n) return(0)
  if (suppressWarnings(is.infinite(factorial(n))) || suppressWarnings(is.infinite(factorial(k)))){
    x <- gmp::factorialZ(n) %/% gmp::factorialZ(k) %/% gmp::factorialZ(n - k) 
    x <- as.numeric(x)
    if(is.infinite(x)) stop("n and/or k is too large for the internal comb function; there are likely too many species", "\n",
                            "in the dataset for calculations to be possible")
  } else{
    x <- factorial(n) / factorial (k) / factorial(n - k)
  }
  return(x)
}



cooc <- function(j, N, N1, N2){
  a <- comb(N1, j)
  b <- comb( (N - N1), (N2 - j))
  d <- comb(N, N2)
  return( (a * b) / d)
}



veech_p<-function(N,N1,N2,obs){
  p<-0
  for (i in obs:min(N1,N2)){
    pc<-cooc(i,N,N1,N2)
    p<-p+pc
  }
  return (list('p'=p))
}

load("Atmar_Patterson_Matrices.RData")
ap_data<-data.Atmar


res_veech<-c()
for (m in 1:length(ap_data)){
  mat<-ap_data[[m]]
  R<-nrow(mat)
  C<-ncol(mat)
  fill<-sum(mat)/(R*C)
  if ((min(R,C)>5)&&((R*C)<500)){
    site_n<-ncol(mat)
    amle_vals<-c()
    prop<-c()
    veech<-c()    
    for (i in 1:nrow(mat)){
      for (j in 1:nrow(mat)){
        if (i<j){
          A<-mat[i,]
          B<-mat[j,]
          C <- A + B
          N1<-sum(A)
          N2<-sum(B)
          X <- sum(C==2) ##total number of sites for both species
          if (N1*N2>0){# && max(N1,N2)/site_n<0.8 && min(N1,N2)/site_n>0.2){
            veech<-c(veech,veech_p(site_n,N1,N2,X)$p)
            prop<-rbind(prop,c(N1,N2,X,site_n))
            alpha_mle_list <- ML.Alpha(x=X, c(sum(A),sum(B),site_n),lev=0.9)
            if (alpha_mle_list != "Degenerate co-occurrence distribution!"){
              alpha_mle <- alpha_mle_list$est
              amle_vals<-c(amle_vals,alpha_mle)} else {amle_vals<-c(amle_vals,NaN)
              }
            }
          }
        }
      }
    res_veech<-rbind(res_veech,cbind(prop,amle_vals,veech))
    }
  print (m)
}


colnames(res_veech)<-c('N1','N2','X','site_n','alphaMLE','Veech_p')
res_veech<-res_veech[is.finite(rowSums(res_veech)),]

write.table(res_veech,'results_veech_comparison.csv',col.names=T,row.names=F,quote=F,sep=',')

res_5<-res_veech[res_veech[,5]>5,]

pdf('veech_vs_aMLE_revision_S1.pdf',width=10,height=5)
par(mfrow=c(1,2))
plot(res_5[,4],res_5[,5],pch=16,col=viridis(1,alpha=0.5)[1],cex=1.5,xlab='site number',ylab='Affinity',cex.lab=1.2,cex.axis=1.2,las=1,main='(a)',cex.main=1.5)
plot(res_5[,4],res_5[,6],pch=16,col=viridis(1,alpha=0.5)[1],cex=1.5,xlab='site number',ylab=expression("p"[V]),cex.lab=1.2,cex.axis=1.2,las=1,main='(b)',cex.main=1.5)
dev.off()




#######compartmented matrix for Fig 1 B and C
mat_file<-'compartmented.csv'
mat<-t(read.csv(mat_file,header=T,row.names=1))
site_n<-ncol(mat)
sp_n<-nrow(mat)

res<-c()
for (i in 1:sp_n){
    for (j in 1:sp_n){
      if (i<j){
        A<-mat[i,]
        B<-mat[j,]
        C <- A + B
        X <- sum(C==2) ##total number of sites for both species
        alpha_mle_list <- ML.Alpha(x=X, c(sum(A),sum(B),site_n),lev=0.9)
        if (alpha_mle_list != "Degenerate co-occurrence distribution!"){alpha_mle <- alpha_mle_list$est} else {alpha_mle<-NaN}
        jaccard0 = X / (X + sum(C==1))
        j_nm<-c()
        for (nm_rep in 1:1000){
          A_<-sample(A,site_n)
          C_ <- A_ + B
          X <- sum(C_==2) ##total number of sites for both species
          j_nm<-c(j_nm, X / (X + sum(C_==1)))
        }
        j_z<-(jaccard0-mean(j_nm))/sd(j_nm)
        res<-rbind(res,c(i,j,alpha_mle,j_z,sum(A),sum(B),sum(C>0),X))
        print(res[nrow(res),])
      }
    }
  }





colnames(res)<-c('sp1','sp2','alpha_mle','j_z','A','B','sp_n','X')

write.table(res,'results_compartmented.csv',col.names=T,row.names=F,quote=F,sep=',')



pdf('fig1_revision.pdf',width=15,height=5)
par(mfrow=c(1,3))
plot(res_veech[,6],res_veech[,5],ylab='Affinity',xlab=expression("p"[V]),cex.lab=1.2,cex.axis=1.2,las=1,main='(a)',cex.main=2,pch=16,ylim=c(-8,8))
plot(res[,4],res[,3],ylab='Affinity',xlab='SES Jaccard',cex.lab=1.2,cex.axis=1.2,las=1,main='(b)',cex.main=2,pch=16,ylim=c(-8,8))
plot(res[,7],res[,3],ylab='Affinity',xlab='number of species',cex.lab=1.2,cex.axis=1.2,las=1,main='(c)',cex.main=2,pch=16,ylim=c(-8,8))
abline(lm(res[,3]~res[,7]))

dev.off()








