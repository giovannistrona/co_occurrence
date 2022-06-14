library(BiasedUrn)
###functions shoudl be obtained from the Science Advances paper SI, dowloaded from https://www.science.org/doi/10.1126/sciadv.abj9204
AlphInt <- dget("functions/AlphInt function.R")
AlphMLE <- dget("functions/AlphMLE.R")
NewAlph <- dget("functions/NewAlph.R")

###the following combination (20 sites, species a present in 7 sites, species b present in 10 sites, 4 sites shared)
###produces the error pasted below.
NewAlph(4, mA=7, mB=10, N=20)$AlphMLE

# Error in pFNCHypergeo(t, mA, N - mA, mB, exp(alp[i])) : 
#   Inconsistency. mean = 16, lower limit = 0, upper limit = 7


###the following code simulate random (possible) values for X, mA, mB and N, and stores
###combinations producing errors; sometimes, my R studio just crashes with no apparent reason
error_pars<-c()
sc<-1
while(T){
  N<-sample(10:500,1)
  X<-sample(1:N,1)
  if (X==N){mA<-X} else {mA<-sample(X:N,1)}
  if (mA==N){mB<-X} else {mB<-sample(X:(X+(N-mA)),1)}
  stop<-tryCatch(NewAlph(X, mA, mB, N)$AlphMLE,error=function(cond){return ('stop')})
  if (!is.null(stop) && stop=='stop'){
    error_pars<-rbind(error_pars,c(X, mA, mB, N))
    print (c(sc,100*nrow(error_pars)/sc))
    write.table(error_pars,'error_parameters.csv',col.names = F,row.names=F,quote=F,sep=',',append = T)}
  sc<-sc+1
}


