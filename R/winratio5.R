winratio.pssm<-function(pssm0,pssm1,cvals,int){
  scaleUp=1/int
  top=ceiling(min(pssm0@intervals/pssm0@rescale,pssm1@intervals/pssm1@rescale))
  range=seq(0,top,1/scaleUp)
  fdistribution=function(pssm)
    data.frame(tdeath=range,F=pssm.survivalcurv(pssm,NULL,NULL)(range)[,3])
  #create fake data for using pssm
  m=length(range)
  xCmatrix=matrix(0,m*(m+1)/2,4)
  ts=0
  for (ideath in 1:length(range)){for (iprog in 1:ideath){
   tprog0=0
   tprog1=range[iprog]
   tdeath=range[ideath]
   ts=ts+1
   xCmatrix[ts,]=c(tprog0,tprog1,tdeath,0)}}
  colnames(xCmatrix)<-c("tprog0","tprog1","tdeath","cdeath")
  xCmatrix=data.frame(xCmatrix)
  xCmatrix=xCmatrix[order(xCmatrix[,3],xCmatrix[,2]),]
  F0=fdistribution(pssm0)
  F1=fdistribution(pssm1)
   gvalues=function(pssm,F){
    lfcn=llikef(cov1=NULL,cov2=NULL,accumulate=FALSE,dat=xCmatrix,rescale=pssm@rescale,m=pssm@intervals)
    pdat=cbind(xCmatrix,loglik=ifelse(xCmatrix$tprog1==0,0,exp(c(lfcn(pssm@estimates)))))
    mdat=merge(F,pdat)
    gdistribution=cbind(mdat,Gdistribution=1-mdat[,"loglik"]/mdat[,"F"])
  return(gdistribution)
    }
  g0=gvalues(pssm0,F0)
  g1=gvalues(pssm1,F1)
  #firstprob=function(f,F) cumsum(f*(F[-(m+1)]+F[-1])/2)
  firstprob=function(F1,F2,order) {
    #F1 is the one that the derivative is take of
    m=dim(F1)[1]
    FF=(F2[-m,2]+F2[-1,2])/2
    dF=(F1[-m,2]-F1[-1,2])
    vs=data.frame(cbind(F1[-1,1],cumsum(FF*dF)))
    names(vs)<-c('time',paste('firstP',order,sep=""))

    return(vs)
  }
  secondprob=function(G1,G2,order){
    names(G1)[c(2,7)]<-c("F1","G1")
    names(G2)[c(2,7)]<-c("F2","G2")
    C=merge(G1[,c("tdeath","tprog1","F1","G1")],
            G2[,c("tdeath","tprog1","F2","G2")],
            by=c("tdeath","tprog1"))
  vs=ddply(C,.(tdeath),.fun=
            function(C){
     m=dim(C)[1]
     Gdist=(C$G2[-1]+C$G2[-m])/2
     Gden=(C$G1[-m]-C$G1[-1])
     min(C$F1)*min(C$F2)*sum(Gdist*Gden)})
   names(vs)<-c("time",paste('secondP',order,sep=""))
   return(vs)
  }
   prob31=firstprob(F1,F0,'10')
   prob21=firstprob(F0,F1,'01')
   prob22=secondprob(g0,g1,'01')
   prob32=secondprob(g1,g0,'10')
   outp=Reduce(function(x,y) merge(x,y,all=TRUE),list(prob31,prob21,prob32,prob22))
   outp$winRatio=(outp$firstP01+outp$secondP01)/(outp$firstP10+outp$secondP10)
   vt=approx(outp$t,1:dim(outp)[1],cvals,method='constant',rule=2,f=1)$y
   return(outp[vt,])
}



