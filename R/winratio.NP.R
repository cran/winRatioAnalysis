winratio.NP<-function(surv0,surv1,
                   cvals,type=c('l')) {
#surv0 is list SurvP,SurvD of survival objects
#list of with survD, followed by SurvP
#survD is exact, SurvP may be interval

maxTime=max(surv0$survP[,1],surv0$survD[,1],surv1$survP[,1],surv1$survD[,1])
mp00=maxTime+100
mp99=maxTime+99

    fdist=function(surv){
    if(dim(surv)[1]>0){
    km=survival::survfit(surv~1)
    dist=function(t) approx(km$time,km$surv,t,method="constant",
                       rule=2,f=0,yleft=1)$y
    times=km$time
    drops=c(1,km$surv)[-(length(km$surv)+1)]-km$surv
    } else {
      dist=function(t) 1
      times=numeric(0)
      drops=numeric(0)
    }
    list(dist=dist,times=times,drops=drops)
    }
   createR=function(surv){
   ss=as.data.frame(cbind(surv$survP[,1:2],surv$survD))
     rect1=matrix(NA,dim(ss)[1],4)
     names(ss)<-c("pTime","pStatus","sTime","sStatus")
     rect1[,3]=ss$sTime
     rect1[,4]=ifelse(ss$sStatus==1,ss$sTime,mp00)
     if(attributes(surv$survP)$type=='right'){
     rect1[,1]=ifelse(ss$pStatus==1,ss$pTime,ifelse(ss$sStatus==1,mp99,ss$pTime))
     rect1[,2]=ifelse(ss$pStatus==1,ss$pTime,ifelse(ss$sStatus==1,mp99,mp00))
   }
     else  {
       names(ss)[1:2]<-c("pTime","pTime1")
       rect1[,1:2]=cbind(ss$pTime,ss$pTime1)
       }
     names(rect1)<-c("X1","X2","Y1","Y2")
     return(rect1)
     }
    firstprob=function(c,F1,F2){
    times=F1$times[F1$times<=c]
    drops=F1$drops[F1$times<=c]
    sum(F2$dist(times)*drops)
  }
  secondprob=function(c,g1,g2,F1,F2){
  if (type=='l') index=c(2,4) else index=c(1,3)
  G=Vectorize(function (t)
    sum(g2$p[(g2$rects[,index[1]]>t|g2$rects[,index[1]]==mp99)&(g2$rects[,index[2]]>c)],na.rm=TRUE)/F2$dist(c),'t')
  g=list()
  pred=g1$rects[,index[2]]>c
  g$times=g1$rects[pred,index[1]]
  g$drops=g1$p[pred]/F1$dist(c)
  pred=(g$times<=c)
  if(any(pred))
  return(F1$dist(c)*F2$dist(c)*sum(G(g$times[pred])*g$drops[pred],na.rm=TRUE))
  else return(0)
  }
  F1=fdist(surv1$survD)
  F0=fdist(surv0$survD)
  R0=createR(surv0)
  R1=createR(surv1)
  g0<-MLEcens::computeMLE(R=R0,c(1,1,1,1),max.inner=100,max.outer=1000,tol=.000001)
  g1=MLEcens::computeMLE(R=R1,c(1,1,1,1),max.inner=100,max.outer=1000,tol=.000001)
  outp1=sapply(1:length(cvals),function(i) {
  c=cvals[i]
   firstP10=firstprob(c,F1,F0)
   firstP01=firstprob(c,F0,F1)
   secondP01=secondprob(c,g0,g1,F0,F1)
   secondP10=secondprob(c,g1,g0,F1,F0)
   ans=c(firstP10,firstP01,secondP01,secondP10)
    probNoInfo=1-sum(ans)
   fractionFromDeath=(firstP10+firstP01)/(1-probNoInfo)
   winRatio=(firstP01+secondP01)/(firstP10+secondP10)
   outp=matrix(c(c,firstP10,firstP01,
                 secondP10,secondP01,winRatio),1,6)
   colnames(outp)=c('time','firstP10','firstP01','secondP10','secondP01',
                    'winRatio')
   return(outp)
} )
outp1=data.frame(t(outp1))
names(outp1)<-c('time','firstP10','firstP01','secondP10','secondP01','winRatio')
return(outp1)
}


