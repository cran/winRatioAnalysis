winratioSim<-function(mod0,mod1,nominalTimes,cvals,mult) {
if(is.null(cvals)) cvals=seq(0,max(nominalTimes),.5)+.1
censorTime=max(nominalTimes)+.15
d0t=data.frame()
for(i in 1:mult) {
  d0t=rbind(d0t,cbind(JM::simulate.jointModel(mod0,nsim=1,times=nominalTimes,
                                  censoring=censorTime)[[1]],study=i))}
d1t=data.frame()
for(i in 1:mult) {
  d1t=rbind(d1t,cbind(JM::simulate.jointModel(mod1,nsim=1,times=nominalTimes,
                                  censoring=censorTime)[[1]],study=i))}
timeVar=colnames(mod0$x$X)[2]
yvar=colnames(d0t)[1]
fp=function(c,x1,x2) kronecker(x1,x2,function(x,y) (y>x)&(x<c))
fc=function(x1,x2)   kronecker(x1,x2,function(x,y) x==1)

probs=function(c,i){
study=NULL
d0=d0t[d0t$study==i,]
d1=d1t[d1t$study==i,]
sd0=plyr::ddply(d0,.(study,id), function(x) utils::tail(x[x[,2]<c,],1))
sd1=plyr::ddply(d1,.(study,id), function(x) utils::tail(x[x[,2]<c,],1))
fp01=fp(c,sd0$Time,sd1$Time)&fc(sd0$Event,sd1$Event)
fp10=fp(c,sd1$Time,sd0$Time)&fc(sd1$Event,sd0$Event)
nfp=(!fp01)&(!fp10)
sp01=fp(Inf,sd0$y,sd1$y)&nfp
sp10=fp(Inf,sd1$y,sd0$y)&nfp
return(c(mean(fp01),mean(fp10),mean(sp01),mean(sp10)))
}

outp1=sapply(cvals,function(c1){
  if(c1 %in% nominalTimes) c=c1+.001
  else c=c1
  vs=matrix(0,mult,4)
   for(i in 1:mult) vs[i,]=probs(c,i)
   firstP10=mean(vs[,2])
   firstP01=mean(vs[,1])
   secondP01=mean(vs[,3])
   secondP10=mean(vs[,4])
   winRatioVector=(vs[,1]+vs[,3])/(vs[,2]+vs[,4])
   winRatio=mean(winRatioVector)
   outp=matrix(c(c,firstP10,firstP01,
                 secondP10,secondP01,winRatio),1,6)
   colnames(outp)=c('time','firstP10','firstP01','secondP10','secondP01','winRatio')
   return(outp)
} )
outp1=data.frame(t(outp1))
names(outp1)<-c('time','firstP10','firstP01','secondP10','secondP01','winRatio')
return(outp1)
}


