winRatio<-function(survivalObject,treatmentVariable,treatmentCodes,data,
                   method=c('pssm','NP'),
                   plotTimeUnit=NULL,
                   secondSurvivalObject=NULL,
                   longitudinalVariable=NULL,
                   timeVar=NULL,
                   subjectId=NULL,
                   plotPoints=NULL,
                   nominalTimes=NULL,
                   pssmIntervals=3,
                   integrationIntervals=1,type='l',
                   mult=100,
                   findValue=function(times,values,c)
                     stats::approx(times,values,xout=c,method="constant",
                            rule=2,f=0)$y){
  survivalFormula=as.formula(paste(deparse(substitute(survivalObject)),'~1'))
  secondSurvivalFormula=as.formula(paste(deparse(substitute(secondSurvivalObject)),'~1'))
  if(as.character(secondSurvivalFormula)[2]=='NULL') secondSurvivalFormula=NULL
  survivalObjectC=as.character(substitute(survivalObject))
  time1=quantile(survival::survfit(Surv(data[,survivalObjectC[2]],
                              data[,survivalObjectC[3]]==0)~1),probs=.5)$quantile

  outpt=cbind(time=time1,cWinRatio(survivalObject,secondSurvivalObject,
                                   longitudinal_outcome=longitudinalVariable,
                                   treatmentVariable=treatmentVariable,
                                   codes=treatmentCodes,timeVar=timeVar,
                                   subjectId=subjectId,
                                   data,nominalTimes=nominalTimes,
                                   findValue=findValue))
  if(is.null(plotPoints)) return(outpt)
  rownames(outpt)[1]<-'Value'
  if(is.null(longitudinalVariable)){
    if(is.null(secondSurvivalFormula)) {print('Needs either secondSurvivalFormula or longitudinalFormula')
      return(outpt)}
  } else {
    method='longitudinal'
    if(is.null(timeVar)|is.null(nominalTimes)) {print('Longitudinal data needs timevar and nominalTimes')
      return(outpt)}
    if(!is.null(secondSurvivalFormula)) print('Longitudinal data secondSurvivalFormula ignored')
    }
  data0=subset(data,data[,treatmentVariable]==treatmentCodes[1])
  data1=subset(data,data[,treatmentVariable]==treatmentCodes[2])
  if(method=='pssm'){
    pssm0=pssm::pssm(secondSurvivalFormula,survivalFormula,data0,intervals=pssmIntervals)
    pssm1=pssm::pssm(secondSurvivalFormula,survivalFormula,data1,intervals=pssmIntervals)
    res1=winratio.pssm(pssm0,pssm1,plotPoints,integrationIntervals)
    }
  else if(method=='NP'){
    surv0=list()
    surv1=list()
    surv0$survP=eval(parse(text=as.character(secondSurvivalFormula)[2]),data0)
    surv0$survD=eval(parse(text=as.character(survivalFormula)[2]),data0)
    surv1$survP=eval(parse(text=as.character(secondSurvivalFormula)[2]),data1)
    surv1$survD=eval(parse(text=as.character(survivalFormula)[2]),data1)
    res1=winratio.NP(surv0,surv1,cvals=plotPoints,type=type)
  }
  else if(method=='longitudinal'){

    fixedLongitudinalFormula=as.formula(paste(longitudinalVariable,'~',timeVar))
    randomLongitudinalFormula=as.formula(paste('~',timeVar,'|',subjectId))
    mod0=createModel(data0,sF=survivalFormula,fLF=fixedLongitudinalFormula,
                     rLF=randomLongitudinalFormula,tV=timeVar)
    mod1=createModel(data1,sF=survivalFormula,fLF=fixedLongitudinalFormula,
                     rLF=randomLongitudinalFormula,tV=timeVar)
    res1=winratioSim(mod0,mod1,nominalTimes,cvals=plotPoints,mult)
  }
  outpt=rbind(outpt,plots=cbind(res1,WinRatioSE=NA))
  #class(outpt)<-c('winratio','data.frame')
  if(!is.null(plotTimeUnit)) {
    x=outpt
    ylimW=max(x$winRatio)
    xlim=max(x$time)
    ylimP=min(2*max(x[,2:5]))
    graphics::plot.default(x$time[-1],x$winRatio[-1],type='l',xlab=plotTimeUnit,
         ylab='win ratio',main=paste('Win Ratio by',plotTimeUnit),
         lty=1,ylim=c(0,1.5*ylimW))
    graphics::abline(a=x$winRatio[1],b=0,lty=2)
    graphics::legend(x=.6*xlim,y=.5*ylimW,lty=1:2,bty='n',pt.cex=0.5,
           legend=c('Win Ratio by Time','Observed Win Ratio'))
    plts=c("firstP10", "firstP01",  "secondP10",  "secondP01" )
    labels=c('Survival:Control>RX','Survival:Rx>Control','Secondary:Control>Rx',
             'Secondary:Rx>Control')
    graphics::plot.default(x$time[-1],x$firstP10[-1],type='l',
         xlab=plotTimeUnit,
         ylab='Probablity',
         lty=1,
         main=paste('Win Ratio by',plotTimeUnit),
         ylim=c(0,ylimP))

    for (i in 2:4){
      graphics::lines(x$time[-1],x[-1,plts[i]],type='l',lty=i)}
    graphics::legend(xlim/5,ylimP,bty='n',pt.cex=.3,lty=1:4,legend=labels)
    }
return(outpt)
}
