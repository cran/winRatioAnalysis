createModel=function(data1,sF,fLF,rLF,tV){
  datas=unique(data1[,all.vars(sF)])
  lmemod=nlme::lme(fixed=fLF,
             random=rLF,
             na.action=na.omit,data=data1,
             control=nlme::lmeControl(msMaxIter=1000,msVerbose=FALSE,tolerance=0.01,opt="optim"))
  survmod=survival::survreg(sF,data=datas,x=TRUE)
  lmemod$call$fixed=force(fLF)
  JM::jointModel(lmemod,survmod,timeVar=tV,
             method=c('weibull-PH-GH'),
             parameterization = 'value')
}

