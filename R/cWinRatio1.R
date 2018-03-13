cWinRatio=function(
  survivalObject,secondSurvivalObject=NULL,longitudinal_outcome=NULL,
  treatmentVariable,codes=c(0,1),timeVar=NULL,
  subjectId=NULL,
  data,
  nominalTimes=NULL,
  findValue=function(times,values,c){
  stats::approx(times,values,xout=c,method="constant",
         rule=2,f=0)$y}) {
  .=first=survtime.y=survtime.x=csurv.x=subject.x=csurv.y=subject.y=
                         progtime.y=progtime.x=cprog.x=cprog.y=NULL
  #kludge to avoid notes
  mm=dim(data)[1]
  rx=data[,treatmentVariable]
  if(is.null(subjectId)) subject=1:mm
      else subject=data[,subjectId]
  secondSurvivalObjectC=as.character(substitute(secondSurvivalObject,env=parent.frame()))
  survivalObjectC=as.character(substitute(survivalObject,env=parent.frame()))
  notThere=rep(NA,mm)
  if(length(secondSurvivalObjectC)==0){cprog=notThere;progtime=notThere}
  else if (length(secondSurvivalObjectC)==4){
    cprog=ifelse( is.na(data[,secondSurvivalObjectC[3]]),0,1)
    progtime=ifelse( cprog==0 ,data[,survivalObjectC[2]], data[,secondSurvivalObjectC[3]])
  }
    else{
    cprog=data[,secondSurvivalObjectC[3]]==1
    progtime=data[,secondSurvivalObjectC[2]]
  }
  if(is.null(longitudinal_outcome)) {time=notThere;y=notThere} else
  {time=data[,timeVar];y=data[,longitudinal_outcome]}

  dat=data.table(
    survtime=data[,survivalObjectC[2]],
    csurv=data[,survivalObjectC[3]]==1,
    progtime=progtime,
    cprog= cprog,
    rx=rx,
    subject=subject,
    time=time,
    y=y,
    keyby=subject)

  dat$first=ifelse(is.null(timeVar),TRUE,dat[,time==min(time),by=.(subject)][,2])
  dat0=data.frame(dat[(dat$rx==codes[1])&(first)])
  dat1=data.frame(dat[(dat$rx==codes[2])&(first)])
  pairs=data.table(merge(dat0,dat1,by=NULL,allow.cartesian=TRUE))
  m=dim(pairs)[1]

  firstP01=dim(pairs[(survtime.y>survtime.x)&(csurv.x),list(subject.x)])[1]/m
  firstP10=dim(pairs[(survtime.x>survtime.y)&(csurv.y),list(subject.y)])[1]/m
  notfirst=pairs[!((survtime.y>survtime.x)&(csurv.x)|
                 (survtime.x>survtime.y)&(csurv.y)),]

  if(is.null(longitudinal_outcome)&!is.null(secondSurvivalObjectC[1])){
    secondP01=dim(notfirst[(progtime.y>progtime.x)&(cprog.x),list(subject.x)])[1]/m
    secondP10=dim(notfirst[(progtime.x>progtime.y)&(cprog.y),list(subject.y)])[1]/m
  }
  else if(!is.null(longitudinal_outcome)) {
    notfirst$followup=pmin(notfirst$survtime.x,notfirst$survtime.y)
    mm=dim(notfirst)[1]

    ts=rowSums(sapply(1:mm,function(i) {
      dat1=dat[subject==notfirst[i,subject.x],list(time,y)]
      x=findValue(dat1$time,dat1$y,notfirst[i,'followup',with=FALSE])
      dat1=dat[subject==notfirst[i,subject.y],list(time,y)]
      y=findValue(dat1$time,dat1$y,notfirst[i,'followup',with=FALSE])
      return(c(secondP01=y>x,secondP10=x>y))
    }))/m
    secondP01=ts[1]
    secondP10=ts[2]
  } else {
   secondP01=0
   secondP10=0
  }

  winRatio=(firstP01+secondP01)/(firstP10+secondP10)
if(is.null(longitudinal_outcome)) se=  WR(data.frame(dat),'rx',codes[2:1],
                                         'survtime','csurv',
                                         'progtime','cprog')
  else
    se=WR(data.frame(dat),'rx',codes[2:1],
     'survtime','csurv',
     'progtime','cprog',longitudinal_outcome="y",
     timeVar='time',subject="subject",findValue=findValue)
if(abs(se$wr-winRatio)>.01) print(paste("Possible Error, computations of WR aren't equal",se$wr,winRatio))

outpt=data.frame(firstP10,
             firstP01,
             secondP10,
             secondP01,
             winRatio=winRatio,WinRatioSE=se$se)

  return(outpt)
}




