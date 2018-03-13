FS_step = function(d)
{
  names(d)<-c( 't1','c1','t2','c2')
  u=ifelse(d$c1 == 1 & d$c2 == 1,sign(d$t1 - d$t2),ifelse(d$c1 == 1 & d$c2 == 0 & d$t1 < d$t2,-1,
                                                  ifelse(d$c1 == 0 & d$c2 == 1 & d$t1 > d$t2, 1,0)))
  return(u)
}

pairwise_comp=function(data1,data2,CVDEATHT,CVDEATH,MIDT,MI,
                       longitudinal_outcome,timeVar,subject,findValue){
 if(missing(longitudinal_outcome))longitudinal_outcome=NULL
  if(is.null(longitudinal_outcome)){
    data1$subject=1:dim(data1)[1]
    data2$subject=1:dim(data2)[1]
    data1=data.table(data1[,c('subject',CVDEATHT,CVDEATH,MIDT,MI)])
    data2=data.table(data2[,c('subject',CVDEATHT,CVDEATH,MIDT,MI)])
    pairs=merge.data.frame(data1,data2,by=NULL,allow.cartesian=TRUE)
    firstc=FS_step(pairs[,c(2:3,7:8)])
    secondc=FS_step(pairs[,c(4:5,9:10)])
  comb=ifelse(firstc!=0,firstc,secondc)
  return(sparseMatrix(i=pairs$subject.x,j=pairs$subject.y,x=comb))
  }  else {
  #Kludge to avoid a name conflict between the subject code from the database and the subject number 1,2,3...
  if(subject=='subject'){
    names(data1)[which(names(data1)==subject)]<-'subjectid'
    names(data2)[which(names(data2)==subject)]<-'subjectid'
    subject='subjectid'
  }

  data1=data1[,c(subject,CVDEATHT,CVDEATH,longitudinal_outcome,timeVar)]
  data2=data2[,c(subject,CVDEATHT,CVDEATH,longitudinal_outcome,timeVar)]
  data1u=unique(data1[,c(subject,CVDEATHT,CVDEATH)])
  data2u=unique(data2[,c(subject,CVDEATHT,CVDEATH)])

  data1u$subject=1:dim(data1u)[1]
  data2u$subject=1:dim(data2u)[1]
  pairs=merge.data.frame(data1u,data2u,by=NULL,allow.cartesian=TRUE)
  vars=c(paste(CVDEATHT,'.x',sep=""),paste(CVDEATH,'.x',sep=""),
         paste(CVDEATHT,'.y',sep=""),paste(CVDEATH,'.y',sep=""))
  pairs$firstc=FS_step(pairs[,vars])
  pairs$notfirst=(pairs$firstc==0)
  pairs$followup=pmin(pairs[,paste(CVDEATHT,'.x',sep='')],
                      pairs[,paste(CVDEATHT,'.y',sep='')])
  mm=dim(pairs)[1]
  pairs$second=sapply(1:mm,function(i) {
    if(pairs$notfirst[i]){
      dat1=data1[data1[,subject]==pairs[i,paste(subject,'.x',sep='')],c(timeVar,longitudinal_outcome)]
      dat2=data2[data2[,subject]==pairs[i,paste(subject,'.y',sep='')],c(timeVar,longitudinal_outcome)]
      dat1=dat1[order(dat1[,timeVar]),]
      dat2=dat2[order(dat2[,timeVar]),]

      x=findValue(unlist(dat1[,timeVar]),unlist(dat1[,longitudinal_outcome]),pairs$followup[i])
      y=findValue(unlist(dat2[,timeVar]),unlist(dat2[,longitudinal_outcome]),pairs$followup[i])
      return(return(sign(x-y)))
    } else return(0)
  })

  return(sparseMatrix(i=pairs$subject.x,j=pairs$subject.y,x=pairs$firstc+pairs$second))
}
}
get_mu_var_twogroups = function(dataa,TX,TXcodes,CVDEATHT,CVDEATH,MIDT,MI,...)
  {
    #note 1 0
    qq = pairwise_comp(dataa[dataa[,TX]==TXcodes[1],],
                                   dataa[dataa[,TX]==TXcodes[2],],
                       CVDEATHT,CVDEATH,MIDT,MI,...)
    qq=as.matrix(qq)
    qq1 <- qq2 <- qq
    qq1[qq1 < 0] = 0
    qq2[qq2 > 0] = 0
    qq2 = - qq2

    n1 = sum(dataa[,TX]==TXcodes[1])
    n2 = sum(dataa[,TX]==TXcodes[2])

    qq1_r = rowSums(qq1)
    qq1_c = colSums(qq1)
    qq2_r = rowSums(qq2)
    qq2_c = colSums(qq2)
    p1 = 1/(n1*n2)*sum(qq1)
    p2 = 1/(n1*n2)*sum(qq2)
    p3 = 1/(n1*n2*(n2-1))*sum(qq1_c*(qq1_c-1))
    p3 = 1/(n1*n2*(n2-1))*sum(qq1_c*(qq1_c-1))
    p4 = 1/(n1*n2*(n1-1))*sum(qq1_r*(qq1_r-1))
    p5 = 1/(n1*n2*(n2-1))*sum(qq2_c*(qq2_c-1))
    p6 = 1/(n1*n2*(n1-1))*sum(qq2_r*(qq2_r-1))
    p7 = 1/(n1*n2*(n2-1))*sum(qq1_c*qq2_c)
    p8 = 1/(n1*n2*(n1-1))*sum(qq1_r*qq2_r)

    xi_10_11 = p3 - p1^2
    xi_01_11 = p4 - p1^2
    xi_10_22 = p5 - p2^2
    xi_01_22 = p6 - p2^2
    xi_10_12 = p7 - p1*p2
    xi_01_12 = p8 - p1*p2

    s11 = (n1+n2)/n1*xi_10_11 + (n1+n2)/n2*xi_01_11
    s22 = (n1+n2)/n1*xi_10_22 + (n1+n2)/n2*xi_01_22
    s12 = (n1+n2)/n1*xi_10_12 + (n1+n2)/n2*xi_01_12
list(mu = c(sum(qq1),sum(qq2))/(n1*n2),Sigma = matrix(c(s11,s12,s12,s22),ncol=2)/(n1+n2))
}
WR_test = function(mu,Sigma)
    {
      out = c(mu[1]/mu[2],
              mu[1]/mu[2] - 1.96*mu[1]/mu[2]*sqrt(Sigma[1,1]/mu[1]^2 + Sigma[2,2]/mu[2]^2-2*Sigma[1,2]/(mu[1]*mu[2])),
              mu[1]/mu[2] + 1.96*mu[1]/mu[2]*sqrt(Sigma[1,1]/mu[1]^2 + Sigma[2,2]/mu[2]^2-2*Sigma[1,2]/(mu[1]*mu[2])),
              (mu[1]/mu[2]-1)/(mu[1]/mu[2]*sqrt(Sigma[1,1]/mu[1]^2 + Sigma[2,2]/mu[2]^2-2*Sigma[1,2]/(mu[1]*mu[2]))))
      names(out) = c("WR","LL","UL","z")
      out
    }
WR=function(dataa,TX,TXcodes,CVDEATHT,CVDEATH,MIDT,MI,...){

  parms=get_mu_var_twogroups(dataa,TX,TXcodes,CVDEATHT,CVDEATH,MIDT,MI,...)
              #longitudinal_outcome=NULL,time=NULL,subject=NULL)
  #findValue=function(times,values,c){
   # stats::approx(times,values,xout=c,method="constant",
  #rule=2,f=0)$y}

wr=parms$mu[1]/parms$mu[2]
se=parms$mu[1]/parms$mu[2]*sqrt(parms$Sigma[1,1]/parms$mu[1]^2 +
                                  parms$Sigma[2,2]/parms$mu[2]^2-2*parms$Sigma[1,2]/(parms$mu[1]*parms$mu[2]))

  list(wr=wr,se=se)}

