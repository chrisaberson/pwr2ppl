#'Compute Precision Analyses for Mean Differences
#'@param m1 Mean of first group
#'@param m2 Mean of second group
#'@param s1 Standard deviation of first group
#'@param s2 Standard deviation of second group
#'@param nlow starting sample size
#'@param nhigh ending sample size
#'@param by Incremental increase in sample (e.g. nlow = 10, nhigh = 24, by = 2, produces estimates of 10, 12, and 14)
#'@param propn1 Proportion in First Group
#'@param ci Type of Confidence Interval (e.g., .95)
#'@return Precision Analyses for Mean Differences
#'@export
#'

md_prec<-function(m1,m2,s1,s2,nlow, nhigh, propn1= .5, ci=.95, by=1)

  {
  for(n in seq(nlow,nhigh, by)){
    n1<-n * propn1
    n2<-n * (1-propn1)
    var1 <- s1*s1
    var2 <- s2*s2
    nxs1 <- (n1-1)*(var1)
    nxs2 <- (n2-1)*(var2)
    s2p<-(nxs1+nxs2)/(n1+n2-2)
    sp <- sqrt(s2p)
    d<-(m1-m2)/sp
    a<-MBESS::ci.smd(smd=d, n.1=n1,n.2=n2, conf.level = .95)
    ll<-a[1]
    ul<-a[3]
    ll<-as.numeric(ll)
    ul<-as.numeric(ul)
    ll_m<-ll*sp
    ul_m<-ul*sp
    ll_m<-round((ll_m),4)
    ul_m<-round((ul_m),4)
    print(paste("n1=",n1,",n2 =",n2,",d = ",d,",LL =  ",ll_m,",UL =  ",ul_m,",precision =",ul_m-ll_m ))}
}

