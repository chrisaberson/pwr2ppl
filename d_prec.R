#'Compute Precision Analyses for Standardized Mean Differences
#'@param d Standardized means difference between groups
#'@param nlow starting sample size
#'@param nhigh ending sample size
#'@param by Incremental increase in sample (e.g. nlow = 10, nhigh = 24, by = 2, produces estimates of 10, 12, and 14)
#'@param propn1 Proportion in First Group
#'@param ci Type of Confidence Interval (e.g., .95)
#'@return Precision Analyses for Standardized Mean Differences
#'@export
#'


d_prec<-function(d,nlow, nhigh, propn1= .5, ci=.95, tails=2, by=1)
{
  for(n in seq(nlow,nhigh, by)){
    n1<-n * propn1
    n2<-n * (1-propn1)
    a<-MBESS::ci.smd(smd=d, n.1=n1,n.2=n2, conf.level = .95)
    ll<-a[1]
    ul<-a[3]
    ll<-round(as.numeric(ll),4)
    ul<-round(as.numeric(ul),4)
    print(paste("n1=",n1,",n2 =",n2,"d = ",d,",LL =  ",ll,",UL =  ",ul,",precision =",ul-ll ))}
}

