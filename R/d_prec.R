#'Compute Precision Analyses for Standardized Mean Differences
#'@param d Standardized means difference between groups
#'@param nlow starting sample size
#'@param nhigh ending sample size
#'@param by Incremental increase in sample (e.g. nlow = 10, nhigh = 24, by = 2, produces estimates of 10, 12, and 14)
#'@param propn1 Proportion in First Group
#'@param ci Type of Confidence Interval (e.g., .95)
#'@param tails number of tails for test (default is 2)
#'@examples
#'d_prec(d=.4,nlow=100, nhigh=2000, propn1=.5, ci=.95, by=100)
#'@return Precision Analyses for Standardized Mean Differences
#'@export
#'


d_prec<-function(d,nlow, nhigh, propn1= .5, ci=.95, tails=2, by=1)
{
    result <- data.frame(matrix(ncol = 6))
    colnames(result) <- c("n1", "n2","d","LL","UL","Precision")
    for(n in seq(nlow,nhigh, by)){
    n1<-n * propn1
    n2<-n * (1-propn1)
    a<-MBESS::ci.smd(smd=d, n.1=n1,n.2=n2, conf.level = .95)
    ll<-a[1]
    ul<-a[3]
    precision<-round((as.numeric(ul)-(as.numeric(ll))),4)
    ll<-round(as.numeric(ll),4)
    ul<-round(as.numeric(ul),4)
    result[n, 1]<-n1
    result[n, 2]<-n2
    result[n, 3]<-d
    result[n, 4]<-ll
    result[n, 5]<-ul
    result[n, 6]<-precision}
    output<-na.omit(result)
    rownames(output)<- c()
    output}


