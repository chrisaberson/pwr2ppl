#'Compute Precision Analyses for R-Squared
#'This approach simply loops a function from MBESS
#'@param R2 R-squared
#'@param pred Number of Predictors
#'@param nlow starting sample size
#'@param nhigh ending sample size
#'@param by Incremental increase in sample (e.g. nlow = 10, nhigh = 24, by = 2, produces estimates of 10, 12, and 14)
#'@param ci Type of Confidence Interval (e.g., .95)
#'@examples
#'R2_prec(R2=.467, nlow=24, nhigh=100, pred=3, by=4)
#'@importFrom MBESS ci.R2 ci.smd ci.cc
#'@return Precision Analyses for R-Squared
#'@export
#'

R2_prec<-function(R2,nlow, nhigh, pred, ci=.95, by=1)
  {
    result <- data.frame(matrix(ncol = 5))
    colnames(result) <- c("n","R Squared","LL","UL","Precision")
    for(n in seq(nlow,nhigh, by)){
    df1<-pred
    df2<-n-pred-1
    a<-MBESS::ci.R2(R2=R2, df.1=df1,df.2=df2, conf.level = .95, Random.Predictors = FALSE)
    ll<-a[1]
    ul<-a[3]
    precision<-round((as.numeric(ul)-(as.numeric(ll))),4)
    ll<-round(as.numeric(ll),4)
    ul<-round(as.numeric(ul),4)
    result[n, 1]<-n
    result[n, 2]<-R2
    result[n, 3]<-ll
    result[n, 4]<-ul
    result[n, 5]<-precision}
    output<-na.omit(result)
    rownames(output)<- c()
    output}
