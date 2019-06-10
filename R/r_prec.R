#'Compute Precision Analyses for Correlations
#'This approach simply loops a function from MBESS
#'@param r Correlation
#'@param nlow starting sample size
#'@param nhigh ending sample size
#'@param by Incremental increase in sample (e.g. nlow = 10, nhigh = 24, by = 2, produces estimates of 10, 12, and 14)
#'@param ci Type of Confidence Interval (e.g., .95)
#'@examples
#'r_prec(r=.3, nlow=80, nhigh=400, by=20, ci=.95)
#'@return Precision Analyses for Correlations
#'@export
#'

r_prec<-function(r,nlow, nhigh, ci=.95, by=1)
{
  result <- data.frame(matrix(ncol = 5))
  colnames(result) <- c("n","r","LL","UL","Precision")
  for(n in seq(nlow,nhigh, by)){
  a<-MBESS::ci.cc(r, n, ci)
  ll<-a[1]
  ul<-a[3]
  precision<-round((as.numeric(ul)-(as.numeric(ll))),4)
  ll<-round(as.numeric(ll),4)
  ul<-round(as.numeric(ul),4)
  result[n, 1]<-n
  result[n, 2]<-r
  result[n, 3]<-ll
  result[n, 4]<-ul
  result[n, 5]<-precision}
  output<-na.omit(result)
  rownames(output)<- c()
  output}

