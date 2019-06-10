#'Compute power for a single sample proportion test
#'Takes phi, degrees of freedom, and a range of sample sizes. Alpha is .05 by default, alternative values may be entered by user
#'@param p1 expected proportion (a.k.a. alternative proportion)
#'@param p0 null proportion
#'@param nlow starting sample size
#'@param nhigh ending sample size
#'@param by Incremental increase in sample (e.g. nlow = 10, nhigh = 24, by = 2, produces estimates of 10, 12, and 14)
#'@param alpha Type I error (default is .05)
#'@param tails number of tails for test (default is 2)
#'@examples
#'prop1(p1=.60, p0=.42,nlow=20,nhigh=100, tails=1, by=10)
#'@return Power for Tests of Single Proportion
#'@export
#'
#'

prop1<-function(p1,p0,nlow, nhigh, alpha=.05, tails=2, by=1)
{
  result <- data.frame(matrix(ncol = 2))
  colnames(result) <- c("n","Power")
  {if(p1<0|p1>1.0|p0<0|p0>1.0){stop("Invalid proportions, must be between 0 and 1.0")
  }
    else
    p1a<-2*asin(p1^.5)
    p0a<-2*asin(p0^.5)
    h = abs(p1a-p0a)
    for(n in seq(nlow,nhigh, by)){
      zlambda<-h*(n^.5)
      prob<-1-(alpha/tails)
      tabled<-abs(stats::qnorm(prob))
      zpower<-tabled-zlambda
      power<-round(1-stats::pnorm(zpower),4)
      result[n, 1]<-n
      result[n, 2]<-power}
      output<-na.omit(result)
      rownames(output)<- c()
      output
  }}


