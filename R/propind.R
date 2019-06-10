#'Compute power for Tests of Two Independent Proportions
#'Takes phi, degrees of freedom, and a range of sample sizes. Alpha is .05 by default, alternative values may be entered by user
#'@param p1 expected proportion Group 1
#'@param p2 expected proportion Group 2
#'@param nlow starting sample size
#'@param nhigh ending sample size
#'@param nratio ratio of sample size of first group to second (default is .5 for equally sized groups)
#'@param by Incremental increase in sample (e.g. nlow = 10, nhigh = 24, by = 2, produces estimates of 10, 12, and 14)
#'@param alpha Type I error (default is .05)
#'@param tails number of tails for test (default is 2)
#'@examples
#'propind(p1=.62, p2=.55,nlow=200,nhigh=2500, by=100,nratio=.2)
#'@return Power for Tests of Two Independent Proportions
#'@export
#'
#'

propind<-function(p1,p2,nlow, nhigh, nratio=0.5, alpha=.05, tails=2, by=1)
{
  result <- data.frame(matrix(ncol = 3))
  colnames(result) <- c("n1","n2","Power")
  {if(p1<0|p1>1.0|p2<0|p2>1.0){stop("Invalid proportions, must be between 0 and 1.0")
  }
    else
    p1a<-2*asin(p1^.5)
    p2a<-2*asin(p2^.5)
    h<- abs(p1a-p2a)
    for(n in seq(nlow,nhigh, by)){
      n1<-n*nratio
      n2<-n-n1
      nharm<-(2*n1*n2)/(n1+n2)
      zlambda<-h*((nharm/2)^.5)
      prob<-1-(alpha/tails)
      tabled<-abs(stats::qnorm(prob))
      zpower<-tabled-zlambda
      power<-round(1-stats::pnorm(zpower),4)
      result[n, 1]<-n1
      result[n, 2]<-n2
      result[n, 3]<-power}
      output<-na.omit(result)
      rownames(output)<- c()
      output

  }}

