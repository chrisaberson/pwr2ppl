#'Compute power for Tests of Two Independent Proportions
#'Takes phi, degrees of freedom, and a range of sameple sizes. Alpha is .05 by default, alterative values may be entered by user
#'@param p1 expected proportion Group 1
#'@param p2 expected proportion Group 2
#'@param nlow starting sample size
#'@param nhigh ending sample size
#'@param by Incremental increase in sample (e.g. nlow = 10, nhigh = 24, by = 2, produces estimates of 10, 12, and 14)
#'@param alpha Type I error (default is .05)
#'@return Power for Tests of Two Independent Proportions
#'@export
#'
#'

propind<-function(p1,p2,nlow, nhigh, nratio=0.5, alpha=.05, tails=2, by=1)
{
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
      tabled<-abs(qnorm(prob))
      zpower<-tabled-zlambda
      power<-round(1-pnorm(zpower),4)
      #print(c(p1a,p2a))}
      print(paste("Power for sample sizes of ", n1, n2, "=", power))}

  }}

