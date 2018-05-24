#'Compute power for Chi Square Based on Effect Size
#'Takes phi, degrees of freedom, and a range of sameple sizes. Alpha is .05 by default, alterative values may be entered by user
#'@param phi phi coefficient (effect size for 2x2)
#'@param df degrees of freedom
#'@param nlow starting sample size
#'@param nhigh ending sample size
#'@param by Incremental increase in sample (e.g. nlow = 10, nhigh = 24, by = 2, produces estimates of 10, 12, and 14)
#'@param alpha Type I error (default is .05)
#'@return Power for Chi Square Based on Effect Size
#'@export
#'
#'

ChiES<-function(phi, df, nlow, nhigh, by = 1, alpha=.05)
{
  {if(phi<0|phi>1.0){stop("Invalid effect size, phi must be between 0 and 1.0")
  }
    else
      for(n in seq(nlow,nhigh, by)){
        lambda<-n*phi^2
        tabled<-qchisq(1-alpha, df=df)
        power<-round(1-pchisq(tabled, df=df, lambda),4)
        print(paste("Power for n of", n, "=", power))}
  }}
