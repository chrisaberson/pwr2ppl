#'Compute power for Chi Square Based on Effect Size
#'Takes phi, degrees of freedom, and a range of sample sizes. Alpha is .05 by default, alternative values may be entered by user
#'@param phi phi coefficient (effect size for 2x2)
#'@param df degrees of freedom
#'@param nlow starting sample size
#'@param nhigh ending sample size
#'@param by Incremental increase in sample (e.g. nlow = 10, nhigh = 24, by = 2, produces estimates of 10, 12, and 14)
#'@param alpha Type I error (default is .05)
#'@examples
#'ChiES(phi=.3,df=1,nlow=10,nhigh=200,by=10, alpha = .01)
#'@return Power for Chi Square Based on Effect Size
#'@export
#'
#'

ChiES<-function(phi, df, nlow, nhigh, by = 1, alpha=.05)
{
  result <- data.frame(matrix(ncol = 2))
  colnames(result) <- c( "n","Power")
  {if(phi<0|phi>1.0){stop("Invalid effect size, phi must be between 0 and 1.0")
  }
   else
      for(n in seq(nlow,nhigh, by)){
    lambda<-n*phi^2
    tabled<-stats::qchisq(1-alpha, df=df)
    power<-round(1-stats::pchisq(tabled, df=df, lambda),4)
    result[n, 1]<-n
    result[n, 2]<-power}
    output<-na.omit(result)
    rownames(output)<- c()
    output
  }}

