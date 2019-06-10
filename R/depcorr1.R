#'Compute Power for Comparing Two Dependent Correlations, One Variable in Common
#'Takes correlations and range of values
#'@param r1y Correlation between the first predictor and the dependent variable
#'@param r2y Correlation between the second predictor and the dependent variable
#'@param r12 Correlation between the first predictor and the second predictor
#'@param nlow Starting sample size
#'@param nhigh Ending sample size
#'@param by Incremental increase in sample size from low to high
#'@param tails one or two-tailed tests (default is 2)
#'@param alpha Type I error (default is .05)
#'@examples
#'depcorr1(r1y=.3,r2y=.04,r12 = .2, nlow=100,nhigh=300,by=10, tails=2)
#'@return Power for Comparing Dependent Correlations, One Variable in Common
#'@export
#'
#'
depcorr1<-function(r1y,r2y,r12, nlow, nhigh, alpha=.05, tails=2, by=1)
{
  result <- data.frame(matrix(ncol = 2))
  colnames(result) <- c("n", "Power")
  for(n in seq(nlow,nhigh, by)){
    df<-n-3
    rdiff<-abs(r1y-r2y)
    rave<-(r1y+r2y)/2
    rdet<-1-(r1y**2)-(r2y**2)-(r12**2)+(2*r1y*r2y*r12)
    numer<-(n-1)*(1+r12)
    denom1<-((2*(n-1))/(n-3))*rdet
    denom2<-(rave**2)*((1-r12)**3)
    denom<-denom1+denom2
    delta<-rdiff*((numer/denom)^.5)
    alphatails<-alpha/tails
    tabled <- stats::qt(1-alphatails, df)
    Power<-round(1-stats::pt(tabled, df,delta),4)
    result[n, 1]<-n
    result[n, 2]<-Power}
    output<-na.omit(result)
    rownames(output)<- c()
    output
    }

