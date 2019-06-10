#'Compute Power for Comparing Two Dependent Correlations, No Variables in Common
#'Takes correlations and range of values. First variable in each pair is termed predictor, second is DV
#'@param r12 Correlation between the predictor and DV (first set of measures)
#'@param rxy Correlation between the predictor and DV (second set of measures)
#'@param r1x Correlation between the predictor (first measure) and the predictor variable (first measure)
#'@param r2x Correlation between the DV (first measure) and the predictor variable (first measure)
#'@param r1y Correlation between the predictor (first measure) and the dependent variable (second measure)
#'@param r2y Correlation between the DV (first measure) and the dependent variable (second measure)
#'@param nlow Starting sample size
#'@param nhigh Ending sample size
#'@param by Incremental increase in sample size from low to high
#'@param tails one or two-tailed tests (default is 2)
#'@param alpha Type I error (default is .05)
#'@examples
#'depcorr0(r12=.4,rxy=.7,r1x=.3,r1y=.1,r2x=.45,r2y=.35,  nlow=20,nhigh=200,by=10, tails=2)
#'@return Power for Comparing Two Dependent Correlations, No Variables in Common
#'@export
#'
#'
depcorr0<-function(r12,rxy,r1x,r1y,r2x,r2y, nlow, nhigh, alpha=.05, tails=2, by=1)
{
  result <- data.frame(matrix(ncol = 2))
  colnames(result) <- c("n", "Power")
  for(n in seq(nlow,nhigh, by)){
    zr12<-0.5*(log((1+r12)/(1-r12)))
    zrxy<-0.5*(log((1+rxy)/(1-rxy)))
    zdiff<-abs(zr12-zrxy)
    rave<-(r12+rxy)/2
    denom <-(1-rave^2)^2
    numer1<-(r1x -(r12*r2x))*(r2y-(r2x*rxy))
    numer2<-(r1y -(r1x*rxy))*(r2x-(r12*r1x))
    numer3<-(r1x -(r1y*rxy))*(r2y-(r12*r1y))
    numer4<-(r1y -(r12*r2y))*(r2x-(r2y*rxy))
    numer<-(numer1 + numer2 +numer3+numer4)/2
    sd<-numer /denom
    z<-(zdiff*((n-3)^.5)) / ((2-(2*sd))^.5)
    alphatails<-alpha/tails
    tabled<-stats::qnorm(1-alphatails)
    zpower<-tabled-z
    Power<-round((1-stats::pnorm(zpower)),4)
    result[n, 1]<-n
    result[n, 2]<-Power}
    output<-na.omit(result)
    rownames(output)<- c()
    output
    }
