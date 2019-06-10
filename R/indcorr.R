#'Compute Power for Comparing Two Independent Correlations
#'Takes correlations and range of values
#'@param r1 Correlation for Group 1
#'@param r2 Correlation for Group 2
#'@param nlow Starting sample size
#'@param nhigh Ending sample size
#'@param by Incremental increase in sample size from low to high
#'@param propn1 Proportion of sample in first group (default is .50 for equally size groups)
#'@param tails one or two-tailed tests (default is 2)
#'@param alpha Type I error (default is .05)
#'@examples
#'indcorr(r1=.3,r2=.1,nlow=200,nhigh=800,by=50, tails=1)
#'@return Power for Comparing Two Independent Correlations
#'@export
#'
#'
indcorr<-function(r1,r2,nlow, nhigh, propn1= .5, alpha=.05, tails=2, by=1)
{
    result <- data.frame(matrix(ncol = 3))
    colnames(result) <- c("n1", "n2","Power")
    for(n in seq(nlow,nhigh, by)){
    n1<-n*propn1
    n2<-n*(1-propn1)
    zr1<-0.5*(log((1+r1)/(1-r1)))
    zr2<-0.5*(log((1+r2)/(1-r2)))
    zdiff<-abs(zr1-zr2)
    sdz<-((1/(n1-3))+(1/(n2-3)))^.5
    z<-zdiff/sdz
    alphatails<-alpha/tails
    tabled<-stats::qnorm(1-alphatails)
    zpower<-tabled-z
    Power<-round((1-stats::pnorm(zpower)),4)
    result[n, 1]<-n1
    result[n, 2]<-n2
    result[n, 3]<-Power}
    output<-na.omit(result)
    rownames(output)<- c()
    output
    }
