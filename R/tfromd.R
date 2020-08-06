#'Compute power for a t test using d statistic
#'Takes d, sample size range, type of test, and tails.
#'@param d standardize mean difference (Cohen's d)
#'@param nlow Starting total sample size
#'@param nhigh Ending total sample size
#'@param by Incremental increase in sample size from low to high
#'@param tails one or two-tailed tests (default is 2)
#'@param test "I" for independent, "P" for paired
#'@param alpha Type I error (default is .05)
#'@examples
#'tfromd(d=.2,nlow=10,nhigh=200,by=10, test="P")
#'tfromd(d=.2,nlow=10,nhigh=200,by=10, test="I")
#'@return Power for the t-test from d statistic
#'@export
#'
#'

tfromd<-function(d,nlow, nhigh, alpha=.05, test="I", tails=2, by=2)
{
  if (test=="I") {
    resultI <- data.frame(matrix(ncol = 2))
    colnames(resultI) <- c("n total","Power")
    d<-abs(d)
    for(n in seq(nlow,nhigh, by)){
      ncalc<-n
      delta<-d*((ncalc/2)^.5)
      lambda<-delta^2
      minusalpha<-1-alpha
      Ft<-stats::qf(minusalpha, 1, n-2)
      Power<-round(1-stats::pf(Ft, 1,n-2,lambda),4)
      resultI[n, 1]<-n*2
      resultI[n, 2]<-Power}
      outputI<-na.omit(resultI)
      rownames(outputI)<- c()
      outputI
      }

  else if (test=="P"){
    resultP <- data.frame(matrix(ncol = 2))
    colnames(resultP) <- c("total n","Power")
    for(n in seq(nlow,nhigh, by)){
      d<-abs(d)
      delta<-d*(n^.5)
      lambda<-delta^2
      minusalpha<-1-alpha
      Ft<-stats::qf(minusalpha, 1, n-1)
      Power<-round(1-stats::pf(Ft, 1,n-2,lambda),3)
      resultP[n, 1]<-n
      resultP[n, 2]<-Power}
      outputP<-na.omit(resultP)
      rownames(outputP)<- c()
      outputP}
      }
