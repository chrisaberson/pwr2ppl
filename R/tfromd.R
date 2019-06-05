#'Compute power for a t test using d statistic
#'Takes d, sample size range, type of test, and tails.
#'@param d standardize mean difference (Cohen's d)
#'@param nlow Starting sample size
#'@param nhigh Ending sample size
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

tfromd<-function(d,nlow, nhigh, alpha=.05, test="I", tails=2, by=1)
{
  if (test=="I") {
    d<-abs(d)
    for(n in seq(nlow,nhigh, by)){
      ncalc<-n/2
      delta<-d*(ncalc^.5)
      lambda<-delta^2
      minusalpha<-1-alpha
      Ft<-stats::qf(minusalpha, 1, n-2)
      Power<-round(1-stats::pf(Ft, 1,n-2,lambda),4)
      print(paste("Power a per group n of (Independent)", n, "=", Power))}
  }
  else if (test=="P")
    for(n in seq(nlow,nhigh, by)){
      d<-abs(d)
      delta<-d*(n^.5)
      lambda<-delta^2
      minusalpha<-1-alpha
      Ft<-stats::qf(minusalpha, 1, n-1)
      Power<-round(1-stats::pf(Ft, 1,n-2,lambda),3)
      print(paste("Power for total n of (Paired)", n, "=", Power))}
  on.exit()}

