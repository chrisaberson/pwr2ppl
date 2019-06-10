#'Compute Power for Regression Interaction (Correlation/Coefficient Approach)
#'@param Group1 Estimates (r or b) for Group 1
#'@param Group2 Estimates (r or b) for Group 2
#'@param Prop_n1 Proportion of Sample in First Group (defaults to equal sample sizes)
#'@param Estimates 1 for Correlations (default), 2 for coefficients
#'@param sx1 Standard deviation of predictor, group 1 (defaults to 1)
#'@param sx2 Standard deviation of predictor, group 2 (defaults to 1)
#'@param sy1 Standard deviation of outcome, group 1 (defaults to 1)
#'@param sy2 Standard deviation of outcome, group 2 (defaults to 1)
#'@param nlow starting sample size
#'@param nhigh ending sample size
#'@param by incremental increase in sample (e.g. nlow = 10, nhigh = 24, by = 2, produces estimates of 10, 12, and 14)
#'@param alpha Type I error (default is .05)
#'@examples
#'regint(Group1=-.26,Group2=.25, alpha=.05,Prop_n1=0.5,nlow=110, nhigh=140,by=2,Estimates=1)
#'@return Power for Regression Interaction (Correlation/Coefficient Approach)
#'@export
#'
#'


regint<-function(Group1,Group2, sx1=1, sx2=1, sy1=1, sy2=1, nlow, nhigh, alpha=.05, Prop_n1=.5, by=2, Estimates=1){
  result <- data.frame(matrix(ncol = 4))
  colnames(result) <- c("n1","n2","R2 Change)","Power")
  for(n in seq(nlow,nhigh, by)){
    n1 <- n * Prop_n1
    n2 <- n * (1-Prop_n1)
    if (Estimates=="1"){
      r1 <- Group1
      r2 <- Group2}
    if (Estimates=="2"){
      r1 <- Group1*(sx1/sy1)
      r2 <- Group2*(sx2/sy2)}
    sx1_sq <- sx1^2
    sx2_sq <- sx2^2
    sy1_sq <- sy1^2
    sy2_sq <- sy2^2
    r1_sq <- r1^2
    r2_sq <- r2^2
    numer1 <- ((n1-1)*r1_sq* sy1_sq) + ((n2-1)*r2_sq* sy2_sq)
    numer2 <- (((n1-1)*r1 * sx1 * sy1)+ ((n2-1)*r2 * sx2 * sy2))^2
    numer3 <- ((n1-1)* sx1_sq)+ ((n2-1)* sx2_sq)
    numer <- numer1 - (numer2 / numer3)
    denom <- ((n1-2)*(1-r1_sq)* sy1_sq) + ((n2-2)*(1-r2_sq)* sy2_sq)
    f2 <- numer/denom
    df1 <- 1
    df2 <- n-4
    lambda <- f2 * df2
    minusalpha<-1-alpha
    Ft<-stats::qf(minusalpha, df1, df2)
    Power<-round(1-stats::pf(Ft, df1,df2,lambda),4)
    R2<-round((f2/(1+f2)),4)
    result[n, 1]<-n1
    result[n, 2]<-n2
    result[n, 3]<-R2
    result[n, 4]<-Power}
    output<-na.omit(result)
    rownames(output)<- c()
    output}
