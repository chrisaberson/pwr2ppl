#'Compute Power for Regression Interaction (Correlation/Coefficient Approach)
#'@param Group 1 Estimates (r or b) for Group 1
#'@param Group 2 Estimates (r or b) for Group 2
#'@param Prop_n1 Proportion of Sample in First Group (defaults to equal sample sizes)
#'@param Estimates 1 for Correlations (default), 2 for coefficients
#'@param nlow starting sample size
#'@param nhigh ending sample size
#'@param by incrimental increase in sample (e.g. nlow = 10, nhigh = 24, by = 2, produces estimates of 10, 12, and 14)
#'@param alpha Type I error (default is .05)
#'@return Power for Regression Interaction (Correlation/Coefficient Approach)
#'@export
#'
#'


regint_pwr<-function(Group1,Group2, sx1=1, sx2=1, sy1=1, sy2=1, nlow, nhigh, alpha=.05, Prop_n1=.5, by=2, Estimates=1){
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
    Ft<-qf(minusalpha, df1, df2)
    Power<-round(1-pf(Ft, df1,df2,lambda),4)
    R2<-round((f2/(1+f2)),4)
    print(paste("Power with n1 = ", n1, "n2 = ", n2, "= ", Power))}
    print(paste("Effect size (R2 Change/Squared Semi Partial) = ", R2, f2, numer,denom))}

