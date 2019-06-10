#'Compute Power for Logistic Regression with Continuous Predictors
#'@param OR Odds Ratio for Predictor of Interest
#'@param r Correlation for Predictor of Interest
#'@param ER Event Ratio Probability of a Desirable Outcome Overall
#'@param R2 How Well Predictor of Interest is Explained by Other Predictors (default is 0)
#'@param alpha Type I error (default is .05)
#'@param power Desired Power
#'@examples
#'LRcont(OR = 4.05, ER = .463,power=.95)
#'@return Power for Logistic Regression with Continuous Predictors
#'@export
#'
#'

LRcont<-function(OR = NA, r = NA, ER=NULL, alpha=.05, power=NULL, R2=.00)
{
  est<-NA
  est[is.na(OR)]<-1 #r
  est[is.na(r)]<-2 #OR
  if(est=="1"){
  tod<-(2*r)/((1-r^2)^.5)
  OR<-exp(((tod)*pi)/(3^.5))}
  lgOR<-log(OR)
  zalpha<-stats::qnorm(1-alpha/2)
  zpower<-stats::qnorm(power)
  num<-zalpha+zpower
  den<-ER*(1-ER)*lgOR^2
  n<-(num^2/den)/(1-R2)
  nprint<-ceiling(n)
  OR<-round((OR),4)
  message("Sample Size = ", nprint, " for Odds Ratio = ", OR)
  result <- data.frame(matrix(ncol = 3))
  colnames(result) <- c("n", "OR","power")
  result[, 1]<-nprint
  result[, 2]<-OR
  result[, 3]<-power
  output<-na.omit(result)
  rownames(output)<- c()
  invisible(output)
}
