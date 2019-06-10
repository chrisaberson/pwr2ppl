#'Compute power for a Paired t-test
#'Takes means, sd, and sample sizes. Alpha is .05 by default, alternative values may be entered by user.
#'correlation (r) defaults to .50.
#'@param m1 Mean for Pre Test
#'@param m2 Mean for Post Test
#'@param s Standard deviation
#'@param r Correlation pre-post measures (default is .50)
#'@param n Sample size
#'@param alpha Type I error (default is .05)
#'@examples
#'pairt(m1=25,m2=20, s = 5, n = 25, r = .5)
#'@return Power for the Paired t-test
#'@export

pairt<-function(m1=NULL,m2=NULL, s=NULL, n=NULL, r = NULL, alpha=.05)
  {
  cov<-s^2
  corr<-r*cov
  data <- MASS::mvrnorm(n, mu = c(m1,m2), Sigma = matrix(c(cov,corr,corr,cov), ncol = 2),
                 empirical = TRUE)
  data<-as.data.frame(data)
  t<-stats::t.test(data$V1,data$V2, paired=TRUE)
  lambda<-abs(t$statistic)^2
  minusalpha<-1-alpha
  Ft<-stats::qf(minusalpha, 1, n-1)
  Power<-round(1-stats::pf(Ft, 1,n-1,lambda),4)
  message("Power for n = ", n, " is ", Power)
  result <- data.frame(matrix(ncol = 2))
  colnames(result) <- c("n", "Power")
  result[, 1]<-n
  result[, 2]<-Power
  output<-na.omit(result)
  rownames(output)<- c()
  invisible(output)
  }

