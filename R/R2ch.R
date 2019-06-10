#'Compute power for R2 change in Multiple Regression (up to three predictors)
#'Requires correlations between all variables as sample size. Means, sds, and alpha are option. Also computes Power(All)
#'Example code below for three predictors. Expand as needed for four or five
#'@param ry1 Correlation between DV (y) and first predictor (1)
#'@param ry2 Correlation between DV (y) and second predictor (2)
#'@param ry3 Correlation between DV (y) and third predictor (3)
#'@param r12 Correlation between first (1) and second predictor (2)
#'@param r13 Correlation between first (1) and third predictor (3)
#'@param r23 Correlation between second (2) and third predictor (3)
#'@param n Sample size
#'@param alpha Type I error (default is .05)
#'@param my Mean of DV (default is 0)
#'@param m1 Mean of first predictor (default is 0)
#'@param m2 Mean of second predictor (default is 0)
#'@param m3 Mean of third predictor (default is 0)
#'@param sy Standard deviation of DV (default is 1)
#'@param s1 Standard deviation of first predictor (default is 1)
#'@param s2 Standard deviation of second predictor (default is 1)
#'@param s3 Standard deviation of third predictor (default is 1)
#'@examples
#'R2ch(ry1=.40,ry2=.40,ry3=-.40, r12=-.15, r13=-.60,r23=.25,n=24)
#'@return Power for R2 change in Multiple Regression (up to three predictors)
#'@export
#'
#'

R2ch<-function(ry1=NULL, ry2=NULL, ry3=NULL, r12=NULL, r13=NULL, r23=NULL,n=NULL, alpha=.05,
                   my=0, m1=0, m2=0,m3=0,s1=1,s2=1,s3=1,sy=1)
{
  pop <- MASS::mvrnorm(n, mu = c(my, m1, m2, m3), Sigma = matrix(c(sy, ry1, ry2, ry3,
                                                             ry1, s1, r12, r13,
                                                             ry2, r12,s2, r23,
                                                             ry3, r13, r23, s3),
                                                           ncol = 4), empirical = TRUE)


  pop2 = data.frame(pop)
  full<-summary(stats::lm(X1~X2+X3+X4, pop2))
  mod1<-summary(stats::lm(X1~X2, pop2))
  mod2<-summary(stats::lm(X1~X3, pop2))
  mod3<-summary(stats::lm(X1~X4, pop2))
  ch23<-round((full$r.squared-mod1$r.squared),4)
  ch13<-round((full$r.squared-mod2$r.squared),4)
  ch12<-round((full$r.squared-mod3$r.squared),4)
  fullR2<-round((full$r.squared),4)
  f2ch23<-ch23/(1-full$r.squared)
  f2ch13<-ch13/(1-full$r.squared)
  f2ch12<-ch12/(1-full$r.squared)
  df1<-full$fstatistic[2]-mod1$fstatistic[2]
  df2<-full$fstatistic[3]
  lambda1<-f2ch23*df2
  lambda2<-f2ch13*df2
  lambda3<-f2ch12*df2
  minusalpha<-1-alpha
  Ft<-stats::qf(minusalpha, df1, df2)
  powerch1<-round(1-stats::pf(Ft, df1,df2,lambda1),4)
  powerch2<-round(1-stats::pf(Ft, df1,df2,lambda2),4)
  powerch3<-round(1-stats::pf(Ft, df1,df2,lambda3),4)
  message("R2 Model = ", fullR2)
  message("R2 Change Vars2 and 3 over Var1 = ", ch23, ", Power = ", powerch1)
  message("R2 Change Vars1 and 3 over Var2 = ", ch13, ", Power = ", powerch2)
  message("R2 Change Vars1 and 2 over Var3 = ", ch12, ", Power = ", powerch3)
  result <- data.frame(matrix(ncol = 8))
  colnames(result) <- c("n", "full R2","R2 Ch V2&3","Power Ch V2&3","R2 Ch V1&3","Power Ch V1&3","R2 Ch V1&2","Power Ch V1&2")
  result[, 1]<-n
  result[, 2]<-fullR2
  result[, 3]<-ch23
  result[, 4]<-powerch1
  result[, 5]<-ch13
  result[, 6]<-powerch2
  result[, 7]<-ch12
  result[, 8]<-powerch3
  output<-na.omit(result)
  rownames(output)<- c()
  invisible(output)
  }
