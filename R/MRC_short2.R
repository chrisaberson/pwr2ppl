#'Compute Multiple Regression shortcuts with three predictors for Ind Coefficients
#'Requires correlations between all variables as sample size. Means and sds are option. Also computes Power(All)
#'@param ry1_1 Correlation between DV (y) and first predictor (1), first group
#'@param ry2_1 Correlation between DV (y) and second predictor (2), first group
#'@param ry3_1 Correlation between DV (y) and third predictor (3), first group
#'@param r12_1 Correlation between first (1) and second predictor (2), first group
#'@param r13_1 Correlation between first (1) and third predictor (3), first group
#'@param r23_1 Correlation between second (2) and third predictor (3), first group
#'@param n1 Sample size, first group
#'@param my_1 Mean of DV (default is 0), first group
#'@param m1_1 Mean of first predictor (default is 0), first group
#'@param m2_1 Mean of second predictor (default is 0), first group
#'@param m3_1 Mean of third predictor (default is 0), first group
#'@param sy_1 Standard deviation of DV (default is 1), first group
#'@param s1_1 Standard deviation of first predictor (default is 1), first group
#'@param s2_1 Standard deviation of second predictor (default is 1), first group
#'@param s3_1 Standard deviation of third predictor (default is 1), first group
#'@param ry1_2 Correlation between DV (y) and first predictor (1), second group
#'@param ry2_2 Correlation between DV (y) and second predictor (2), second group
#'@param ry3_2 Correlation between DV (y) and third predictor (3), second group
#'@param r12_2 Correlation between first (1) and second predictor (2), second group
#'@param r13_2 Correlation between first (1) and third predictor (3), second group
#'@param r23_2 Correlation between second (2) and third predictor (3), second group
#'@param n2 Sample size, second group
#'@param my_2 Mean of DV (default is 0), second group
#'@param m1_2 Mean of first predictor (default is 0), second group
#'@param m2_2 Mean of second predictor (default is 0), second group
#'@param m3_2 Mean of third predictor (default is 0), second group
#'@param sy_2 Standard deviation of DV (default is 1), second group
#'@param s1_2 Standard deviation of first predictor (default is 1), second group
#'@param s2_2 Standard deviation of second predictor (default is 1), second group
#'@param s3_2 Standard deviation of third predictor (default is 1), second group
#'@param alpha Type I error (default is .05)
#'@examples
#'MRC_short2(ry1_1=.40, ry2_1=.40, ry3_1 =-.40, r12_1=-.15,r13_1=-.60, r23_1=.25,
#'ry1_2=.40, ry2_2=.10, ry3_2 =-.40, r12_2=-.15,r13_2=-.60, r23_2=.25,
#'n1=50,n2=50,alpha=.05,my_1=1,m1_1=1,m2_1=1,m3_1=1,
#'sy_1=7,s1_1=1,s2_1=1,s3_1=2,
#'my_2=1,m1_2=1,m2_2=1,m3_2=1,sy_2=7,s1_2=1,s2_2=1,s3_2=2)
#'@return Multiple Regression shortcuts with three predictors for Ind Coefficients
#'@export
#'

MRC_short2<-function(ry1_1, ry2_1, ry3_1=NULL, r12_1, r13_1=NULL, r23_1=NULL,n1,
                     ry1_2, ry2_2, ry3_2=NULL, r12_2, r13_2=NULL, r23_2=NULL,n2, alpha=.05,
                     my_1=0, m1_1=0, m2_1=0,m3_1=0,s1_1=1,s2_1=1,s3_1=1,sy_1=1,
                     my_2=0, m1_2=0, m2_2=0,m3_2=0,s1_2=1,s2_2=1,s3_2=1,sy_2=1)
{
  pred<-NA
  pred[is.null(r23_1)]<-2
  pred[!is.null(r23_1)]<-3


  vary_1<-sy_1^2; var1_1<-s1_1^2;var2_1<-s2_1^2;var3_1<-s3_1^2
  covy1_1<-ry1_1*sy_1*s1_1
  covy2_1<-ry2_1*sy_1*s2_1
  covy3_1<-ry3_1*sy_1*s3_1
  cov12_1<-r12_1*s1_1*s2_1
  cov13_1<-r13_1*s1_1*s3_1
  cov23_1<-r23_1*s2_1*s3_1

  pop1 <- MASS::mvrnorm(n1, mu = c(my_1, m1_1, m2_1, m3_1), Sigma = matrix(c(vary_1, covy1_1, covy2_1, covy3_1,
                                                                       covy1_1, var1_1, r12_1, cov13_1,
                                                                       covy2_1, cov12_1,var2_1,cov23_1,
                                                                       covy3_1, cov13_1, cov23_1, var3_1),
                                                                     ncol = 4), empirical = TRUE)

  vary_2<-sy_2^2; var1_2<-s1_2^2;var2_2<-s2_2^2;var3_2<-s3_2^2
  covy1_2<-ry1_2*sy_2*s1_2
  covy2_2<-ry2_2*sy_2*s2_2
  covy3_2<-ry3_2*sy_2*s3_2
  cov12_2<-r12_2*s1_2*s2_2
  cov13_2<-r13_2*s1_2*s3_2
  cov23_2<-r23_2*s2_2*s3_2

  pop2 <- MASS::mvrnorm(n2, mu = c(my_2, m1_2, m2_2, m3_2), Sigma = matrix(c(vary_2, covy1_2, covy2_2, covy3_2,
                                                                       covy1_2, var1_2, r12_2, cov13_2,
                                                                       covy2_2, cov12_2,var2_2,cov23_2,
                                                                       covy3_2, cov13_2, cov23_2, var3_2),
                                                                     ncol = 4), empirical = TRUE)

  pop1<-data.frame(pop1)
  pop2<-data.frame(pop2)
  values1<-stats::lm(X1~X2+X3+X4, pop1)
  values1<-summary(values1)
  values1b<-stats::lm(X3~X2+X4, pop1)
  values1b<-summary(values1b)
  values2<-stats::lm(X1~X2+X3+X4, pop2)
  values2<-summary(values2)
  values2b<-stats::lm(X3~X2+X4,pop2)
  values2b<-summary(values2b)

output<-list(values1, values1b, values2, values2b)
names(output)<-c("Overall Analyses for R2 Full model and coefficients, First Group",
                 "Analyses for R2i (how well predictor is explained by other predictors, First Group",
                 "Overall Analyses for R2 Full model and coefficients, Second Group",
                 "Analyses for R2i (how well predictor is explained by other predictors, Second Group")
output
  }
