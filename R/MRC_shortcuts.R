#'Compute Multiple Regression shortcuts with three predictors (will expand to handle two to five)
#'Requires correlations between all variables as sample size. Means and sds are option. Also computes Power(All)
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
#'MRC_shortcuts(ry1=.40,ry2=.40,ry3=-.40, r12=-.15, r13=-.60,r23=.25,
#'n=110, my=1,m1=1,m2=1,m3=1,sy=7,s1=1,s2=1,s3=2)
#'@return Multiple Regression shortcuts with three predictors
#'@export
#'
MRC_shortcuts<-function(ry1=NULL, ry2=NULL, ry3=NULL, r12=NULL, r13=NULL, r23=NULL,n=100, alpha=.05,
                   my=0, m1=0, m2=0,m3=0,s1=1,s2=1,s3=1,sy=1){
vary<-sy^2; var1<-s1^2;var2<-s2^2;var3<-s3^2
covy1<-ry1*sy*s1
covy2<-ry2*sy*s2
covy3<-ry3*sy*s3
cov12<-r12*s1*s2
cov13<-r13*s1*s3
cov23<-r23*s2*s3
pop <- MASS::mvrnorm(n, mu = c(my, m1, m2, m3), Sigma = matrix(c(vary, covy1, covy2, covy3,
                                                           covy1, var1, cov12, cov13,
                                                           covy2, cov12,var2,cov23,
                                                           covy3, cov13, cov23, var3),
                                                           ncol = 4), empirical = TRUE)


pop2 = data.frame(pop)
pred<-NA
pred[is.null(r23)]<-2
pred[!is.null(r23)]<-3

full<-summary(stats::lm(X1~X2+X3+X4, pop2))
full
reduced<-summary(stats::lm(X1~X2+X3, pop2))
reduced
}
