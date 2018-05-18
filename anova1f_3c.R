#'Compute power for a One Factor ANOVA with three levels and contrasts.
#'Takes means, sds, and sample sizes for each group. Alpha is .05 by default, alterative values may be entered by user
#'@param m1 Mean of first group
#'@param m2 Mean of second group
#'@param m3 Mean of third group
#'@param s1 Standard deviation of first group
#'@param s2 Standard deviation of second group
#'@param s3 Standard deviation of third group
#'@param n1 Sample size for first group
#'@param n2 Sample size for second group
#'@param n3 Sample size for third group
#'@param alpha Type I error (default is .05)
#'@param c1 Weight for Contrast 1 (default is 0)
#'@param c2 Weight for Contrast 2 (default is 0)
#'@param c3 Weight for Contrast 3 (default is 0)
#'@return Power for the One Factor ANOVA
#'@export
#'
#'
anova1f_3c<-function(m1=NULL,m2=NULL,m3=NULL,s1=NULL,s2=NULL,s3=NULL,n1=NULL,n2=NULL,n3=NULL,alpha=.05, c1 =0, c2=0, c3=0)
  {
  x<-rnorm(n1,m1,s1)
  X<-x
  MEAN<-m1
  SD<-s1
  Z <- (((X - mean(X, na.rm = TRUE))/sd(X, na.rm = TRUE))) * SD
  y<-MEAN + Z
  group<-rep("A1",n1)
  l1<-data.frame(y, group)
  x<-rnorm(n2,m2,s2)
  X<-x
  MEAN<-m2
  SD<-s2
  Z <- (((X - mean(X, na.rm = TRUE))/sd(X, na.rm = TRUE))) * SD
  y<-MEAN + Z
  group<-rep("A2",n2)
  l2<-data.frame(y, group)
  x<-rnorm(n3,m3,s3)
  X<-x
  MEAN<-m3
  SD<-s3
  Z <- (((X - mean(X, na.rm = TRUE))/sd(X, na.rm = TRUE))) * SD
  y<-MEAN + Z
  group<-rep("A3",n3)
  l3<-data.frame(y, group)
  simdat<-rbind(l1,l2,l3)
  anova<-aov(y~group, data=simdat)
  anova<-Anova(anova, type="III")
  SSA<-anova[2,1] #column, row
  SSwin<-anova[3,1]
  dfwin<-anova[3,2]
  mswin<-SSwin/dfwin
  dfbg<-anova[2,2]
  eta2<-SSA/(SSA+SSwin)
  f2<-eta2/(1-eta2)
  lambda<-f2*dfwin
  minusalpha<-1-alpha
  Ft<-qf(minusalpha, dfbg, dfwin)
  power<-1-pf(Ft, dfbg,dfwin,lambda)
  delta=((c1*m1)+(c2*m2)+(c3*m3))/((mswin*((c1^2/n1)+(c2^2/n2)+(c3^2/n3))))^.5
  lambda.c=delta^2
  Ft.c<-qf(minusalpha, 1, dfwin)
  power.contrast<-1-pf(Ft.c, 1,dfwin,lambda.c)
  list(Power.for.Contrast = power.contrast)}
