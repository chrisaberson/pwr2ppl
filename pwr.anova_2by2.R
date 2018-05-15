#'Compute power for a Two by Two Between Subjects ANOVA with four levels.
#'Takes means, sds, and sample sizes for each group. Alpha is .05 by default, alterative values may be entered by user
#'@param m1.1 Cell mean for First level of Factor A, First level of Factor B
#'@param m1.2 Cell mean for First level of Factor A, Second level of Factor B
#'@param m2.1 Cell mean for Second level of Factor A, First level of Factor B
#'@param m2.2 Cell mean for Second level of Factor A, Second level of Factor B
#'@param s1.1 Cell standard deviation for First level of Factor A, First level of Factor B
#'@param s1.2 Cell standard deviation for First level of Factor A, Second level of Factor B
#'@param s2.1 Cell standard deviation for Second level of Factor A, First level of Factor B
#'@param s2.2 Cell standard deviation for Second level of Factor A, Second level of Factor B
#'@param n1.1 Cell sample size for First level of Factor A, First level of Factor B
#'@param n1.2 Cell sample size for First level of Factor A, Second level of Factor B
#'@param n2.1 Cell sample size for Second level of Factor A, First level of Factor B
#'@param n2.2 Cell sample size for Second level of Factor A, Second level of Factor B
#'@param alpha Type I error (default is .05)
#'@return Power for the One Factor ANOVA
#'@export

pwr.anova_2by2<-function(m1.1=NULL,m1.2=NULL,m2.1=NULL,m2.2=NULL, s1.1=NULL,s1.2=NULL,s2.1=NULL,s2.2=NULL,
                         n1.1=NULL,n1.2=NULL,n2.1=NULL,n2.2=NULL, alpha=.05, all="OFF"){
  x<-rnorm(n1.1,m1.1,s1.1)
  X<-x
  MEAN<-m1.1
  SD<-s1.1
  Z <- (((X - mean(X, na.rm = TRUE))/sd(X, na.rm = TRUE))) * SD
  y<-MEAN + Z
  A<-rep("A1",n1.1)
  B<-rep("B1",n1.1)
  l1.1<-data.frame(y, A, B)
  x<-rnorm(n1.2,m1.2,s1.2)
  X<-x
  MEAN<-m1.2
  SD<-s1.2
  Z <- (((X - mean(X, na.rm = TRUE))/sd(X, na.rm = TRUE))) * SD
  y<-MEAN + Z
  A<-rep("A1",n1.2)
  B<-rep("B2",n1.2)
  l1.2<-data.frame(y, A, B)
  x<-rnorm(n2.1,m2.1,s2.1)
  X<-x
  MEAN<-m2.1
  SD<-s2.1
  Z <- (((X - mean(X, na.rm = TRUE))/sd(X, na.rm = TRUE))) * SD
  y<-MEAN + Z
  A<-rep("A2",n2.1)
  B<-rep("B1",n2.1)
  l2.1<-data.frame(y, A, B)
  x<-rnorm(n2.2,m2.2,s2.2)
  X<-x
  MEAN<-m2.2
  SD<-s2.2
  Z <- (((X - mean(X, na.rm = TRUE))/sd(X, na.rm = TRUE))) * SD
  y<-MEAN + Z
  A<-rep("A2",n2.2)
  B<-rep("B2",n2.2)
  l2.2<-data.frame(y, A, B)
  simdat<-rbind(l1.1,l1.2,l2.1,l2.2)
  options(contrasts=c("contr.sum", "contr.poly"))
  anova<-aov(y~A*B, data=simdat)
  anova<-Anova(anova, type="III")
  SSA<-anova[2,1] #column, row
  SSB<-anova[3,1]
  SSAB<-anova[4,1]
  SSwin<-anova[5,1]
  dfA<-anova[2,2]
  dfB<-anova[3,2]
  dfAB<-anova[4,2]
  dfwin<-anova[5,2]
  MSwin<-SSwin/dfwin
  eta2A<-SSA/(SSA+SSwin)
  f2A<-eta2A/(1-eta2A)
  lambdaA<-f2A*dfwin
  minusalpha<-1-alpha
  FtA<-qf(minusalpha, dfA, dfwin)
  power.A<-round(1-pf(FtA, dfA,dfwin,lambdaA),3)
  eta2B<-SSB/(SSB+SSwin)
  f2B<-eta2B/(1-eta2B)
  lambdaB<-f2B*dfwin
  FtB<-qf(minusalpha, dfB, dfwin)
  power.B<-round(1-pf(FtB, dfB,dfwin,lambdaB),3)
  eta2AB<-SSAB/(SSAB+SSwin)
  f2AB<-eta2AB/(1-eta2AB)
  lambdaAB<-f2AB*dfwin
  FtAB<-qf(minusalpha,dfAB, dfwin)
  power.AB<-round(1-pf(FtAB,dfAB,dfwin,lambdaAB),3)
  power.All<-power.A*power.B*power.AB
  if (all=="OFF")
  {print(paste("Power for Main Effect Factor A=", power.A))
  print(paste("Power for Main Effect Factor B=", power.B))
  print(paste("Power for Interaction AB=", power.AB))}
  else
  {print(paste("Power for Main Effect Factor A=", power.A))
  print(paste("Power for Main Effect Factor B=", power.B))
  print(paste("Power for Interaction AB=", power.AB))
  print(paste("Power(All)=", power.All))}
#    list(Power.Factor.A = power.A, Power.Factor.B = power.B, Power.Factor.AB = power.AB, Power.All=power.All)
}
