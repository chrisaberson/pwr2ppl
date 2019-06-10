#'Compute power for Simple Effects in a Two by Two Between Subjects ANOVA with two levels for each factor.
#'Takes means, sds, and sample sizes for each group. Alpha is .05 by default, alternative values may be entered by user
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
#'examples
#'anova2x2_se(m1.1=0.85, m1.2=0.85, m2.1=0.00, m2.2=0.60,
#'s1.1=1.7, s1.2=1.7, s2.1=1.7, s2.2=1.7,
#'n1.1=250, n1.2=250, n2.1=250, n2.2=250, alpha=.05)
#'@return Power for Simple Effects Tests in a Two By Two ANOVA
#'@export
#'
#'

anova2x2_se<-function(m1.1=NULL,m1.2=NULL,m2.1=NULL,m2.2=NULL, s1.1=NULL,s1.2=NULL,s2.1=NULL,s2.2=NULL,
                          n1.1=NULL,n1.2=NULL,n2.1=NULL,n2.2=NULL, alpha=.05){
  oldoption<-options(contrasts=c("contr.helmert", "contr.poly"))
  oldoption
  on.exit(options(oldoption))
  x<-stats::rnorm(n1.1,m1.1,s1.1)
  X<-x
  MEAN<-m1.1
  SD<-s1.1
  Z <- (((X - mean(X, na.rm = TRUE))/stats::sd(X, na.rm = TRUE))) * SD
  y<-MEAN + Z
  A<-rep("A1",n1.1)
  B<-rep("B1",n1.1)
  l1.1<-data.frame(y, A, B)
  x<-stats::rnorm(n1.2,m1.2,s1.2)
  X<-x
  MEAN<-m1.2
  SD<-s1.2
  Z <- (((X - mean(X, na.rm = TRUE))/stats::sd(X, na.rm = TRUE))) * SD
  y<-MEAN + Z
  A<-rep("A1",n1.2)
  B<-rep("B2",n1.2)
  l1.2<-data.frame(y, A, B)
  x<-stats::rnorm(n2.1,m2.1,s2.1)
  X<-x
  MEAN<-m2.1
  SD<-s2.1
  Z <- (((X - mean(X, na.rm = TRUE))/stats::sd(X, na.rm = TRUE))) * SD
  y<-MEAN + Z
  A<-rep("A2",n2.1)
  B<-rep("B1",n2.1)
  l2.1<-data.frame(y, A, B)
  x<-stats::rnorm(n2.2,m2.2,s2.2)
  X<-x
  MEAN<-m2.2
  SD<-s2.2
  Z <- (((X - mean(X, na.rm = TRUE))/stats::sd(X, na.rm = TRUE))) * SD
  y<-MEAN + Z
  A<-rep("A2",n2.2)
  B<-rep("B2",n2.2)
  l2.2<-data.frame(y, A, B)
  simdat<-rbind(l1.1,l1.2,l2.1,l2.2)
  dataA1<-subset(simdat, A=="A1")
  dataA2<-subset(simdat, A=="A2")
  dataB1<-subset(simdat, B=="B1")
  dataB2<-subset(simdat, B=="B2")
  options(contrasts=c("contr.sum", "contr.poly"))
  anova<-stats::aov(y~A*B, data=simdat)
  anova<-car::Anova(anova, type="III")
  SSwin<-anova[5,1] #row column
  dfwin<-anova[5,2]
  SSA<-anova[2,1] #column, row
  SSB<-anova[3,1]
  SSAB<-anova[4,1]
  SST<-SSA+SSB+SSAB+SSwin
  MSwin<-SSwin/dfwin

  options(contrasts=c("contr.sum", "contr.poly"))
  anoAatB1<-stats::aov(y~A, data=dataB1)
  anoAatB1<-car::Anova(anoAatB1, type="III")
  options(contrasts=c("contr.sum", "contr.poly"))
  anoAatB2<-stats::aov(y~A, data=dataB2)
  anoAatB2<-car::Anova(anoAatB2, type="III")

  options(contrasts=c("contr.sum", "contr.poly"))
  anoBatA1<-stats::aov(y~B, data=dataA1)
  anoBatA1<-car::Anova(anoBatA1,type="III")
  options(contrasts=c("contr.sum", "contr.poly"))
  anoBatA2<-stats::aov(y~B, data=dataA2)
  anoBatA2<-car::Anova(anoBatA2, type="III")


  dfwinSE<-dfwin+2
  SSBatA1<-anoBatA1[2,1]
  dfBatA1<-anoBatA1[2,2]
  eta2BatA1<-SSBatA1/SST
  f2BatA1<-eta2BatA1/(1-eta2BatA1)
  lambdaBatA1<-f2BatA1*dfwinSE
  minusalpha<-1-alpha
  FtBatA1<-stats::qf(minusalpha, dfBatA1, dfwinSE)
  power.BatA1<-round(1-stats::pf(FtBatA1, dfBatA1,dfwinSE,lambdaBatA1),3)

  SSBatA2<-anoBatA2[2,1]
  dfBatA2<-anoBatA2[2,2]
  eta2BatA2<-SSBatA2/SST
  f2BatA2<-eta2BatA2/(1-eta2BatA2)
  lambdaBatA2<-f2BatA2*dfwinSE
  FtBatA2<-stats::qf(minusalpha, dfBatA2, dfwinSE)
  power.BatA2<-round(1-stats::pf(FtBatA2, dfBatA2,dfwinSE,lambdaBatA2),3)

  SSAatB1<-anoAatB1[2,1]
  dfAatB1<-anoAatB1[2,2]
  dfwinAat<-anoAatB1[3,2]
  eta2AatB1<-SSAatB1/SST
  f2AatB1<-eta2AatB1/(1-eta2AatB1)
  lambdaAatB1<-f2AatB1*dfwinSE
  FtAatB1<-stats::qf(minusalpha, dfAatB1, dfwinSE)
  power.AatB1<-round(1-stats::pf(FtAatB1, dfAatB1,dfwinSE,lambdaAatB1),3)

  SSAatB2<-anoAatB2[2,1]
  dfAatB2<-anoAatB2[2,2]
  eta2AatB2<-SSAatB2/SST
  f2AatB2<-eta2AatB2/(1-eta2AatB2)
  lambdaAatB2<-f2AatB2*dfwinSE
  FtAatB2<-stats::qf(minusalpha, dfAatB2, dfwinSE)
  power.AatB2<-round(1-stats::pf(FtAatB2, dfAatB2,dfwinSE,lambdaAatB2),3)

  message("Simple Effect Comparing M = ",m1.1, " and M = ", m2.1,". Power = ", power.AatB1)
  message("Simple Effect Comparing M= ",m1.2, " and M = ", m2.2,". Power = ", power.AatB2)
  message("Simple Effect Comparing M = ",m1.1, " and M = ", m1.2,". Power = ", power.BatA1)
  message("Simple Effect Comparing M = ",m2.1, " and M = ", m2.2,". Power = ", power.BatA2)

  result <- data.frame(matrix(ncol = 8))
  colnames(result) <- c("Eta-squared A at B1","Power A at B1","Eta-squared A at B2","Power A at B2","Eta-squared B at A1","Power B at A1","Eta-squared B at A2","Power B at A2")
  result[, 1]<-eta2AatB1
  result[, 2]<-power.AatB1
  result[, 3]<-eta2AatB2
  result[, 4]<-power.AatB2
  result[, 5]<-eta2BatA1
  result[, 6]<-power.BatA1
  result[, 7]<-eta2BatA2
  result[, 8]<-power.BatA2
  output<-na.omit(result)
  rownames(output)<- c()
  invisible(output)
  }
