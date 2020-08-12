#'Compute Power for One or Two Factor ANCOVA with a single covariate
#'Takes means, sds, and sample sizes for each group. Alpha is .05 by default, alternative values may be entered by user
#'Factor A can have up to four levels, Factor B, if used, can only be two
#'@param m1.1 Cell mean for First level of Factor A, First level of Factor B
#'@param m2.1 Cell mean for Second level of Factor A, First level of Factor B
#'@param m3.1 Cell mean for Third level of Factor A, First level of Factor B
#'@param m4.1 Cell mean for Fourth level of Factor A, First level of Factor B
#'@param m1.2 Cell mean for First level of Factor A, Second level of Factor B
#'@param m2.2 Cell mean for Second level of Factor A, Second level of Factor B
#'@param m3.2 Cell mean for Third level of Factor A, Second level of Factor B
#'@param m4.2 Cell mean for Fourth level of Factor A, Second level of Factor B
#'@param s1.1 Cell standard deviation for First level of Factor A, First level of Factor B
#'@param s2.1 Cell standard deviation for Second level of Factor A, First level of Factor B
#'@param s3.1 Cell standard deviation for Third level of Factor A, First level of Factor B
#'@param s4.1 Cell standard deviation for Fourth level of Factor A, First level of Factor B
#'@param s1.2 Cell standard deviation for First level of Factor A, Second level of Factor B
#'@param s2.2 Cell standard deviation for Second level of Factor A, Second level of Factor B
#'@param s3.2 Cell standard deviation for Third level of Factor A, Second level of Factor B
#'@param s4.2 Cell standard deviation for Fourth level of Factor A, Second level of Factor B
#'@param s Overall standard deviation. Sets all cell sds equal
#'@param r Correlation between covariate and dependent variable.
#'@param n Sample Size per cell
#'@param factors Number of factors (1 or 2)
#'@param alpha Type I error (default is .05)
#'@param levelsA levels for factor A (up to four)
#'@importFrom car Anova
#'@examples
#' anc(m1.1=.85,m2.1=2.5, s1.1 = 1.7, s2.1=1,
#' m1.2=0.85, m2.2= 2.5, s1.2 = 1.7, s2.2=1,
#' m3.1=0.0,m3.2=2.5, s3.1 = 1.7, s3.2=1,
#' m4.1=0.6, m4.2 = 2.5, s4.1 = 1.7, s4.2=1, r= 0.4,
#' n=251, factors =2,levelsA = 4)
#'@return Power for One or Two Factor ANCOVA with a single covariate
#'@export

anc=function(m1.1,m2.1,m1.2,m2.2,m3.1=NULL,m3.2=NULL,m4.1=NULL,m4.2=NULL,
             s1.1=NULL,s2.1=NULL,s1.2=NULL,s2.2=NULL,s3.1=NULL,s3.2=NULL,s4.1=NULL,s4.2=NULL,
             r,s=NULL,alpha=.05,factors,levelsA=NULL,n){
{
V1<-V2<-ivbg<-iv1<-iv2<-NULL
oldoption<-options(contrasts=c("contr.helmert", "contr.poly"))
oldoption
on.exit(options(oldoption))

if (factors=="1" & levelsA=="2")
{
var1<-s1.1^2; var2<-s2.1^2;varcov<-1
#Note - all s for covariate are set to 1
cov1<-r*s1.1
cov2<-r*s2.1
mcov<-0

out1 <- MASS::mvrnorm(n, mu = c(m1.1,mcov),
                Sigma = matrix(c(var1,cov1,
                                 cov1,varcov), ncol = 2),
                empirical = TRUE)
out1<-as.data.frame(out1)
out1$ivbg<-NA
out1$ivbg<-1 #identifies group

out2 <- MASS::mvrnorm(n, mu = c(m2.1,mcov),
                Sigma = matrix(c(var2,cov2,
                                 cov2,var2), ncol = 2),
                empirical = TRUE)
out2<-as.data.frame(out2)
out2$ivbg<-NA
out2$ivbg<-2

out<-rbind(out1,out2)
out<-as.data.frame(out)
out<-dplyr::rename(out, y1 = V1, cov = V2, ivbg = ivbg)
out$ivbg<-as.factor(out$ivbg)

anc<-stats::aov(y1~cov+ivbg, data=out)
sum<-car::Anova(anc, type="III")
SSA<-sum[3,1] #column, row
SSwin<-sum[4,1]
dfA<-sum[3,2]
dfwin<-sum[4,2]
MSwin<-SSwin/dfwin
eta2A<-SSA/(SSA+SSwin)
f2A<-eta2A/(1-eta2A)
lambdaA<-f2A*dfwin
minusalpha<-1-alpha
FtA<-stats::qf(minusalpha, dfA, dfwin)
power.A<-round(1-stats::pf(FtA, dfA,dfwin,lambdaA),4)
nall<-n*2
eta2A<-round((eta2A),3)
message("Sample size per cell = ",n)

message("Power IV1 = ", power.A, " for partial eta-squared = ", eta2A)
result <- data.frame(matrix(ncol = 3))
colnames(result) <- c("n", "Eta-squared IV1","Power IV1")
result[n, 1]<-n
result[n, 2]<-eta2A
result[n, 3]<-power.A
output<-na.omit(result)
rownames(output)<- c()
}

if (factors=="1" & levelsA=="3")
{
  cov1<-r*s1.1
  cov2<-r*s2.1
  mcov<-0
  var1<-s1.1^2; var2<-s2.1^2;var3<-s3.1^2;varcov<-1

  out1 <- MASS::mvrnorm(n, mu = c(m1.1,mcov),
                        Sigma = matrix(c(var1,cov1,
                                         cov1,varcov), ncol = 2),
                        empirical = TRUE)
  out1<-as.data.frame(out1)
  out1$ivbg<-NA
  out1$ivbg<-1 #identifies group

  out2 <- MASS::mvrnorm(n, mu = c(m2.1,mcov),
                        Sigma = matrix(c(var2,cov2,
                                         cov2,var2), ncol = 2),
                        empirical = TRUE)
  out2<-as.data.frame(out2)
  out2$ivbg<-NA
  out2$ivbg<-2

  out3 <- MASS::mvrnorm(n, mu = c(m3.1,mcov),
                        Sigma = matrix(c(var3,cov3,
                                         cov3,var3), ncol = 2),
                        empirical = TRUE)
  out3<-as.data.frame(out2)
  out3$ivbg<-NA
  out3$ivbg<-3


  out<-rbind(out1,out2,out3)
  out<-as.data.frame(out)
  out<-dplyr::rename(out, y1 = V1, cov = V2, ivbg = ivbg)
  out$ivbg<-as.factor(out$ivbg)

  anc<-stats::aov(y1~cov+ivbg, data=out)
  sum<-car::Anova(anc, type="III")
  SSA<-sum[3,1] #column, row
  SSwin<-sum[4,1]
  dfA<-sum[3,2]
  dfwin<-sum[4,2]
  MSwin<-SSwin/dfwin
  eta2A<-SSA/(SSA+SSwin)
  f2A<-eta2A/(1-eta2A)
  lambdaA<-f2A*dfwin
  minusalpha<-1-alpha
  FtA<-stats::qf(minusalpha, dfA, dfwin)
  power.A<-round(1-stats::pf(FtA, dfA,dfwin,lambdaA),4)
  nall<-n*2
  eta2A<-round((eta2A),3)
  message("Sample size per cell = ",n)

  message("Power IV1 = ", power.A, " for partial eta-squared = ", eta2A)
  result <- data.frame(matrix(ncol = 3))
  colnames(result) <- c("n cell", "Eta-squared IV1","Power IV1")
  result[n, 1]<-n
  result[n, 2]<-eta2A
  result[n, 3]<-power.A
  output<-na.omit(result)
  rownames(output)<- c()

}
if (factors=="1" & levelsA=="4")
{
  cov1<-r*s1.1
  cov2<-r*s2.1
  cov3<-r*s3.1
  cov4<-r*s4.1
  mcov<-0
  var1<-s1.1^2; var2<-s2.1^2;var3<-s3.1^2;var4<-s4.1^2;varcov<-1


  out1 <- MASS::mvrnorm(n, mu = c(m1.1,mcov),
                        Sigma = matrix(c(var1,cov1,
                                         cov1,varcov), ncol = 2),
                        empirical = TRUE)
  out1<-as.data.frame(out1)
  out1$ivbg<-NA
  out1$ivbg<-1 #identifies group

  out2 <- MASS::mvrnorm(n, mu = c(m2.1,mcov),
                        Sigma = matrix(c(var2,cov2,
                                         cov2,var2), ncol = 2),
                        empirical = TRUE)
  out2<-as.data.frame(out2)
  out2$ivbg<-NA
  out2$ivbg<-2

  out3 <- MASS::mvrnorm(n, mu = c(m3.1,mcov),
                        Sigma = matrix(c(var3,cov3,
                                         cov3,var3), ncol = 2),
                        empirical = TRUE)
  out3<-as.data.frame(out2)
  out3$ivbg<-NA
  out3$ivbg<-3


  out4 <- MASS::mvrnorm(n, mu = c(m4.1,mcov),
                        Sigma = matrix(c(var4,cov4,
                                         cov4,var4), ncol = 2),
                        empirical = TRUE)
  out4<-as.data.frame(out2)
  out4$ivbg<-NA
  out4$ivbg<-4

  out<-rbind(out1,out2,out3,out4)
  out<-as.data.frame(out)
  out<-dplyr::rename(out, y1 = V1, cov = V2, ivbg = ivbg)
  out$ivbg<-as.factor(out$ivbg)

  anc<-stats::aov(y1~cov+ivbg, data=out)
  sum<-car::Anova(anc, type="III")
  SSA<-sum[3,1] #column, row
  SSwin<-sum[4,1]
  dfA<-sum[3,2]
  dfwin<-sum[4,2]
  MSwin<-SSwin/dfwin
  eta2A<-SSA/(SSA+SSwin)
  f2A<-eta2A/(1-eta2A)
  lambdaA<-f2A*dfwin
  minusalpha<-1-alpha
  FtA<-stats::qf(minusalpha, dfA, dfwin)
  power.A<-round(1-stats::pf(FtA, dfA,dfwin,lambdaA),4)
  nall<-n*2
  eta2A<-round((eta2A),3)
  message("Sample size per cell = ",n)

  message("Power IV1 = ", power.A, " for partial eta-squared = ", eta2A)
  result <- data.frame(matrix(ncol = 3))
  colnames(result) <- c("n cell", "Eta-squared IV1","Power IV1")
  result[n, 1]<-n
  result[n, 2]<-eta2A
  result[n, 3]<-power.A
  output<-na.omit(result)
  rownames(output)<- c()

}

if (factors==2 & levelsA==2)
{
var1<-s1.1^2; var2<-s2.1^2;var3<-s1.2^2; var4<-s2.2^2;varcov<-1
#Note - all s for covariate are set to 1
cov1<-r*s1.1
cov2<-r*s2.1
cov3<-r*s1.2
cov4<-r*s2.2
mcov<-0

out1 <- MASS::mvrnorm(n, mu = c(m1.1,mcov),
                      Sigma = matrix(c(var1,cov1,
                                       cov1,varcov), ncol = 2),
                      empirical = TRUE)
out1<-as.data.frame(out1)
out1$iv1<-NA
out1$iv2<-NA
out1$iv1<-1
out1$iv2<-1 #identifies group

out2 <- MASS::mvrnorm(n, mu = c(m1.2,mcov),
                      Sigma = matrix(c(var2,cov2,
                                       cov2,varcov), ncol = 2),
                      empirical = TRUE)
out2<-as.data.frame(out2)
out2$iv1<-NA
out2$iv2<-NA
out2$iv1<-2
out2$iv2<-1

out3 <- MASS::mvrnorm(n, mu = c(m2.1,mcov),
                      Sigma = matrix(c(var3,cov3,
                                       cov3,varcov), ncol = 2),
                      empirical = TRUE)
out3<-as.data.frame(out3)
out3$iv1<-NA
out3$iv2<-NA
out3$iv1<-1
out3$iv2<-2

out4 <- MASS::mvrnorm(n, mu = c(m2.2,mcov),
                      Sigma = matrix(c(var4,cov4,
                                       cov4,varcov), ncol = 2),
                      empirical = TRUE)
out4<-as.data.frame(out4)
out4$iv1<-NA
out4$iv2<-NA
out4$iv1<-2
out4$iv2<-2

out<-rbind(out1,out2,out3,out4)
out<-as.data.frame(out)
out<-dplyr::rename(out, y1 = V1, cov = V2, iv1 = iv1,iv2=iv2)
out$iv1<-as.factor(out$iv1)
out$iv2<-as.factor(out$iv2)
anc<-stats::aov(y1~cov+iv1*iv2, data=out)
sum<-car::Anova(anc, type="III")
SSA<-sum[3,1] #column, row
SSB<-sum[4,1]
SSAB<-sum[5,1]
SSwin<-sum[6,1]
dfA<-sum[3,2]
dfB<-sum[4,2]
dfAB<-sum[5,2]
dfwin<-sum[6,2]
MSwin<-SSwin/dfwin
eta2A<-SSA/(SSA+SSwin)
eta2B<-SSB/(SSB+SSwin)
eta2AB<-SSAB/(SSAB+SSwin)
f2A<-eta2A/(1-eta2A)
f2B<-eta2B/(1-eta2B)
f2AB<-eta2AB/(1-eta2AB)
lambdaA<-f2A*dfwin
lambdaB<-f2B*dfwin
lambdaAB<-f2AB*dfwin
minusalpha<-1-alpha
FtA<-stats::qf(minusalpha, dfA, dfwin)
FtB<-stats::qf(minusalpha, dfB, dfwin)
FtAB<-stats::qf(minusalpha, dfAB, dfwin)
power.A<-round(1-stats::pf(FtA, dfA,dfwin,lambdaA),4)
power.B<-round(1-stats::pf(FtB, dfB,dfwin,lambdaB),4)
power.AB<-round(1-stats::pf(FtAB, dfAB,dfwin,lambdaAB),4)
eta2A<-round((eta2A),3)
eta2B<-round((eta2B),3)
eta2AB<-round((eta2AB),3)
message("Sample size per cell = ",n)
message("Power Factor A = ", power.A, " for partial eta-squared = ", eta2A)
message("Power Factor B = ", power.B, " for partial eta-squared = ", eta2B)
message("Power Factor AxB = ", power.AB, " for partial eta-squared = ", eta2AB)
result <- data.frame(matrix(ncol = 7))
colnames(result) <- c("n", "Eta-squared IV1","Power IV1","Eta-squared IV2","Power IV2","Power Interaction","Eta-squared Interaction")
result[n, 1]<-n
result[n, 2]<-eta2A
result[n, 3]<-power.A
result[n, 4]<-eta2B
result[n, 5]<-power.B
result[n, 6]<-eta2AB
result[n, 7]<-power.AB

output<-na.omit(result)
rownames(output)<- c()


}

if (factors==2 & levelsA==3)
{
  var1<-s1.1^2; var2<-s2.1^2;var3<-s1.2^2;
  var4<-s2.2^2;var5<-s3.1^2;var6<-s3.2^2;varcov<-1
  #Note - all s for covariate are set to 1
  cov1<-r*s1.1
  cov2<-r*s2.1
  cov3<-r*s1.2
  cov4<-r*s2.2
  cov5<-r*s3.1
  cov6<-r*s3.2
  mcov<-0

  out1 <- MASS::mvrnorm(n, mu = c(m1.1,mcov),
                        Sigma = matrix(c(var1,cov1,
                                         cov1,varcov), ncol = 2),
                        empirical = TRUE)
  out1<-as.data.frame(out1)
  out1$iv1<-NA
  out1$iv2<-NA
  out1$iv1<-1
  out1$iv2<-1 #identifies group

  out2 <- MASS::mvrnorm(n, mu = c(m2.1,mcov),
                        Sigma = matrix(c(var2,cov2,
                                         cov2,varcov), ncol = 2),
                        empirical = TRUE)
  out2<-as.data.frame(out2)
  out2$iv1<-NA
  out2$iv2<-NA
  out2$iv1<-2
  out2$iv2<-1

  out3 <- MASS::mvrnorm(n, mu = c(m1.2,mcov),
                        Sigma = matrix(c(var3,cov3,
                                         cov3,varcov), ncol = 2),
                        empirical = TRUE)
  out3<-as.data.frame(out3)
  out3$iv1<-NA
  out3$iv2<-NA
  out3$iv1<-1
  out3$iv2<-2

  out4 <- MASS::mvrnorm(n, mu = c(m2.2,mcov),
                        Sigma = matrix(c(var4,cov4,
                                         cov4,varcov), ncol = 2),
                        empirical = TRUE)
  out4<-as.data.frame(out4)
  out4$iv1<-NA
  out4$iv2<-NA
  out4$iv1<-2
  out4$iv2<-2

  out5 <- MASS::mvrnorm(n, mu = c(m3.1,mcov),
                        Sigma = matrix(c(var5,cov5,
                                         cov5,varcov), ncol = 2),
                        empirical = TRUE)
  out5<-as.data.frame(out5)
  out5$iv1<-NA
  out5$iv2<-NA
  out5$iv1<-3
  out5$iv2<-1

  out6 <- MASS::mvrnorm(n, mu = c(m3.2,mcov),
                        Sigma = matrix(c(var6,cov6,
                                         cov6,varcov), ncol = 2),
                        empirical = TRUE)
  out6<-as.data.frame(out6)
  out6$iv1<-NA
  out6$iv2<-NA
  out6$iv1<-3
  out6$iv2<-2


  out<-rbind(out1,out2,out3,out4,out5,out6)
  out<-as.data.frame(out)
  out<-dplyr::rename(out, y1 = V1, cov = V2, iv1 = iv1,iv2=iv2)
  out$iv1<-as.factor(out$iv1)
  out$iv2<-as.factor(out$iv2)
  anc<-stats::aov(y1~cov+iv1*iv2, data=out)
  sum<-car::Anova(anc, type="III")
  SSA<-sum[3,1] #column, row
  SSB<-sum[4,1]
  SSAB<-sum[5,1]
  SSwin<-sum[6,1]
  dfA<-sum[3,2]
  dfB<-sum[4,2]
  dfAB<-sum[5,2]
  dfwin<-sum[6,2]
  MSwin<-SSwin/dfwin
  eta2A<-SSA/(SSA+SSwin)
  eta2B<-SSB/(SSB+SSwin)
  eta2AB<-SSAB/(SSAB+SSwin)
  f2A<-eta2A/(1-eta2A)
  f2B<-eta2B/(1-eta2B)
  f2AB<-eta2AB/(1-eta2AB)
  lambdaA<-f2A*dfwin
  lambdaB<-f2B*dfwin
  lambdaAB<-f2AB*dfwin
  minusalpha<-1-alpha
  FtA<-stats::qf(minusalpha, dfA, dfwin)
  FtB<-stats::qf(minusalpha, dfB, dfwin)
  FtAB<-stats::qf(minusalpha, dfAB, dfwin)
  power.A<-round(1-stats::pf(FtA, dfA,dfwin,lambdaA),4)
  power.B<-round(1-stats::pf(FtB, dfB,dfwin,lambdaB),4)
  power.AB<-round(1-stats::pf(FtAB, dfAB,dfwin,lambdaAB),4)
  eta2A<-round((eta2A),3)
  eta2B<-round((eta2B),3)
  eta2AB<-round((eta2AB),3)
  message("Sample size per cell = ",n)
  message("Power Factor A = ", power.A, " for partial eta-squared = ", eta2A)
  message("Power Factor B = ", power.B, " for partial eta-squared = ", eta2B)
  message("Power Factor AxB = ", power.AB, " for partial eta-squared = ", eta2AB)
  result <- data.frame(matrix(ncol = 7))
  colnames(result) <- c("n", "Eta-squared IV1","Power IV1","Eta-squared IV2","Power IV2","Power Interaction","Eta-squared Interaction")
  result[n, 1]<-n
  result[n, 2]<-eta2A
  result[n, 3]<-power.A
  result[n, 4]<-eta2B
  result[n, 5]<-power.B
  result[n, 6]<-eta2AB
  result[n, 7]<-power.AB

  output<-na.omit(result)
  rownames(output)<- c()


}

if (factors==2 & levelsA==4)
{
  var1<-s1.1^2; var2<-s2.1^2;var3<-s1.2^2;
  var4<-s2.2^2;var5<-s3.1^2;var6<-s3.2^2;
  var7<-s4.1^2; var8<-s4.2^2; varcov<-1
  #Note - all s for covariate are set to 1
  cov1<-r*s1.1
  cov2<-r*s2.1
  cov3<-r*s1.2
  cov4<-r*s2.2
  cov5<-r*s3.1
  cov6<-r*s3.2
  cov7<-r*s4.1
  cov8<-r*s4.2
  mcov<-0

  out1 <- MASS::mvrnorm(n, mu = c(m1.1,mcov),
                        Sigma = matrix(c(var1,cov1,
                                         cov1,varcov), ncol = 2),
                        empirical = TRUE)
  out1<-as.data.frame(out1)
  out1$iv1<-NA
  out1$iv2<-NA
  out1$iv1<-1
  out1$iv2<-1 #identifies group

  out2 <- MASS::mvrnorm(n, mu = c(m2.1,mcov),
                        Sigma = matrix(c(var2,cov2,
                                         cov2,varcov), ncol = 2),
                        empirical = TRUE)
  out2<-as.data.frame(out2)
  out2$iv1<-NA
  out2$iv2<-NA
  out2$iv1<-2
  out2$iv2<-1

  out3 <- MASS::mvrnorm(n, mu = c(m1.2,mcov),
                        Sigma = matrix(c(var3,cov3,
                                         cov3,varcov), ncol = 2),
                        empirical = TRUE)
  out3<-as.data.frame(out3)
  out3$iv1<-NA
  out3$iv2<-NA
  out3$iv1<-1
  out3$iv2<-2

  out4 <- MASS::mvrnorm(n, mu = c(m2.2,mcov),
                        Sigma = matrix(c(var4,cov4,
                                         cov4,varcov), ncol = 2),
                        empirical = TRUE)
  out4<-as.data.frame(out4)
  out4$iv1<-NA
  out4$iv2<-NA
  out4$iv1<-2
  out4$iv2<-2

  out5 <- MASS::mvrnorm(n, mu = c(m3.1,mcov),
                        Sigma = matrix(c(var5,cov5,
                                         cov5,varcov), ncol = 2),
                        empirical = TRUE)
  out5<-as.data.frame(out5)
  out5$iv1<-NA
  out5$iv2<-NA
  out5$iv1<-3
  out5$iv2<-1

  out6 <- MASS::mvrnorm(n, mu = c(m3.2,mcov),
                        Sigma = matrix(c(var6,cov6,
                                         cov6,varcov), ncol = 2),
                        empirical = TRUE)
  out6<-as.data.frame(out6)
  out6$iv1<-NA
  out6$iv2<-NA
  out6$iv1<-3
  out6$iv2<-2


  out7 <- MASS::mvrnorm(n, mu = c(m4.1,mcov),
                        Sigma = matrix(c(var7,cov7,
                                         cov7,varcov), ncol = 2),
                        empirical = TRUE)
  out7<-as.data.frame(out7)
  out7$iv1<-NA
  out7$iv2<-NA
  out7$iv1<-4
  out7$iv2<-1

  out8 <- MASS::mvrnorm(n, mu = c(m4.2,mcov),
                        Sigma = matrix(c(var8,cov8,
                                         cov8,varcov), ncol = 2),
                        empirical = TRUE)
  out8<-as.data.frame(out8)
  out8$iv1<-NA
  out8$iv2<-NA
  out8$iv1<-4
  out8$iv2<-2

  out<-rbind(out1,out2,out3,out4,out5,out6,out7,out8)
  out<-as.data.frame(out)
  out<-dplyr::rename(out, y1 = V1, cov = V2, iv1 = iv1,iv2=iv2)
  out$iv1<-as.factor(out$iv1)
  out$iv2<-as.factor(out$iv2)
  anc<-stats::aov(y1~cov+iv1*iv2, data=out)
  sum<-car::Anova(anc, type="III")
  SSA<-sum[3,1] #column, row
  SSB<-sum[4,1]
  SSAB<-sum[5,1]
  SSwin<-sum[6,1]
  dfA<-sum[3,2]
  dfB<-sum[4,2]
  dfAB<-sum[5,2]
  dfwin<-sum[6,2]
  MSwin<-SSwin/dfwin
  eta2A<-SSA/(SSA+SSwin)
  eta2B<-SSB/(SSB+SSwin)
  eta2AB<-SSAB/(SSAB+SSwin)
  f2A<-eta2A/(1-eta2A)
  f2B<-eta2B/(1-eta2B)
  f2AB<-eta2AB/(1-eta2AB)
  lambdaA<-f2A*dfwin
  lambdaB<-f2B*dfwin
  lambdaAB<-f2AB*dfwin
  minusalpha<-1-alpha
  FtA<-stats::qf(minusalpha, dfA, dfwin)
  FtB<-stats::qf(minusalpha, dfB, dfwin)
  FtAB<-stats::qf(minusalpha, dfAB, dfwin)
  power.A<-round(1-stats::pf(FtA, dfA,dfwin,lambdaA),4)
  power.B<-round(1-stats::pf(FtB, dfB,dfwin,lambdaB),4)
  power.AB<-round(1-stats::pf(FtAB, dfAB,dfwin,lambdaAB),4)
  eta2A<-round((eta2A),3)
  eta2B<-round((eta2B),3)
  eta2AB<-round((eta2AB),3)
  message("Sample size per cell = ",n)
  message("Power Factor A = ", power.A, " for partial eta-squared = ", eta2A)
  message("Power Factor B = ", power.B, " for partial eta-squared = ", eta2B)
  message("Power Factor AxB = ", power.AB, " for partial eta-squared = ", eta2AB)
  result <- data.frame(matrix(ncol = 7))
  colnames(result) <- c("n", "Eta-squared IV1","Power IV1","Eta-squared IV2","Power IV2","Power Interaction","Eta-squared Interaction")
  result[n, 1]<-n
  result[n, 2]<-eta2A
  result[n, 3]<-power.A
  result[n, 4]<-eta2B
  result[n, 5]<-power.B
  result[n, 6]<-eta2AB
  result[n, 7]<-power.AB

  output<-na.omit(result)
  rownames(output)<- c()


  }

{
invisible(output)
}}}
