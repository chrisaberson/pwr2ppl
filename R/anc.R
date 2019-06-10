#'Compute Power for One or Two Factor ANCOVA with a single covariate
#'Takes means, sds, and sample sizes for each group. Alpha is .05 by default, alternative values may be entered by user
#'@param m1.1 Cell mean for First level of Factor A, First level of Factor B
#'@param m1.2 Cell mean for First level of Factor A, Second level of Factor B
#'@param m1.3 Cell mean for First level of Factor A, Third level of Factor B
#'@param m1.4 Cell mean for First level of Factor A, Fourth level of Factor B
#'@param m2.1 Cell mean for Second level of Factor A, First level of Factor B
#'@param m2.2 Cell mean for Second level of Factor A, Second level of Factor B
#'@param m2.3 Cell mean for Second level of Factor A, Third level of Factor B
#'@param m2.4 Cell mean for Second level of Factor A, Fourth level of Factor B
#'@param s1.1 Cell standard deviation for First level of Factor A, First level of Factor B
#'@param s1.2 Cell standard deviation for First level of Factor A, Second level of Factor B
#'@param s1.3 Cell standard deviation for First level of Factor A, Third level of Factor B
#'@param s1.4 Cell standard deviation for First level of Factor A, Fourth level of Factor B
#'@param s2.1 Cell standard deviation for Second level of Factor A, First level of Factor B
#'@param s2.2 Cell standard deviation for Second level of Factor A, Second level of Factor B
#'@param s2.3 Cell standard deviation for Second level of Factor A, Third level of Factor B
#'@param s2.4 Cell standard deviation for Second level of Factor A, Fourth level of Factor B
#'@param s Overall standard deviation. Sets all cell sds equal
#'@param r Correlation between covariate and dependent variable.
#'@param n Sample Size per cell
#'@param factors Number of factors (1 or 2)
#'@param alpha Type I error (default is .05)
#'@importFrom car Anova
#'@examples
#' anc(m1.1=.85,m2.1=2.5, s1.1 = 1.7, s2.1=1,
#' m1.2=0.85, m2.2= 2.5, s1.2 = 1.7, s2.2=1,
#' m1.3=0.0,m2.3=2.5, s1.3 = 1.7, s2.3=1,
#' m1.4=0.6, m2.4 = 2.5, s1.4 = 1.7, s2.4=1, r= 0.4,
#' n=251, factors =2)
#'@return Power for One or Two Factor ANCOVA with a single covariate
#'@export

anc=function(m1.1,m2.1,m1.2,m2.2,m1.3=NULL,m2.3=NULL,m1.4=NULL,m2.4=NULL,
             s1.1=NULL,s2.1=NULL,s1.2=NULL,s2.2=NULL,s1.3=NULL,s2.3=NULL,s1.4=NULL,s2.4=NULL,
             r,s=NULL,alpha=.05,factors,n){
{
V1<-V2<-ivbg<-iv1<-iv2<-NULL
oldoption<-options(contrasts=c("contr.helmert", "contr.poly"))
oldoption
on.exit(options(oldoption))

if (factors=="1")
{
var1<-s1.1^2; var2<-s2.1^2;var3<-s1.2^2;var4<-s2.2^2
cov12<-r*s1.1*s2.1
cov34<-r*s1.2*s2.2

out1 <- MASS::mvrnorm(n, mu = c(m1.1,m2.1),
                Sigma = matrix(c(var1,cov12,
                                 cov12,var2), ncol = 2),
                empirical = TRUE)
out1<-as.data.frame(out1)
out1$ivbg<-NA
out1$ivbg<-1 #identifies group

out2 <- MASS::mvrnorm(n, mu = c(m1.2,m2.2),
                Sigma = matrix(c(var3,cov34,
                                 cov34,var4), ncol = 2),
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
message("Sample size overall = ",nall)
message("Power IV1 = ", power.A, " for eta-squared = ", eta2A)
result <- data.frame(matrix(ncol = 3))
colnames(result) <- c("n", "Eta-squared IV1","Power IV1")
result[n, 1]<-nall
result[n, 2]<-eta2A
result[n, 3]<-power.A
output<-na.omit(result)
rownames(output)<- c()


}

if (factors==2)
{

  var1<-s1.1^2; var2<-s2.1^2;var3<-s1.2^2;var4<-s2.2^2
  var5<-s1.3^2; var6<-s2.3^2;var7<-s1.4^2;var8<-s2.4^2
  cov12<-r*s1.1*s2.1
  cov34<-r*s1.2*s2.2
  cov56<-r*s1.3*s2.3
  cov78<-r*s1.4*s2.4

out1 <- MASS::mvrnorm(n, mu = c(m1.1,m2.1),
                Sigma = matrix(c(var1,cov12,
                                 cov12,var2), ncol = 2),
                empirical = TRUE)
out1<-as.data.frame(out1)
out1$iv1<-NA
out1$iv1<-1 #identifies group
out1$iv2<-NA
out1$iv2<-1 #identifies group

out2 <- MASS::mvrnorm(n, mu = c(m1.2,m2.2),
                Sigma = matrix(c(var3,cov34,
                                 cov34,var4), ncol = 2),
                empirical = TRUE)
out2<-as.data.frame(out2)
out2$iv1<-NA
out2$iv1<-1
out2$iv2<-NA
out2$iv2<-2 #identifies group

out3 <- MASS::mvrnorm(n, mu = c(m1.3,m2.3),
                Sigma = matrix(c(var5,cov56,
                                 cov56,var6), ncol = 2),
                empirical = TRUE)
out3<-as.data.frame(out3)
out3$iv1<-NA
out3$iv1<-2 #identifies group
out3$iv2<-NA
out3$iv2<-1 #identifies group

out4 <- MASS::mvrnorm(n, mu = c(m1.4,m2.4),
                Sigma = matrix(c(var7,cov78,
                                 cov78,var8), ncol = 2),
                empirical = TRUE)
out4<-as.data.frame(out4)
out4$iv1<-NA
out4$iv1<-2 #identifies group
out4$iv2<-NA
out4$iv2<-2 #identifies group


out<-rbind(out1,out2,out3,out4)
out<-as.data.frame(out)
out<-dplyr::rename(out, y1 = V1, cov = V2, iv1 = iv1, iv2=iv2)
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
f2A<-eta2A/(1-eta2A)
lambdaA<-f2A*dfwin
minusalpha<-1-alpha
FtA<-stats::qf(minusalpha, dfA, dfwin)
power.A<-round(1-stats::pf(FtA, dfA,dfwin,lambdaA),4)
eta2B<-SSB/(SSB+SSwin)
f2B<-eta2B/(1-eta2B)
lambdaB<-f2B*dfwin
FtB<-stats::qf(minusalpha, dfB, dfwin)
power.B<-round(1-stats::pf(FtB, dfB,dfwin,lambdaB),4)
eta2AB<-SSAB/(SSAB+SSwin)
f2AB<-eta2AB/(1-eta2AB)
lambdaAB<-f2AB*dfwin
FtAB<-stats::qf(minusalpha,dfAB, dfwin)
power.AB<-round(1-stats::pf(FtAB,dfAB,dfwin,lambdaAB),4)
eta2A<-round((eta2A),3)
eta2B<-round((eta2B),3)
eta2AB<-round((eta2AB),3)
nall<-n*4

message("Sample size per cell = ",n)
message("Sample size overall = ",nall)
message("Power IV1 = ", power.A, " for eta-squared = ", eta2A)
message("Power IV2 = ", power.B, " for eta-squared = ", eta2B)
message("Power IV1*IV2 = ", power.AB, " for eta-squared = ", eta2AB)


#Return results silently. Available if written to object
result <- data.frame(matrix(ncol = 7))
colnames(result) <- c("n", "Eta-squared IV1","Power IV1", "Eta-squared IV2", "Power IV2","Eta-squared IV1*IV2", "Power IV1*IV2")
result[n, 1]<-nall
result[n, 2]<-eta2A
result[n, 3]<-power.A
result[n, 4]<-eta2B
result[n, 5]<-power.B
result[n, 6]<-eta2AB
result[n, 7]<-power.AB
output<-na.omit(result)
rownames(output)<- c()
#invisible(output2)

}
invisible(output)
}}
