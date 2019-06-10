#'Power for Comparing Independent R2 in Multiple Regression with Two or Three Predictors
#'Requires correlations between all variables as sample size. Means, sds, and alpha are option. Also computes Power(All)
#'@param ry1_1 Correlation between DV (y) and first predictor (1), first test
#'@param ry2_1 Correlation between DV (y) and second predictor (2), first test
#'@param ry3_1 Correlation between DV (y) and third predictor (3), first test
#'@param r12_1 Correlation between first (1) and second predictor (2), first test
#'@param r13_1 Correlation between first (1) and third predictor (3), first test
#'@param r23_1 Correlation between second (2) and third predictor (3), first test
#'@param ry1_2 Correlation between DV (y) and first predictor (1), second test
#'@param ry2_2 Correlation between DV (y) and second predictor (2), second test
#'@param ry3_2 Correlation between DV (y) and third predictor (3), second test
#'@param r12_2 Correlation between first (1) and second predictor (2), second test
#'@param r13_2 Correlation between first (1) and third predictor (3), second test
#'@param r23_2 Correlation between second (2) and third predictor (3), second test
#'@param n1 Sample size first test
#'@param n2 Sample size second test
#'@param alpha Type I error (default is .05)
#'@param tails number of tails for test (default is 2)
#'@examples
#'indR2(ry1_1=.40, ry2_1=.40, ry3_1 =-.40, r12_1=-.15,r13_1=-.60, r23_1=.25,
#'ry1_2=.40, ry2_2=.10, ry3_2 =-.40, r12_2=-.15,r13_2=-.60, r23_2=.25,
#'n1=115,n2=115, alpha=.05)
#'@return Power for Comparing R2 Coefficients in Multiple Regression
#'@export
#'
#'
indR2<-function(ry1_1, ry2_1, ry3_1=NULL, r12_1, r13_1=NULL, r23_1=NULL,n1,
                   ry1_2, ry2_2, ry3_2=NULL, r12_2, r13_2=NULL, r23_2=NULL,n2, alpha=.05, tails=2)
{
  pred<-NA
  pred[is.null(r23_1)]<-2
  pred[!is.null(r23_1)]<-3

  if (pred=="2")
  {pop1 <- MASS::mvrnorm(n1, mu<-c(0,0,0), Sigma<-matrix(c(1, ry1_1, ry2_1,
                                                   ry1_1, 1, r12_1,
                                                   ry2_1, r12_1, 1),
                                                 ncol=3), empirical=TRUE)
  pop1<-data.frame(pop1)

  pop2 <- MASS::mvrnorm(n2, mu<-c(0,0,0), Sigma<-matrix(c(1, ry1_2, ry2_2,
                                                   ry1_2, 1, r12_2,
                                                   ry2_2, r12_2, 1),
                                                 ncol=3), empirical=TRUE)
  pop2<-data.frame(pop2)


values1<-stats::lm(X1~X2+X3, pop1)
values1<-summary(values1)
R2_1<-(values1$r.squared)
values2<-stats::lm(X1~X2+X3, pop2)
values2<-summary(values2)
R2_2<-(values2$r.squared)
SER2_1<-((4*R2_1)*(1-R2_1)^2)*((n1-pred-1)^2) / ((n1^2 - 1)* (n1+3))
SER2_2<-((4*R2_2)*(1-R2_2)^2)*((n2-pred-1)^2) / ((n2^2 - 1)* (n2+3))
SER2<-(SER2_1 + SER2_2)^.5
diff<-abs(R2_1-R2_2)
alphatails<-alpha/tails
df<-n1+n2-pred-pred-1
delta<-diff/SER2
tabled<-stats::qt(1-alphatails, df)
Power<-round(1-stats::pt(tabled,df,delta),3)
LL_diff<-round(diff - (tabled*SER2),3)
UL_diff<-round(diff + (tabled*SER2),3)


message("Power = ", Power, ", n1 = ", n1,", n2 = ", n2, ", LLdiff = ", LL_diff, ", ULdiff = ", UL_diff)
result <- data.frame(matrix(ncol = 5))
colnames(result) <- c("n1", "n2","LLdiff","ULdiff","Power")
result[, 1]<-n1
result[, 2]<-n2
result[, 3]<-LL_diff
result[, 4]<-UL_diff
result[, 5]<-Power
output<-na.omit(result)
rownames(output)<- c()}

if (pred=="3")
  {
    pop1 <- MASS::mvrnorm(n1, mu<-c(0, 0, 0, 0), Sigma<-matrix(c(1, ry1_1, ry2_1, ry3_1,
                                                         ry1_1, 1, r12_1, r13_1,
                                                         ry2_1, r12_1,1, r23_1,
                                                         ry3_1, r13_1, r23_1, 1),
                                                       ncol=4), empirical=TRUE)

    pop2 <- MASS::mvrnorm(n2, mu<-c(0, 0, 0, 0), Sigma<-matrix(c(1, ry1_2, ry2_2, ry3_2,
                                                          ry1_2, 1, r12_2, r13_2,
                                                          ry2_2, r12_2,1, r23_2,
                                                          ry3_2, r13_2, r23_2, 1),
                                                        ncol=4), empirical=TRUE)
pop1<-data.frame(pop1)
pop2<-data.frame(pop2)
values1<-stats::lm(X1~X2+X3+X4, pop1)
values1<-summary(values1)
R2_1<-(values1$r.squared)
values2<-stats::lm(X1~X2+X3+X4, pop2)
values2<-summary(values2)
R2_2<-(values2$r.squared)
SER2_1<-((4*R2_1)*(1-R2_1)^2)*((n1-pred-1)^2) / ((n1^2 - 1)* (n1+3))
SER2_2<-((4*R2_2)*(1-R2_2)^2)*((n2-pred-1)^2) / ((n2^2 - 1)* (n2+3))
SER2<-(SER2_1 + SER2_2)^.5
diff<-abs(R2_1-R2_2)
alphatails<-alpha/tails
df<-n1+n2-pred-pred-1
delta<-diff/SER2
tabled<-stats::qt(1-alphatails, df)
Power<-round(1-stats::pt(tabled,df,delta),3)
LL_diff<-round(diff - (tabled*SER2),3)
UL_diff<-round(diff + (tabled*SER2),3)
message("Power = ", Power, ", n1 = ", n1,", n2 = ", n2, ", LLdiff = ", LL_diff, ", ULdiff = ", UL_diff)
result <- data.frame(matrix(ncol = 5))
colnames(result) <- c("n1", "n2","LLdiff","ULdiff","Power")
result[, 1]<-n1
result[, 2]<-n2
result[, 3]<-LL_diff
result[, 4]<-UL_diff
result[, 5]<-Power
output<-na.omit(result)
rownames(output)<- c()}
invisible(output)
}

