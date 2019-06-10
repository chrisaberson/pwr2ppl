#'Power for Comparing Independent Coefficients in Multiple Regression with Two or Three Predictors
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
#'@examples
#'indb(ry1_1=.40, ry2_1=.40, ry3_1 =-.40, r12_1=-.15,r13_1=-.60, r23_1=.25,
#'ry1_2=.40, ry2_2=.10, ry3_2 =-.40, r12_2=-.15,r13_2=-.60, r23_2=.25,
#'n1=50,n2=50, alpha=.05)
#'@return Power for Comparing Independent Coefficients in Multiple Regression
#'@export
#'
#'
indb<-function(ry1_1, ry2_1, ry3_1=NULL, r12_1, r13_1=NULL, r23_1=NULL,n1,
                   ry1_2, ry2_2, ry3_2=NULL, r12_2, r13_2=NULL, r23_2=NULL,n2, alpha=.05)
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
  b1_1<-(values1$coefficients)[2,1] #grabs b from each analysis
  b2_1<-(values1$coefficients)[3,1]
  seb1_1<-(values1$coefficients)[2,2]
  seb2_1<-(values1$coefficients)[3,2]
  sebb1<-((seb1_1^2)+(seb2_1^2))^.5
  values2<-stats::lm(X1~X2+X3, pop2)
  values2<-summary(values2)
  b1_2<-(values2$coefficients)[2,1] #grabs b from each analysis
  b2_2<-(values2$coefficients)[3,1]
  seb1_2<-(values2$coefficients)[2,2]
  seb2_2<-(values2$coefficients)[3,2]
  sebb2<-((seb2_1^2)+(seb2_2^2))^.5
  t1<-(abs(b1_1-b1_2))/sebb1
  t2<-(abs(b2_1-b2_2))/sebb2
  lambda1<-t1^2
  lambda2<-t2^2
  df<-n1+n2-pred-pred-2
  minusalpha<-1-alpha
  Fb<-stats::qf(minusalpha, 1, df)
  power1<-round(1-stats::pf(Fb, 1,df,lambda1),3)
  power2<-round(1-stats::pf(Fb, 1,df,lambda2),3)
  message("Sample size Group 1 = ",n1, ", Group 2 = ", n2)
  message("Power comparing b1 across samples = ", power1)
  message("Power comparing b2 across samples = ", power2)
  result <- data.frame(matrix(ncol = 4))
  colnames(result) <- c("n1","n2","Power b1","Power b2")
  result[, 1]<-n1
  result[, 2]<-n2
  result[, 3]<-power1
  result[, 4]<-power2
  output<-na.omit(result)
  rownames(output)<- c()

}

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
    values1<-stats::lm(X1~X2+X3+X4, pop1)
    values1<-summary(values1)
    b1_1<-(values1$coefficients)[2,1] #grabs b from each analysis
    b2_1<-(values1$coefficients)[3,1]
    b3_1<-(values1$coefficients)[4,1]
    seb1_1<-(values1$coefficients)[2,2]
    seb2_1<-(values1$coefficients)[3,2]
    seb3_1<-(values1$coefficients)[4,2]
    pop2<-data.frame(pop2)
    values2<-stats::lm(X1~X2+X3+X4, pop2)
    values2<-summary(values2)
    b1_2<-(values2$coefficients)[2,1] #grabs b from each analysis
    b2_2<-(values2$coefficients)[3,1]
    b3_2<-(values2$coefficients)[4,1]
    seb1_2<-(values2$coefficients)[2,2]
    seb2_2<-(values2$coefficients)[3,2]
    seb3_2<-(values2$coefficients)[4,2]
    sebb1<-((seb1_1^2)+(seb1_2^2))^.5
    sebb2<-((seb2_1^2)+(seb2_2^2))^.5
    sebb3<-((seb3_1^2)+(seb3_2^2))^.5
    t1<-(abs(b1_1-b1_2))/sebb1
    t2<-(abs(b2_1-b2_2))/sebb2
    t3<-(abs(b3_1-b3_2))/sebb3
    lambda1<-t1^2
    lambda2<-t2^2
    lambda3<-t3^2
    df<-n1+n2-pred-pred-2
    minusalpha<-1-alpha
    Fb<-stats::qf(minusalpha, 1, df)
    power1<-round(1-stats::pf(Fb, 1,df,lambda1),3)
    power2<-round(1-stats::pf(Fb, 1,df,lambda2),3)
    power3<-round(1-stats::pf(Fb, 1,df,lambda3),3)
    message("Sample size Group 1 = ",n1, ", Group 2 = ", n2)
    message("Power comparing b1 across samples = ", power1)
    message("Power comparing b2 across samples = ", power2)
    message("Power comparing b3 across samples = ", power3)
    result <- data.frame(matrix(ncol = 5))
    colnames(result) <- c("n1","n2","Power b1","Power b2","Power b3")
    result[, 1]<-n1
    result[, 2]<-n2
    result[, 3]<-power1
    result[, 4]<-power2
    result[, 5]<-power3
    output<-na.omit(result)
    rownames(output)<- c()

  }
  invisible(output)
  }
