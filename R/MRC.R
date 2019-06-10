#'Compute power for Multiple Regression with up to Five Predictors
#'Example code below for three predictors. Expand as needed for four or five
#'@param ry1 Correlation between DV (y) and first predictor (1)
#'@param ry2 Correlation between DV (y) and second predictor (2)
#'@param ry3 Correlation between DV (y) and third predictor (3)
#'@param ry4 Correlation between DV (y) and fourth predictor (4)
#'@param ry5 Correlation between DV (y) and fifth predictor (5)
#'@param r12 Correlation between first (1) and second predictor (2)
#'@param r13 Correlation between first (1) and third predictor (3)
#'@param r14 Correlation between first (1) and fourth predictor (4)
#'@param r15 Correlation between first (1) and fifth predictor (5)
#'@param r23 Correlation between second (2) and third predictor (3)
#'@param r24 Correlation between second (2) and fourth predictor (4)
#'@param r25 Correlation between second (2) and fifth predictor (5)
#'@param r34 Correlation between third (3) and fourth predictor (4)
#'@param r35 Correlation between third (3) and fifth predictor (5)
#'@param r45 Correlation between fourth (4) and fifth predictor (5)
#'@param n Sample size
#'@param alpha Type I error (default is .05)
#'@examples
#'MRC(ry1=.40,ry2=.40, r12=-.15,n=30)
#'MRC(ry1=.40,ry2=.40,ry3=-.40, r12=-.15, r13=-.60,r23=.25,n=24)
#'
#'@return Power for Multiple Regression with Two to Five Predictors
#'@export
#'
#'

MRC<-function(ry1=NULL, ry2=NULL, ry3=NULL, ry4=NULL, ry5=NULL,
                  r12=NULL, r13=NULL,r14=NULL,r15=NULL,
                  r23=NULL, r24=NULL, r25=NULL,
                  r34=NULL, r35=NULL,
                  r45=NULL,
                  n=NULL, alpha=.05)
  {

pred<-NA
pred[!is.null(ry2)]<-2
pred[!is.null(ry3)]<-3
pred[!is.null(ry4)]<-4
pred[!is.null(ry5)]<-5
vary<-NA
vary<-1;var1<-1;var2<-1; var3<-1;var4<-1;var5<-1

  if (pred=="2")
  {pop <- MASS::mvrnorm(n, mu = c(0, 0, 0), Sigma = matrix(c(vary, ry1, ry2,
                                                          ry1, var1, r12,
                                                          ry2, r12, var2),
                                                          ncol = 3), empirical = TRUE)
    pop2 = data.frame(pop)

    values<-stats::lm(X1~X2+X3, pop2)
    values<-summary(values)

    int<-(values$coefficients)[1,3]
    tb1<-(values$coefficients)[2,3] #grabs t from each analysis
    tb2<-(values$coefficients)[3,3]
    R2<-values$r.squared
    F<-values$fstatistic[1]
    df1<-values$fstatistic[2]
    df2<-values$fstatistic[3]

    f2<-R2/(1-R2)
    lambdaR2<-f2*df2
    minusalpha<-1-alpha
    FtR2<-stats::qf(minusalpha, df1, df2)
    powerR2<-round(1-stats::pf(FtR2, df1,df2,lambdaR2),3)

    lambdab1<-tb1^2
    lambdab2<-tb2^2
    Fb<-stats::qf(minusalpha, 1, df2)
    powerb1<-round(1-stats::pf(Fb, 1,df2,lambdab1),3)
    powerb2<-round(1-stats::pf(Fb, 1,df2,lambdab2),3)

    message("Sample size is ",n)
    message("Power R2 = ", powerR2)
    message("Power b1 = ", powerb1)
    message("Power b2 = ", powerb2)
    result <- data.frame(matrix(ncol = 4))
    colnames(result) <- c( "n","Power R2", "Power b1", "Power b2")
    result[, 1]<-n
    result[, 2]<-powerR2
    result[, 3]<-powerb1
    result[, 4]<-powerb2
    output<-na.omit(result)
    rownames(output)<- c()
     }

  if (pred=="3")
        {
  pop <- MASS::mvrnorm(n, mu = c(0, 0, 0, 0), Sigma = matrix(c(vary, ry1, ry2, ry3,
                                                             ry1, var1, r12, r13,
                                                             ry2, r12, var2, r23,
                                                             ry3, r13, r23, var3),
                 ncol = 4), empirical = TRUE)
  pop2 = data.frame(pop)

  values<-stats::lm(X1~X2+X3+X4, pop2)
  values<-summary(values)

  int<-(values$coefficients)[1,3]
  tb1<-(values$coefficients)[2,3] #grabs t from each analysis
  tb2<-(values$coefficients)[3,3]
  tb3<-(values$coefficients)[4,3]
  R2<-values$r.squared
  F<-values$fstatistic[1]
  df1<-values$fstatistic[2]
  df2<-values$fstatistic[3]

  f2<-R2/(1-R2)
  lambdaR2<-f2*df2
  minusalpha<-1-alpha
  FtR2<-stats::qf(minusalpha, df1, df2)
  powerR2<-round(1-stats::pf(FtR2, df1,df2,lambdaR2),3)

  lambdab1<-tb1^2
  lambdab2<-tb2^2
  lambdab3<-tb3^2
  Fb<-stats::qf(minusalpha, 1, df2)
  powerb1<-round(1-stats::pf(Fb, 1,df2,lambdab1),3)
  powerb2<-round(1-stats::pf(Fb, 1,df2,lambdab2),3)
  powerb3<-round(1-stats::pf(Fb, 1,df2,lambdab3),3)


  message("Sample size is ",n)
  message("Power R2 = ", powerR2)
  message("Power b1 = ", powerb1)
  message("Power b2 = ", powerb2)
  message("Power b3 = ", powerb3)
  result <- data.frame(matrix(ncol = 5))
  colnames(result) <- c( "n","Power R2", "Power b1", "Power b2",
                         "Power b3")
  result[, 1]<-n
  result[, 2]<-powerR2
  result[, 3]<-powerb1
  result[, 4]<-powerb2
  result[, 5]<-powerb3
  output<-na.omit(result)
  rownames(output)<- c()
  }

  if (pred=="4")
  {
    pop <- MASS::mvrnorm(n, mu = c(0, 0, 0, 0,0), Sigma = matrix(c(vary, ry1, ry2, ry3, ry4,
                                                                  ry1, var1, r12, r13, r14,
                                                                  ry2, r12, var2, r23, r24,
                                                                  ry3, r13, r23, var3, r34,
                                                                  ry4, r14, r24, r34, var4),
                                                             ncol = 5), empirical = TRUE)
    pop2 = data.frame(pop)

    values<-stats::lm(X1~X2+X3+X4+X5, pop2)
    values<-summary(values)

    int<-(values$coefficients)[1,3]
    tb1<-(values$coefficients)[2,3] #grabs t from each analysis
    tb2<-(values$coefficients)[3,3]
    tb3<-(values$coefficients)[4,3]
    tb4<-(values$coefficients)[5,3]
    R2<-values$r.squared
    F<-values$fstatistic[1]
    df1<-values$fstatistic[2]
    df2<-values$fstatistic[3]

    f2<-R2/(1-R2)
    lambdaR2<-f2*df2
    minusalpha<-1-alpha
    FtR2<-stats::qf(minusalpha, df1, df2)
    powerR2<-round(1-stats::pf(FtR2, df1,df2,lambdaR2),3)

    lambdab1<-tb1^2
    lambdab2<-tb2^2
    lambdab3<-tb3^2
    lambdab4<-tb4^2
    Fb<-stats::qf(minusalpha, 1, df2)
    powerb1<-round(1-stats::pf(Fb, 1,df2,lambdab1),3)
    powerb2<-round(1-stats::pf(Fb, 1,df2,lambdab2),3)
    powerb3<-round(1-stats::pf(Fb, 1,df2,lambdab3),3)
    powerb4<-round(1-stats::pf(Fb, 1,df2,lambdab4),3)

    message("Sample size is ",n)
    message("Power R2 = ", powerR2)
    message("Power b1 = ", powerb1)
    message("Power b2 = ", powerb2)
    message("Power b3 = ", powerb3)
    message("Power b4 = ", powerb4)
    result <- data.frame(matrix(ncol = 6))
    colnames(result) <- c( "n","Power R2", "Power b1", "Power b2",
                           "Power b3", "Power b4")
    result[, 1]<-n
    result[, 2]<-powerR2
    result[, 3]<-powerb1
    result[, 4]<-powerb2
    result[, 5]<-powerb3
    result[, 6]<-powerb4
    output<-na.omit(result)
    rownames(output)<- c()
    }

    if (pred=="5")
    {

      pop <- MASS::mvrnorm(n, mu = c(0, 0, 0, 0,0,0), Sigma = matrix(c(vary, ry1, ry2, ry3, ry4, ry5,
                                                                    ry1, var1, r12, r13, r14,r15,
                                                                    ry2, r12, var2, r23, r24,r25,
                                                                    ry3, r13, r23, var3, r34,r35,
                                                                    ry4, r14, r24, r34, var4,r45,
                                                                    ry5,r15,r25,r34,r45,var5),
                                                                  ncol = 6), empirical = TRUE)
      pop2 = data.frame(pop)

      values<-stats::lm(X1~X2+X3+X4+X5+X6, pop2)
      values<-summary(values)

      int<-(values$coefficients)[1,3]
      tb1<-(values$coefficients)[2,3] #grabs t from each analysis
      tb2<-(values$coefficients)[3,3]
      tb3<-(values$coefficients)[4,3]
      tb4<-(values$coefficients)[5,3]
      tb5<-(values$coefficients)[6,3]
      R2<-values$r.squared
      F<-values$fstatistic[1]
      df1<-values$fstatistic[2]
      df2<-values$fstatistic[3]

      f2<-R2/(1-R2)
      lambdaR2<-f2*df2
      minusalpha<-1-alpha
      FtR2<-stats::qf(minusalpha, df1, df2)
      powerR2<-round(1-stats::pf(FtR2, df1,df2,lambdaR2),3)

      lambdab1<-tb1^2
      lambdab2<-tb2^2
      lambdab3<-tb3^2
      lambdab4<-tb4^2
      lambdab5<-tb5^2
      Fb<-stats::qf(minusalpha, 1, df2)
      powerb1<-round(1-stats::pf(Fb, 1,df2,lambdab1),3)
      powerb2<-round(1-stats::pf(Fb, 1,df2,lambdab2),3)
      powerb3<-round(1-stats::pf(Fb, 1,df2,lambdab3),3)
      powerb4<-round(1-stats::pf(Fb, 1,df2,lambdab4),3)
      powerb5<-round(1-stats::pf(Fb, 1,df2,lambdab5),3)

      message("Sample size is ",n)
      message("Power R2 = ", powerR2)
      message("Power b1 = ", powerb1)
      message("Power b2 = ", powerb2)
      message("Power b3 = ", powerb3)
      message("Power b4 = ", powerb4)
      message("Power b5 = ", powerb5)
      result <- data.frame(matrix(ncol = 7))
      colnames(result) <- c( "n","Power R2", "Power b1", "Power b2",
                             "Power b3", "Power b4", "Power b5")
      result[, 1]<-n
      result[, 2]<-powerR2
      result[, 3]<-powerb1
      result[, 4]<-powerb2
      result[, 5]<-powerb3
      result[, 6]<-powerb4
      result[, 7]<-powerb5
      output<-na.omit(result)
      rownames(output)<- c()
      }
      invisible(output)
}
