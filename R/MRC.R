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
#'@param rep number of replications (default is 10000)
#'@param my Mean of DV (default is 0)
#'@param m1 Mean of first predictor (default is 0)
#'@param m2 Mean of second redictor (default is 0)
#'@param m3 Mean of third predictor (default is 0)
#'@param m4 Mean of fourth predictor (default is 0)
#'@param m5 Mean of fifth predictor (default is 0)
#'@param sy Standard deviation of DV (default is 1)
#'@param s1 Standard deviation of first predictor (default is 1)
#'@param s2 Standard deviation of second predictor (default is 1)
#'@param s3 Standard deviation of third predictor (default is 1)
#'@param s4 Standard deviation of fourth predictor (default is 1)
#'@param s5 Standard deviation of fifth predictor (default is 1)
#'@return Power for Multiple Regression with Two to Five Predictors
#'@export
#'
#'

MRC<-function(ry1=NULL, ry2=NULL, ry3=NULL, ry4=NULL, ry5=NULL,
                  r12=NULL, r13=NULL,r14=NULL,r15=NULL,
                  r23=NULL, r24=NULL, r25=NULL,
                  r34=NULL, r35=NULL,
                  r45=NULL,
                  n=NULL, alpha=.05, rep = 10000, my=0,m1=0,m2=0,m3=0, m4=0, m5=0,
                  sy=1,s1=1,s2=1,s3=1, s4=1, s5=1)
  {

pred<-NA
pred[is.null(r23)]<-2
pred[!is.null(r23)]<-3
pred[!is.null(ry4)]<-4
vary<-NA
vary<-sy^2;var1<-s1^2;var2<-s2^2; var3<-s3^2;var4<-s4^2;var5<-s5^2

  if (pred=="2")
  {pop <- MASS::mvrnorm(n, mu = c(my, m1, m2), Sigma = matrix(c(vary, ry1, ry2,
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

    {print(paste("Sample size is ",n))}
    {print(paste("Power R2 = ", powerR2))}
    {print(paste("Power b1 = ", powerb1))}
    {print(paste("Power b2 = ", powerb2))}
     }

  if (pred=="3")
        {
  pop <- MASS::mvrnorm(n, mu = c(my, m1, m2, m3), Sigma = matrix(c(vary, ry1, ry2, ry3,
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


  {print(paste("Sample size is ",n))}
  {print(paste("Power R2 = ", powerR2))}
  {print(paste("Power b1 = ", powerb1))}
  {print(paste("Power b2 = ", powerb2))}
  {print(paste("Power b3 = ", powerb3))}
  }

  if (pred=="4")
  {
    pop <- MASS::mvrnorm(n, mu = c(my, m1, m2, m3,m4), Sigma = matrix(c(vary, ry1, ry2, ry3, ry4,
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

    print(paste("Sample size is ",n))
    print(paste("Power R2 = ", powerR2))
    print(paste("Power b1 = ", powerb1))
    print(paste("Power b2 = ", powerb2))
    print(paste("Power b3 = ", powerb3))
    print(paste("Power b4 = ", powerb4))
    }

    if (pred=="5")
    {

      pop <- MASS::mvrnorm(n, mu = c(my, m1, m2, m3,m4,m5), Sigma = matrix(c(vary, ry1, ry2, ry3, ry4, ry5,
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

      print(paste("Sample size is ",n))
      print(paste("Power R2 = ", powerR2))
      print(paste("Power b1 = ", powerb1))
      print(paste("Power b2 = ", powerb2))
      print(paste("Power b3 = ", powerb3))
      print(paste("Power b4 = ", powerb4))
      print(paste("Power b5 = ", powerb5))
    }
      }
