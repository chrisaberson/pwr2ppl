#'Compute power for Multiple Regression with Up to Five Predictors
#'Requires correlations between all variables as sample size. Means, sds, and alpha are option. Also computes Power(All)
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
#'@examples
#'\donttest{MRC_all(ry1=.50,ry2=.50,ry3=.50, r12=.2, r13=.3,r23=.4,n=82, rep=10000)}
#'@return Power for Multiple Regression (ALL)
#'@export
#'
#'


MRC_all<-function(ry1=NULL, ry2=NULL, ry3=NULL, ry4=NULL, ry5=NULL,
                  r12=NULL, r13=NULL,r14=NULL,r15=NULL,
                  r23=NULL, r24=NULL, r25=NULL,
                  r34=NULL, r35=NULL,
                  r45=NULL, n=NULL, alpha=.05, rep = 10000)
  {

pred<-NA
pred[!is.null(ry2)]<-2
pred[!is.null(ry3)]<-3
pred[!is.null(ry4)]<-4
pred[!is.null(ry5)]<-5

vary<-NA
vary<-1;var1<-1;var2<-1; var3<-1; var4<-1; var5<-1

samp = data.frame(MASS::mvrnorm(n, mu = c(0, 0, 0),
                                Sigma = matrix(c(vary, ry1, ry2,
                                                 ry1, var1, r12,
                                                 ry2, r12, var2),
                                               ncol = 3), empirical = FALSE))

  if (pred=="2")

    {
    nruns = rep
    int = numeric(nruns)
    b1 = numeric(nruns)
    b2 = numeric(nruns)
    R2 = numeric(nruns)
    F = numeric(nruns)
    df1 = numeric(nruns)
    df2 = numeric(nruns)
    for (i in 1:nruns)
    {samp <- data.frame(MASS::mvrnorm(n, mu = c(0, 0, 0), Sigma = matrix(c(vary, ry1, ry2,
                                             ry1, var1, r12,
                                             ry2, r12, var2),
                                             ncol = 3), empirical = FALSE))
    test <- stats::lm(formula = X1 ~ X2+ X3, data = samp)
    c<-summary(test)
    int[i] = stats::coef(summary(test))[1,4]
    b1[i] = stats::coef(summary(test))[2,4] #grabs p from each analysis
    b2[i] = stats::coef(summary(test))[3,4]
    R2[i] = c$r.squared
    F[i]<-c$fstatistic[1]
    df1[i]<-c$fstatistic[2]
    df2[i]<-c$fstatistic[3]}
    Powerall = data.frame(int = int, b1 = b1, b2 = b2)
    Powerall[4:5, "rejectb1"]<-NA
    Powerall$rejectb1 [ b1 < alpha] <- 1
    Powerall$rejectb1 [ b1 >= alpha] <- 0
    Powerall[4:5, "rejectb2"]<-NA
    Powerall$rejectb2 [ b2 < alpha] <- 1
    Powerall$rejectb2 [ b2 >= alpha] <- 0
    Powerall[4:5, "rejecttotal"]<-NA
    Powerall$rejectall <- (Powerall$rejectb1+ Powerall$rejectb2)

    Reject.None <-NA
    Reject.None [Powerall$rejectall == 0]<-1
    Reject.None [Powerall$rejectall > 0]<-0
    Reject.One <-NA
    Reject.One [Powerall$rejectall == 1]<-1
    Reject.One [Powerall$rejectall != 1]<-0
    Reject.All <-NA
    Reject.All [Powerall$rejectall == 2]<-1
    Reject.All [Powerall$rejectall != 2]<-0
    is.numeric(Reject.None)
    is.numeric(Reject.One)
    is.numeric(Reject.All)

    Power_b1<-mean(Powerall$rejectb1)
    Power_b2<-mean(Powerall$rejectb2)
    pR2<-1-stats::pf(F,df1, df2)
    Powerall$rejectR2 [pR2 < alpha] <- 1
    Powerall$rejectR2 [pR2 >= alpha] <- 0
    Power_R2<-mean(Powerall$rejectR2)
    PowerAll_R0<-mean(Reject.None)
    PowerAll_R1<-mean(Reject.One)
    PowerAll_R2<-mean(Reject.All)


   message("Sample size is ",n)
   message("Power R2 = ", Power_R2)
   message("Power b1 = ", Power_b1)
   message("Power b2 = ", Power_b2)
   message("Proportion Rejecting None = ", PowerAll_R0)
   message("Proportion Rejecting One = ", PowerAll_R1)
   message("Power ALL (Proportion Rejecting All) = ", PowerAll_R2)
   result <- data.frame(matrix(ncol = 7))
   colnames(result) <- c( "n","Power R2", "Power b1", "Power b2",
                          "Power Reject None",
                          "Power Reject One","Power Reject All")
   result[, 1]<-n
   result[, 2]<-Power_R2
   result[, 3]<-Power_b1
   result[, 4]<-Power_b2
   result[, 5]<-PowerAll_R0
   result[, 6]<-PowerAll_R1
   result[, 7]<-PowerAll_R2


   output<-na.omit(result)
   rownames(output)<- c()
     }

  if (pred=="3")
    {
  nruns = rep
  int = numeric(nruns)
  b1 = numeric(nruns)
  b2 = numeric(nruns)
  b3 = numeric(nruns)
  R2 = numeric(nruns)
  F = numeric(nruns)
  df1 = numeric(nruns)
  df2 = numeric(nruns)
  for (i in 1:nruns)
  {samp <- data.frame(MASS::mvrnorm(n, mu = c(0, 0, 0, 0),
                                        Sigma = matrix(c(vary, ry1, ry2, ry3,
                                                         ry1, var1, r12, r13,
                                                         ry2, r12, var2, r23,
                                                         ry3, r13, r23, var3),
                                                       ncol = 4), empirical = FALSE))

  test <- stats::lm(formula = X1 ~ X2+ X3+ X4, data = samp)
  c<-summary(test)
  int[i] = stats::coef(summary(test))[1,4]
  b1[i] = stats::coef(summary(test))[2,4] #grabs p from each analysis
  b2[i] = stats::coef(summary(test))[3,4]
  b3[i] = stats::coef(summary(test))[4,4]
  R2[i] = c$r.squared
  F[i]<-c$fstatistic[1]
  df1[i]<-c$fstatistic[2]
  df2[i]<-c$fstatistic[3]}
  Powerall = data.frame(int = int, b1 = b1, b2 = b2, b3 = b3)
  Powerall[4:5, "rejectb1"]<-NA
  Powerall$rejectb1 [ b1 < alpha] <- 1
  Powerall$rejectb1 [ b1 >= alpha] <- 0
  Powerall[4:5, "rejectb2"]<-NA
  Powerall$rejectb2 [ b2 < alpha] <- 1
  Powerall$rejectb2 [ b2 >= alpha] <- 0
  Powerall[4:5, "rejectb3"]<-NA
  Powerall$rejectb3 [ b3 < alpha] <- 1
  Powerall$rejectb3 [ b3 >= alpha] <- 0
  Powerall[4:5, "rejecttotal"]<-NA
  Powerall$rejectall <- (Powerall$rejectb1+ Powerall$rejectb2+ Powerall$rejectb3)

  Reject.None <-NA
  Reject.None [Powerall$rejectall == 0]<-1
  Reject.None [Powerall$rejectall > 0]<-0
  Reject.One <-NA
  Reject.One [Powerall$rejectall == 1]<-1
  Reject.One [Powerall$rejectall != 1]<-0
  Reject.Two <-NA
  Reject.Two [Powerall$rejectall == 2]<-1
  Reject.Two [Powerall$rejectall != 2]<-0
  Reject.All <-NA
  Reject.All [Powerall$rejectall == 3]<-1
  Reject.All [Powerall$rejectall != 3]<-0
  is.numeric(Reject.None)
  is.numeric(Reject.One)
  is.numeric(Reject.Two)
  is.numeric(Reject.All)

  Power_b1<-mean(Powerall$rejectb1)
  Power_b2<-mean(Powerall$rejectb2)
  Power_b3<-mean (Powerall$rejectb3)
  pR2<-1-stats::pf(F,df1, df2)
  Powerall$rejectR2 [pR2 < alpha] <- 1
  Powerall$rejectR2 [pR2 >= alpha] <- 0
  Power_R2<-mean(Powerall$rejectR2)
  PowerAll_R0<-mean(Reject.None)
  PowerAll_R1<-mean(Reject.One)
  PowerAll_R2<-mean(Reject.Two)
  PowerAll_R3<-mean(Reject.All)


 message("Sample size is ",n)
 message("Power R2 = ", Power_R2)
 message("Power b1 = ", Power_b1)
 message("Power b2 = ", Power_b2)
 message("Power b3 = ", Power_b3)
 message("Proportion Rejecting None = ", PowerAll_R0)
 message("Proportion Rejecting One = ", PowerAll_R1)
 message("Proportion Rejecting Two = ", PowerAll_R2)
 message("Power ALL (Proportion Rejecting All) = ", PowerAll_R3)
 result <- data.frame(matrix(ncol = 9))
 colnames(result) <- c( "n","Power R2", "Power b1", "Power b2",
                        "Power b3", "Power Reject None",
                        "Power Reject One","Power Reject Two","Power Reject All")
 result[, 1]<-n
 result[, 2]<-Power_R2
 result[, 3]<-Power_b1
 result[, 4]<-Power_b2
 result[, 5]<-Power_b3
 result[, 6]<-PowerAll_R0
 result[, 7]<-PowerAll_R1
 result[, 8]<-PowerAll_R2
 result[, 9]<-PowerAll_R3

 output<-na.omit(result)
 rownames(output)<- c()
  }

if (pred=="4")
{
  nruns = rep
  int = numeric(nruns)
  b1 = numeric(nruns)
  b2 = numeric(nruns)
  b3 = numeric(nruns)
  b4 = numeric(nruns)
  R2 = numeric(nruns)
  F = numeric(nruns)
  df1 = numeric(nruns)
  df2 = numeric(nruns)
  for (i in 1:nruns)
  { samp <- data.frame(MASS::mvrnorm(n, mu = c(0, 0, 0, 0,0), Sigma = matrix(c(vary, ry1, ry2, ry3, ry4,
                                              ry1, var1, r12, r13, r14,
                                              ry2, r12, var2, r23, r24,
                                              ry3, r13, r23, var3, r34,
                                              ry4, r14, r24, r34, var4),
                                            ncol = 5), empirical = FALSE))
  test <- stats::lm(formula = X1 ~ X2+ X3+ X4 + X5, data = samp)
  c<-summary(test)
  int[i] = stats::coef(summary(test))[1,4]
  b1[i] = stats::coef(summary(test))[2,4] #grabs p from each analysis
  b2[i] = stats::coef(summary(test))[3,4]
  b3[i] = stats::coef(summary(test))[4,4]
  b4[i] = stats::coef(summary(test))[5,4]
  R2[i] = c$r.squared
  F[i]<-c$fstatistic[1]
  df1[i]<-c$fstatistic[2]
  df2[i]<-c$fstatistic[3]}
  Powerall = data.frame(int = int, b1 = b1, b2 = b2, b3 = b3, b4 = b4)
  Powerall[4:5, "rejectb1"]<-NA
  Powerall$rejectb1 [ b1 < alpha] <- 1
  Powerall$rejectb1 [ b1 >= alpha] <- 0
  Powerall[4:5, "rejectb2"]<-NA
  Powerall$rejectb2 [ b2 < alpha] <- 1
  Powerall$rejectb2 [ b2 >= alpha] <- 0
  Powerall[4:5, "rejectb3"]<-NA
  Powerall$rejectb3 [ b3 < alpha] <- 1
  Powerall$rejectb3 [ b3 >= alpha] <- 0
  Powerall[4:5, "rejectb4"]<-NA
  Powerall$rejectb4 [ b4 < alpha] <- 1
  Powerall$rejectb4 [ b4 >= alpha] <- 0
  Powerall[4:5, "rejecttotal"]<-NA
  Powerall$rejectall <- (Powerall$rejectb1+ Powerall$rejectb2+ Powerall$rejectb3+ Powerall$rejectb4)

  Reject.None <-NA
  Reject.None [Powerall$rejectall == 0]<-1
  Reject.None [Powerall$rejectall > 0]<-0
  Reject.One <-NA
  Reject.One [Powerall$rejectall == 1]<-1
  Reject.One [Powerall$rejectall != 1]<-0
  Reject.Two <-NA
  Reject.Two [Powerall$rejectall == 2]<-1
  Reject.Two [Powerall$rejectall != 2]<-0
  Reject.Three <-NA
  Reject.Three [Powerall$rejectall == 3]<-1
  Reject.Three [Powerall$rejectall != 3]<-0
  Reject.All <-NA
  Reject.All [Powerall$rejectall == 4]<-1
  Reject.All [Powerall$rejectall != 4]<-0
  is.numeric(Reject.None)
  is.numeric(Reject.One)
  is.numeric(Reject.Two)
  is.numeric(Reject.Three)
  is.numeric(Reject.All)

  Power_b1<-mean(Powerall$rejectb1)
  Power_b2<-mean(Powerall$rejectb2)
  Power_b3<-mean (Powerall$rejectb3)
  Power_b4<-mean (Powerall$rejectb4)
  pR2<-1-stats::pf(F,df1, df2)
  Powerall$rejectR2 [pR2 < alpha] <- 1
  Powerall$rejectR2 [pR2 >= alpha] <- 0
  Power_R2<-mean(Powerall$rejectR2)
  PowerAll_R0<-mean(Reject.None)
  PowerAll_R1<-mean(Reject.One)
  PowerAll_R2<-mean(Reject.Two)
  PowerAll_R3<-mean(Reject.Three)
  PowerAll_R4<-mean(Reject.All)


 message("Sample size is ",n)
 message("Power R2 = ", Power_R2)
 message("Power b1 = ", Power_b1)
 message("Power b2 = ", Power_b2)
 message("Power b3 = ", Power_b3)
 message("Power b4 = ", Power_b4)
 message("Proportion Rejecting None = ", PowerAll_R0)
 message("Proportion Rejecting One = ", PowerAll_R1)
 message("Proportion Rejecting Two = ", PowerAll_R2)
 message("Proportion Rejecting Three = ", PowerAll_R3)
 message("Power ALL (Proportion Rejecting All) = ", PowerAll_R4)
 result <- data.frame(matrix(ncol = 11))
 colnames(result) <- c( "n","Power R2", "Power b1", "Power b2",
                        "Power b3", "Power b4", "Power Reject None",
                        "Power Reject One","Power Reject Two","Power Reject Three",
                        "Power Reject All")
 result[, 1]<-n
 result[, 2]<-Power_R2
 result[, 3]<-Power_b1
 result[, 4]<-Power_b2
 result[, 5]<-Power_b3
 result[, 6]<-Power_b4
 result[, 7]<-PowerAll_R0
 result[, 8]<-PowerAll_R1
 result[, 9]<-PowerAll_R2
 result[, 10]<-PowerAll_R3
 result[, 11]<-PowerAll_R4
 output<-na.omit(result)
 rownames(output)<- c()
}

  if (pred=="5")
  {
    nruns = rep
    int = numeric(nruns)
    b1 = numeric(nruns)
    b2 = numeric(nruns)
    b3 = numeric(nruns)
    b4 = numeric(nruns)
    b5 = numeric(nruns)
    R2 = numeric(nruns)
    F = numeric(nruns)
    df1 = numeric(nruns)
    df2 = numeric(nruns)
    for (i in 1:nruns)
    {samp <- data.frame(MASS::mvrnorm(n, mu = c(0, 0, 0, 0,0,0), Sigma = matrix(c(vary, ry1, ry2, ry3, ry4, ry5,
                                             ry1, var1, r12, r13, r14,r15,
                                             ry2, r12, var2, r23, r24,r25,
                                             ry3, r13, r23, var3, r34,r35,
                                             ry4, r14, r24, r34, var4,r45,
                                             ry5,r15,r25,r35,r45,var5),
                                             ncol = 6), empirical = FALSE))

    test <- stats::lm(formula = X1 ~ X2+ X3+ X4 + X5+ X6, data = samp)
    c<-summary(test)
    int[i] = stats::coef(summary(test))[1,4]
    b1[i] = stats::coef(summary(test))[2,4] #grabs p from each analysis
    b2[i] = stats::coef(summary(test))[3,4]
    b3[i] = stats::coef(summary(test))[4,4]
    b4[i] = stats::coef(summary(test))[5,4]
    b5[i] = stats::coef(summary(test))[6,4]
    R2[i] = c$r.squared
    F[i]<-c$fstatistic[1]
    df1[i]<-c$fstatistic[2]
    df2[i]<-c$fstatistic[3]}
    Powerall = data.frame(int = int, b1 = b1, b2 = b2, b3 = b3, b4 = b4, b5 = b5)
    Powerall[4:5, "rejectb1"]<-NA
    Powerall$rejectb1 [ b1 < alpha] <- 1
    Powerall$rejectb1 [ b1 >= alpha] <- 0
    Powerall[4:5, "rejectb2"]<-NA
    Powerall$rejectb2 [ b2 < alpha] <- 1
    Powerall$rejectb2 [ b2 >= alpha] <- 0
    Powerall[4:5, "rejectb3"]<-NA
    Powerall$rejectb3 [ b3 < alpha] <- 1
    Powerall$rejectb3 [ b3 >= alpha] <- 0
    Powerall[4:5, "rejectb4"]<-NA
    Powerall$rejectb4 [ b4 < alpha] <- 1
    Powerall$rejectb4 [ b4 >= alpha] <- 0
    Powerall[4:5, "rejectb5"]<-NA
    Powerall$rejectb5 [ b5 < alpha] <- 1
    Powerall$rejectb5 [ b5 >= alpha] <- 0
    Powerall[4:5, "rejecttotal"]<-NA
    Powerall$rejectall <- (Powerall$rejectb1+ Powerall$rejectb2+ Powerall$rejectb3+ Powerall$rejectb4+ Powerall$rejectb5)

    Reject.None <-NA
    Reject.None [Powerall$rejectall == 0]<-1
    Reject.None [Powerall$rejectall > 0]<-0
    Reject.One <-NA
    Reject.One [Powerall$rejectall == 1]<-1
    Reject.One [Powerall$rejectall != 1]<-0
    Reject.Two <-NA
    Reject.Two [Powerall$rejectall == 2]<-1
    Reject.Two [Powerall$rejectall != 2]<-0
    Reject.Three <-NA
    Reject.Three [Powerall$rejectall == 3]<-1
    Reject.Three [Powerall$rejectall != 3]<-0
    Reject.Four <-NA
    Reject.Four [Powerall$rejectall == 4]<-1
    Reject.Four [Powerall$rejectall != 4]<-0
    Reject.All <-NA
    Reject.All [Powerall$rejectall == 5]<-1
    Reject.All [Powerall$rejectall != 5]<-0
    is.numeric(Reject.None)
    is.numeric(Reject.One)
    is.numeric(Reject.Two)
    is.numeric(Reject.Three)
    is.numeric(Reject.Four)
    is.numeric(Reject.All)

    Power_b1<-mean(Powerall$rejectb1)
    Power_b2<-mean(Powerall$rejectb2)
    Power_b3<-mean (Powerall$rejectb3)
    Power_b4<-mean (Powerall$rejectb4)
    Power_b5<-mean (Powerall$rejectb5)
    pR2<-1-stats::pf(F,df1, df2)
    Powerall$rejectR2 [pR2 < alpha] <- 1
    Powerall$rejectR2 [pR2 >= alpha] <- 0
    Power_R2<-mean(Powerall$rejectR2)
    PowerAll_R0<-mean(Reject.None)
    PowerAll_R1<-mean(Reject.One)
    PowerAll_R2<-mean(Reject.Two)
    PowerAll_R3<-mean(Reject.Three)
    PowerAll_R4<-mean(Reject.Four)
    PowerAll_R5<-mean(Reject.All)


   message("Sample size is ",n)
   message("Power R2 = ", Power_R2)
   message("Power b1 = ", Power_b1)
   message("Power b2 = ", Power_b2)
   message("Power b3 = ", Power_b3)
   message("Power b4 = ", Power_b4)
   message("Power b5 = ", Power_b5)
   message("Proportion Rejecting None = ", PowerAll_R0)
   message("Proportion Rejecting One = ", PowerAll_R1)
   message("Proportion Rejecting Two = ", PowerAll_R2)
   message("Proportion Rejecting Three = ", PowerAll_R3)
   message("Proportion Rejecting Four = ", PowerAll_R4)
   message("Power ALL (Proportion Rejecting All) = ", PowerAll_R5)
   result <- data.frame(matrix(ncol = 13))
   colnames(result) <- c( "n","Power R2", "Power b1", "Power b2",
                          "Power b3", "Power b4", "Power b5", "Power Reject None",
                          "Power Reject One","Power Reject Two","Power Reject Three",
                          "Power Reject Four","Power Reject All")
   result[, 1]<-n
   result[, 2]<-Power_R2
   result[, 3]<-Power_b1
   result[, 4]<-Power_b2
   result[, 5]<-Power_b3
   result[, 6]<-Power_b4
   result[, 7]<-Power_b5
   result[, 8]<-PowerAll_R0
   result[, 9]<-PowerAll_R1
   result[, 10]<-PowerAll_R2
   result[, 11]<-PowerAll_R3
   result[, 12]<-PowerAll_R4
   result[, 13]<-PowerAll_R5
   output<-na.omit(result)
   rownames(output)<- c()

   }
   invisible(output)
    }





