#'Compute power for Mutliple Regression with Three Predictors
#'Requires correlatiosn between all variables as sample size. Means, sds, and alpha are option. Also computes Power(All)
#'@param ry1 Correlation between DV (y) and first predictor (1)
#'@param ry2 Correlation between DV (y) and second predictor (2)
#'@param ry3 Correlation between DV (y) and third predictor (3)
#'@param r12 Correlation between first (1) and second predictor (2)
#'@param r13 Correlation between first (1) and third predictor (3)
#'@param r23 Correlation between second (2) and third predictor (3)
#'@param n Sample size
#'@param alpha Type I error (default is .05)
#'@param rep number of replications (default is 10000)
#'@param my Mean of DV (default is 0)
#'@param m1 Mean of first predictor (default is 0)
#'@param m2 Mean of second redictor (default is 0)
#'@param m3 Mean of third predictor (default is 0)
#'@param sy Standard deviation of DV (default is 1)
#'@param s1 Standard deviation of first predictor (default is 1)
#'@param s2 Standard deviation of second predictor (default is 1)
#'@param s3 Standard deviation of third predictor (default is 1)
#'@param all default is OFF, ON returns Power(All)
#'@return Power for Multiple Regression (ALL)
#'@export
#'
#'

MRC_all<-function(ry1=NULL, ry2=NULL, ry3=NULL, r12=NULL, r13=NULL, r23=NULL,n=NULL, alpha=.05, rep = 10000,
                        my=0,m1=0,m2=0,m3=0, sy=1,s1=1,s2=1,s3=1)
  {

pred<-NA
pred[is.null(r23)]<-2
pred[!is.null(r23)]<-3
vary<-NA
vary<-sy^2;var1<-s1^2;var2<-s2^2; var3<-s3^2

  if (pred=="2")

    {pop <- mvrnorm(100000, mu = c(my, m1, m2), Sigma = matrix(c(vary, ry1, ry2,
                                                                     ry1, var1, r12,
                                                                     ry2, r12, var2),
                                                                   ncol = 3), empirical = TRUE)
    pop2 = data.frame(pop)
    nruns = rep
    int = numeric(nruns)
    b1 = numeric(nruns)
    b2 = numeric(nruns)
    R2 = numeric(nruns)
    F = numeric(nruns)
    df1 = numeric(nruns)
    df2 = numeric(nruns)
    for (i in 1:nruns)
    {samp <- pop2[ sample(nrow(pop2), n), ]
    test <- lm(formula = X1 ~ X2+ X3, data = samp)
    c<-summary(test)
    int[i] = coef(summary(test))[1,4]
    b1[i] = coef(summary(test))[2,4] #grabs p from each analysis
    b2[i] = coef(summary(test))[3,4]
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
    pR2<-1-pf(F,df1, df2)
    Powerall$rejectR2 [pR2 < alpha] <- 1
    Powerall$rejectR2 [pR2 >= alpha] <- 0
    Power_R2<-mean(Powerall$rejectR2)
    PowerAll_R0<-mean(Reject.None)
    PowerAll_R1<-mean(Reject.One)
    PowerAll_R2<-mean(Reject.All)


    {print(paste("Sample size is ",n))}
    {print(paste("Power R2 = ", Power_R2))}
    {print(paste("Power b1 = ", Power_b1))}
    {print(paste("Power b2 = ", Power_b2))}
    {print(paste("Proportion Rejecting None = ", PowerAll_R0))}
    {print(paste("Proportion Rejecting One = ", PowerAll_R1))}
    {print(paste("Power ALL (Proportion Rejecting All) = ", PowerAll_R2))}
     }

  if (pred=="3")
    {
  pop <- mvrnorm(100000, mu = c(my, m1, m2, m3),
                 Sigma = matrix(c(vary, ry1, ry2, ry3,
                                  ry1, var1, r12, r13,
                                  ry2, r12, var2, r23,
                                  ry3, r13, r23, var3),
                 ncol = 4), empirical = TRUE)
  pop2 = data.frame(pop)
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
  {samp <- pop2[ sample(nrow(pop2), n), ]
  test <- lm(formula = X1 ~ X2+ X3+ X4, data = samp)
  c<-summary(test)
  int[i] = coef(summary(test))[1,4]
  b1[i] = coef(summary(test))[2,4] #grabs p from each analysis
  b2[i] = coef(summary(test))[3,4]
  b3[i] = coef(summary(test))[4,4]
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
  pR2<-1-pf(F,df1, df2)
  Powerall$rejectR2 [pR2 < alpha] <- 1
  Powerall$rejectR2 [pR2 >= alpha] <- 0
  Power_R2<-mean(Powerall$rejectR2)
  PowerAll_R0<-mean(Reject.None)
  PowerAll_R1<-mean(Reject.One)
  PowerAll_R2<-mean(Reject.Two)
  PowerAll_R3<-mean(Reject.All)


  {print(paste("Sample size is ",n))}
  {print(paste("Power R2 = ", Power_R2))}
  {print(paste("Power b1 = ", Power_b1))}
  {print(paste("Power b2 = ", Power_b2))}
  {print(paste("Power b3 = ", Power_b3))}
  {print(paste("Proportion Rejecting None = ", PowerAll_R0))}
  {print(paste("Proportion Rejecting One = ", PowerAll_R1))}
  {print(paste("Proportion Rejecting Two = ", PowerAll_R2))}
  {print(paste("Power ALL (Proportion Rejecting All) = ", PowerAll_R3))}
   }}
