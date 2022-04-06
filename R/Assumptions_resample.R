#'Compute power for Multiple Regression with Violated assumptions (Beta)
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
#'@param sy Skew of outcome variable
#'@param s1 Skew of first predictor
#'@param s2 Skew of second predictor
#'@param s3 Skew of third predictor
#'@param s4 Skew of fourth predictor
#'@param s5 Skew of fifth predictor
#'@param ky Kurtosis of outcome variable
#'@param k1 Kurtosis of first predictor
#'@param k2 Kurtosis of second predictor
#'@param k3 Kurtosis of third predictor
#'@param k4 Kurtosis of fourth predictor
#'@param k5 Kurtosis of fifth predictor
#'@param test type of test (boot,jack,perm)
#'@param nruns number of runs, default is 5000
#'
#'@examples
#'Assumptions_resample(ry1=.0, ry2=.3, ry3=.3, ry4=.1, r12 = .0,
#'r13=.0, r14=.0, r23=.0, r24=.0,r34=0,
#'sy=1,s1=2,s2=2, s3=1,s4=1, ky=1,k1=1,k2=1,
#'k3=1,k4=1, n=100, test="sqrt")
#'
#'@return Power for Multiple Regression with Non Normal Variables
#'@export
#'
#'

Assumptions_resample<-function(ry1=NULL, ry2=NULL, ry3=NULL, ry4=NULL, ry5=NULL,
                  r12=NULL, r13=NULL,r14=NULL,r15=NULL,
                  r23=NULL, r24=NULL, r25=NULL,
                  r34=NULL, r35=NULL,
                  r45=NULL, sy=NULL, s1=NULL, s2=NULL, s3=NULL, s4=NULL,s5=NULL,
                  ky=NULL,k1=NULL,k2=NULL,k3=NULL,k4=NULL,k5=NULL,
                  n=NULL, alpha=.05, test=NULL, reps=200, boots=500)
  {




set.seed(8675309)

pred<-NA
pred[!is.null(ry2)]<-2
pred[!is.null(ry3)]<-3
pred[!is.null(ry4)]<-4
pred[!is.null(ry5)]<-5
vary<-NA
vary<-1;var1<-1;var2<-1; var3<-1;var4<-1;var5<-1
b1<-NA;b2<-NA;b3<-NA;b4<-NA;b5<-NA
ll1<-NA;ll2<-NA;ll3<-NA;ll4<-NA;ll5<-NA
ul1<-NA;ul2<-NA;ul3<-NA;ul4<-NA;ul5<-NA
bb1<-NA;bb2<-NA;bb3<-NA;bb4<-NA;bb5<-NA


nruns = reps*boots
b1 = numeric(nruns)
b2 = numeric(nruns)
b3 = numeric(nruns)
b4 = numeric(nruns)
b5 = numeric(nruns)

rejectb1<-NA
rejectb2<-NA
rejectb3<-NA
rejectb4<-NA
rejectb5<-NA
sample<-n


# Bootstrap

if (pred==2 && test=="boot")
 {
    for (i in 1:nruns)
    {samp <- semTools::mvrnonnorm(sample, mu = c(0, 0, 0),
                                  Sigma = matrix(c(vary, ry1, ry2,
                                                   ry1, var1, r12,
                                                   ry2, r12, var2),
                                                 ncol = 3),
                                  skewness=c(sy,s1,s2) ,
                                  kurtosis=c(ky,k1,k2))
    samp = data.frame(samp)
    values <- stats::lm(formula = X1 ~ X2+ X3, data = samp)
    testx<-summary(values)
    b1[i] = testx$coefficients[2,4] #grabs p from each analysis
    b2[i] = testx$coefficients[3,4]}
    message("Adjustment: ", test)
    rejectb1<-NA
    rejectb1 [ b1 < alpha] <- 1
    rejectb1 [ b1 >= alpha] <- 0
    rejectb2<-NA
    rejectb2 [ b2 < alpha] <- 1
    rejectb2 [ b2 >= alpha] <- 0
    message("Sample size is ",n)
    message("Power b1 = ", mean(rejectb1))
    message("Power b2 = ", mean(rejectb2))
     }

if (pred==3 && test=="boot")
  {
  for (i in 1:nruns)
  {samp <- semTools::mvrnonnorm(sample, mu = c(0, 0, 0, 0),
                Sigma = matrix(c(vary, ry1, ry2, ry3,
                                ry1, var1, r12, r13,
                                ry2, r12, var2, r23,
                                ry3, r13, r23, var3),
                                ncol = 4),
                                skewness=c(sy,s1,s2,s3) ,
                                kurtosis=c(ky,k1,k2,k3))
                                samp = data.frame(samp)

  values<-stats::lm(X1~X2+X3+X4, samp)
  testx<-summary(values)

  b1[i] = testx$coefficients[2,4] #grabs p from each analysis
  b2[i] = testx$coefficients[3,4]
  b3[i] = testx$coefficients[4,4]}

  message("Adjustment: ", test)
  rejectb1<-NA
  rejectb1 [ b1 < alpha] <- 1
  rejectb1 [ b1 >= alpha] <- 0
  rejectb2<-NA
  rejectb2 [ b2 < alpha] <- 1
  rejectb2 [ b2 >= alpha] <- 0
  rejectb3<-NA
  rejectb3 [ b3 < alpha] <- 1
  rejectb3 [ b3 >= alpha] <- 0

  message("Sample size is ",n)
  message("Power b1 = ", mean(rejectb1))
  message("Power b2 = ", mean(rejectb2))
  message("Power b3 = ", mean(rejectb3))
  }

  if (pred==4 && test=="boot")
  {
    for (i in 1:nruns)
    {
    samp <- semTools::mvrnonnorm(sample, mu = c(0, 0, 0, 0, 0),
                                  Sigma = matrix(c(vary, ry1, ry2, ry3, ry4,
                                       ry1, var1, r12, r13, r14,
                                       ry2, r12, var2, r23, r24,
                                       ry3, r13, r23, var3, r34,
                                       ry4, r14, r24, r34, var4),
                                      ncol = 5),
                                  skewness=c(sy,s1,s2,s3,s4) ,
                                  kurtosis=c(ky,k1,k2,k3,k4))
    samp = data.frame(samp)

    values<-stats::lm(X1~X2+X3+X4+X5,samp)
    testx<-summary(values)

    b1[i] = testx$coefficients[2,4] #grabs p from each analysis
    b2[i] = testx$coefficients[3,4]
    b3[i] = testx$coefficients[4,4]
    b4[i] = testx$coefficients[5,4]}
    message("Adjustment: ", test)
    rejectb1<-NA
    rejectb1 [ b1 < alpha] <- 1
    rejectb1 [ b1 >= alpha] <- 0
    rejectb2<-NA
    rejectb2 [ b2 < alpha] <- 1
    rejectb2 [ b2 >= alpha] <- 0
    rejectb3<-NA
    rejectb3 [ b3 < alpha] <- 1
    rejectb3 [ b3 >= alpha] <- 0
    rejectb4<-NA
    rejectb4 [ b4 < alpha] <- 1
    rejectb4 [ b4 >= alpha] <- 0

    message("Sample size is ",n)
    message("Power b1 = ", mean(rejectb1))
    message("Power b2 = ", mean(rejectb2))
    message("Power b3 = ", mean(rejectb3))
    message("Power b4 = ", mean(rejectb4))

  }

    if (pred==5 && test=="boot")
    {
      for (i in 1:nruns)
      {
      samp <- semTools::mvrnonnorm(sample, mu = c(0, 0, 0,0,0,0),
                                     Sigma = matrix(c(vary, ry1, ry2, ry3, ry4,ry5,
                                     ry1, var1, r12, r13, r14,r15,
                                     ry2, r12, var2, r23, r24,r25,
                                     ry3, r13, r23, var3, r34,r35,
                                     ry4, r14, r24, r34, var4,r45,
                                     ry5,r15,r25,r35,r45,var5),
                                          ncol = 6),
                                    skewness=c(sy,s1,s2,s3,s4,s5) ,
                                    kurtosis=c(ky,k1,k2,k3,k4,k5))
        samp = data.frame(samp)
        values<-stats::lm(X1~X2+X3+X4+X5+X6,samp)
        testx<-summary(values)

        b1[i] = testx$coefficients[2,4] #grabs p from each analysis
        b2[i] = testx$coefficients[3,4]
        b3[i] = testx$coefficients[4,4]
        b4[i] = testx$coefficients[5,4]
        b5[i] = testx$coefficients[6,4]}

      message("Adjustment: ", test)
      rejectb1<-NA
      rejectb1 [ b1 < alpha] <- 1
      rejectb1 [ b1 >= alpha] <- 0
      rejectb2<-NA
      rejectb2 [ b2 < alpha] <- 1
      rejectb2 [ b2 >= alpha] <- 0
      rejectb3<-NA
      rejectb3 [ b3 < alpha] <- 1
      rejectb3 [ b3 >= alpha] <- 0
      rejectb4<-NA
      rejectb4 [ b4 < alpha] <- 1
      rejectb4 [ b4 >= alpha] <- 0
      rejectb5<-NA
      rejectb5 [ b5 < alpha] <- 1
      rejectb5 [ b5 >= alpha] <- 0

      message("Sample size is ",n)
      message("Power b1 = ", mean(rejectb1))
      message("Power b2 = ", mean(rejectb2))
      message("Power b3 = ", mean(rejectb3))
      message("Power b4 = ", mean(rejectb4))
      message("Power b5 = ", mean(rejectb5))}

# Permutation

if (pred==2 && test=="perm")
{
  for (i in 1:nruns)
  {samp <- semTools::mvrnonnorm(sample, mu = c(0, 0, 0),
                                Sigma = matrix(c(vary, ry1, ry2,
                                                 ry1, var1, r12,
                                                 ry2, r12, var2),
                                               ncol = 3),
                                skewness=c(sy,s1,s2) ,
                                kurtosis=c(ky,k1,k2))

  samp = data.frame(samp)
  invisible(capture.output(testx<-(summary(lmPerm::lmp(X1~ X2+ X3, data = samp)))))
  b1[i] = testx$coefficients[2,3] #grabs p from each analysis
  b2[i] = testx$coefficients[3,3]}
  message("Adjustment: ", test)
  rejectb1<-NA
  rejectb1 [ b1 < alpha] <- 1
  rejectb1 [ b1 >= alpha] <- 0
  rejectb2<-NA
  rejectb2 [ b2 < alpha] <- 1
  rejectb2 [ b2 >= alpha] <- 0
  message("Sample size is ",n)
  message("Power b1 = ", mean(rejectb1))
  message("Power b2 = ", mean(rejectb2))
}

if (pred==3 && test=="perm")
{
  for (i in 1:nruns)
  {samp <- semTools::mvrnonnorm(sample, mu = c(0, 0, 0, 0),
                                Sigma = matrix(c(vary, ry1, ry2, ry3,
                                                 ry1, var1, r12, r13,
                                                 ry2, r12, var2, r23,
                                                 ry3, r13, r23, var3),
                                               ncol = 4),
                                skewness=c(sy,s1,s2,s3) ,
                                kurtosis=c(ky,k1,k2,k3))
  samp = data.frame(samp)

  invisible(capture.output(testx<-summary(lmPerm::lmp(X1~ X2+ X3+X4, data = samp))))

  b1[i] = testx$coefficients[2,3] #grabs p from each analysis
  b2[i] = testx$coefficients[3,3]
  b3[i] = testx$coefficients[4,3]}

  message("Adjustment: ", test)
  rejectb1<-NA
  rejectb1 [ b1 < alpha] <- 1
  rejectb1 [ b1 >= alpha] <- 0
  rejectb2<-NA
  rejectb2 [ b2 < alpha] <- 1
  rejectb2 [ b2 >= alpha] <- 0
  rejectb3<-NA
  rejectb3 [ b3 < alpha] <- 1
  rejectb3 [ b3 >= alpha] <- 0

  message("Sample size is ",n)
  message("Power b1 = ", mean(rejectb1))
  message("Power b2 = ", mean(rejectb2))
  message("Power b3 = ", mean(rejectb3))
}

if (pred==4 && test=="perm")
{
  for (i in 1:nruns)
  {
    samp <- semTools::mvrnonnorm(sample, mu = c(0, 0, 0, 0, 0),
                                 Sigma = matrix(c(vary, ry1, ry2, ry3, ry4,
                                                  ry1, var1, r12, r13, r14,
                                                  ry2, r12, var2, r23, r24,
                                                  ry3, r13, r23, var3, r34,
                                                  ry4, r14, r24, r34, var4),
                                                ncol = 5),
                                 skewness=c(sy,s1,s2,s3,s4) ,
                                 kurtosis=c(ky,k1,k2,k3,k4))
    samp = data.frame(samp)

    invisible(capture.output(testx<-summary(lmPerm::lmp(X1~ X2+ X3+X4+X5, data = samp))))

    b1[i] = testx$coefficients[2,3] #grabs p from each analysis
    b2[i] = testx$coefficients[3,3]
    b3[i] = testx$coefficients[4,3]
    b4[i] = testx$coefficients[5,3]}
  message("Adjustment: ", test)
  rejectb1<-NA
  rejectb1 [ b1 < alpha] <- 1
  rejectb1 [ b1 >= alpha] <- 0
  rejectb2<-NA
  rejectb2 [ b2 < alpha] <- 1
  rejectb2 [ b2 >= alpha] <- 0
  rejectb3<-NA
  rejectb3 [ b3 < alpha] <- 1
  rejectb3 [ b3 >= alpha] <- 0
  rejectb4<-NA
  rejectb4 [ b4 < alpha] <- 1
  rejectb4 [ b4 >= alpha] <- 0

  message("Sample size is ",n)
  message("Power b1 = ", mean(rejectb1))
  message("Power b2 = ", mean(rejectb2))
  message("Power b3 = ", mean(rejectb3))
  message("Power b4 = ", mean(rejectb4))

}

if (pred==5 && test=="perm")
{
  for (i in 1:nruns)
  {
    samp <- semTools::mvrnonnorm(sample, mu = c(0, 0, 0,0,0,0),
                                 Sigma = matrix(c(vary, ry1, ry2, ry3, ry4,ry5,
                                                  ry1, var1, r12, r13, r14,r15,
                                                  ry2, r12, var2, r23, r24,r25,
                                                  ry3, r13, r23, var3, r34,r35,
                                                  ry4, r14, r24, r34, var4,r45,
                                                  ry5,r15,r25,r35,r45,var5),
                                                ncol = 6),
                                 skewness=c(sy,s1,s2,s3,s4,s5) ,
                                 kurtosis=c(ky,k1,k2,k3,k4,k5))
    samp = data.frame(samp)
    invisible(capture.output(testx<-summary(lmPerm::lmp(X1~ X2+ X3+X4+X5+X6, data = samp))))

    b1[i] = testx$coefficients[2,3] #grabs p from each analysis
    b2[i] = testx$coefficients[3,3]
    b3[i] = testx$coefficients[4,3]
    b4[i] = testx$coefficients[5,3]
    b5[i] = testx$coefficients[6,3]}

  message("Adjustment: ", test)
  rejectb1<-NA
  rejectb1 [ b1 < alpha] <- 1
  rejectb1 [ b1 >= alpha] <- 0
  rejectb2<-NA
  rejectb2 [ b2 < alpha] <- 1
  rejectb2 [ b2 >= alpha] <- 0
  rejectb3<-NA
  rejectb3 [ b3 < alpha] <- 1
  rejectb3 [ b3 >= alpha] <- 0
  rejectb4<-NA
  rejectb4 [ b4 < alpha] <- 1
  rejectb4 [ b4 >= alpha] <- 0
  rejectb5<-NA
  rejectb5 [ b5 < alpha] <- 1
  rejectb5 [ b5 >= alpha] <- 0

  message("Sample size is ",n)
  message("Power b1 = ", mean(rejectb1))
  message("Power b2 = ", mean(rejectb2))
  message("Power b3 = ", mean(rejectb3))
  message("Power b4 = ", mean(rejectb4))
  message("Power b5 = ", mean(rejectb5))}

# Jackknife

if (pred==2 && test=="knife")
{
  for (i in 1:nruns)
  {samp <- semTools::mvrnonnorm(sample, mu = c(0, 0, 0),
                                Sigma = matrix(c(vary, ry1, ry2,
                                                 ry1, var1, r12,
                                                 ry2, r12, var2),
                                               ncol = 3),
                                skewness=c(sy,s1,s2) ,
                                kurtosis=c(ky,k1,k2))
  samp = data.frame(samp)
  testx <- pls:::plsr(X1~ X2+ X3, data = samp, ncomp=2,validation = "LOO", jackknife = TRUE)
  test2<-pls::jack.test(testx, ncomp = 2)
  b1[i]<-test2$pvalues[1]
  b2[i]<-test2$pvalues[2]}
  message("Adjustment: ", test)
  rejectb1<-NA
  rejectb1 [ b1 < alpha] <- 1
  rejectb1 [ b1 >= alpha] <- 0
  rejectb2<-NA
  rejectb2 [ b2 < alpha] <- 1
  rejectb2 [ b2 >= alpha] <- 0
  message("Sample size is ",n)
  message("Power b1 = ", mean(rejectb1))
  message("Power b2 = ", mean(rejectb2))
}

if (pred==3 && test=="knife")
{
  for (i in 1:nruns)
  {samp <- semTools::mvrnonnorm(sample, mu = c(0, 0, 0, 0),
                                Sigma = matrix(c(vary, ry1, ry2, ry3,
                                                 ry1, var1, r12, r13,
                                                 ry2, r12, var2, r23,
                                                 ry3, r13, r23, var3),
                                               ncol = 4),
                                skewness=c(sy,s1,s2,s3) ,
                                kurtosis=c(ky,k1,k2,k3))
  samp = data.frame(samp)

  testx <- pls:::plsr(X1~ X2+ X3+ X4, data = samp, ncomp=2,validation = "LOO", jackknife = TRUE)
  test2<-pls::jack.test(testx, ncomp = 2)
  b1[i]<-test2$pvalues[1]
  b2[i]<-test2$pvalues[2]
  b3[i]<-test2$pvalues[3]
}

  message("Adjustment: ", test)
  rejectb1<-NA
  rejectb1 [ b1 < alpha] <- 1
  rejectb1 [ b1 >= alpha] <- 0
  rejectb2<-NA
  rejectb2 [ b2 < alpha] <- 1
  rejectb2 [ b2 >= alpha] <- 0
  rejectb3<-NA
  rejectb3 [ b3 < alpha] <- 1
  rejectb3 [ b3 >= alpha] <- 0

  message("Sample size is ",n)
  message("Power b1 = ", mean(rejectb1))
  message("Power b2 = ", mean(rejectb2))
  message("Power b3 = ", mean(rejectb3))
}

if (pred==4 && test=="knife")
{
  for (i in 1:nruns)
  {
    samp <- semTools::mvrnonnorm(sample, mu = c(0, 0, 0, 0, 0),
                                 Sigma = matrix(c(vary, ry1, ry2, ry3, ry4,
                                                  ry1, var1, r12, r13, r14,
                                                  ry2, r12, var2, r23, r24,
                                                  ry3, r13, r23, var3, r34,
                                                  ry4, r14, r24, r34, var4),
                                                ncol = 5),
                                 skewness=c(sy,s1,s2,s3,s4) ,
                                 kurtosis=c(ky,k1,k2,k3,k4))
    samp = data.frame(samp)

    testx <- pls:::plsr(X1~ X2+ X3+ X4+ X5, data = samp, ncomp=2,validation = "LOO", jackknife = TRUE)
    test2<-pls::jack.test(testx, ncomp = 2)
    b1[i]<-test2$pvalues[1]
    b2[i]<-test2$pvalues[2]
    b3[i]<-test2$pvalues[3]
    b4[i]<-test2$pvalues[4]
  }
  message("Adjustment: ", test)
  rejectb1<-NA
  rejectb1 [ b1 < alpha] <- 1
  rejectb1 [ b1 >= alpha] <- 0
  rejectb2<-NA
  rejectb2 [ b2 < alpha] <- 1
  rejectb2 [ b2 >= alpha] <- 0
  rejectb3<-NA
  rejectb3 [ b3 < alpha] <- 1
  rejectb3 [ b3 >= alpha] <- 0
  rejectb4<-NA
  rejectb4 [ b4 < alpha] <- 1
  rejectb4 [ b4 >= alpha] <- 0

  message("Sample size is ",n)
  message("Power b1 = ", mean(rejectb1))
  message("Power b2 = ", mean(rejectb2))
  message("Power b3 = ", mean(rejectb3))
  message("Power b4 = ", mean(rejectb4))

}

if (pred==5 && test=="knife")
{
  for (i in 1:nruns)
  {
    samp <- semTools::mvrnonnorm(sample, mu = c(0, 0, 0,0,0,0),
                                 Sigma = matrix(c(vary, ry1, ry2, ry3, ry4,ry5,
                                                  ry1, var1, r12, r13, r14,r15,
                                                  ry2, r12, var2, r23, r24,r25,
                                                  ry3, r13, r23, var3, r34,r35,
                                                  ry4, r14, r24, r34, var4,r45,
                                                  ry5,r15,r25,r35,r45,var5),
                                                ncol = 6),
                                 skewness=c(sy,s1,s2,s3,s4,s5) ,
                                 kurtosis=c(ky,k1,k2,k3,k4,k5))
    samp = data.frame(samp)
    testx <- pls:::plsr(X1~ X2+ X3+ X4+ X5+ X6, data = samp, ncomp=2,validation = "LOO", jackknife = TRUE)
    test2<-pls::jack.test(testx, ncomp = 2)
    b1[i]<-test2$pvalues[1]
    b2[i]<-test2$pvalues[2]
    b3[i]<-test2$pvalues[3]
    b4[i]<-test2$pvalues[4]
    b5[i]<-test2$pvalues[5]}

  message("Adjustment: ", test)
  rejectb1<-NA
  rejectb1 [ b1 < alpha] <- 1
  rejectb1 [ b1 >= alpha] <- 0
  rejectb2<-NA
  rejectb2 [ b2 < alpha] <- 1
  rejectb2 [ b2 >= alpha] <- 0
  rejectb3<-NA
  rejectb3 [ b3 < alpha] <- 1
  rejectb3 [ b3 >= alpha] <- 0
  rejectb4<-NA
  rejectb4 [ b4 < alpha] <- 1
  rejectb4 [ b4 >= alpha] <- 0
  rejectb5<-NA
  rejectb5 [ b5 < alpha] <- 1
  rejectb5 [ b5 >= alpha] <- 0

  message("Sample size is ",n)
  message("Power b1 = ", mean(rejectb1))
  message("Power b2 = ", mean(rejectb2))
  message("Power b3 = ", mean(rejectb3))
  message("Power b4 = ", mean(rejectb4))
  message("Power b5 = ", mean(rejectb5))}


}

