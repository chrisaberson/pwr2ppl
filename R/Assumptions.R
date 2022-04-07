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
#'@param test type of test (none, sqrt, log, inv, robust, boot, quantile, hc0, hc1, hc2, hc3)
#'@param nruns number of runs, default is 500
#'
#'@examples
#'\donttest{Assumptions(ry1=.0,ry2=.3,r12=.3,sy=1,s1=2,s2=2,ky=1,k1=1,k2=1,n=100,nruns=20,test="sqrt")}
#'
#'@return Power for Resampled Multiple Regression with Non Normal Variables
#'@export
#'
#'

Assumptions<-function(ry1=NULL, ry2=NULL, ry3=NULL, ry4=NULL, ry5=NULL,
                  r12=NULL, r13=NULL,r14=NULL,r15=NULL,
                  r23=NULL, r24=NULL, r25=NULL,
                  r34=NULL, r35=NULL,
                  r45=NULL, sy=NULL, s1=NULL, s2=NULL, s3=NULL, s4=NULL,s5=NULL,
                  ky=NULL,k1=NULL,k2=NULL,k3=NULL,k4=NULL,k5=NULL,
                  n=NULL, alpha=.05, test=NULL, nruns=500)
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

b1 = numeric(n)
b2 = numeric(n)
b3 = numeric(n)
b4 = numeric(n)
b5 = numeric(n)
ll1 = numeric(n)
ll2 = numeric(n)
ll3 = numeric(n)
ll4 = numeric(n)
ll5 = numeric(n)

ul1 = numeric(n)
ul2 = numeric(n)
ul3 = numeric(n)
ul4 = numeric(n)
ul5 = numeric(n)

rejectb1<-NA
rejectb2<-NA
rejectb3<-NA
rejectb4<-NA
rejectb5<-NA


#For boots

reps = 200
runs = reps*nruns
bb1 = numeric(runs)
bb2 = numeric(runs)
bb3 = numeric(runs)
bb4 = numeric(runs)
bb5 = numeric(runs)

# Build population

if (pred==2)
{pop <- semTools::mvrnonnorm(100000, mu = c(0, 0, 0), Sigma = matrix(c(vary, ry1, ry2,
                                                                       ry1, var1, r12,
                                                                       ry2, r12, var2),
                                                                     ncol = 3),
                             skewness=c(sy,s1,s2) ,
                             kurtosis=c(ky,k1,k2))
pop2 = data.frame(pop)
}

if (pred==3)
{
  pop <- semTools::mvrnonnorm(100000, mu = c(0, 0, 0, 0), Sigma = matrix(c(vary, ry1, ry2, ry3,
                                                                           ry1, var1, r12, r13,
                                                                           ry2, r12, var2, r23,
                                                                           ry3, r13, r23, var3),
                                                                         ncol = 4),
                              skewness=c(sy,s1,s2,s3) ,
                              kurtosis=c(ky,k1,k2,k3))
  pop2 = data.frame(pop)
}

  if (pred==4)
  {
    pop <- semTools::mvrnonnorm(100000, mu = c(0, 0, 0, 0,0), Sigma = matrix(c(vary, ry1, ry2, ry3, ry4,
                                                                               ry1, var1, r12, r13, r14,
                                                                               ry2, r12, var2, r23, r24,
                                                                               ry3, r13, r23, var3, r34,
                                                                               ry4, r14, r24, r34, var4),
                                                                             ncol = 5),
                                skewness=c(sy,s1,s2,s3,s4) ,
                                kurtosis=c(ky,k1,k2,k3,k4))
    pop2 = data.frame(pop)
  }

if (pred==5)
{

  pop <- semTools::mvrnonnorm(100000, mu = c(0, 0, 0, 0,0,0), Sigma = matrix(c(vary, ry1, ry2, ry3, ry4,ry5,
                                                                               ry1, var1, r12, r13, r14,r15,
                                                                               ry2, r12, var2, r23, r24,r25,
                                                                               ry3, r13, r23, var3, r34,r35,
                                                                               ry4, r14, r24, r34, var4,r45,
                                                                               ry5,r15,r25,r35,r45,var5),
                                                                             ncol = 6),
                              skewness=c(sy,s1,s2,s3,s4,s5) ,
                              kurtosis=c(ky,k1,k2,k3,k4,k5))
  pop2 = data.frame(pop)
}


# No adjustment

if (pred==2 && test=="none")
 {
    for (i in 1:nruns)
    {samp <- pop2[ sample(nrow(pop2), n), ]
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

if (pred==3 && test=="none")
        {
  for (i in 1:nruns)
  {samp <- pop2[ sample(nrow(pop2), n), ]
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

  if (pred==4 && test=="none")
  {
    for (i in 1:nruns)
    {samp <- pop2[ sample(nrow(pop2), n), ]
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

    if (pred==5 && test=="none")
    {
      for (i in 1:nruns)
      {samp <- pop2[ sample(nrow(pop2), n), ]
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

# Square Root Transform

if (pred==2 && test=="sqrt")
{
pop2$X1<-sqrt(pop2$X1+10)
pop2$X2<-sqrt(pop2$X2+10)
pop2$X3<-sqrt(pop2$X3+10)

for (i in 1:nruns)
{samp <- pop2[ sample(nrow(pop2), n), ]
values <- stats::lm(formula = X1 ~ X2+ X3, data = samp)
testx<-summary(values)
c<-summary(testx$coefficients)
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

if (pred==3 && test=="sqrt")
{
  pop2$X1<-sqrt(pop2$X1+10)
  pop2$X2<-sqrt(pop2$X2+10)
  pop2$X3<-sqrt(pop2$X3+10)
  pop2$X4<-sqrt(pop2$X4+10)
  for (i in 1:nruns)
  {samp <- pop2[ sample(nrow(pop2), n), ]
  values<-stats::lm(X1~X2+X3+X4,samp)
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

if (pred==4 && test=="sqrt")
{
  pop2$X1<-sqrt(pop2$X1+10)
  pop2$X2<-sqrt(pop2$X2+10)
  pop2$X3<-sqrt(pop2$X3+10)
  pop2$X4<-sqrt(pop2$X4+10)
  pop2$X5<-sqrt(pop2$X5+10)
  for (i in 1:nruns)
  {samp <- pop2[ sample(nrow(pop2), n), ]
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
  message("Power b4 = ", mean(rejectb4))}


if (pred==5 && test=="sqrt")
{
  pop2 = data.frame(pop)
  pop2$X1<-sqrt(pop2$X1+10)
  pop2$X2<-sqrt(pop2$X2+10)
  pop2$X3<-sqrt(pop2$X3+10)
  pop2$X4<-sqrt(pop2$X4+10)
  pop2$X5<-sqrt(pop2$X5+10)
  pop2$X6<-sqrt(pop2$X6+10)

  for (i in 1:nruns)
  {samp <- pop2[ sample(nrow(pop2), n), ]
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

# Log Transform

if (pred==2 && test=="log")
{
pop2$X1<-log10(pop2$X1+10)
pop2$X2<-log10(pop2$X2+10)
pop2$X3<-log10(pop2$X3+10)
for (i in 1:nruns)
{samp <- pop2[ sample(nrow(pop2), n), ]
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

if (pred==3 && test=="log")
{
  pop2$X1<-log10(pop2$X1+10)
  pop2$X2<-log10(pop2$X2+10)
  pop2$X3<-log10(pop2$X3+10)
  pop2$X4<-log10(pop2$X4+10)
  for (i in 1:nruns)
  {samp <- pop2[ sample(nrow(pop2), n), ]
  values<-stats::lm(X1~X2+X3+X4,samp)
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
  message("Power b3 = ", mean(rejectb3)) }

if (pred==4 && test=="log")
{
  pop2$X1<-log10(pop2$X1+10)
  pop2$X2<-log10(pop2$X2+10)
  pop2$X3<-log10(pop2$X3+10)
  pop2$X4<-log10(pop2$X4+10)
  pop2$X5<-log10(pop2$X5+10)
  for (i in 1:nruns)
  {samp <- pop2[ sample(nrow(pop2), n), ]
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
  message("Power b4 = ", mean(rejectb4))}


if (pred==5 && test=="log")
{
  pop2$X1<-log10(pop2$X1+10)
  pop2$X2<-log10(pop2$X2+10)
  pop2$X3<-log10(pop2$X3+10)
  pop2$X4<-log10(pop2$X4+10)
  pop2$X5<-log10(pop2$X5+10)
  pop2$X6<-log10(pop2$X6+10)

  for (i in 1:nruns)
  {samp <- pop2[ sample(nrow(pop2), n), ]
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

# Inverse Transform

if (pred==2 && test=="inv")
{
pop2$X1<-1/(pop2$X1+10)
pop2$X2<-1/(pop2$X2+10)
pop2$X3<-1/(pop2$X3+10)
for (i in 1:nruns)
{samp <- pop2[ sample(nrow(pop2), n), ]
values <- stats::lm(formula = X1 ~ X2+ X3, data = samp)
testx<-summary(values)
c<-summary(testx$coefficients)
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

if (pred==3 && test=="inv")
{
  pop2$X1<-1/(pop2$X1+10)
  pop2$X2<-1/(pop2$X2+10)
  pop2$X3<-1/(pop2$X3+10)
  pop2$X4<-1/(pop2$X4+10)
  for (i in 1:nruns)
  {samp <- pop2[ sample(nrow(pop2), n), ]
  values<-stats::lm(X1~X2+X3+X4,samp)
  testx<-summary(values)

  b1[i] = testx$coefficients[2,4] #grabs p from each analysis
  b2[i] = testx$coefficients[3,4]
  b3[i] = testx$coefficients[4,4]}
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
  message("Power b3 = ", mean(rejectb3)) }

if (pred==4 && test=="inv")
  {
  pop2$X1<-1/(pop2$X1+10)
  pop2$X2<-1/(pop2$X2+10)
  pop2$X3<-1/(pop2$X3+10)
  pop2$X4<-1/(pop2$X4+10)
  pop2$X5<-1/(pop2$X5+10)
  for (i in 1:nruns)
  {samp <- pop2[ sample(nrow(pop2), n), ]
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
  message("Power b4 = ", mean(rejectb4))}


if (pred==5 && test=="inv")
{
  pop2$X1<-1/(pop2$X1+10)
  pop2$X2<-1/(pop2$X2+10)
  pop2$X3<-1/(pop2$X3+10)
  pop2$X4<-1/(pop2$X4+10)
  pop2$X5<-1/(pop2$X5+10)
  pop2$X6<-1/(pop2$X6+10)
  for (i in 1:nruns)
  {samp <- pop2[ sample(nrow(pop2), n), ]
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

# Robust

if (pred==2 && test=="robust")
{
for (i in 1:nruns)
{samp <- pop2[ sample(nrow(pop2),n), ]
values <- MASS::rlm(formula = X1~ X2+ X3, data = samp)
testx<-summary(values)
b1[i] = 2*pt(abs(testx$coefficients[2,3]),n-pred,lower.tail=FALSE)
b2[i] = 2*pt(abs(testx$coefficients[3,3]),n-pred,lower.tail=FALSE)
}
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

if (pred==3 && test=="robust")
{
  for (i in 1:nruns)
  {samp <- pop2[ sample(nrow(pop2),n), ]
  values <- MASS::rlm(formula = X1~ X2+ X3 +X4, data = samp)
  testx<-summary(values)
  b1[i] = 2*pt(abs(testx$coefficients[2,3]),n-pred,lower.tail=FALSE)
  b2[i] = 2*pt(abs(testx$coefficients[3,3]),n-pred,lower.tail=FALSE)
  b3[i] = 2*pt(abs(testx$coefficients[4,3]),n-pred,lower.tail=FALSE)
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
  message("Power b3 = ", mean(rejectb3)) }

if (pred==4 && test=="robust")
{
  for (i in 1:nruns)
  {samp <- pop2[ sample(nrow(pop2),n), ]
  values <- MASS::rlm(formula = X1~ X2+ X3+X4+X5, data = samp)
  testx<-summary(values)
  b1[i] = 2*pt(abs(testx$coefficients[2,3]),n-pred,lower.tail=FALSE)
  b2[i] = 2*pt(abs(testx$coefficients[3,3]),n-pred,lower.tail=FALSE)
  b3[i] = 2*pt(abs(testx$coefficients[4,3]),n-pred,lower.tail=FALSE)
  b4[i] = 2*pt(abs(testx$coefficients[5,3]),n-pred,lower.tail=FALSE)
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
  message("Power b4 = ", mean(rejectb4))}


if (pred==5 && test=="robust")
{
for (i in 1:nruns)
{samp <- pop2[ sample(nrow(pop2),n), ]
values <- MASS::rlm(formula = X1~ X2+ X3 +X4+X5+X6, data = samp)
testx<-summary(values)
b1[i] = 2*pt(abs(testx$coefficients[2,3]),n-pred,lower.tail=FALSE)
b2[i] = 2*pt(abs(testx$coefficients[3,3]),n-pred,lower.tail=FALSE)
b3[i] = 2*pt(abs(testx$coefficients[4,3]),n-pred,lower.tail=FALSE)
b4[i] = 2*pt(abs(testx$coefficients[5,3]),n-pred,lower.tail=FALSE)
b5[i] = 2*pt(abs(testx$coefficients[6,3]),n-pred,lower.tail=FALSE)
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
  rejectb5<-NA
  rejectb5 [ b5 < alpha] <- 1
  rejectb5 [ b5 >= alpha] <- 0

  message("Sample size is ",n)
  message("Power b1 = ", mean(rejectb1))
  message("Power b2 = ", mean(rejectb2))
  message("Power b3 = ", mean(rejectb3))
  message("Power b4 = ", mean(rejectb4))
  message("Power b5 = ", mean(rejectb5))}

# Quantile

if (pred==2 && test=="quantile")
{
for (i in 1:nruns)
{samp <- pop2[ sample(nrow(pop2),n), ]
testx <- quantreg::rq(X1~X2+X3, data=samp)
xx<-broom::tidy(testx)
ll1[i] = xx[2,3]
ll2[i] = xx[3,3]
ul1[i] = xx[2,4]
ul2[i] = xx[3,4]
}
  message("Adjustment: ", test)
  rejectb1<-NA
  rejectb1 [(ll1< 0 & ul1<0)|(ll1> 0 & ul1>0)] <- 1
  rejectb1 [(ll1< 0 & ul1>0)|(ll1> 0 & ul1<0)] <- 0
  rejectb2<-NA
  rejectb2 [(ll2< 0 & ul2<0)|(ll2> 0 & ul2>0)] <- 1
  rejectb2 [(ll2< 0 & ul2>0)|(ll2> 0 & ul2<0)] <- 0
  message("Sample size is ",n)
  message("Power b1 = ", mean(rejectb1))
  message("Power b2 = ", mean(rejectb2))
}

if (pred==3 && test=="quantile")
{
for (i in 1:nruns)
{samp <- pop2[ sample(nrow(pop2),n), ]
testx <- quantreg::rq(X1~X2+X3+X4, data=samp)
xx<-broom::tidy(testx)
ll1[i] = xx[2,3]
ll2[i] = xx[3,3]
ll3[i] = xx[4,3]
ul1[i] = xx[2,4]
ul2[i] = xx[3,4]
ul3[i] = xx[4,4]
}
  message("Adjustment: ", test)
  rejectb1<-NA
  rejectb1 [(ll1< 0 & ul1<0)|(ll1> 0 & ul1>0)] <- 1
  rejectb1 [(ll1< 0 & ul1>0)|(ll1> 0 & ul1<0)] <- 0
  rejectb2<-NA
  rejectb2 [(ll2< 0 & ul2<0)|(ll2> 0 & ul2>0)] <- 1
  rejectb2 [(ll2< 0 & ul2>0)|(ll2> 0 & ul2<0)] <- 0
  rejectb3<-NA
  rejectb3 [(ll3< 0 & ul3<0)|(ll3> 0 & ul3>0)] <- 1
  rejectb3 [(ll3< 0 & ul3>0)|(ll3> 0 & ul3<0)] <- 0

  message("Sample size is ",n)
  message("Power b1 = ", mean(rejectb1))
  message("Power b2 = ", mean(rejectb2))
  message("Power b3 = ", mean(rejectb3)) }

if (pred==4 && test=="quantile")
{
  for (i in 1:nruns)
  {samp <- pop2[ sample(nrow(pop2),n), ]
  testx <- quantreg::rq(X1~X2+X3+X4+X5, data=samp)
  xx<-broom::tidy(testx)
  ll1[i] = xx[2,3]
  ll2[i] = xx[3,3]
  ll3[i] = xx[4,3]
  ll4[i] = xx[5,3]
  ul1[i] = xx[2,4]
  ul2[i] = xx[3,4]
  ul3[i] = xx[4,4]
  ul4[i] = xx[5,4]
  }
  message("Adjustment: ", test)
  rejectb1<-NA
  rejectb1 [(ll1< 0 & ul1<0)|(ll1> 0 & ul1>0)] <- 1
  rejectb1 [(ll1< 0 & ul1>0)|(ll1> 0 & ul1<0)] <- 0
  rejectb2<-NA
  rejectb2 [(ll2< 0 & ul2<0)|(ll2> 0 & ul2>0)] <- 1
  rejectb2 [(ll2< 0 & ul2>0)|(ll2> 0 & ul2<0)] <- 0
  rejectb3<-NA
  rejectb3 [(ll3< 0 & ul3<0)|(ll3> 0 & ul3>0)] <- 1
  rejectb3 [(ll3< 0 & ul3>0)|(ll3> 0 & ul3<0)] <- 0
  rejectb4<-NA
  rejectb4 [(ll4< 0 & ul4<0)|(ll4> 0 & ul4>0)] <- 1
  rejectb4 [(ll4< 0 & ul4>0)|(ll4> 0 & ul4<0)] <- 0

  message("Sample size is ",n)
  message("Power b1 = ", mean(rejectb1))
  message("Power b2 = ", mean(rejectb2))
  message("Power b3 = ", mean(rejectb3))
  message("Power b4 = ", mean(rejectb4))}


if (pred==5 && test=="quantile")
{
  for (i in 1:nruns)
  {samp <- pop2[ sample(nrow(pop2),n), ]
  testx <- quantreg::rq(X1~X2+X3+X4+X5+X6, data=samp)
  xx<-broom::tidy(testx)
  ll1[i] = xx[2,3]
  ll2[i] = xx[3,3]
  ll3[i] = xx[4,3]
  ll4[i] = xx[5,3]
  ll5[i] = xx[6,3]
  ul1[i] = xx[2,4]
  ul2[i] = xx[3,4]
  ul3[i] = xx[4,4]
  ul4[i] = xx[5,4]
  ul5[i] = xx[6,4]
  }
  message("Adjustment: ", test)
  rejectb1<-NA
  rejectb1 [(ll1< 0 & ul1<0)|(ll1> 0 & ul1>0)] <- 1
  rejectb1 [(ll1< 0 & ul1>0)|(ll1> 0 & ul1<0)] <- 0
  rejectb2<-NA
  rejectb2 [(ll2< 0 & ul2<0)|(ll2> 0 & ul2>0)] <- 1
  rejectb2 [(ll2< 0 & ul2>0)|(ll2> 0 & ul2<0)] <- 0
  rejectb3<-NA
  rejectb3 [(ll3< 0 & ul3<0)|(ll3> 0 & ul3>0)] <- 1
  rejectb3 [(ll3< 0 & ul3>0)|(ll3> 0 & ul3<0)] <- 0
  rejectb4<-NA
  rejectb4 [(ll4< 0 & ul4<0)|(ll4> 0 & ul4>0)] <- 1
  rejectb4 [(ll4< 0 & ul4>0)|(ll4> 0 & ul4<0)] <- 0
  rejectb5<-NA
  rejectb5 [(ll5< 0 & ul5<0)|(ll5> 0 & ul5>0)] <- 1
  rejectb5 [(ll5< 0 & ul5>0)|(ll5> 0 & ul5<0)] <- 0

  message("Sample size is ",n)
  message("Power b1 = ", mean(rejectb1))
  message("Power b2 = ", mean(rejectb2))
  message("Power b3 = ", mean(rejectb3))
  message("Power b4 = ", mean(rejectb4))
  message("Power b5 = ", mean(rejectb5))}


# Bootstrap

if (pred==2 && test=="boot")
{for (i in 1:runs)
{samp <- pop2[ sample(nrow(pop2), n, replace=TRUE), ]
values <- stats::lm(formula = X1~ X2+ X3, data = samp)
testx<-summary(values)
bb1[i] = testx$coefficients[2,4] #grabs p from each analysis
bb2[i] = testx$coefficients[3,4]
}
  message("Adjustment: ", test)
  rejectb1<-NA
  rejectb1 [ bb1 < alpha] <- 1
  rejectb1 [ bb1 >= alpha] <- 0
  rejectb2<-NA
  rejectb2 [ bb2 < alpha] <- 1
  rejectb2 [ bb2 >= alpha] <- 0
  message("Sample size is ",n)
  message("Power b1 = ", mean(rejectb1))
  message("Power b2 = ", mean(rejectb2))
}

if (pred==3 && test=="boot")
{
  for (i in 1:runs)
  {samp <- pop2[ sample(nrow(pop2), n, replace=TRUE), ]
  values <- stats::lm(formula = X1~ X2+ X3+X4, data = samp)
  testx<-summary(values)
  bb1[i] = testx$coefficients[2,4] #grabs p from each analysis
  bb2[i] = testx$coefficients[3,4]
  bb3[i] = testx$coefficients[4,4]
  }
  message("Adjustment: ", test)
  rejectb1<-NA
  rejectb1 [ bb1 < alpha] <- 1
  rejectb1 [ bb1 >= alpha] <- 0
  rejectb2<-NA
  rejectb2 [ bb2 < alpha] <- 1
  rejectb2 [ bb2 >= alpha] <- 0
  rejectb3<-NA
  rejectb3 [ bb3 < alpha] <- 1
  rejectb3 [ bb3 >= alpha] <- 0

  message("Sample size is ",n)
  message("Power b1 = ", mean(rejectb1))
  message("Power b2 = ", mean(rejectb2))
  message("Power b3 = ", mean(rejectb3)) }


if (pred==4 && test=="boot")
{
  for (i in 1:runs)
  {samp <- pop2[ sample(nrow(pop2), n, replace=TRUE), ]
  values <- stats::lm(formula = X1~ X2+X3+X4+X5, data = samp)
  testx<-summary(values)

  bb1[i] = testx$coefficients[2,4] #grabs p from each analysis
  bb2[i] = testx$coefficients[3,4]
  bb3[i] = testx$coefficients[4,4]
  bb4[i] = testx$coefficients[5,4]
  }
  message("Adjustment: ", test)
  rejectb1<-NA
  rejectb1 [ bb1 < alpha] <- 1
  rejectb1 [ bb1 >= alpha] <- 0
  rejectb2<-NA
  rejectb2 [ bb2 < alpha] <- 1
  rejectb2 [ bb2 >= alpha] <- 0
  rejectb3<-NA
  rejectb3 [ bb3 < alpha] <- 1
  rejectb3 [ bb3 >= alpha] <- 0
  rejectb4<-NA
  rejectb4 [ bb4 < alpha] <- 1
  rejectb4 [ bb4 >= alpha] <- 0

  message("Sample size is ",n)
  message("Power b1 = ", mean(rejectb1))
  message("Power b2 = ", mean(rejectb2))
  message("Power b3 = ", mean(rejectb3))
  message("Power b4 = ", mean(rejectb4))}


if (pred==5 && test=="boot")
{
  for (i in 1:runs)
  {samp <- pop2[ sample(nrow(pop2), n, replace=TRUE), ]
  values <- stats::lm(formula = X1~ X2+ X3+X4+X5+X6, data = samp)
  testx<-summary(values)

  bb1[i] = testx$coefficients[2,4] #grabs p from each analysis
  bb2[i] = testx$coefficients[3,4]
  bb3[i] = testx$coefficients[4,4]
  bb4[i] = testx$coefficients[5,4]
  bb5[i] = testx$coefficients[6,4]
  }
  message("Adjustment: ", test)
  rejectb1<-NA
  rejectb1 [ bb1 < alpha] <- 1
  rejectb1 [ bb1 >= alpha] <- 0
  rejectb2<-NA
  rejectb2 [ bb2 < alpha] <- 1
  rejectb2 [ bb2 >= alpha] <- 0
  rejectb3<-NA
  rejectb3 [ bb3 < alpha] <- 1
  rejectb3 [ bb3 >= alpha] <- 0
  rejectb4<-NA
  rejectb4 [ bb4 < alpha] <- 1
  rejectb4 [ bb4 >= alpha] <- 0
  rejectb5<-NA
  rejectb5 [ bb5 < alpha] <- 1
  rejectb5 [ bb5 >= alpha] <- 0

  message("Sample size is ",n)
  message("Power b1 = ", mean(rejectb1))
  message("Power b2 = ", mean(rejectb2))
  message("Power b3 = ", mean(rejectb3))
  message("Power b4 = ", mean(rejectb4))
  message("Power b5 = ", mean(rejectb5))}

# Heteroscedasticity Adjusted hc1

if (pred==2 && test=="hc1")
{
for (i in 1:nruns)
{samp <- pop2[ sample(nrow(pop2),n), ]
reg<-lm(X1~ X2+ X3, data = samp)
hccm<-car::hccm(reg, type="hc1")
testx<-lmtest::coeftest(reg,vcov.=hccm)
b1[i] = testx[2,4]
b2[i] = testx[3,4]
}
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

if (pred==3 && test=="hc1")
{
for (i in 1:nruns)
{samp <- pop2[ sample(nrow(pop2),n), ]
reg<-lm(X1~ X2+ X3+X4, data = samp)
hccm<-car::hccm(reg, type="hc1")
testx<-lmtest::coeftest(reg,vcov.=hccm)
b1[i] = testx[2,4]
b2[i] = testx[3,4]
b3[i] = testx[4,4]
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

if (pred==4 && test=="hc1")
{for (i in 1:nruns)
{samp <- pop2[ sample(nrow(pop2),n), ]
reg<-lm(X1~ X2+ X3+X4+X5, data = samp)
hccm<-car::hccm(reg, type="hc1")
testx<-lmtest::coeftest(reg,vcov.=hccm)
b1[i] = testx[2,4]
b2[i] = testx[3,4]
b3[i] = testx[4,4]
b4[i] = testx[5,4]
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
  message("Power b4 = ", mean(rejectb4))}


if (pred==5 && test=="hc1")
{
  for (i in 1:nruns)
  {samp <- pop2[ sample(nrow(pop2),n), ]
  reg<-lm(X1~ X2+ X3+X4+X5+X6, data = samp)
  hccm<-car::hccm(reg, type="hc1")
  testx<-lmtest::coeftest(reg,vcov.=hccm)
  b1[i] = testx[2,4]
  b2[i] = testx[3,4]
  b3[i] = testx[4,4]
  b4[i] = testx[5,4]
  b5[i] = testx[6,4]
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
  rejectb5<-NA
  rejectb5 [ b5 < alpha] <- 1
  rejectb5 [ b5 >= alpha] <- 0

  message("Sample size is ",n)
  message("Power b1 = ", mean(rejectb1))
  message("Power b2 = ", mean(rejectb2))
  message("Power b3 = ", mean(rejectb3))
  message("Power b4 = ", mean(rejectb4))
  message("Power b5 = ", mean(rejectb5))}

# Heteroscedasticity Adjusted hc0

if (pred==2 && test=="hc0")
{
for (i in 1:nruns)
{samp <- pop2[ sample(nrow(pop2),n), ]
reg<-lm(X1~ X2+ X3, data = samp)
hccm<-car::hccm(reg, type="hc0")
testx<-lmtest::coeftest(reg,vcov.=hccm)
b1[i] = testx[2,4]
b2[i] = testx[3,4]
}
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

if (pred==3 && test=="hc0")
{
for (i in 1:nruns)
{samp <- pop2[ sample(nrow(pop2),n), ]
reg<-lm(X1~ X2+ X3+X4, data = samp)
hccm<-car::hccm(reg, type="hc0")
testx<-lmtest::coeftest(reg,vcov.=hccm)
b1[i] = testx[2,4]
b2[i] = testx[3,4]
b3[i] = testx[4,4]
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
  message("Power b3 = ", mean(rejectb3)) }

if (pred==4 && test=="hc0")
{
for (i in 1:nruns)
{samp <- pop2[ sample(nrow(pop2),n), ]
reg<-lm(X1~ X2+ X3+X4+X5, data = samp)
hccm<-car::hccm(reg, type="hc0")
testx<-lmtest::coeftest(reg,vcov.=hccm)
b1[i] = testx[2,4]
b2[i] = testx[3,4]
b3[i] = testx[4,4]
b4[i] = testx[5,4]
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
  message("Power b4 = ", mean(rejectb4))}



if (pred==5 && test=="hc0")
{
  for (i in 1:nruns)
  {samp <- pop2[ sample(nrow(pop2),n), ]
  reg<-lm(X1~ X2+ X3+X4+X5+X6, data = samp)
  hccm<-car::hccm(reg, type="hc0")
  b1[i] = testx[2,4]
  b2[i] = testx[3,4]
  b3[i] = testx[4,4]
  b4[i] = testx[5,4]
  b5[i] = testx[6,4]
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
  rejectb5<-NA
  rejectb5 [ b5 < alpha] <- 1
  rejectb5 [ b5 >= alpha] <- 0

  message("Sample size is ",n)
  message("Power b1 = ", mean(rejectb1))
  message("Power b2 = ", mean(rejectb2))
  message("Power b3 = ", mean(rejectb3))
  message("Power b4 = ", mean(rejectb4))
  message("Power b5 = ", mean(rejectb5))}

# Heteroscedasticity Adjusted hc2

if (pred==2 && test=="hc2")
{for (i in 1:nruns)
{samp <- pop2[ sample(nrow(pop2),n), ]
reg<-lm(X1~ X2+ X3, data = samp)
hccm<-car::hccm(reg, type="hc2")
testx<-lmtest::coeftest(reg,vcov.=hccm)
b1[i] = testx[2,4]
b2[i] = testx[3,4]
}
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

if (pred==3 && test=="hc2")
{for (i in 1:nruns)
{samp <- pop2[ sample(nrow(pop2),n), ]
reg<-lm(X1~ X2+ X3+X4, data = samp)
hccm<-car::hccm(reg, type="hc2")
testx<-lmtest::coeftest(reg,vcov.=hccm)
b1[i] = testx[2,4]
b2[i] = testx[3,4]
b3[i] = testx[4,4]
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
  message("Power b3 = ", mean(rejectb3)) }

if (pred==4 && test=="hc2")
{
for (i in 1:nruns)
{samp <- pop2[ sample(nrow(pop2),n), ]
reg<-lm(X1~ X2+ X3+X4+X5, data = samp)
hccm<-car::hccm(reg, type="hc2")
testx<-lmtest::coeftest(reg,vcov.=hccm)
b1[i] = testx[2,4]
b2[i] = testx[3,4]
b3[i] = testx[4,4]
b4[i] = testx[5,4]
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
  message("Power b4 = ", mean(rejectb4))}

if (pred==5 && test=="hc2")
{
  for (i in 1:nruns)
  {samp <- pop2[ sample(nrow(pop2),n), ]
  reg<-lm(X1~ X2+ X3+X4+X5+X6, data = samp)
  hccm<-car::hccm(reg, type="hc2")
  testx<-lmtest::coeftest(reg,vcov.=hccm)
  xx<-broom::tidy(testx$coefficients)
  b1[i] = testx[2,4]
  b2[i] = testx[3,4]
  b3[i] = testx[4,4]
  b4[i] = testx[5,4]
  b5[i] = testx[6,4]
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
  rejectb5<-NA
  rejectb5 [ b5 < alpha] <- 1
  rejectb5 [ b5 >= alpha] <- 0

  message("Sample size is ",n)
  message("Power b1 = ", mean(rejectb1))
  message("Power b2 = ", mean(rejectb2))
  message("Power b3 = ", mean(rejectb3))
  message("Power b4 = ", mean(rejectb4))
  message("Power b5 = ", mean(rejectb5))}

# Heteroscedasticity Adjusted hc3

if (pred==2 && test=="hc3")
{
for (i in 1:nruns)
{samp <- pop2[ sample(nrow(pop2),n), ]
reg<-lm(X1~ X2+ X3, data = samp)
hccm<-car::hccm(reg, type="hc3")
testx<-lmtest::coeftest(reg,vcov.=hccm)
b1[i] = testx[2,4]
b2[i] = testx[3,4]
}
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

if (pred==3 && test=="hc3")
{
for (i in 1:nruns)
{samp <- pop2[ sample(nrow(pop2),n), ]
reg<-lm(X1~ X2+ X3+X4, data = samp)
hccm<-car::hccm(reg, type="hc3")
testx<-lmtest::coeftest(reg,vcov.=hccm)
b1[i] = testx[2,4]
b2[i] = testx[3,4]
b3[i] = testx[4,4]
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
  message("Power b3 = ", mean(rejectb3)) }

if (pred==4 && test=="hc3")
{
for (i in 1:nruns)
{samp <- pop2[ sample(nrow(pop2),n), ]
reg<-lm(X1~ X2+ X3+X4+X5, data = samp)
hccm<-car::hccm(reg, type="hc3")
testx<-lmtest::coeftest(reg,vcov.=hccm)
b1[i] = testx[2,4]
b2[i] = testx[3,4]
b3[i] = testx[4,4]
b4[i] = testx[5,4]
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
  message("Power b4 = ", mean(rejectb4))}


if (pred==5 && test=="hc3")
{
  for (i in 1:nruns)
  {samp <- pop2[ sample(nrow(pop2),n), ]
  reg<-lm(X1~ X2+ X3+X4+X5+X6, data = samp)
  hccm<-car::hccm(reg, type="hc3")
  testx<-lmtest::coeftest(reg,vcov.=hccm)
  b1[i] = testx[2,4]
  b2[i] = testx[3,4]
  b3[i] = testx[4,4]
  b4[i] = testx[5,4]
  b5[i] = testx[6,4]
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
