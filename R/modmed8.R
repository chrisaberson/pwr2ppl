#'Compute Power for Power for Model 8 Conditional Processes Using Joint Significance
#'Requires correlations between all variables as sample size.
#'This is the recommended approach for determining power
#'@param rxy Correlation between DV (y) and predictor (x)
#'@param rxm Correlation between predictor (x) and mediator (m)
#'@param rxw Correlation between predictor (x) and moderator (w)
#'@param rxxw Correlation between predictor (x) and interaction term (xw) - defaults to 0
#'@param rmy Correlation between DV (y) and mediator (m)
#'@param rwy Correlation between DV (y) and moderator (w)
#'@param rwm Correlation between moderator (w) and mediator (m)
#'@param rxwy Correlation between DV (y) and interaction (xw) - defaults to 0
#'@param rwxw Correlation between moderator (w) and interaction (xw) - defaults to 0
#'@param rxwm Correlation between mediator (m) and interaction (xw) - Key value
#'@param n Sample size
#'@param alpha Type I error (default is .05)
#'@param rep Number of samples drawn (defaults to 5000)
#'@examples \donttest{modmed8(rxw<-.21, rxm<-.31, rxxw=0, rxy=.32,rwm=.40,
#'rmy=.19,rwy=.22,rwxw=.23,rxwm=.24,rxwy=.18,alpha=.05,rep=1000,n=400)}
#'@return Power for Model 8 Conditional Processes
#'@export

modmed8<-function(rxw, rxm, rxxw, rxy,
                   rwm=0, rwy,
                   rxwm, rxwy,rwxw,
                   rmy, n,alpha=.05,rep=5000)
{
V1<-NA;V2<-NA;V3<-NA;V4<-NA;V5<-NA
set.seed(1235)
out <- MASS::mvrnorm(100000, mu = c(0,0,0,0,0),
                     Sigma = matrix(c(1.0,rxw,rxm,rxxw,rxy,
                                      rxw,1.0,rwm,rwxw,rwy,
                                      rxm,rwm,1.0,rxwm,rmy,
                                      rxxw,rwxw,rxwm,1.0,rxwy,
                                      rxy,rwy,rmy,rxwy,1.0),
                                    ncol = 5),
                     empirical = TRUE)
out<-as.data.frame(out)

out<-dplyr::rename(out, x = V1, w= V2, m= V3, xw = V4, y=V5)


nruns = rep
a = numeric(nruns)
b = numeric(nruns)
pa1<-NA
pa2<-NA
pb1<-NA

for (i in 1:nruns)
{
  samp <- out[ sample(nrow(out), n), ]
  test_a <- stats::lm(m ~ x + w + xw , data = samp)
  test_b <- stats::lm(y ~ x + m + w + xw, data = samp)

  apath<-summary(test_a)
  bpath<-summary(test_b)
  pa1[i]<-apath$coefficients[2,4]
  pa2[i]<-bpath$coefficients[4,4]
  pb1[i]<-bpath$coefficients[3,4]
  power = data.frame(pa1=pa1,pa2=pa2,pb1=pb1)

  power$jsa1["JoxwSiga1"]<-NA
  power$jsa1[pa1 < alpha] <- 1
  power$jsa1[pa1 >= alpha] <- 0
  power$jsa2["JoxwSiga2"]<-NA
  power$jsa2[pa2 < alpha] <- 1
  power$jsa2[pa2 >= alpha] <- 0
  power$jsb1["JoxwSigb1"]<-NA
  power$jsb1[pb1 < alpha] <- 1
  power$jsb1[pb1 >= alpha] <- 0

  power$joxw_med<-power$jsa1*power$jsb1
  JSa<-mean(power$joxw_med)
  power$joxw_dv<-power$jsa2*power$jsb1
  JSb<-mean(power$joxw_dv)
}

message("Sample size is ",n)
message("Power for Conditional Indirect Effect (On Mediator) ", JSa)
message("Power for Conditional Indirect Effect (On DV) ", JSb)
}
