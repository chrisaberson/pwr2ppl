#'Compute Power for Conditional Process Model 15 Joint Significance
#'Requires correlations between all variables as sample size.
#'This is the recommended approach for determining power
#'@param rxy Correlation between DV (y) and predictor (x)
#'@param rxm Correlation between predictor (x) and mediator (m)
#'@param rxw Correlation between predictor (x) and moderator (w)
#'@param rxwm Correlation between predictor (x) and interaction term (xw) - defaults to 0
#'@param rym Correlation between DV (y) and mediator (m)
#'@param ryw Correlation between DV (y) and moderator (w)
#'@param rwm Correlation between moderator (w) and mediator (m)
#'@param rmwy Correlation between DV (y) and interaction (xw) - defaults to 0
#'@param rmww Correlation between moderator (w) and interaction (xw) - defaults to 0
#'@param rmwm Correlation between mediator (m) and interaction (xw) - Key value
#'@param n Sample size
#'@param alpha Type I error (default is .05)
#'@param alpha rep Number of samples drawn (defaults to 5000)
#'@examples modmed15(rxw=.40, rxm=.4, rxwm=.0, rxy=.5, rwm=.45, rmww=.0, rwy=.2,
#'rmwm=.45,rmy=.30,rmwy=.0,rep=5000, alpha=.05, n=400)
#'@return Power for Model 7 Conditional Processes
#'@export
#'
#'

modmed15<-function(rxw, rxm, rxwm, rxy,
                  rwm=0, rwwm, rwxw,rxwwx, rwy,
                  rmwm=0, rmwx, rmwy,rwxy,rwmwx,
                  rmy,rxwx, rxwy,rmwwx,
                  n,alpha=.05,rep=5000)
{
set.seed(1235)
out <- MASS::mvrnorm(100000, mu = c(0,0,0,0,0,0),
                     Sigma = matrix(c(1.0,rxw,rxm,rxwm,rxwx,rxy,
                                      rxw,1.0,rwm,rwwm,rmwx,rwy,
                                      rxm,rwm,1.0,rmwm,rwxw,rmy,
                                      rxwm,rwwm,rmwm,1.0,rwmwx,rmwy,
                                      rxwx,rwxw,rmwx,rwmwx,1.0,rwxy,
                                      rxy,rwy,rmy,rmwy,rwxy,1.0),
                                    ncol = 6),
                     empirical = TRUE)
out<-as.data.frame(out)

out<-dplyr::rename(out, x = V1, w= V2, m= V3, xw = V4, mw = V5, y=V6)


nruns = rep
a = numeric(nruns)
b = numeric(nruns)
pa1<-NA
pb1<-NA
pb2<-NA

for (i in 1:nruns)
{
  samp <- out[ sample(nrow(out), n), ]
  test_a <- stats::lm(m ~ x , data = samp)
  test_b <- stats::lm(y ~ x + m + w + xw + mw, data = samp)
  apath<-summary(test_a)
  bpath<-summary(test_b)
  pa1[i]<-apath$coefficients[2,4]
  pb1[i]<-bpath$coefficients[5,4]
  pb2[i]<-bpath$coefficients[6,4]
  power = data.frame(pa1=pa1,pb1=pb1,pb2=pb2)

  power$jsa1["JointSiga1"]<-NA
  power$jsa1[pa1 < alpha] <- 1
  power$jsa1[pa1 >= alpha] <- 0
  power$jsb1["JointSigb1"]<-NA
  power$jsb1[pb1 < alpha] <- 1
  power$jsb1[pb1 >= alpha] <- 0
  power$jsb2["JointSigb2"]<-NA
  power$jsb2[pb2 < alpha] <- 1
  power$jsb2[pb2 >= alpha] <- 0

  power$joint_med<-power$jsa1*power$jsb1
  JSa<-mean(power$joint_med)
  power$joint_dv<-power$jsa1*power$jsb2
  JSb<-mean(power$joint_dv)
  }
message("Sample size is ",n)
message("Power for Conditional Indirect Effect (On Mediator) ", JSa)
message("Power for Conditional Indirect Effect (On DV) ", JSb)
}
