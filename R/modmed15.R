#'Compute Power for Conditional Process Model 15 Joint Significance
#'Requires correlations between all variables as sample size.
#'This is the recommended approach for determining power
#'@param rxy Correlation between DV (y) and predictor (x)
#'@param rxm Correlation between predictor (x) and mediator (m)
#'@param rxw Correlation between predictor (x) and moderator (w)
#'@param rmy Correlation between DV (y) and mediator (m)
#'@param rwy Correlation between DV (y) and moderator (w)
#'@param rwm Correlation between moderator (w) and mediator (m)
#'@param rmww Correlation between mediator (m) and interaction (mw) - defaults to 0
#'@param rmxw Correlation between mediator (m) and interaction (xw) - defaults to 0
#'@param rxwy Correlation between DV (y) and interaction (xw) - defaults to 0
#'@param rxwx Correlation between moderator (w) and interaction (xw) - defaults to 0
#'@param rmwxw Correlation between inteaction (xw) and interaction (mw) - Key value
#'@param rxmw Correlation between predictor (x) and interaction (mw) - Key value
#'@param rmwm Correlation between mediator (m) and interaction (xmw) - Key value
#'@param rmwy Correlation between dv (y) and interaction (mw) - Key value
#'@param rwxw Correlation between moderator (w) and interaction (xw) - Key value
#'@param rxwmw Correlation between interaction (mw) and interaction (mw) - Key value
#'@param n Sample size
#'@param alpha Type I error (default is .05)
#'@param rep Number of samples drawn (defaults to 5000)
#'@examples \donttest{modmed15(rxw=.40, rxm=.42, rxy=.5, rwm=.45, rmxw=.0,rmww=.01, rwy=.2,
#'rmwm=.46,rwxw=.21,rxwy=.31,rmy=.30,rxwx=.1, rmwy=.02,rxmw=.21,rmwxw=.22,rep=5000, alpha=.05, n=400)}
#'@return Power for Model 15 Conditional Processes
#'@export
#'
#'

modmed15<-function(rxw, rxm, rxmw, rxy,
                  rwm=0, rmww, rwxw, rwy,
                  rmwm=0, rmxw, rmwy,rxwy,rxwmw,
                  rmy,rxwx, rmwxw,
                  n,alpha=.05,rep=5000)
{

  V1<-NA;V2<-NA;V3<-NA;V4<-NA;V5<-NA;V6<-NA

set.seed(1235)
out <- MASS::mvrnorm(100000, mu = c(0,0,0,0,0,0),
                     Sigma = matrix(c(1.0,rxw,rxm,rxmw,rxwx,rxy,
                                      rxw,1.0,rwm,rmww,rmxw,rwy,
                                      rxm,rwm,1.0,rmwm,rwxw,rmy,
                                      rxmw,rmww,rmwm,1.0,rmwxw,rmwy,
                                      rxwx,rwxw,rmxw,rmwxw,1.0,rxwy,
                                      rxy,rwy,rmy,rmwy,rxwy,1.0),
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
