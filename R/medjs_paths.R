#'Compute Power for Mediated (Indirect) Effects Using Joint Significance
#'Requires paths for all effects (and if 2 mediators, correlation)
#'Standard deviations/variances set to 1.0 so paths are technically standardized
#'@param a1 path between predictor and first mediator
#'@param a2 path between predictor and first mediator
#'@param b1 Path between first mediator and dependent variable
#'@param b2 Path between first mediator and dependent variable
#'@param cprime Path between predictor and dependent variable
#'@param rm1m2 Correlation first mediator (m1) and second mediator (m2)
#'@param n Sample size
#'@param mvars Number of Mediators
#'@param alpha Type I error (default is .05)
#'@param rep number of repetitions (1000 is default)
#'@examples
#'\donttest{medjs_paths(a1=.25, b1=-.5,cprime=.2,mvars=1, n=150)}
#'\donttest{medjs_paths(a1=.25, a2=.1, b1=-.5,b2=-.2,cprime=.2,mvars=1, n=150)}
#'@return Power for Mediated (Indirect) Effects using Paths Coefficients
#'@export
#'
#'

medjs_paths<-function(a1, a2=NULL,b1,b2=NULL,rm1m2=NULL,cprime,n,
                              alpha=.05,mvars,rep=1000)
{
  V1<-V2<-V3<-V4<-V5<-V6<-V7<-NA
  if(mvars==1)  {
   {rx1m1<-a1; rym1<-b1+(a1*cprime); rx1y<-a1*b1+cprime}
      out <- MASS::mvrnorm(100000, mu = c(0,0,0),
                   Sigma = matrix(c(1.0,rx1m1,rx1y,
                                    rx1m1,1.0,rym1,
                                    rx1y,rym1,1.0),
                                    ncol = 3),
                   empirical = TRUE)
    pop<-as.data.frame(out)
    pop<-dplyr::rename(pop, x = V1, m1 = V2, y = V3)
    set.seed(1234)
    nruns = rep
    a = numeric(nruns)
    b = numeric(nruns)
    pa<-NA
    pb<-NA
    for (i in 1:nruns)
    {
    samp <- pop[ sample(nrow(pop), n), ]
    test_a <- stats::lm(m1 ~ x, data = samp)
    test_b <- stats::lm(y ~ x + m1, data = samp)
    apath<-summary(test_a)
    bpath<-summary(test_b)
    pa[i]<-apath$coefficients[2,4]
    pb[i]<-bpath$coefficients[3,4]
    power = data.frame(pa=pa, pb=pb)
    power$jsa["JointSiga"]<-NA
    power$jsa[pa < alpha] <- 1
    power$jsa[pa >= alpha] <- 0
    power$jsb["JointSigb"]<-NA
    power$jsb[pb < alpha] <- 1
    power$jsb[pb >= alpha] <- 0
    power$joint<-power$jsa*power$jsb
    JS<-mean(power$joint)}
    message("Sample size is ",n)
    message("Power is ", JS)
  }
  if(mvars==2){
   {rx1m1<-a1; rx1m2<-a2;
    rx1y<-cprime + a1*b1 + a2*b2;rym1<-a1*cprime + b1 + b2*rx1m1;
    rym2<-a2*cprime+b2+b1* rx1m1}
    out <- MASS::mvrnorm(100000, mu = c(0,0,0,0),
                         Sigma = matrix(c(1.0,rx1m1,rx1m2,rx1y,
                                          rx1m1,1.0,rm1m2,rym1,
                                          rx1m2,rm1m2,1,rym2,
                                          rx1y,rym1,rym2,1.0),
                                        ncol = 4),
                         empirical = TRUE)
    pop<-as.data.frame(out)
    pop<-dplyr::rename(pop, x = V1, m1 = V2, m2=V3,y = V4)
    set.seed(1234)
    nruns = rep
    a = numeric(nruns)
    b = numeric(nruns)
    pa1<-NA
    pa2<-NA
    pb1<-NA
    pb2<-NA
    for (i in 1:nruns)
    {
      samp <- pop[ sample(nrow(pop), n), ]
      test_a1 <- stats::lm(m1 ~ x , data = samp)
      test_a2 <- stats::lm(m2 ~ x , data = samp)
      test_b <- stats::lm(y ~ x + m1 + m2, data = samp)
      apath1<-summary(test_a1)
      apath2<-summary(test_a2)
      bpath<-summary(test_b)
      pa1[i]<-apath1$coefficients[2,4]
      pa2[i]<-apath2$coefficients[2,4]
      pb1[i]<-bpath$coefficients[3,4]
      pb2[i]<-bpath$coefficients[4,4]
      power = data.frame(pa1=pa1,pa2=pa2, pb1=pb1,pb2=pb2)
      power$jsa1["JointSiga1"]<-NA
      power$jsa1[pa1 < alpha] <- 1
      power$jsa1[pa1 >= alpha] <- 0
      power$jsa2["JointSiga2"]<-NA
      power$jsa2[pa2 < alpha] <- 1
      power$jsa2[pa2 >= alpha] <- 0
      power$jsb1["JointSigb1"]<-NA
      power$jsb1[pb1 < alpha] <- 1
      power$jsb1[pb1 >= alpha] <- 0
      power$jsb2["JointSigb2"]<-NA
      power$jsb2[pb2 < alpha] <- 1
      power$jsb2[pb2 >= alpha] <- 0
      power$jointb1<-power$jsa1*power$jsb1
      power$jointb2<-power$jsa2*power$jsb2
      JSm1<-mean(power$jointb1)
      JSm2<-mean(power$jointb2)}
    message("Sample size is ",n)
    message("Power for M1 ", JSm1)
    message("Power for M2 ", JSm2)
  }
  }


