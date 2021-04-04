#'Compute Power for Mediated (Indirect) Effects Using Joint Significance
#'Requires correlations between all variables as sample size.
#'This is the recommended approach for determining power
#'@param rx1y Correlation between DV (y) and first predictor (x1)
#'@param rx1m1 Correlation between first predictor (x1) and first mediator (m1)
#'@param rx1m2 Correlation between first predictor (x1) and second mediator (m2)
#'@param rx1m3 Correlation between first predictor (x1) and third mediator (m3)
#'@param rx1m4 Correlation between first predictor (x1) and fourth mediator (m4)
#'@param rx1x2 Correlation between first predictor (x1) and second predictor (x2)
#'@param rx2y Correlation between DV (y) and second predictor (x2)
#'@param rx2m1 Correlation between second predictor (x2) and first mediator (m1)
#'@param rx2m2 Correlation between second predictor (x2) and second mediator (m2)
#'@param rx2m3 Correlation between second predictor (x2) and third mediator (m3)
#'@param rx2m4 Correlation between second predictor (x2) and fourth mediator (m4)
#'@param rym1 Correlation between DV (y) and first mediator (m1)
#'@param rym2 Correlation between DV (y) and second mediator (m2)
#'@param rym3 Correlation DV (y) and third mediator (m3)
#'@param rym4 Correlation DV (y) and fourth mediator (m4)
#'@param rm1m2 Correlation first mediator (m1) and second mediator (m2)
#'@param rm1m3 Correlation first mediator (m1) and third mediator (m3)
#'@param rm1m4 Correlation first mediator (m1) and fourth mediator (m4)
#'@param rm2m3 Correlation second mediator (m2) and third mediator (m3)
#'@param rm2m4 Correlation second mediator (m2) and fourth mediator (m4)
#'@param rm3m4 Correlation third mediator (m3) and fourth mediator (m4)
#'@param n Sample size
#'@param mvars Number of Mediators
#'@param alpha Type I error (default is .05)
#'@param rep number of repetitions (1000 is default)
#'@param pred number of predictors (default is one)
#'@examples \donttest{medjs(rx1m1=.3, rx1m2=.3, rx1m3=.25, rx1y=-.35, rym1=-.5,rym2=-.5, rym3 = -.5,
#'rm1m2=.7, rm1m3=.4,rm2m3=.4, mvars=3, n=150)}
#'@return Power for Mediated (Indirect) Effects
#'@export
#'
#'

medjs<-function(rx1x2=NULL, rx1m1, rx1m2=NULL, rx1m3=NULL, rx1m4=NULL,rx1y,
                              rx2m1=NULL, rx2m2=NULL, rx2m3=NULL, rx2m4=NULL,rx2y,
                              rym1, rym2=NULL,rym3=NULL, rym4=NULL, rm1m2=NULL,rm1m3=NULL,
                              rm1m4=NULL, rm2m3=NULL, rm2m4=NULL, rm3m4=NULL,n,
                              alpha=.05,mvars,rep=1000, pred=1)
{
  V1<-V2<-V3<-V4<-V5<-V6<-V7<-NA
  if(mvars==1 & pred==1)  {
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
  if(mvars==2 & pred==1){
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
  if(mvars==3 & pred==1){
    out <- MASS::mvrnorm(100000, mu = c(0,0,0,0,0),
                         Sigma = matrix(c(1.0,rx1m1,rx1m2,rx1m3,rx1y,
                                          rx1m1,1.0,rm1m2,rm1m3,rym1,
                                          rx1m2,rm1m2,1,rm2m3,rym2,
                                          rx1m3,rm1m3,rm2m3,1,rym3,
                                          rx1y,rym1,rym2,rym3,1.0),
                                        ncol = 5),
                         empirical = TRUE)
    pop<-as.data.frame(out)
    pop<-dplyr::rename(pop, x = V1, m1 = V2, m2=V3,m3= V4,y = V5)
    set.seed(1234)
    nruns = rep
    a = numeric(nruns)
    b = numeric(nruns)
    pa1<-NA
    pa2<-NA
    pa3<-NA
    pb1<-NA
    pb2<-NA
    pb3<-NA
    for (i in 1:nruns)
    {
      samp <- pop[ sample(nrow(pop), n), ]
      test_a1 <- stats::lm(m1 ~ x, data = samp)
      test_a2 <- stats::lm(m2 ~ x, data = samp)
      test_a3 <- stats::lm(m3 ~ x, data = samp)
      test_b <- stats::lm(y ~ x + m1 + m2 + m3, data = samp)
      apath1<-summary(test_a1)
      apath2<-summary(test_a2)
      apath3<-summary(test_a3)
      bpath<-summary(test_b)
      pa1[i]<-apath1$coefficients[2,4]
      pa2[i]<-apath2$coefficients[2,4]
      pa3[i]<-apath3$coefficients[2,4]
      pb1[i]<-bpath$coefficients[3,4]
      pb2[i]<-bpath$coefficients[4,4]
      pb3[i]<-bpath$coefficients[5,4]
      power = data.frame(pa1=pa1,pa2=pa2,pa3=pa3,pb1=pb1,pb2=pb2,pb3=pb3)
      power$jsa1["JointSiga1"]<-NA
      power$jsa1[pa1 < alpha] <- 1
      power$jsa1[pa1 >= alpha] <- 0
      power$jsa2["JointSiga2"]<-NA
      power$jsa2[pa2 < alpha] <- 1
      power$jsa2[pa2 >= alpha] <- 0
      power$jsa3["JointSiga3"]<-NA
      power$jsa3[pa3 < alpha] <- 1
      power$jsa3[pa3 >= alpha] <- 0

      power$jsb1["JointSigb1"]<-NA
      power$jsb1[pb1 < alpha] <- 1
      power$jsb1[pb1 >= alpha] <- 0
      power$jsb2["JointSigb2"]<-NA
      power$jsb2[pb2 < alpha] <- 1
      power$jsb2[pb2 >= alpha] <- 0
      power$jsb3["JointSigb3"]<-NA
      power$jsb3[pb3 < alpha] <- 1
      power$jsb3[pb3 >= alpha] <- 0
      power$jointb1<-power$jsa1*power$jsb1
      power$jointb2<-power$jsa2*power$jsb2
      power$jointb3<-power$jsa3*power$jsb3
      JSm1<-mean(power$jointb1)
      JSm2<-mean(power$jointb2)
      JSm3<-mean(power$jointb3)}
    message("Sample size is ",n)
    message("Power for M1 ", JSm1)
    message("Power for M2 ", JSm2)
    message("Power for M3 ", JSm3)
  }

  if(mvars==4 & pred==1){
    out <- MASS::mvrnorm(100000, mu = c(0,0,0,0,0,0),
                         Sigma = matrix(c(1.0,rx1m1,rx1m2,rx1m3,rx1m4,rx1y,
                                          rx1m1,1.0,rm1m2,rm1m3,rm1m4,rym1,
                                          rx1m2,rm1m2,1.0,rm2m3,rm2m4,rym2,
                                          rx1m3,rm1m3,rm2m3,1.0,rm3m4,rym3,
                                          rx1m4,rm1m4,rm2m4,rm3m4,1.0,rym4,
                                          rx1y,rym1,rym2,rym3,rym4,1.0),
                                        ncol = 6),
                         empirical = TRUE)
    pop<-as.data.frame(out)
    pop<-dplyr::rename(pop, x = V1, m1 = V2, m2=V3,m3= V4,m4=V5,y = V6)
    set.seed(1234)
    nruns = rep
    a = numeric(nruns)
    b = numeric(nruns)
    pa1<-NA
    pa2<-NA
    pa3<-NA
    pa4<-NA
    pb1<-NA
    pb2<-NA
    pb3<-NA
    pb4<-NA
    for (i in 1:nruns)
    {
      samp <- pop[ sample(nrow(pop), n), ]
      test_a1 <- stats::lm(m1 ~ x, data = samp)
      test_a2 <- stats::lm(m2 ~ x, data = samp)
      test_a3 <- stats::lm(m3 ~ x, data = samp)
      test_a4 <- stats::lm(m4 ~ x, data = samp)
      test_b <- stats::lm(y ~ x + m1 + m2 + m3 + m4, data = samp)
      apath1<-summary(test_a1)
      apath2<-summary(test_a2)
      apath3<-summary(test_a3)
      apath4<-summary(test_a4)
      bpath<-summary(test_b)
      pa1[i]<-apath1$coefficients[2,4]
      pa2[i]<-apath2$coefficients[2,4]
      pa3[i]<-apath3$coefficients[2,4]
      pa4[i]<-apath4$coefficients[2,4]
      pb1[i]<-bpath$coefficients[3,4]
      pb2[i]<-bpath$coefficients[4,4]
      pb3[i]<-bpath$coefficients[5,4]
      pb4[i]<-bpath$coefficients[6,4]
      power = data.frame(pa1=pa1,pa2=pa2,pa3=pa3,pa4=pa4,pb1=pb1,pb2=pb2,pb3=pb3,pb4=pb4)
      power$jsa1["JointSiga1"]<-NA
      power$jsa1[pa1 < alpha] <- 1
      power$jsa1[pa1 >= alpha] <- 0
      power$jsa2["JointSiga2"]<-NA
      power$jsa2[pa2 < alpha] <- 1
      power$jsa2[pa2 >= alpha] <- 0
      power$jsa3["JointSiga3"]<-NA
      power$jsa3[pa3 < alpha] <- 1
      power$jsa3[pa3 >= alpha] <- 0
      power$jsa4["JointSiga4"]<-NA
      power$jsa4[pa3 < alpha] <- 1
      power$jsa4[pa3 >= alpha] <- 0

      power$jsb1["JointSigb1"]<-NA
      power$jsb1[pb1 < alpha] <- 1
      power$jsb1[pb1 >= alpha] <- 0
      power$jsb2["JointSigb2"]<-NA
      power$jsb2[pb2 < alpha] <- 1
      power$jsb2[pb2 >= alpha] <- 0
      power$jsb3["JointSigb3"]<-NA
      power$jsb3[pb3 < alpha] <- 1
      power$jsb3[pb3 >= alpha] <- 0
      power$jsb4["JointSigb4"]<-NA
      power$jsb4[pb4 < alpha] <- 1
      power$jsb4[pb4 >= alpha] <- 0

      power$jointb1<-power$jsa1*power$jsb1
      power$jointb2<-power$jsa2*power$jsb2
      power$jointb3<-power$jsa3*power$jsb3
      power$jointb4<-power$jsa4*power$jsb4
      JSm1<-mean(power$jointb1)
      JSm2<-mean(power$jointb2)
      JSm3<-mean(power$jointb3)
      JSm4<-mean(power$jointb4)}
    message("Power using Joint Significance Test")
    message("Sample size is ",n)
    message("Power for M1 ", JSm1)
    message("Power for M2 ", JSm2)
    message("Power for M3 ", JSm3)
    message("Power for M4 ", JSm4)
  }


  if(mvars==1 & pred==2)  {
    out <- MASS::mvrnorm(100000, mu = c(0,0,0,0),
                         Sigma = matrix(c(1.0,rx1x2,rx1m1,rx1y,
                                          rx1x2, 1.0,rx2m1,rx2y,
                                          rx1m1, rx2m1, 1.0,rym1,
                                          rx1y,rx2y,rym1,1.0),
                                        ncol = 4),
                         empirical = TRUE)
    pop<-as.data.frame(out)
    pop<-dplyr::rename(pop, x1 = V1, x2 = V2,m1 = V3, y = V4,)
    set.seed(1234)
    nruns = rep
    a = numeric(nruns)
    b = numeric(nruns)
    pa1<-NA
    pa2<-NA
    pb<-NA
    for (i in 1:nruns)
    {
      samp <- pop[ sample(nrow(pop), n), ]
      test_a <- stats::lm(m1 ~ x1 + x2, data = samp)
      test_b <- stats::lm(y ~ x1 + x2+m1, data = samp)
      apath<-summary(test_a)
      bpath<-summary(test_b)
      pa1[i]<-apath$coefficients[2,4]
      pa2[i]<-apath$coefficients[3,4]
      pb[i]<-bpath$coefficients[4,4]
      power = data.frame(pa1=pa1,pa2=pa2,pb=pb)
      power$jsa1["JointSiga1"]<-NA
      power$jsa1[pa1 < alpha] <- 1
      power$jsa1[pa1 >= alpha] <- 0
      power$jsa2["JointSiga2"]<-NA
      power$jsa2[pa2 < alpha] <- 1
      power$jsa2[pa2 >= alpha] <- 0
      power$jsb["JointSigb"]<-NA
      power$jsb[pb < alpha] <- 1
      power$jsb[pb >= alpha] <- 0
      power$joint1<-power$jsa1*power$jsb
      JS1<-mean(power$joint1)
      power$joint2<-power$jsa2*power$jsb
      JS2<-mean(power$joint2)
    }
    message("Sample size is ",n)
    message("Power for predictor 1, mediator 1 is ", JS1)
    message("Power for predictor 2, mediator 1 is ", JS2)
  }
  if(mvars==2 & pred==2){
    out <- MASS::mvrnorm(100000, mu = c(0,0,0,0,0),
                         Sigma = matrix(c(1.0,rx1x2,rx1m1,rx1m2,rx1y,
                                          rx1x2,1.0,rx2m1,rx2m2,rx2y,
                                          rx1m1,rx2m1,1.0,rm1m2,rym1,
                                          rx1m2,rx2m2,rm1m2,1,rym2,
                                          rx1y,rx2y,rym1,rym2,1.0),
                                        ncol = 5),
                         empirical = TRUE)
    pop<-as.data.frame(out)
    pop<-dplyr::rename(pop, x1 = V1, x2 = V2, m1 = V3, m2=V4,y = V5)
    set.seed(1234)
    nruns = rep
    a = numeric(nruns)
    b = numeric(nruns)
    pa1m1<-NA
    pa1m2<-NA
    pa2m1<-NA
    pa2m2<-NA
    pb1<-NA
    pb2<-NA
    for (i in 1:nruns)
    {
      samp <- pop[ sample(nrow(pop), n), ]
      test_a1 <- stats::lm(m1 ~ x1+x2, data = samp)
      test_a2 <- stats::lm(m2 ~ x1+x2, data = samp)
      test_b <- stats::lm(y ~ x1 + x2 + m1 + m2, data = samp)
      apath1<-summary(test_a1)
      apath2<-summary(test_a2)
      bpath<-summary(test_b)
      pa1m1[i]<-apath1$coefficients[2,4]
      pa1m2[i]<-apath1$coefficients[3,4]
      pa2m1[i]<-apath2$coefficients[2,4]
      pa2m2[i]<-apath2$coefficients[3,4]
      pb1[i]<-bpath$coefficients[4,4]
      pb2[i]<-bpath$coefficients[5,4]
      power = data.frame(pa1m1=pa1m1,pa1m2=pa1m2,pa2m1=pa2m1,pa2m2=pa2m2, pb1=pb1,pb2=pb2)
      power$jsa1m1["JointSiga1m1"]<-NA
      power$jsa1m1[pa1m1 < alpha] <- 1
      power$jsa1m1[pa1m1 >= alpha] <- 0
      power$jsa2m1["JointSiga2m1"]<-NA
      power$jsa2m1[pa2m1 < alpha] <- 1
      power$jsa2m1[pa2m1 >= alpha] <- 0
      power$jsa1m2["JointSiga1m2"]<-NA
      power$jsa1m2[pa1m2 < alpha] <- 1
      power$jsa1m2[pa1m2 >= alpha] <- 0
      power$jsa2m2["JointSiga2m2"]<-NA
      power$jsa2m2[pa2m2 < alpha] <- 1
      power$jsa2m2[pa2m2 >= alpha] <- 0
      power$jsb1["JointSigb1"]<-NA
      power$jsb1[pb1 < alpha] <- 1
      power$jsb1[pb1 >= alpha] <- 0
      power$jsb2["JointSigb2"]<-NA
      power$jsb2[pb2 < alpha] <- 1
      power$jsb2[pb2 >= alpha] <- 0
      power$jointa1b1m1<-power$jsa1m1*power$jsb1
      power$jointa2b1m1<-power$jsa2m1*power$jsb1
      power$jointa1b2m2<-power$jsa1m2*power$jsb2
      power$jointa2b2m2<-power$jsa2m2*power$jsb2
      JSx1m1<-mean(power$jointa1b1m1)
      JSx2m1<-mean(power$jointa2b1m1)
      JSx1m2<-mean(power$jointa1b2m2)
      JSx2m2<-mean(power$jointa2b2m2)
    }
    message("Sample size is ",n)
    message("Power for X1-M1 ", JSx1m1)
    message("Power for X2-M1 ", JSx2m1)
    message("Power for X1-M2 ", JSx1m2)
    message("Power for X2-M2 ", JSx2m2)
  }
  if(mvars==3 & pred==2){
    out <- MASS::mvrnorm(100000, mu = c(0,0,0,0,0,0),
                         Sigma = matrix(c(1.0,rx1x2,rx1m1,rx1m2,rx1m3,rx1y,
                                          rx1x2,1.0,rx2m1,rx2m2,rx2m3,rx2y,
                                          rx1m1,rx2m1,1.0,rm1m2,rm1m3,rym1,
                                          rx1m2,rx2m2,rm1m2,1.0,rm2m3,rym2,
                                          rx1m3,rx2m3,rm1m3,rm2m3,1.0,rym3,
                                          rx1y,rx2y,rym1,rym2,rym3,1.0),
                                        ncol = 6),
                         empirical = TRUE)
    pop<-as.data.frame(out)
    pop<-dplyr::rename(pop, x1 = V1, x2 = V2, m1 = V3, m2=V4,m3=V5,y = V6)
    set.seed(1234)
    nruns = rep
    a = numeric(nruns)
    b = numeric(nruns)
    pa1m1<-NA
    pa2m1<-NA
    pa1m2<-NA
    pa2m2<-NA
    pa1m3<-NA
    pa2m3<-NA
    pb1<-NA
    pb2<-NA
    pb3<-NA
    for (i in 1:nruns)
    {
      samp <- pop[ sample(nrow(pop), n), ]
      test_a1 <- stats::lm(m1 ~ x1+x2, data = samp)
      test_a2 <- stats::lm(m2 ~ x1+x2, data = samp)
      test_a3 <- stats::lm(m3 ~ x1+x2, data = samp)
      test_b <- stats::lm(y ~ x1 + x2 + m1 + m2 + m3, data = samp)
      apath1<-summary(test_a1)
      apath2<-summary(test_a2)
      apath3<-summary(test_a3)
      bpath<-summary(test_b)
      pa1m1[i]<-apath1$coefficients[2,4]
      pa1m2[i]<-apath1$coefficients[3,4]
      pa2m1[i]<-apath2$coefficients[2,4]
      pa2m2[i]<-apath2$coefficients[3,4]
      pa1m3[i]<-apath3$coefficients[2,4]
      pa2m3[i]<-apath3$coefficients[3,4]
      pb1[i]<-bpath$coefficients[4,4]
      pb2[i]<-bpath$coefficients[5,4]
      pb3[i]<-bpath$coefficients[6,4]
      power = data.frame(pa1m1=pa1m1,pa1m2=pa1m2,pa1m3=pa1m3,pa2m1=pa2m1,pa2m2=pa2m2,pa2m3=pa2m3,pb1=pb1,pb2=pb2,pb3=pb3)
      power$jsa1m1["JointSiga1m1"]<-NA
      power$jsa1m1[pa1m1 < alpha] <- 1
      power$jsa1m1[pa1m1 >= alpha] <- 0
      power$jsa2m1["JointSiga2m1"]<-NA
      power$jsa2m1[pa2m1 < alpha] <- 1
      power$jsa2m1[pa2m1 >= alpha] <- 0
      power$jsa1m2["JointSiga1m2"]<-NA
      power$jsa1m2[pa1m2 < alpha] <- 1
      power$jsa1m2[pa1m2 >= alpha] <- 0
      power$jsa2m2["JointSiga2m2"]<-NA
      power$jsa2m2[pa2m2 < alpha] <- 1
      power$jsa2m2[pa2m2 >= alpha] <- 0
      power$jsa1m3["JointSiga1m3"]<-NA
      power$jsa1m3[pa1m3 < alpha] <- 1
      power$jsa1m3[pa1m3 >= alpha] <- 0
      power$jsa2m3["JointSiga2m2"]<-NA
      power$jsa2m3[pa2m3 < alpha] <- 1
      power$jsa2m3[pa2m3 >= alpha] <- 0
      power$jsb1["JointSigb1"]<-NA
      power$jsb1[pb1 < alpha] <- 1
      power$jsb1[pb1 >= alpha] <- 0
      power$jsb2["JointSigb2"]<-NA
      power$jsb2[pb2 < alpha] <- 1
      power$jsb2[pb2 >= alpha] <- 0
      power$jsb3["JointSigb3"]<-NA
      power$jsb3[pb3 < alpha] <- 1
      power$jsb3[pb3 >= alpha] <- 0
      power$jointa1b1m1<-power$jsa1m1*power$jsb1
      power$jointa2b1m1<-power$jsa2m1*power$jsb1
      power$jointa1b2m2<-power$jsa1m2*power$jsb2
      power$jointa2b2m2<-power$jsa2m2*power$jsb2
      power$jointa1b3m3<-power$jsa1m3*power$jsb3
      power$jointa2b3m3<-power$jsa2m3*power$jsb3
      JSx1m1<-mean(power$jointa1b1)
      JSx2m1<-mean(power$jointa2b1)
      JSx1m2<-mean(power$jointa1b2)
      JSx2m2<-mean(power$jointa2b2)
      JSx1m3<-mean(power$jointa1b3)
      JSx2m3<-mean(power$jointa2b3)
    }
    message("Sample size is ",n)
    message("Power for X1-M1 ", JSx1m1)
    message("Power for X2-M1 ", JSx2m1)
    message("Power for X1-M2 ", JSx1m2)
    message("Power for X2-M2 ", JSx2m2)
    message("Power for X1-M3 ", JSx1m3)
    message("Power for X2-M3 ", JSx2m3)
  }

  if(mvars==4 & pred==2){
    out <- MASS::mvrnorm(100000, mu = c(0,0,0,0,0,0,0),
                         Sigma = matrix(c(1.0,rx1x2,rx1m1,rx1m2,rx1m3,rx1m4,rx1y,
                                          rx1x2,1.0,rx2m1,rx2m2,rx2m3,rx2m4,rx2y,
                                          rx1m1,rx2m1,1.0,rm1m2,rm1m3,rm1m4,rym1,
                                          rx1m2,rx2m2,rm1m2,1.0,rm2m3,rm2m4,rym2,
                                          rx1m3,rx2m3,rm1m3,rm2m3,1.0,rm3m4,rym3,
                                          rx1m4,rx2m4,rm1m4,rm2m4,rm3m4,1.0,rym4,
                                          rx1y,rx2y,rym1,rym2,rym3,rym4,1.0),
                                        ncol = 7),
                         empirical = TRUE)
    pop<-as.data.frame(out)
    pop<-dplyr::rename(pop, x1 = V1, x2 = V2, m1 = V3, m2=V4,m3=V5,m4=V6,y = V7)
    set.seed(1234)
    nruns = rep
    a = numeric(nruns)
    b = numeric(nruns)
    pa1m1<-NA
    pa2m1<-NA
    pa1m2<-NA
    pa2m2<-NA
    pa1m3<-NA
    pa2m3<-NA
    pa1m4<-NA
    pa2m4<-NA
    pb1<-NA
    pb2<-NA
    pb3<-NA
    pb4<-NA
    for (i in 1:nruns)
    {
      samp <- pop[ sample(nrow(pop), n), ]
      test_a1 <- stats::lm(m1 ~ x1+x2, data = samp)
      test_a2 <- stats::lm(m2 ~ x1+x2, data = samp)
      test_a3 <- stats::lm(m3 ~ x1+x2, data = samp)
      test_a4 <- stats::lm(m4 ~ x1+x2, data = samp)
      test_b <- stats::lm(y ~ x1 + x2 + m1 + m2 + m3+m4, data = samp)
      apath1<-summary(test_a1)
      apath2<-summary(test_a2)
      apath3<-summary(test_a3)
      apath4<-summary(test_a4)
      bpath<-summary(test_b)
      pa1m1[i]<-apath1$coefficients[2,4]
      pa1m2[i]<-apath1$coefficients[3,4]
      pa2m1[i]<-apath2$coefficients[2,4]
      pa2m2[i]<-apath2$coefficients[3,4]
      pa1m3[i]<-apath3$coefficients[2,4]
      pa2m3[i]<-apath3$coefficients[3,4]
      pa1m4[i]<-apath4$coefficients[2,4]
      pa2m4[i]<-apath4$coefficients[3,4]
      pb1[i]<-bpath$coefficients[4,4]
      pb2[i]<-bpath$coefficients[5,4]
      pb3[i]<-bpath$coefficients[6,4]
      pb4[i]<-bpath$coefficients[7,4]
      power = data.frame(pa1m1=pa1m1,pa1m2=pa1m2,pa1m3=pa1m3,pa1m4=pa1m4,pa2m1=pa2m1,pa2m2=pa2m2,pa2m3=pa2m3,pa2m4=pa2m4,pb1=pb1,pb2=pb2,pb3=pb3)
      power$jsa1m1["JointSiga1m1"]<-NA
      power$jsa1m1[pa1m1 < alpha] <- 1
      power$jsa1m1[pa1m1 >= alpha] <- 0
      power$jsa2m1["JointSiga2m1"]<-NA
      power$jsa2m1[pa2m1 < alpha] <- 1
      power$jsa2m1[pa2m1 >= alpha] <- 0
      power$jsa1m2["JointSiga1m2"]<-NA
      power$jsa1m2[pa1m2 < alpha] <- 1
      power$jsa1m2[pa1m2 >= alpha] <- 0
      power$jsa2m2["JointSiga2m2"]<-NA
      power$jsa2m2[pa2m2 < alpha] <- 1
      power$jsa2m2[pa2m2 >= alpha] <- 0
      power$jsa1m3["JointSiga1m3"]<-NA
      power$jsa1m3[pa1m3 < alpha] <- 1
      power$jsa1m3[pa1m3 >= alpha] <- 0
      power$jsa2m3["JointSiga2m3"]<-NA
      power$jsa2m3[pa2m3 < alpha] <- 1
      power$jsa2m3[pa2m3 >= alpha] <- 0
      power$jsa1m4["JointSiga1m4"]<-NA
      power$jsa1m4[pa1m4 < alpha] <- 1
      power$jsa1m4[pa1m4 >= alpha] <- 0
      power$jsa2m4["JointSiga2m4"]<-NA
      power$jsa2m4[pa2m4 < alpha] <- 1
      power$jsa2m4[pa2m4 >= alpha] <- 0
      power$jsb1["JointSigb1"]<-NA
      power$jsb1[pb1 < alpha] <- 1
      power$jsb1[pb1 >= alpha] <- 0
      power$jsb2["JointSigb2"]<-NA
      power$jsb2[pb2 < alpha] <- 1
      power$jsb2[pb2 >= alpha] <- 0
      power$jsb3["JointSigb3"]<-NA
      power$jsb3[pb3 < alpha] <- 1
      power$jsb3[pb3 >= alpha] <- 0
      power$jsb4["JointSigb4"]<-NA
      power$jsb4[pb3 < alpha] <- 1
      power$jsb4[pb3 >= alpha] <- 0
      power$jointa1b1m1<-power$jsa1m1*power$jsb1
      power$jointa2b1m1<-power$jsa2m1*power$jsb1
      power$jointa1b2m2<-power$jsa1m2*power$jsb2
      power$jointa2b2m2<-power$jsa2m2*power$jsb2
      power$jointa1b3m3<-power$jsa1m3*power$jsb3
      power$jointa2b3m3<-power$jsa2m3*power$jsb3
      power$jointa1b4m4<-power$jsa1m4*power$jsb4
      power$jointa2b4m4<-power$jsa2m4*power$jsb4
      JSx1m1<-mean(power$jointa1b1)
      JSx2m1<-mean(power$jointa2b1)
      JSx1m2<-mean(power$jointa1b2)
      JSx2m2<-mean(power$jointa2b2)
      JSx1m3<-mean(power$jointa1b3)
      JSx2m3<-mean(power$jointa2b3)
      JSx1m4<-mean(power$jointa1b4)
      JSx2m4<-mean(power$jointa2b4)

    }
    message("Sample size is ",n)
    message("Power for X1-M1 ", JSx1m1)
    message("Power for X2-M1 ", JSx2m1)
    message("Power for X1-M2 ", JSx1m2)
    message("Power for X2-M2 ", JSx2m2)
    message("Power for X1-M3 ", JSx1m3)
    message("Power for X2-M3 ", JSx2m3)
    message("Power for X1-M4 ", JSx1m4)
    message("Power for X2-M4 ", JSx2m4)
  }

  }


