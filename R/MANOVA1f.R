#'Compute power for a One Factor MANOVA with up to two levels and up to four measures.
#'Takes means, sds, and sample sizes for each group. Alpha is .05 by default, alternative values may be entered by user
#'@param m1.1 Mean of first DV, 1st level Between Factor
#'@param m2.1 Mean of second DV, 1st level Between Factor
#'@param m3.1 Mean of third DV, 1st level Between Factor
#'@param m4.1 Mean of fourth DV, 1st level Between Factor
#'@param m1.2 Mean of first DV, 2nd level Between Factor
#'@param m2.2 Mean of second DV, 2nd level Between Factor
#'@param m3.2 Mean of third DV, 2nd level Between Factor
#'@param m4.2 Mean of fourth DV, 2nd level Between Factor
#'@param s1.1 Standard deviation of first DV, 1st level Between Factor
#'@param s2.1 Standard deviation of second DV, 1st level Between Factor
#'@param s3.1 Standard deviation of third DV, 1st level Between Factor
#'@param s4.1 Standard deviation of forth DV, 1st level Between Factor
#'@param s1.2 Standard deviation of first DV, 2nd level Between Factor
#'@param s2.2 Standard deviation of second DV, 2nd level Between Factor
#'@param s3.2 Standard deviation of third DV, 2nd level Between Factor
#'@param s4.2 Standard deviation of forth DV, 2nd level Between Factor
#'@param r1.2_1 correlation DV 1 and DV 2, 1st level Between
#'@param r1.3_1 correlation DV 1 and DV 3, 1st level Between
#'@param r1.4_1 correlation DV 1 and DV 4, 1st level Between
#'@param r2.3_1 correlation DV 1 and DV 3, 1st level Between
#'@param r2.4_1 correlation DV 1 and DV 4, 1st level Between
#'@param r3.4_1 correlation DV 1 and DV 4, 1st level Between
#'@param r1.2_2 correlation DV 1 and DV 2, 2nd level Between
#'@param r1.3_2 correlation DV 1 and DV 3, 2nd level Between
#'@param r1.4_2 correlation DV 1 and DV 4, 2nd level Between
#'@param r2.3_2 correlation DV 1 and DV 3, 2nd level Between
#'@param r2.4_2 correlation DV 1 and DV 4, 2nd level Between
#'@param r3.4_2 correlation DV 1 and DV 4, 2nd level Between
#'@param r sets same correlations between DVs on all factor levels (seriously, just use this)
#'@param s sets same standard deviation for factor levels (see comment above)
#'@param n Sample size for first group
#'@param alpha Type I error (default is .05)
#'@examples
#'MANOVA1f(n=40,m1.1=0,m2.1=1,m3.1=2.4,m4.1=-0.7,
#'m1.2=-0.25,m2.2=-2,m3.2=2,m4.2=-1,
#'s1.1=.4,s2.1=5,s3.1=1.6,s4.1=1.2,
#'s1.2=.4,s2.2=5,s3.2=1.6,s4.2=1.2,
#'r1.2_1=.1,r1.3_1=.1,r1.4_1=.1,
#'r2.3_1=.35,r2.4_1=.45,r3.4_1=.40,
#'r1.2_2=.1,r1.3_2=.1,r1.4_2=.1,
#'r2.3_2=.35,r2.4_2=.45,r3.4_2=.40,alpha=.05)
#'MANOVA1f(n=40,m1.1=0,m2.1=1,m3.1=2.4,m4.1=-0.7,
#'m1.2=-0.25,m2.2=-2,m3.2=2,m4.2=-1,
#'s=.4,r=.5,alpha=.05)
#'@return Power for the One Factor Within Subjects and One Factor Between ANOVA
#'@export
MANOVA1f<-function(m1.1,m2.1,m3.1=NA,m4.1=NA,m1.2,m2.2,m3.2=NA,m4.2=NA,
                   s1.1=NA,s2.1=NA,s3.1=NA,s4.1=NA,s1.2=NA,s2.2=NA,s3.2=NA,s4.2=NA,
                   r1.2_1=NULL, r1.3_1=NULL, r1.4_1=NULL,
                   r2.3_1=NULL, r2.4_1=NULL,
                   r3.4_1=NULL,
                   r1.2_2=NULL, r1.3_2=NULL, r1.4_2=NULL,
                   r2.3_2=NULL, r2.4_2=NULL,
                   r3.4_2=NULL, r=NULL, s = NULL, n, alpha=.05)

{
  V1<-V2<-ivbg<-V3<-V4<-NULL

  levels<-NA
  levels[is.na(m4.1) & is.na(m3.1)]<-2
  levels[is.na(m4.1) & !is.na(m3.1)]<-3
  levels[!is.na(m4.1)]<-4
  oldoption<-options(contrasts=c("contr.helmert", "contr.poly"))
  oldoption
  on.exit(options(oldoption))

  if (levels=="2"){
    if (!is.null(s)){
      s1.1<-s; s2.1<-s;s1.2<-s;s2.2<-s
      var1<-s^2; var2<-s^2;var3<-s^2;var4<-s^2}
    if (is.null(s)){var1<-s1.1^2; var2<-s2.1^2;var3<-s1.2^2;var4<-s2.2^2}
    if (!is.null(r)){r1.2_1<-r;r1.2_2<-r}
    cov12<-r1.2_1*s1.1*s2.1
    cov34<-r1.2_2*s1.2*s2.2;
    out1 <- MASS::mvrnorm(n, mu = c(m1.1,m2.1),
                    Sigma = matrix(c(var1,cov12,
                                     cov12,var2), ncol = 2),
                    empirical = TRUE)
    out1<-as.data.frame(out1)
    out1$ivbg<-NA
    out1$ivbg<-1 #identifies group

    out2 <- MASS::mvrnorm(n, mu = c(m1.2,m2.2),
                    Sigma = matrix(c(var3,cov34,
                                     cov34,var4), ncol = 2),
                    empirical = TRUE)
    out2<-as.data.frame(out2)
    out2$ivbg<-NA
    out2$ivbg<-2

    out<-rbind(out1,out2)
    out<-as.data.frame(out)
    out<-dplyr::rename(out, y1 = V1, y2 = V2, ivbg = ivbg)
    out$ivbg<-as.factor(out$ivbg)
    dvs<-cbind(out$y1,out$y2)

    mmodel<-stats::manova(dvs~ivbg, data=out)
    values<-summary(mmodel, intercept=TRUE)
    pill<-values$stats[2,2]
    f2<-pill/(1-pill)
    dfb<-values$stats[2,4]
    dfw<-values$stats[2,5]
    lambda<-f2*dfw
    minusalpha<-1-alpha
    Ft<-stats::qf(minusalpha, dfb, dfw)
    power<-round(1-stats::pf(Ft, dfb,dfw,lambda),4)
    message("Power for n = ",n," is ", power)
    result <- data.frame(matrix(ncol = 2))
    colnames(result) <- c("n","power")
    result[, 1]<-n
    result[, 2]<-power
    output<-na.omit(result)
    rownames(output)<- c()
  }

  if (levels=="3"){
    if (!is.null(s)){
      s1.1<-s; s2.1<-s;s3.1<-s;s1.2<-s;s2.2<-s;s3.2<-s
      var1<-s^2; var2<-s^2;var3<-s^2;var4<-s^2;var5<-s^2;var6<-s^2}
    if (is.null(s)){var1<-s1.1^2; var2<-s2.1^2;var3<-s3.1^2;var4<-s1.2^2;var5<-s2.2^2; var6<-s3.2^2}
    if (!is.null(r)){r1.2_1<-r;r1.3_1<-r;
    r2.3_1<-r;
    r1.2_2<-r;r1.3_2<-r;
    r2.3_2<-r}
    cov12<-r1.2_1*s1.1*s2.1;cov13<-r1.3_1*s1.1*s3.1;
    cov23<-r2.3_1*s2.1*s3.1;
    cov45<-r1.2_2*s1.2*s2.2;cov46<-r1.3_2*s1.2*s3.2;
    cov56<-r2.3_2*s2.2*s3.2

    out1 <- MASS::mvrnorm(n, mu = c(m1.1,m2.1, m3.1),
                    Sigma = matrix(c(var1,cov12,cov13,
                                     cov12,var2,cov23,
                                     cov13,cov23,var3), ncol = 3),
                    empirical = TRUE)
    out1<-as.data.frame(out1)
    out1$ivbg<-NA
    out1$ivbg<-1 #identifies group

    out2 <- MASS::mvrnorm(n, mu = c(m1.2,m2.2,m3.2),
                    Sigma = matrix(c(var4,cov45,cov46,
                                     cov45,var5,cov56,
                                     cov46,cov56,var6), ncol = 3),
                    empirical = TRUE)
    out2<-as.data.frame(out2)
    out2$ivbg<-NA
    out2$ivbg<-2


    out<-rbind(out1,out2)
    out<-as.data.frame(out)
    out<-dplyr::rename(out, y1 = V1, y2 = V2,  y3=V3, ivbg = ivbg)
    out$ivbg<-as.factor(out$ivbg)

    dvs<-cbind(out$y1,out$y2,out$y3)

    mmodel<-stats::manova(dvs~ivbg, data=out)
    values<-summary(mmodel, intercept=TRUE)
    pill<-values$stats[2,2]
    f2<-pill/(1-pill)
    dfb<-values$stats[2,4]
    dfw<-values$stats[2,5]
    lambda<-f2*dfw
    minusalpha<-1-alpha
    Ft<-stats::qf(minusalpha, dfb, dfw)
    power<-round(1-stats::pf(Ft, dfb,dfw,lambda),4)
    message("Power for n = ",n," is ", power)
    result <- data.frame(matrix(ncol = 2))
    colnames(result) <- c("n","power")
    result[, 1]<-n
    result[, 2]<-power
    output<-na.omit(result)
    rownames(output)<- c()
  }

  if (levels=="4"){
    if (!is.null(s)){
      s1.1<-s; s2.1<-s;s3.1<-s;s4.1<-s;s1.2<-s;s2.2<-s;s3.2<-s;s4.2<-s
      var1<-s^2; var2<-s^2;var3<-s^2;var4<-s^2;var5<-s^2;var6<-s^2;var7<-s^2;var8<-s^2}
    if (is.null(s)){var1<-s1.1^2; var2<-s2.1^2;var3<-s3.1^2;var4<-s4.1^2;var5<-s1.2^2;var6<-s2.2^2;var7<-s3.2^2;var8<-s4.2^2}
    if (!is.null(r)){r1.2_1<-r;r1.3_1<-r;r1.4_1<-r;r2.3_1<-r;r2.4_1<-r;
    r3.4_1<-r;
    r1.2_2<-r;r1.3_2<-r;r1.4_2<-r
    r2.3_2<-r;r2.4_2<-r;
    r3.4_2<-r}
    cov12<-r1.2_1*s1.1*s2.1;cov13<-r1.3_1*s1.1*s3.1;cov14<-r1.4_1*s1.1*s4.1
    cov23<-r2.3_1*s2.1*s3.1;cov24<-r2.4_1*s2.1*s4.1;
    cov34<-r3.4_1*s3.1*s4.1;
    cov56<-r1.2_2*s1.2*s2.2;cov57<-r1.3_2*s1.2*s3.2;cov58<-r1.4_2*s1.2*s4.2
    cov67<-r2.3_2*s2.2*s3.2;cov68<-r2.4_2*s2.2*s4.2
    cov78<-r3.4_2*s3.2*s4.2


    out1 <- MASS::mvrnorm(n, mu = c(m1.1,m2.1, m3.1, m4.1),
                    Sigma = matrix(c(var1,cov12,cov13,cov14,
                                     cov12,var2,cov23,cov24,
                                     cov13,cov23,var3, cov34,
                                     cov14,cov24,cov34,var4), ncol = 4),
                    empirical = TRUE)
    out1<-as.data.frame(out1)
    out1$ivbg<-NA
    out1$ivbg<-1 #identifies group

    out2 <- MASS::mvrnorm(n, mu = c(m1.2,m2.2,m3.2,m4.2),
                    Sigma = matrix(c(var5,cov56,cov57,cov58,
                                     cov56,var6,cov67,cov68,
                                     cov57,cov67,var7,cov78,
                                     cov58,cov68,cov78,var8), ncol = 4),
                    empirical = TRUE)
    out2<-as.data.frame(out2)
    out2$ivbg<-NA
    out2$ivbg<-2

    out<-rbind(out1,out2)
    out<-as.data.frame(out)
    out<-dplyr::rename(out, y1 = V1, y2 = V2, y3=V3, y4=V4,ivbg = ivbg)
    out$ivbg<-as.factor(out$ivbg)
    dvs<-cbind(out$y1,out$y2,out$y3,out$y4)

    mmodel<-stats::manova(dvs~ivbg, data=out)
    values<-summary(mmodel, intercept=TRUE)
    pill<-values$stats[2,2]
    f2<-pill/(1-pill)
    dfb<-values$stats[2,4]
    dfw<-values$stats[2,5]
    lambda<-f2*dfw
    minusalpha<-1-alpha
    Ft<-stats::qf(minusalpha, dfb, dfw)
    power<-round(1-stats::pf(Ft, dfb,dfw,lambda),4)
    message("Power for n = ",n," is ", power)
    result <- data.frame(matrix(ncol = 2))
    colnames(result) <- c("n","power")
    result[, 1]<-n
    result[, 2]<-power
    output<-na.omit(result)
    rownames(output)<- c()
  }
    invisible(output)
}
