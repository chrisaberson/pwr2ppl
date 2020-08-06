#'Compute power for a One Factor Within Subjects and One Factor Between LMM with up to two by four levels (within).
#'Takes means, sds, and sample sizes for each group. Alpha is .05 by default, alternative values may be entered by user
#'@param m1.1 Mean of first level Within Factor, 1st level Between Factor
#'@param m2.1 Mean of second level Within Factor, 1st level Between Factor
#'@param m3.1 Mean of third level Within Factor, 1st level Between Factor
#'@param m4.1 Mean of fourth level Within Factor, 1st level Between Factor
#'@param m1.2 Mean of first level Within Factor, 2nd level Between Factor
#'@param m2.2 Mean of second level Within Factor, 2nd level Between Factor
#'@param m3.2 Mean of third level Within Factor, 2nd level Between Factor
#'@param m4.2 Mean of fourth level Within Factor, 2nd level Between Factor
#'@param s1.1 Standard deviation of first level Within Factor, 1st level Between Factor
#'@param s2.1 Standard deviation of second level Within Factor, 1st level Between Factor
#'@param s3.1 Standard deviation of third level Within Factor, 1st level Between Factor
#'@param s4.1 Standard deviation of forth level Within Factor, 1st level Between Factor
#'@param s1.2 Standard deviation of first level Within Factor, 2nd level Between Factor
#'@param s2.2 Standard deviation of second level Within Factor, 2nd level Between Factor
#'@param s3.2 Standard deviation of third level Within Factor, 2nd level Between Factor
#'@param s4.2 Standard deviation of forth level Within Factor, 2nd level Between Factor
#'@param r1.2_1 correlation Within Factor Level 1 and Within Factor, Level 2, 1st level Between
#'@param r1.3_1 correlation Within Factor Level 1 and Within Factor, Level 3, 1st level Between
#'@param r1.4_1 correlation Within Factor Level 1 and Within Factor, Level 4, 1st level Between
#'@param r2.3_1 correlation Within Factor Level 1 and Within Factor, Level 3, 1st level Between
#'@param r2.4_1 correlation Within Factor Level 1 and Within Factor, Level 4, 1st level Between
#'@param r3.4_1 correlation Within Factor Level 1 and Within Factor, Level 4, 1st level Between
#'@param r1.2_2 correlation Within Factor Level 1 and Within Factor, Level 2, 2nd level Between
#'@param r1.3_2 correlation Within Factor Level 1 and Within Factor, Level 3, 2nd level Between
#'@param r1.4_2 correlation Within Factor Level 1 and Within Factor, Level 4, 2nd level Between
#'@param r2.3_2 correlation Within Factor Level 1 and Within Factor, Level 3, 2nd level Between
#'@param r2.4_2 correlation Within Factor Level 1 and Within Factor, Level 4, 2nd level Between
#'@param r3.4_2 correlation Within Factor Level 1 and Within Factor, Level 4, 2nd level Between
#'@param r sets same correlations between DVs on all factor levels (seriously, just use this)
#'@param s sets same standard deviation for factor levels (see comment above)
#'@param n n for each between group level
#'@param alpha Type I error (default is .05)
#'@examples
#'lmm1w1b(m1.1 = -.25, m2.1=0, m3.1=0.10, m4.1=.15,
#'m1.2=-.25,m2.2=-.25,m3.2=-.25, m4.2=-.25,
#'s1.1 = .4, s2.1=.5, s3.1=0.6, s4.1=.7,
#'s1.2=.4,s2.2=.5,s3.2=.6, s4.2=.7,n = 50,
#'r1.2_1=.5,r1.3_1=.3,r1.4_1=.15,r2.3_1=.5,r2.4_1=.3,r3.4_1=.5,
#'r1.2_2=.5,r1.3_2=.3,r1.4_2=.15, r2.3_2=.5,r2.4_2=.3,r3.4_2=.5)
#'lmm1w1b(m1.1 = -.25, m2.1=0, m3.1=0.10, m4.1=.15,
#'m1.2=-.25,m2.2=-.25,m3.2=-.25, m4.2=-.25, s=.4, r = .5, n=100)
#'@return Power for the One Factor Within Subjects and One Factor Between LMM
#'@export
lmm1w1b<-function(m1.1,m2.1,m3.1=NA,m4.1=NA,m1.2,m2.2,m3.2=NA,m4.2=NA,
                  s1.1=NA,s2.1=NA,s3.1=NA,s4.1=NA,s1.2=NA,s2.2=NA,s3.2=NA,s4.2=NA,
                  r1.2_1=NULL, r1.3_1=NULL, r1.4_1=NULL,
                  r2.3_1=NULL, r2.4_1=NULL,
                  r3.4_1=NULL,
                  r1.2_2=NULL, r1.3_2=NULL, r1.4_2=NULL,
                  r2.3_2=NULL, r2.4_2=NULL,
                  r3.4_2=NULL, r=NULL, s = NULL, n, alpha=.05)

{
  V1<-V2<-V3<-V4<-id<-ivbg<-NULL
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
    out$id <- rep(1:nrow(out))
    out$id<-as.factor(out$id)
    out<-tidyr::gather(out,key="ivw",value="dv",-id, -ivbg)
    out$ivw<-as.factor(out$ivw)
    out$ivbg<-as.factor(out$ivbg)
    base<-nlme::lme(dv~1, random = ~1|id/ivw, data=out,method="ML")
    model1<-nlme::lme(dv~ivw, random = ~1|id/ivw, data=out,method="ML") #factor A
    model2<-nlme::lme(dv~ivw+ivbg, random = ~1|id/ivw, data=out,method="ML") #factor B
    model3<-nlme::lme(dv~ivw+ivbg+ivw*ivbg, random = ~1|id/ivw, data=out,method="ML") #AxB
    lm<-stats::anova(base,model1,model2,model3)
    df1<-lm$df[2]-lm$df[1]
    df2<-lm$df[3]-lm$df[2]
    df3<-lm$df[4]-lm$df[3]
    lambdalm1<-lm$L.Ratio[2]
    lambdalm2<-lm$L.Ratio[3]
    lambdalm3<-lm$L.Ratio[4]
    tabledlm1<-stats::qchisq(.95, df1)
    tabledlm2<-stats::qchisq(.95, df2)
    tabledlm3<-stats::qchisq(.95, df3)
    powerlm1<-round(1-stats::pchisq(tabledlm1, df1, lambdalm1),3)
    powerlm2<-round(1-stats::pchisq(tabledlm2, df2, lambdalm2),3)
    powerlm3<-round(1-stats::pchisq(tabledlm3, df3, lambdalm3),3)
    message("Power Factor A (Between) for n = ",n," is ", powerlm2)
    message("Power Factor B (Within) for n = ",n," is ", powerlm1)
    message("Power AxB for n = ",n," is ", powerlm3)
    result <- data.frame(matrix(ncol = 4))
    colnames(result) <- c("n", "Power A (Between)", "Power B (Within)", "Power AxB")
    result[, 1]<-n
    result[, 2]<-powerlm2
    result[, 3]<-powerlm1
    result[, 4]<-powerlm3
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
    out$id <- rep(1:nrow(out))
    out$id<-as.factor(out$id)
    out<-tidyr::gather(out,key="ivw",value="dv",-id, -ivbg)
    out$ivw<-as.factor(out$ivw)
    out$ivbg<-as.factor(out$ivbg)
    base<-nlme::lme(dv~1, random = ~1|id/ivw, data=out,method="ML")
    model1<-nlme::lme(dv~ivw, random = ~1|id/ivw, data=out,method="ML") #factor A
    model2<-nlme::lme(dv~ivw+ivbg, random = ~1|id/ivw, data=out,method="ML") #factor B
    model3<-nlme::lme(dv~ivw+ivbg+ivw*ivbg, random = ~1|id/ivw, data=out,method="ML") #AxB
    lm<-stats::anova(base,model1,model2,model3)
    df1<-lm$df[2]-lm$df[1]
    df2<-lm$df[3]-lm$df[2]
    df3<-lm$df[4]-lm$df[3]
    lambdalm1<-lm$L.Ratio[2]
    lambdalm2<-lm$L.Ratio[3]
    lambdalm3<-lm$L.Ratio[4]
    tabledlm1<-stats::qchisq(.95, df1)
    tabledlm2<-stats::qchisq(.95, df2)
    tabledlm3<-stats::qchisq(.95, df3)
    powerlm1<-round(1-stats::pchisq(tabledlm1, df1, lambdalm1),3)
    powerlm2<-round(1-stats::pchisq(tabledlm2, df2, lambdalm2),3)
    powerlm3<-round(1-stats::pchisq(tabledlm3, df3, lambdalm3),3)
    message("Power Factor A (Between) for n = ",n," is ", powerlm2)
    message("Power Factor B (Within) for n = ",n," is ", powerlm1)
    message("Power AxB for n = ",n," is ", powerlm3)
    result <- data.frame(matrix(ncol = 4))
    colnames(result) <- c("n", "Power A (Between)", "Power B (Within)", "Power AxB")
    result[, 1]<-n
    result[, 2]<-powerlm2
    result[, 3]<-powerlm1
    result[, 4]<-powerlm3
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
    out$id <- rep(1:nrow(out))
    out$id<-as.factor(out$id)
    out<-tidyr::gather(out,key="ivw",value="dv",-id, -ivbg)
    out$ivw<-as.factor(out$ivw)
    out$ivbg<-as.factor(out$ivbg)

    base<-nlme::lme(dv~1, random = ~1|id/ivw, data=out,method="ML")
    model1<-nlme::lme(dv~ivw, random = ~1|id/ivw, data=out,method="ML") #factor A
    model2<-nlme::lme(dv~ivw+ivbg, random = ~1|id/ivw, data=out,method="ML") #factor B
    model3<-nlme::lme(dv~ivw+ivbg+ivw*ivbg, random = ~1|id/ivw, data=out,method="ML") #AxB
    lm<-stats::anova(base,model1,model2,model3)
    df1<-lm$df[2]-lm$df[1]
    df2<-lm$df[3]-lm$df[2]
    df3<-lm$df[4]-lm$df[3]
    lambdalm1<-lm$L.Ratio[2]
    lambdalm2<-lm$L.Ratio[3]
    lambdalm3<-lm$L.Ratio[4]
    tabledlm1<-stats::qchisq(.95, df1)
    tabledlm2<-stats::qchisq(.95, df2)
    tabledlm3<-stats::qchisq(.95, df3)
    powerlm1<-round(1-stats::pchisq(tabledlm1, df1, lambdalm1),3)
    powerlm2<-round(1-stats::pchisq(tabledlm2, df2, lambdalm2),3)
    powerlm3<-round(1-stats::pchisq(tabledlm3, df3, lambdalm3),3)
    message("Power Factor A (Between) for n = ",n," is ", powerlm2)
    message("Power Factor B (Within) for n = ",n," is ", powerlm1)
    message("Power AxB for n = ",n," is ", powerlm3)
    result <- data.frame(matrix(ncol = 4))
    colnames(result) <- c("n", "Power A (Between)", "Power B (Within)", "Power AxB")
    result[, 1]<-n
    result[, 2]<-powerlm2
    result[, 3]<-powerlm1
    result[, 4]<-powerlm3
    output<-na.omit(result)
    rownames(output)<- c()
  }
  invisible(output)
}


