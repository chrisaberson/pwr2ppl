#'Compute power for a One Factor Within Subjects LMM Trends with up to four levels.
#'Takes means, sds, and sample sizes for each group. Alpha is .05 by default, alternative values may be entered by user
#'@param m1 Mean of first time point
#'@param m2 Mean of second time point
#'@param m3 Mean of third time point
#'@param m4 Mean of fourth time point
#'@param s1 Standard deviation of first time point
#'@param s2 Standard deviation of second time point
#'@param s3 Standard deviation of third time point
#'@param s4 Standard deviation of forth time point
#'@param r12 correlation Time 1 and Time 2
#'@param r13 correlation Time 1 and Time 3
#'@param r14 correlation Time 1 and Time 4
#'@param r23 correlation Time 2 and Time 3
#'@param r24 correlation Time 2 and Time 4
#'@param r34 correlation Time 3 and Time 4
#'@param n Sample size for first group
#'@param alpha Type I error (default is .05)
#'@examples
#'lmm1Ftrends(m1=-.25,m2=-.15,m3=-.05,m4=.05,s1=.4,s2=.5,s3=.6,s4=.7,
#'r12=.50, r13=.30, r14=.15, r23=.5, r24=.30, r34=.50, n=25)
#'@return Power for the One Factor Within Subjects LMM Trends
#'@export

lmm1Ftrends<-function(m1,m2,m3=NA,m4=NA, s1, s2, s3=NULL,s4=NULL,
                r12, r13=NULL, r14=NULL, r23=NULL, r24=NULL, r34=NULL,
                n, alpha=.05)
{
  V1<-V2<-id<-anova<-V3<-V4<-NULL
  levels<-NA
  levels[is.na(m4) & is.na(m3)]<-2
  levels[is.na(m4) & !is.na(m3)]<-3
  levels[!is.na(m4)]<-4

  if(levels<2|levels>4){stop("Function requires 2 to 4 levels")}
  if(levels=="2"){
    var1<-s1^2
    var2<-s2^2
    cov12<-r12*s1*s2
    out <- MASS::mvrnorm(n, mu = c(m1,m2), Sigma = matrix(c(var1,cov12,
                                                      cov12,var2)
                                                    , ncol = 2),
                   empirical = TRUE)
    out<-as.data.frame(out)
    out<-dplyr::rename(out, y1 = V1, y2 = V2)
    out$id <- rep(1:nrow(out))
    out$id<-as.factor(out$id)
    out<-tidyr::gather(out,key="iv",value="dv",-id)
    out$iv<-as.ordered(out$iv)
    options(contrasts=c("contr.henlme::lmert", "contr.poly"))
    base<-nlme::lme(dv~1, random = ~1|id/iv, data=out,method="ML")
    model1<-nlme::lme(dv~iv, random = ~1|id/iv, data=out,method="ML")
    lm<-stats::anova(base,model1)
    df1<-lm$df[2]-lm$df[1]
    lambdalm<-lm$L.Ratio[2]
    minusalpha<-1-alpha
    tabledlm<-stats::qchisq(minusalpha, df1)
    powerlm<-round(1-stats::pchisq(tabledlm, df1, lambdalm),3)
    {print(paste("Power for n =",n,"=", powerlm))}
    {print(paste("Tests use df =", df))}
    {print(paste("This test is useless for two levels"))}
    }

  if(levels==3){
    var1<-s1^2
    var2<-s2^2
    var3<-s3^2
    cov12<-r12*s1*s2
    cov13<-r13*s1*s3
    cov23<-r23*s2*s3
    out <- MASS::mvrnorm(n, mu = c(m1,m2,m3), Sigma = matrix(c(var1,cov12,cov13,
                                                         cov12,var2,cov23,
                                                         cov13, cov23,var3), ncol = 3),
                   empirical = TRUE)
    out<-as.data.frame(out)
    out<-dplyr::rename(out, y1 = V1, y2 = V2, y3 = V3)
    out$id <- rep(1:nrow(out))
    out$id<-as.factor(out$id)
    out<-tidyr::gather(out,key="iv",value="dv",-id)
    out$iv<-as.ordered(out$iv)
    options(contrasts=c("contr.henlme::lmert", "contr.poly"))
    model1<-nlme::lme(dv~iv, random = ~1|id/iv, data=out,method="ML")
    trends<-summary(model1)
    row<-trends$tTable[1:17]
    lambdaLML<-(row[11])^2
    lambdaLMQ<-(row[12])^2
    df<-(row[7])
    minusalpha<-1-alpha
    tabled<-stats::qf(minusalpha,1,df)
    powerLM.L<-round(1-stats::pf(tabled, 1, df, lambdaLML),3)
    powerLM.Q<-round(1-stats::pf(tabled, 1, df, lambdaLMQ),3)
    {print(paste("Power Linear Trend for n =",n,"=", powerLM.L))}
    {print(paste("Power Quadratic Trend for n =",n,"=", powerLM.Q))}
    {print(paste("Tests use df =", df))}
    }



  if (levels=="4"){
    var1<-s1^2
    var2<-s2^2
    var3<-s3^2
    var4<-s4^2
    cov12<-r12*s1*s2
    cov13<-r13*s1*s3
    cov14<-r14*s1*s4
    cov23<-r23*s2*s3
    cov24<-r24*s2*s4
    cov34<-r34*s3*s4
    out <- MASS::mvrnorm(n, mu = c(m1,m2,m3,m4), Sigma = matrix(c(var1,cov12,cov13, cov14,
                                                            cov12,var2,cov23, cov24,
                                                            cov13, cov23,var3, cov34,
                                                            cov14, cov24, cov34, var4), ncol = 4),
                   empirical = TRUE)
    out<-as.data.frame(out)
    out<-dplyr::rename(out, y1 = V1, y2 = V2, y3 = V3, y4 = V4)
    out$id <- rep(1:nrow(out))
    out$id<-as.factor(out$id)
    out<-tidyr::gather(out,key="iv",value="dv",-id)
    out$iv<-as.ordered(out$iv)
    options(contrasts=c("contr.henlme::lmert", "contr.poly"))
    base<-nlme::lme(dv~1, random = ~1|id/iv, data=out,method="ML")
    model1<-nlme::lme(dv~iv, random = ~1|id/iv, data=out,method="ML")
    trends<-summary(model1)
    row<-trends$tTable[1:18]
    lambdaLML<-(row[14])^2
    lambdaLMQ<-(row[15])^2
    lambdaLMC<-(row[16])^2
    df<-(row[9])
    minusalpha<-1-alpha
    tabled<-stats::qf(minusalpha,1,df)
    powerLM.L<-round(1-stats::pf(tabled, 1, df, lambdaLML),3)
    powerLM.Q<-round(1-stats::pf(tabled, 1, df, lambdaLMQ),3)
    powerLM.C<-round(1-stats::pf(tabled, 1, df, lambdaLMC),3)
    {print(paste("Power Linear Trend for n =",n,"=", powerLM.L))}
    {print(paste("Power Quadratic Trend for n =",n,"=", powerLM.Q))}
    {print(paste("Power Cubic Trend for n =",n,"=", powerLM.C))}
    {print(paste("Tests use df =", df))}
    }
  on.exit()}



