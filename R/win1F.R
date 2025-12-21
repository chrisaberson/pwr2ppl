#'Compute power for a One Factor Within Subjects ANOVA with up to four levels.
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
#'win1F(m1=-.25,m2=.00,m3=.10,m4=.15,s1=.4,s2=.5,s3=.6,s4=.7,
#'r12=.50, r13=.30, r14=.15, r23=.5, r24=.30, r34=.50, n=25)
#'win1F(m1=-.25,m2=.00,m3=.10,m4=.15,s1=.4,s2=.5,s3=2.5,s4=2.0,
#'r12=.50, r13=.30, r14=.10, r23=.5, r24=.30, r34=.40, n=100)
#'@return Power for the One Factor Within Subjects ANOVA
#'@export

win1F<-function(m1,m2,m3=NA,m4=NA, s1, s2, s3=NULL,s4=NULL,
                r12, r13=NULL, r14=NULL, r23=NULL, r24=NULL, r34=NULL,
                n, alpha=.05)
{
  V1<-V2<-V3<-V4<-dv<-iv<-id<-NULL
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
    options(contrasts=c("contr.helmert", "contr.poly"))
    model<-afex::aov_ez("id", "dv" , out, within = c("iv"), type=3, detailed=TRUE)
    xx<-summary(model)
    df1<-xx$univariate.tests[2,2]
    df2<-xx$univariate.tests[2,4]
    SSB<-xx$univariate.tests[2,1]
    SSW<-xx$univariate.tests[2,3]
    eta2<-SSB/(SSB+SSW)
    f2<-eta2/(1-eta2)
    lambda<-f2*df2
    minusalpha<-1-alpha
    Ft<-qf(minusalpha, df1, df2)
    power<-round(1-pf(Ft, df1,df2,lambda),3)
    gge<-1
    hfe<-1
    ggdf1<-gge*df1
    ggdf2<-gge*df2
    hfdf1<-hfe*df1
    hfdf2<-hfe*df2
    lambdagg<-f2*ggdf2
    lambdahf<-f2*hfdf2
    Ftgg<-stats::qf(minusalpha, ggdf1, ggdf2)
    Fthf<-stats::qf(minusalpha, hfdf1, hfdf2)
    powergg<-round(1-stats::pf(Ftgg, ggdf1,ggdf2,lambdagg),3)
    powerhf<-round(1-stats::pf(Fthf, hfdf1,hfdf2,lambdahf),3)
    eta2<-round((eta2),3)
    message("partial eta-squared = ", eta2)
    message("Power (Unadjusted) for n = ",n," is ", power)
    message("Adjusted Power not relevant with two levels")
    result <- data.frame(matrix(ncol = 3))
    colnames(result) <- c("n", "eta2","Power")
    result[, 1]<-n
    result[, 2]<-eta2
    result[, 3]<-power
    output<-na.omit(result)
    rownames(output)<- c()
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
    options(contrasts=c("contr.helmert", "contr.poly"))
    model<-afex::aov_ez("id", "dv" , out, within = c("iv"), type=3, detailed=TRUE)
    xx<-summary(model)
    df1<-xx$univariate.tests[2,2]
    df2<-xx$univariate.tests[2,4]
    SSB<-xx$univariate.tests[2,1]
    SSW<-xx$univariate.tests[2,3]
    eta2<-SSB/(SSB+SSW)
    f2<-eta2/(1-eta2)
    lambda<-f2*df2
    minusalpha<-1-alpha
    Ft<-qf(minusalpha, df1, df2)
    power<-round(1-pf(Ft, df1,df2,lambda),3)
    gge<-round(xx$pval.adjustments[1],3)
    hfe<-round(xx$pval.adjustments[3],3)
    ggdf1<-gge*df1
    ggdf2<-gge*df2
    hfdf1<-hfe*df1
    hfdf2<-hfe*df2
    lambdagg<-f2*ggdf2
    lambdahf<-f2*hfdf2
    Ftgg<-stats::qf(minusalpha, ggdf1, ggdf2)
    Fthf<-stats::qf(minusalpha, hfdf1, hfdf2)
    powergg<-round(1-stats::pf(Ftgg, ggdf1,ggdf2,lambdagg),3)
    powerhf<-round(1-stats::pf(Fthf, hfdf1,hfdf2,lambdahf),3)
    eta2<-round((eta2),3)
    message("partial eta-squared = ", eta2)
    message("Power (Unadjusted) for n = ",n," is ", power)
    message("Power H-F Adjusted (Epsilon = ",hfe ,") for n = ",n, " is ", powerhf)
    message("Power G-G Adjusted (Epsilon = ", gge,") for n = ",n, " is ", powergg)
    result <- data.frame(matrix(ncol = 7))
    colnames(result) <- c("n", "eta2","Power (Unadjusted) ","H-F epsilon",
                          "Power (H-F)","GG Epsilon","Power (G-G)")
    result[, 1]<-n
    result[, 2]<-eta2
    result[, 3]<-power
    result[, 4]<-hfe
    result[, 5]<-powerhf
    result[, 6]<-gge
    result[, 7]<-powergg
    output<-na.omit(result)
    rownames(output)<- c()
    }
  if (levels==4){
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
    options(contrasts=c("contr.helmert", "contr.poly"))
    model<-afex::aov_ez("id", "dv" , out, within = c("iv"), type=3, detailed=TRUE)
    xx<-summary(model)
    df1<-xx$univariate.tests[2,2]
    df2<-xx$univariate.tests[2,4]
    SSB<-xx$univariate.tests[2,1]
    SSW<-xx$univariate.tests[2,3]
    eta2<-SSB/(SSB+SSW)
    f2<-eta2/(1-eta2)
    lambda<-f2*df2
    minusalpha<-1-alpha
    Ft<-qf(minusalpha, df1, df2)
    power<-round(1-pf(Ft, df1,df2,lambda),3)
    gge<-round(xx$pval.adjustments[1],3)
    hfe<-round(xx$pval.adjustments[3],3)
    ggdf1<-gge*df1
    ggdf2<-gge*df2
    hfdf1<-hfe*df1
    hfdf2<-hfe*df2
    lambdagg<-f2*ggdf2
    lambdahf<-f2*hfdf2
    Ftgg<-stats::qf(minusalpha, ggdf1, ggdf2)
    Fthf<-stats::qf(minusalpha, hfdf1, hfdf2)
    powergg<-round(1-stats::pf(Ftgg, ggdf1,ggdf2,lambdagg),3)
    powerhf<-round(1-stats::pf(Fthf, hfdf1,hfdf2,lambdahf),3)
    eta2<-round((eta2),3)
    message("partial eta-squared = ", eta2)
    message("Power (Unadjusted) for n = ",n," is ", power)
    message("Power H-F Adjusted (Epsilon = ",hfe ,") for n = ",n, " is ", powerhf)
    message("Power G-G Adjusted (Epsilon = ", gge,") for n = ",n, " is ", powergg)
    result <- data.frame(matrix(ncol = 7))
    colnames(result) <- c("n", "eta2","Power (Unadjusted) ","H-F epsilon",
                          "Power (H-F)","GG Epsilon","Power (G-G)")
    result[, 1]<-n
    result[, 2]<-eta2
    result[, 3]<-power
    result[, 4]<-hfe
    result[, 5]<-powerhf
    result[, 6]<-gge
    result[, 7]<-powergg
    output<-na.omit(result)
    rownames(output)<- c()
    }
 invisible(output)
  }
