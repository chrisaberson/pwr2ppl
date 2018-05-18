#'Compute power for a One Factor Within Subjects ANOVA with up to four levels.
#'Takes means, sds, and sample sizes for each group. Alpha is .05 by default, alterative values may be entered by user
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
#'@return Power for the One Factor Within Subjects ANOVA
#'@export

win1F<-function(m1,m2,m3=NA,m4=NA, s1, s2, s3=NULL,s4=NULL,
                    r12, r13=NULL, r14=NULL, r23=NULL, r24=NULL, r34=NULL,
                    n, alpha=.05)
{
levels<-NA
levels[is.na(m4) & is.na(m3)]<-2
levels[is.na(m4) & !is.na(m3)]<-3
levels[!is.na(m4)]<-4

if(levels<2|levels>4){stop("Function requires 2 to 4 levels")}
if(levels=="2"){
  var1<-s1^2
  var2<-s2^2
  cov12<-r12*s1*s2
  out <- mvrnorm(n, mu = c(m1,m2), Sigma = matrix(c(var1,cov12,
                                                      cov12,var2)
                                                    , ncol = 2),
                   empirical = TRUE)
    out<-as.data.frame(out)
    out<-rename(out, y1 = V1, y2 = V2)
    out$id <- rep(1:nrow(out))
    out$id<-as.factor(out$id)
    out<-gather(out,key="iv",value="dv",-id)
    out$iv<-as.ordered(out$iv)
    options(contrasts=c("contr.helmert", "contr.poly"))
    model<-ezANOVA(data=out, dv=.(dv), wid=.(id), within = .(iv), type=3, detailed=TRUE)
    df1<-model$ANOVA$DFn[2]
    df2<-model$ANOVA$DFd[2]
    SSB<-model$ANOVA$SSn[2]
    SSW<-model$ANOVA$SSd[2]
    eta2<-SSB/(SSB+SSW)
    f2<-eta2/(1-eta2)
    lambda<-f2*df2
    minusalpha<-1-alpha
    Ft<-qf(minusalpha, df1, df2)
    power<-round(1-pf(Ft, df1,df2,lambda),3)
    gge<-model$`Sphericity Corrections`$GGe
    hfe<-model$`Sphericity Corrections`$HFe
    ggdf1<-gge*df1
    ggdf2<-gge*df2
    hfdf1<-hfe*df1
    hfdf2<-hfe*df2
    lambdagg<-f2*ggdf2
    lambdahf<-f2*hfdf2
    Ftgg<-qf(minusalpha, ggdf1, ggdf2)
    Fthf<-qf(minusalpha, hfdf1, hfdf2)
    powergg<-round(1-pf(Ftgg, ggdf1,ggdf2,lambdagg),3)
    powerhf<-round(1-pf(Fthf, hfdf1,hfdf2,lambdahf),3)
    {print(paste("Power (Unadjusted) for n =",n,"=", power))}
    {print(paste("Adjusted Power not relevant with two levels"))}}

    if(levels==3){
    var1<-s1^2
    var2<-s2^2
    var3<-s3^2
    cov12<-r12*s1*s2
    cov13<-r13*s1*s3
    cov23<-r23*s2*s3
    out <- mvrnorm(n, mu = c(m1,m2,m3), Sigma = matrix(c(var1,cov12,cov13,
                                                         cov12,var2,cov23,
                                                         cov13, cov23,var3), ncol = 3),
                                                         empirical = TRUE)
    out<-as.data.frame(out)
    out<-rename(out, y1 = V1, y2 = V2, y3 = V3)
    out$id <- rep(1:nrow(out))
    out$id<-as.factor(out$id)
    out<-gather(out,key="iv",value="dv",-id)
    out$iv<-as.ordered(out$iv)
    options(contrasts=c("contr.helmert", "contr.poly"))
    model<-ezANOVA(data=out, dv=.(dv), wid=.(id), within = .(iv), type=3, detailed=TRUE)
    df1<-model$ANOVA$DFn[2]
    df2<-model$ANOVA$DFd[2]
    SSB<-model$ANOVA$SSn[2]
    SSW<-model$ANOVA$SSd[2]
    eta2<-SSB/(SSB+SSW)
    f2<-eta2/(1-eta2)
    lambda<-f2*df2
    minusalpha<-1-alpha
    Ft<-qf(minusalpha, df1, df2)
    power<-round(1-pf(Ft, df1,df2,lambda),3)
    gge<-round(model$`Sphericity Corrections`$GGe,3)
    hfe<-round(model$`Sphericity Corrections`$HFe,3)
    ggdf1<-gge*df1
    ggdf2<-gge*df2
    hfdf1<-hfe*df1
    hfdf2<-hfe*df2
    lambdagg<-f2*ggdf2
    lambdahf<-f2*hfdf2
    Ftgg<-qf(minusalpha, ggdf1, ggdf2)
    Fthf<-qf(minusalpha, hfdf1, hfdf2)
    powergg<-round(1-pf(Ftgg, ggdf1,ggdf2,lambdagg),3)
    powerhf<-round(1-pf(Fthf, hfdf1,hfdf2,lambdahf),3)
    {print(paste("Power (Unadjusted) for n =",n,"=", power))}
    {print(paste("Power H-F Adjusted (Epsilon = ",hfe ,") for n =",n, "=", powerhf))}
    {print(paste("Power G-G Adjusted (Epsilon = ", gge,") for n =",n, "=", powergg))}}
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
      out <- mvrnorm(n, mu = c(m1,m2,m3,m4), Sigma = matrix(c(var1,cov12,cov13, cov14,
                                                              cov12,var2,cov23, cov24,
                                                              cov13, cov23,var3, cov34,
                                                              cov14, cov24, cov34, var4), ncol = 4),
                     empirical = TRUE)
      out<-as.data.frame(out)
      out<-rename(out, y1 = V1, y2 = V2, y3 = V3, y4 = V4)
      out$id <- rep(1:nrow(out))
      out$id<-as.factor(out$id)
      out<-gather(out,key="iv",value="dv",-id)
      out$iv<-as.ordered(out$iv)
      options(contrasts=c("contr.helmert", "contr.poly"))
      model<-ezANOVA(data=out, dv=.(dv), wid=.(id), within = .(iv), type=3, detailed=TRUE)
      df1<-model$ANOVA$DFn[2]
      df2<-model$ANOVA$DFd[2]
      SSB<-model$ANOVA$SSn[2]
      SSW<-model$ANOVA$SSd[2]
      eta2<-SSB/(SSB+SSW)
      f2<-eta2/(1-eta2)
      lambda<-f2*df2
      minusalpha<-1-alpha
      Ft<-qf(minusalpha, df1, df2)
      power<-round(1-pf(Ft, df1,df2,lambda),3)
      gge<-round(model$`Sphericity Corrections`$GGe,3)
      hfe<-round(model$`Sphericity Corrections`$HFe,3)
      ggdf1<-gge*df1
      ggdf2<-gge*df2
      hfdf1<-hfe*df1
      hfdf2<-hfe*df2
      lambdagg<-f2*ggdf2
      lambdahf<-f2*hfdf2
      Ftgg<-qf(minusalpha, ggdf1, ggdf2)
      Fthf<-qf(minusalpha, hfdf1, hfdf2)
      powergg<-round(1-pf(Ftgg, ggdf1,ggdf2,lambdagg),3)
      powerhf<-round(1-pf(Fthf, hfdf1,hfdf2,lambdahf),3)
      {print(paste("Power (Unadjusted) for n =",n,"=", power))}
      {print(paste("Power H-F Adjusted (Epsilon = ",hfe ,") for n =",n, "=", powerhf))}
      {print(paste("Power G-G Adjusted (Epsilon = ", gge,") for n =",n, "=", powergg))}}
}


