#'Compute power for a Two Factor Within Subjects Using Linear Mixed Models with up to two by four levels.
#'Takes means, sds, and sample sizes for each group. Alpha is .05 by default, alternative values may be entered by user
#'@param m1.1 Mean of first level factor 1, 1st level factor two
#'@param m2.1 Mean of second level factor 1, 1st level factor two
#'@param m3.1 Mean of third level factor 1, 1st level factor two
#'@param m4.1 Mean of fourth level factor 1, 1st level factor two
#'@param m1.2 Mean of first level factor 1, 2nd level factor two
#'@param m2.2 Mean of second level factor 1, 2nd level factor two
#'@param m3.2 Mean of third level factor 1, 2nd level factor two
#'@param m4.2 Mean of fourth level factor 1, 2nd level factor two
#'@param s1.1 Standard deviation of first level factor 1, 1st level factor two
#'@param s2.1 Standard deviation of second level factor 1, 1st level factor two
#'@param s3.1 Standard deviation of third level factor 1, 1st level factor two
#'@param s4.1 Standard deviation of forth level factor 1, 1st level factor two
#'@param s1.2 Standard deviation of first level factor 1, 2nd level factor two
#'@param s2.2 Standard deviation of second level factor 1, 2nd level factor two
#'@param s3.2 Standard deviation of third level factor 1, 2nd level factor two
#'@param s4.2 Standard deviation of forth level factor 1, 2nd level factor two
#'@param r12 correlation Factor 1, Level 1 and Factor 1, Level 2
#'@param r13 correlation Factor 1, Level 1 and Factor 1, Level 3
#'@param r14 correlation Factor 1, Level 1 and Factor 1, Level 4
#'@param r15 correlation Factor 1, Level 1 and Factor 2, Level 1
#'@param r16 correlation Factor 1, Level 1 and Factor 2, Level 2
#'@param r17 correlation Factor 1, Level 1 and Factor 2, Level 3
#'@param r18 correlation Factor 1, Level 1 and Factor 2, Level 4
#'@param r23 correlation Factor 1, Level 2 and Factor 1, Level 3
#'@param r24 correlation Factor 1, Level 2 and Factor 1, Level 4
#'@param r25 correlation Factor 1, Level 2 and Factor 2, Level 1
#'@param r26 correlation Factor 1, Level 2 and Factor 2, Level 2
#'@param r27 correlation Factor 1, Level 2 and Factor 2, Level 3
#'@param r28 correlation Factor 1, Level 2 and Factor 2, Level 4
#'@param r34 correlation Factor 1, Level 3 and Factor 1, Level 4
#'@param r35 correlation Factor 1, Level 3 and Factor 2, Level 1
#'@param r36 correlation Factor 1, Level 3 and Factor 2, Level 2
#'@param r37 correlation Factor 1, Level 3 and Factor 2, Level 3
#'@param r38 correlation Factor 1, Level 3 and Factor 2, Level 4
#'@param r45 correlation Factor 1, Level 4 and Factor 2, Level 1
#'@param r46 correlation Factor 1, Level 4 and Factor 2, Level 2
#'@param r47 correlation Factor 1, Level 4 and Factor 2, Level 3
#'@param r48 correlation Factor 1, Level 4 and Factor 2, Level 4
#'@param r56 correlation Factor 2, Level 1 and Factor 2, Level 2
#'@param r57 correlation Factor 2, Level 1 and Factor 2, Level 3
#'@param r58 correlation Factor 2, Level 1 and Factor 2, Level 4
#'@param r67 correlation Factor 2, Level 2 and Factor 2, Level 3
#'@param r68 correlation Factor 2, Level 2 and Factor 2, Level 4
#'@param r78 correlation Factor 2, Level 3 and Factor 2, Level 4
#'@param r sets same correlations between DVs on all factor levels (seriously, just use this)
#'@param s sets same standard deviation for factor levels (see comment above)
#'@param n Sample size for first group
#'@param alpha Type I error (default is .05)
#'@examples
#'lmm2Fse(m1.1=-.25,m2.1=0,m3.1=.10,m4.1=.15,m1.2=-.25,m2.2=.10,m3.2=.30,m4.2=.35,
#'s1.1=.4,s2.1=.5,s3.1=2.5,s4.1=2.0,s1.2=.4,s2.2=.5,s3.2=2.5,s4.2=2.0,r=.5,n=220)
#'@return Power for Simple Effects in Two Factor Within Subjects LMM
#'@export

lmm2Fse<-function(m1.1,m2.1,m3.1=NA,m4.1=NA,m1.2,m2.2,m3.2=NA,m4.2=NA,
                s1.1=NA,s2.1=NA,s3.1=NA,s4.1=NA,s1.2=NA,s2.2=NA,s3.2=NA,s4.2=NA,
                r12=NULL, r13=NULL, r14=NULL, r15=NULL, r16=NULL, r17=NULL, r18=NULL,
                r23=NULL, r24=NULL, r25=NULL, r26=NULL, r27=NULL, r28=NULL,
                r34=NULL, r35=NULL, r36=NULL, r37=NULL, r38=NULL,
                r45=NULL, r46=NULL, r47=NULL, r48=NULL,
                r56=NULL, r57=NULL, r58=NULL,
                r67=NULL, r68=NULL,
                r78=NULL, r=NULL, s = NULL, n, alpha=.05)
{
  V1<-V2<-V3<-V4<-V5<-V6<-V7<-V8<-iv1<-iv2<-id<-NULL
  levels<-NA
  levels[is.na(m4.1) & is.na(m4.2)]<-2
  levels[!is.na(m3.1) & !is.na(m3.2)]<-3
  levels[!is.na(m4.1)&!is.na(m4.2)]<-4
  oldoption<-options(contrasts=c("contr.helmert", "contr.poly"))
  oldoption
  on.exit(options(oldoption))

  if (levels=="2"){
    if (!is.null(s)){
      s1.1<-s; s2.1<-s;s1.2<-s;s2.2<-s
      var1<-s^2; var2<-s^2;var3<-s^2;var4<-s^2}
    if (is.null(s)){var1<-s1.1^2; var2<-s2.1^2;var3<-s1.2^2;var4<-s2.2^2}
    if (!is.null(r)){r12<-r;r13<-r;r14<-r;
    r23<-r;r24<-r;
    r34<-r;}
    cov12<-r12*s1.1*s2.1;cov13<-r13*s1.1*s1.2;cov14<-r14*s1.1*s2.2;
    cov23<-r23*s2.1*s1.2;cov24<-r24*s2.1*s2.2;
    cov34<-r34*s2.1*s2.2;
    out <- MASS::mvrnorm(n, mu = c(m1.1,m2.1,m1.2,m2.2),
                   Sigma = matrix(c(var1,cov12,cov13, cov14,
                                    cov12,var2,cov23, cov24,
                                    cov13, cov23,var3,cov34,
                                    cov14, cov24, cov34, var4), ncol = 4),
                   empirical = TRUE)
    out<-as.data.frame(out)
    out<-dplyr::rename(out, y1 = V1, y2 = V2, y3 = V3, y4 = V4)
    out$id <- rep(1:nrow(out))
    out$id<-as.factor(out$id)
    out<-tidyr::gather(out,key="time",value="dv",-id)
    out$time<-as.factor(out$time)
    out$time<-as.numeric(out$time)
    out$iv1<-NA
    out$iv1[out$time==1|out$time==3]<-1
    out$iv1[out$time==2|out$time==4]<-2
    out$iv2<-NA
    out$iv2[out$time==1|out$time==2]<-1
    out$iv2[out$time==3|out$time==4]<-2
    out$iv1<-as.ordered(out$iv1)
    out$iv2<-as.ordered(out$iv2)
    options(contrasts=c("contr.helmert", "contr.poly"))
    data.ab1<-subset(out, iv2==1)
    base<- nlme::lme(dv~1, random = ~1|id/iv1, data=data.ab1,method="ML")
    modelab1<- nlme::lme(dv~iv1, random = ~1|id/iv1, data=data.ab1,method="ML") #A at B1
    lmab1<-stats::anova(base,modelab1)
    dfab1<-lmab1$df[2]-lmab1$df[1]
    lambdalmab1<-lmab1$L.Ratio[2]
    tabledlab1<-stats::qchisq(.95, dfab1)
    powerlab1<-round(1-stats::pchisq(tabledlab1, dfab1, lambdalmab1),3)

    data.ab2<-subset(out, iv2==2)
    base<- nlme::lme(dv~1, random = ~1|id/iv1, data=data.ab2,method="ML")
    modelab2<- nlme::lme(dv~iv1, random = ~1|id/iv1, data=data.ab2,method="ML") #A at B1
    lmab2<-stats::anova(base,modelab2)
    dfab2<-lmab2$df[2]-lmab2$df[1]
    lambdalmab2<-lmab2$L.Ratio[2]
    tabledlab2<-stats::qchisq(.95, dfab2)
    powerlab2<-round(1-stats::pchisq(tabledlab2, dfab2, lambdalmab2),3)


    data.ba1<-subset(out, iv1==1)
    base<- nlme::lme(dv~1, random = ~1|id/iv2, data=data.ba1,method="ML")
    modelba1<- nlme::lme(dv~iv2, random = ~1|id/iv2, data=data.ba1,method="ML") #A at B1
    lmba1<-stats::anova(base,modelba1)
    dfba1<-lmba1$df[2]-lmba1$df[1]
    lambdalmba1<-lmba1$L.Ratio[2]
    tabledlba1<-stats::qchisq(.95, dfba1)
    powerlba1<-round(1-stats::pchisq(tabledlba1, dfba1, lambdalmba1),3)

    data.ba2<-subset(out, iv1==2)
    base<- nlme::lme(dv~1, random = ~1|id/iv2, data=data.ba2,method="ML")
    modelba2<- nlme::lme(dv~iv2, random = ~1|id/iv2, data=data.ba2,method="ML") #A at B1
    lmba2<-stats::anova(base,modelba2)
    dfba2<-lmba2$df[2]-lmba2$df[1]
    lambdalmba2<-lmba2$L.Ratio[2]
    tabledlba2<-stats::qchisq(.95, dfba2)
    powerlba2<-round(1-stats::pchisq(tabledlba2, dfba2, lambdalmba2),3)

    message("Power A at B1 for n = ",n," is ", powerlab1)
    message("Power A at B2 for n = ",n," is ", powerlab2)
    message("Power B at A1 for n = ",n," is ", powerlba1)
    message("Power B at A2 for n = ",n," is ", powerlba2)
    result <- data.frame(matrix(ncol = 5))
    colnames(result) <- c("n", "Power A at B1", "Power A at B2", "Power B at A1", "Power B at A2")
    result[, 1]<-n
    result[, 2]<-powerlab1
    result[, 3]<-powerlab2
    result[, 4]<-powerlba1
    result[, 5]<-powerlba2
    output<-na.omit(result)
    rownames(output)<- c()
  }


  if (levels=="3"){
    if (!is.null(s)){
      s1.1<-s; s2.1<-s;s3.1<-s;s1.2<-s;s2.2<-s;s3.2<-s
      var1<-s^2; var2<-s^2;var3<-s^2;var4<-s^2;var5<-s^2;var6<-s^2}
    if (is.null(s)){var1<-s1.1^2; var2<-s2.1^2;var3<-s3.1^2;var4<-s1.2^2;var5<-s2.2^2; var6<-s3.2^2}
    if (!is.null(r)){r12<-r;r13<-r;r14<-r;r15<-r;r16<-r;
    r23<-r;r24<-r;r25<-r;r26<-r;
    r34<-r;r35<-r;r36<-r;
    r45<-r;r46<-r;
    r56<-r}
    cov12<-r12*s1.1*s2.1;cov13<-r13*s1.1*s3.1;cov14<-r14*s1.1*s1.2;cov15<-r15*s1.1*s2.2;cov16<-r16*s1.1*s3.2;
    cov23<-r23*s2.1*s3.1;cov24<-r24*s2.1*s1.2;cov25<-r25*s2.1*s2.2;cov26<-r26*s2.1*s3.2;
    cov34<-r34*s3.1*s1.2;cov35<-r35*s3.1*s2.2;cov36<-r36*s3.1*s3.2;
    cov45<-r45*s1.2*s2.2;cov46<-r46*s1.2*s3.2;
    cov56<-r56*s2.2*s3.2
    out <- MASS::mvrnorm(n, mu = c(m1.1,m2.1,m3.1,m1.2,m2.2,m3.2),
                   Sigma = matrix(c(var1,cov12,cov13, cov14, cov15, cov16,
                                    cov12,var2,cov23, cov24, cov25, cov26,
                                    cov13, cov23,var3,cov34, cov35, cov36,
                                    cov14, cov24, cov34, var4, cov45, cov46,
                                    cov15, cov25, cov35, cov45, var5, cov56,
                                    cov16, cov26, cov36, cov46, cov56, var6), ncol = 6),
                   empirical = TRUE)
    out<-as.data.frame(out)
    out<-dplyr::rename(out, y1 = V1, y2 = V2, y3 = V3, y4 = V4, y5 = V5, y6 = V6)
    out$id <- rep(1:nrow(out))
    out$id<-as.factor(out$id)
    out<-tidyr::gather(out,key="time",value="dv",-id)
    out$time<-as.factor(out$time)
    out$time<-as.numeric(out$time)
    out$iv1<-NA
    out$iv1[out$time==1|out$time==4]<-1
    out$iv1[out$time==2|out$time==5]<-2
    out$iv1[out$time==3|out$time==6]<-3
    out$iv2<-NA
    out$iv2[out$time==1|out$time==2|out$time==3]<-1
    out$iv2[out$time==4|out$time==5|out$time==6]<-2
    out$iv1<-as.ordered(out$iv1)
    out$iv2<-as.ordered(out$iv2)
    options(contrasts=c("contr.helmert", "contr.poly"))

    data.ab1<-subset(out, iv2==1)
    base<- nlme::lme(dv~1, random = ~1|id/iv1, data=data.ab1,method="ML")
    modelab1<- nlme::lme(dv~iv1, random = ~1|id/iv1, data=data.ab1,method="ML") #A at B1
    lmab1<-stats::anova(base,modelab1)
    dfab1<-lmab1$df[2]-lmab1$df[1]
    lambdalmab1<-lmab1$L.Ratio[2]
    tabledlab1<-stats::qchisq(.95, dfab1)
    powerlab1<-round(1-stats::pchisq(tabledlab1, dfab1, lambdalmab1),3)

    data.ab2<-subset(out, iv2==2)
    base<- nlme::lme(dv~1, random = ~1|id/iv1, data=data.ab2,method="ML")
    modelab2<- nlme::lme(dv~iv1, random = ~1|id/iv1, data=data.ab2,method="ML") #A at B1
    lmab2<-stats::anova(base,modelab2)
    dfab2<-lmab2$df[2]-lmab2$df[1]
    lambdalmab2<-lmab2$L.Ratio[2]
    tabledlab2<-stats::qchisq(.95, dfab2)
    powerlab2<-round(1-stats::pchisq(tabledlab2, dfab2, lambdalmab2),3)


    data.ba1<-subset(out, iv1==1)
    base<- nlme::lme(dv~1, random = ~1|id/iv2, data=data.ba1,method="ML")
    modelba1<- nlme::lme(dv~iv2, random = ~1|id/iv2, data=data.ba1,method="ML") #A at B1
    lmba1<-stats::anova(base,modelba1)
    dfba1<-lmba1$df[2]-lmba1$df[1]
    lambdalmba1<-lmba1$L.Ratio[2]
    tabledlba1<-stats::qchisq(.95, dfba1)
    powerlba1<-round(1-stats::pchisq(tabledlba1, dfba1, lambdalmba1),3)

    data.ba2<-subset(out, iv1==2)
    base<- nlme::lme(dv~1, random = ~1|id/iv2, data=data.ba2,method="ML")
    modelba2<- nlme::lme(dv~iv2, random = ~1|id/iv2, data=data.ba2,method="ML") #A at B1
    lmba2<-stats::anova(base,modelba2)
    dfba2<-lmba2$df[2]-lmba2$df[1]
    lambdalmba2<-lmba2$L.Ratio[2]
    tabledlba2<-stats::qchisq(.95, dfba2)
    powerlba2<-round(1-stats::pchisq(tabledlba2, dfba2, lambdalmba2),3)

    data.ba3<-subset(out, iv1==3)
    base<- nlme::lme(dv~1, random = ~1|id/iv2, data=data.ba3,method="ML")
    modelba3<- nlme::lme(dv~iv2, random = ~1|id/iv2, data=data.ba3,method="ML") #A at B1
    lmba3<-stats::anova(base,modelba3)
    dfba3<-lmba3$df[2]-lmba3$df[1]
    lambdalmba3<-lmba3$L.Ratio[2]
    tabledlba3<-stats::qchisq(.95, dfba3)
    powerlba3<-round(1-stats::pchisq(tabledlba3, dfba3, lambdalmba3),3)

    message("Power A at B1 for n = ",n," is ", powerlab1)
    message("Power A at B2 for n = ",n," is ", powerlab2)
    message("Power B at A1 for n = ",n," is ", powerlba1)
    message("Power B at A2 for n = ",n," is ", powerlba2)
    message("Power B at A3 for n = ",n," is ", powerlba3)
    result <- data.frame(matrix(ncol = 6))
    colnames(result) <- c("n", "Power A at B1", "Power A at B2", "Power B at A1", "Power B at A2", "Power B at A3")
    result[, 1]<-n
    result[, 2]<-powerlab1
    result[, 3]<-powerlab2
    result[, 4]<-powerlba1
    result[, 5]<-powerlba2
    result[, 6]<-powerlba3
    output<-na.omit(result)
    rownames(output)<- c()
    }

  if (levels=="4"){
    if (!is.null(s)){
      s1.1<-s; s2.1<-s;s3.1<-s;s4.1<-s;s1.2<-s;s2.2<-s;s3.2<-s;s4.2<-s
      var1<-s^2; var2<-s^2;var3<-s^2;var4<-s^2;var5<-s^2;var6<-s^2;var7<-s^2;var8<-s^2}
    if (is.null(s)){var1<-s1.1^2; var2<-s2.1^2;var3<-s3.1^2;var4<-s4.1^2;var5<-s1.2^2;var6<-s2.2^2;var7<-s3.2^2;var8<-s4.2^2}
    if (!is.null(r)){r12<-r;r13<-r;r14<-r;r15<-r;r16<-r;r17<-r;r18<-r;r23<-r;r24<-r;r25<-r;r26<-r;r27<-r;r28<-r
    r34<-r;r35<-r;r36<-r;r37<-r;r38<-r;r45<-r;r46<-r;r47<-r;r48<-r;r56<-r;r57<-r;r58<-r
    r67<-r;r68<-r;r78<-r}
    cov12<-r12*s1.1*s2.1;cov13<-r13*s1.1*s3.1;cov14<-r14*s1.1*s4.1;cov15<-r15*s1.1*s1.2;cov16<-r16*s1.1*s2.2;cov17<-r17*s1.1*s3.2;cov18<-r18*s1.1*s4.2
    cov23<-r23*s2.1*s3.1;cov24<-r24*s2.1*s4.1;cov25<-r25*s2.1*s1.2;cov26<-r26*s2.1*s2.2;cov27<-r27*s2.1*s3.2;cov28<-r28*s2.1*s4.2
    cov34<-r34*s3.1*s4.1;cov35<-r35*s3.1*s1.2;cov36<-r36*s3.1*s2.2;cov37<-r37*s3.1*s3.2;cov38<-r38*s3.1*s4.2
    cov45<-r45*s4.1*s1.2;cov46<-r46*s4.1*s2.2;cov47<-r47*s4.1*s3.2;cov48<-r48*s4.1*s4.2
    cov56<-r56*s1.2*s2.2;cov57<-r57*s1.2*s3.2;cov58<-r58*s1.2*s4.2
    cov67<-r67*s2.2*s3.2;cov68<-r68*s2.2*s4.2
    cov78<-r78*s3.2*s4.2
    out <- MASS::mvrnorm(n, mu = c(m1.1,m2.1,m3.1,m4.1,m1.2,m2.2,m3.2,m4.2),
                   Sigma = matrix(c(var1,cov12,cov13, cov14, cov15, cov16, cov17, cov18,
                                    cov12,var2,cov23, cov24, cov25, cov26, cov27, cov28,
                                    cov13, cov23,var3,cov34, cov35, cov36, cov37, cov38,
                                    cov14, cov24, cov34, var4, cov45, cov46, cov47, cov48,
                                    cov15, cov25, cov35, cov45, var5, cov56, cov57, cov58,
                                    cov16, cov26, cov36, cov46, cov56, var6, cov67, cov68,
                                    cov17, cov27, cov37, cov47, cov57, cov67, var7, cov78,
                                    cov18, cov28, cov38, cov48, cov58, cov68, cov78, var8), ncol = 8),
                   empirical = TRUE)
    out<-as.data.frame(out)
    out<-dplyr::rename(out, y1 = V1, y2 = V2, y3 = V3, y4 = V4, y5 = V5, y6 = V6, y7 = V7, y8 = V8)
    out$id <- rep(1:nrow(out))
    out$id<-as.factor(out$id)
    out<-tidyr::gather(out,key="time",value="dv",-id)
    out$time<-as.factor(out$time)
    out$time<-as.numeric(out$time)
    out$iv1<-NA
    out$iv1[out$time==1|out$time==5]<-1
    out$iv1[out$time==2|out$time==6]<-2
    out$iv1[out$time==3|out$time==7]<-3
    out$iv1[out$time==4|out$time==8]<-4
    out$iv2<-NA
    out$iv2[out$time==1|out$time==2|out$time==3|out$time==4]<-1
    out$iv2[out$time==5|out$time==6|out$time==7|out$time==8]<-2
    out$iv1<-as.ordered(out$iv1)
    out$iv2<-as.ordered(out$iv2)
    options(contrasts=c("contr.helmert", "contr.poly"))

    data.ab1<-subset(out, iv2==1)
    base<- nlme::lme(dv~1, random = ~1|id/iv1, data=data.ab1,method="ML")
    modelab1<- nlme::lme(dv~iv1, random = ~1|id/iv1, data=data.ab1,method="ML") #A at B1
    lmab1<-stats::anova(base,modelab1)
    dfab1<-lmab1$df[2]-lmab1$df[1]
    lambdalmab1<-lmab1$L.Ratio[2]
    tabledlab1<-stats::qchisq(.95, dfab1)
    powerlab1<-round(1-stats::pchisq(tabledlab1, dfab1, lambdalmab1),3)

    data.ab2<-subset(out, iv2==2)
    base<- nlme::lme(dv~1, random = ~1|id/iv1, data=data.ab2,method="ML")
    modelab2<- nlme::lme(dv~iv1, random = ~1|id/iv1, data=data.ab2,method="ML") #A at B1
    lmab2<-stats::anova(base,modelab2)
    dfab2<-lmab2$df[2]-lmab2$df[1]
    lambdalmab2<-lmab2$L.Ratio[2]
    tabledlab2<-stats::qchisq(.95, dfab2)
    powerlab2<-round(1-stats::pchisq(tabledlab2, dfab2, lambdalmab2),3)


    data.ba1<-subset(out, iv1==1)
    base<- nlme::lme(dv~1, random = ~1|id/iv2, data=data.ba1,method="ML")
    modelba1<- nlme::lme(dv~iv2, random = ~1|id/iv2, data=data.ba1,method="ML") #A at B1
    lmba1<-stats::anova(base,modelba1)
    dfba1<-lmba1$df[2]-lmba1$df[1]
    lambdalmba1<-lmba1$L.Ratio[2]
    tabledlba1<-stats::qchisq(.95, dfba1)
    powerlba1<-round(1-stats::pchisq(tabledlba1, dfba1, lambdalmba1),3)

    data.ba2<-subset(out, iv1==2)
    base<- nlme::lme(dv~1, random = ~1|id/iv2, data=data.ba2,method="ML")
    modelba2<- nlme::lme(dv~iv2, random = ~1|id/iv2, data=data.ba2,method="ML") #A at B1
    lmba2<-stats::anova(base,modelba2)
    dfba2<-lmba2$df[2]-lmba2$df[1]
    lambdalmba2<-lmba2$L.Ratio[2]
    tabledlba2<-stats::qchisq(.95, dfba2)
    powerlba2<-round(1-stats::pchisq(tabledlba2, dfba2, lambdalmba2),3)

    data.ba3<-subset(out, iv1==3)
    base<- nlme::lme(dv~1, random = ~1|id/iv2, data=data.ba3,method="ML")
    modelba3<- nlme::lme(dv~iv2, random = ~1|id/iv2, data=data.ba3,method="ML") #A at B1
    lmba3<-stats::anova(base,modelba3)
    dfba3<-lmba3$df[2]-lmba3$df[1]
    lambdalmba3<-lmba3$L.Ratio[2]
    tabledlba3<-stats::qchisq(.95, dfba3)
    powerlba3<-round(1-stats::pchisq(tabledlba3, dfba3, lambdalmba3),3)

    data.ba4<-subset(out, iv1==4)
    base<- nlme::lme(dv~1, random = ~1|id/iv2, data=data.ba4,method="ML")
    modelba4<- nlme::lme(dv~iv2, random = ~1|id/iv2, data=data.ba4,method="ML") #A at B1
    lmba4<-stats::anova(base,modelba4)
    dfba4<-lmba4$df[2]-lmba4$df[1]
    lambdalmba4<-lmba4$L.Ratio[2]
    tabledlba4<-stats::qchisq(.95, dfba4)
    powerlba4<-round(1-stats::pchisq(tabledlba4, dfba4, lambdalmba4),3)

    message("Power A at B1 for n = ",n," is ", powerlab1)
    message("Power A at B2 for n = ",n," is ", powerlab2)
    message("Power B at A1 for n = ",n," is ", powerlba1)
    message("Power B at A2 for n = ",n," is ", powerlba2)
    message("Power B at A3 for n = ",n," is ", powerlba3)
    message("Power B at A4 for n = ",n," is ", powerlba4)
    result <- data.frame(matrix(ncol = 7))
    colnames(result) <- c("n", "Power A at B1", "Power A at B2", "Power B at A1", "Power B at A2", "Power B at A3", "Power B at A4")
    result[, 1]<-n
    result[, 2]<-powerlab1
    result[, 3]<-powerlab2
    result[, 4]<-powerlba1
    result[, 5]<-powerlba2
    result[, 6]<-powerlba3
    result[, 7]<-powerlba4
    output<-na.omit(result)
    rownames(output)<- c()
    }
    invisible(output)
    }

