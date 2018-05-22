#'Compute power for a Two Factor Within Subjects ANOVA with up to two by four levels.
#'Takes means, sds, and sample sizes for each group. Alpha is .05 by default, alterative values may be entered by user
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
#'@return Power for the Two Factor Within Subjects ANOVA
#'@export

win2F<-function(m1.1,m2.1,m3.1=NA,m4.1=NA,m1.2,m2.2,m3.2=NA,m4.2=NA,
                s1.1=NA,s2.1=NA,s3.1=NA,s4.1=NA,s1.2=NA,s2.2=NA,s3.2=NA,s4.2=NA,
                r12=NULL, r13=NULL, r14=NULL, r15=NULL, r16=NULL, r17=NULL, r18=NULL,
                r23=NULL, r24=NULL, r25=NULL, r26=NULL, r27=NULL, r28=NULL,
                r34=NULL, r35=NULL, r36=NULL, r37=NULL, r38=NULL,
                r45=NULL, r46=NULL, r47=NULL, r48=NULL,
                r56=NULL, r57=NULL, r58=NULL,
                r67=NULL, r68=NULL,
                r78=NULL, r=NULL, s = NULL, n, alpha=.05)

{

  levels<-NA
  levels[is.na(m4.1) & is.na(m3.1)]<-2
  levels[is.na(m4.1) & !is.na(m3.1)]<-3
  levels[!is.na(m4.1)]<-4

  if (levels=="2"){
    if (!is.null(s)){
      s1.1<-s; s2.1<-s;s1.2<-s;s2.2<-s
      var1<-s^2; var2<-s^2;var3<-s^2;var4<-s^2}
    if (is.null(s)){var1<-s1.1^2; var2<-s2.1^2;var3<-s1.2^2;var4<-s2.2^2}
    r12<-r;r13<-r;r14<-r;
    r23<-r;r24<-r;
    r34<-r;
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

    model<-ez::ezANOVA(data=out, dv=.(dv), wid=.(id), within = .(iv1,iv2), type=3, detailed=TRUE)
    dfA<-model$ANOVA$DFn[2]
    dfB<-model$ANOVA$DFn[3]
    dfAB<-model$ANOVA$DFn[4]
    dfWA<-model$ANOVA$DFd[2]
    dfWB<-model$ANOVA$DFd[3]
    dfWAB<-model$ANOVA$DFd[4]
    SSA<-model$ANOVA$SSn[2]
    SSB<-model$ANOVA$SSn[3]
    SSAB<-model$ANOVA$SSn[4]
    SSWA<-model$ANOVA$SSd[2]
    SSWB<-model$ANOVA$SSd[3]
    SSWAB<-model$ANOVA$SSd[4]
    eta2A<-SSA/(SSA+SSWA)
    eta2B<-SSB/(SSB+SSWB)
    eta2AB<-SSAB/(SSAB+SSWAB)
    f2A<-eta2A/(1-eta2A)
    f2B<-eta2B/(1-eta2B)
    f2AB<-eta2AB/(1-eta2AB)
    lambdaA<-f2A*dfWA
    lambdaB<-f2B*dfWB
    lambdaAB<-f2AB*dfWAB
    minusalpha<-1-alpha
    FtA<-qf(minusalpha, dfA, dfWA)
    FtB<-qf(minusalpha, dfB, dfWB)
    FtAB<-qf(minusalpha, dfAB, dfWAB)
    powerA<-round(1-pf(FtA, dfA,dfWA,lambdaA),3)
    powerB<-round(1-pf(FtB, dfB,dfWB,lambdaB),3)
    powerAB<-round(1-pf(FtAB, dfAB,dfWAB,lambdaAB),3)
    {print(paste("Power Factor A (Unadjusted) for n =",n,"=", powerA))}
    {print(paste("Power Factor B (Unadjusted) for n =",n,"=", powerB))}
    {print(paste("Power Factor AB (Unadjusted) for n =",n,"=", powerAB))}
    {print(paste("Both Factors Have 2 levels - There is no adjustment when levels = 2"))}

  }

    if (levels=="3"){
      if (!is.null(s)){
        s1.1<-s; s2.1<-s;s3.1<-s;s1.2<-s;s2.2<-s;s3.2<-s
        var1<-s^2; var2<-s^2;var3<-s^2;var4<-s^2;var5<-s^2;var6<-s^2}
      if (is.null(s)){var1<-s1.1^2; var2<-s2.1^2;var3<-s3.1^2;var4<-s1.2^2;var5<-s2.2^2; var6<-s3.2^2}
      r12<-r;r13<-r;r14<-r;r15<-r;r16<-r;
      r23<-r;r24<-r;r25<-r;r26<-r;
      r34<-r;r35<-r;r36<-r;
      r45<-r;r46<-r;
      r56<-r
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
      model<-ez::ezANOVA(data=out, dv=.(dv), wid=.(id), within = .(iv1,iv2), type=3, detailed=TRUE)
      dfA<-model$ANOVA$DFn[2]
      dfB<-model$ANOVA$DFn[3]
      dfAB<-model$ANOVA$DFn[4]
      dfWA<-model$ANOVA$DFd[2]
      dfWB<-model$ANOVA$DFd[3]
      dfWAB<-model$ANOVA$DFd[4]
      SSA<-model$ANOVA$SSn[2]
      SSB<-model$ANOVA$SSn[3]
      SSAB<-model$ANOVA$SSn[4]
      SSWA<-model$ANOVA$SSd[2]
      SSWB<-model$ANOVA$SSd[3]
      SSWAB<-model$ANOVA$SSd[4]
      eta2A<-SSA/(SSA+SSWA)
      eta2B<-SSB/(SSB+SSWB)
      eta2AB<-SSAB/(SSAB+SSWAB)
      f2A<-eta2A/(1-eta2A)
      f2B<-eta2B/(1-eta2B)
      f2AB<-eta2AB/(1-eta2AB)
      lambdaA<-f2A*dfWA
      lambdaB<-f2B*dfWB
      lambdaAB<-f2AB*dfWAB
      minusalpha<-1-alpha
      FtA<-qf(minusalpha, dfA, dfWA)
      FtB<-qf(minusalpha, dfB, dfWB)
      FtAB<-qf(minusalpha, dfAB, dfWAB)
      powerA<-round(1-pf(FtA, dfA,dfWA,lambdaA),3)
      powerB<-round(1-pf(FtB, dfB,dfWB,lambdaB),3)
      powerAB<-round(1-pf(FtAB, dfAB,dfWAB,lambdaAB),3)
      ggeA<-round(model$`Sphericity Corrections`$GGe[1],3)
      ggeAB<-round(model$`Sphericity Corrections`$GGe[2],3)
      hfeA<-round(model$`Sphericity Corrections`$HFe[1],3)
      hfeAB<-round(model$`Sphericity Corrections`$HFe[2],3)
      hfeA[hfeA>1]<-1
      hfeAB[hfeAB>1]<-1
      ggdfA<-ggeA*dfA
      ggdfAB<-ggeAB*dfAB
      ggdfWA<-ggeA*dfWA
      ggdfWAB<-ggeAB*dfWAB
      hfdfA<-hfeA*dfA
      hfdfAB<-hfeAB*dfAB
      hfdfWA<-hfeA*dfWA
      hfdfWAB<-hfeAB*dfWAB
      lambdaggA<-f2A*ggdfWA
      lambdaggAB<-f2AB*ggdfWAB
      lambdahfA<-f2A*hfdfWA
      lambdahfAB<-f2AB*hfdfWAB
      FtggA<-qf(minusalpha, ggdfA, ggdfWA)
      FtggAB<-qf(minusalpha, ggdfAB, ggdfWAB)
      FthfA<-qf(minusalpha, hfdfA, hfdfWA)
      FthfAB<-qf(minusalpha, hfdfAB, hfdfWAB)
      powerggA<-round(1-pf(FtggA, ggdfA,ggdfWA,lambdaggA),3)
      powerggAB<-round(1-pf(FtggAB, ggdfAB,ggdfWAB,lambdaggAB),3)
      powerhfA<-round(1-pf(FthfA, hfdfA,hfdfWA,lambdahfA),3)
      powerhfAB<-round(1-pf(FthfAB, hfdfAB,hfdfWAB,lambdahfAB),3)
      {print(paste("Power Factor A (Unadjusted) for n =",n,"=", powerA))}
      {print(paste("Power Factor A H-F Adjusted (Epsilon = ",hfeA ,") for n =",n, "=", powerhfA))}
      {print(paste("Power Factor A G-G Adjusted (Epsilon = ", ggeA,") for n =",n, "=", powerggA))}
      {print(paste("Power Factor B (Unadjusted) for n =",n,"=", powerB))}
      {print(paste("Power Factor B Adjusted - There is no adjustment when levels = 2"))}
      {print(paste("Power Factor AB (Unadjusted) for n =",n,"=", powerAB))}
      {print(paste("Power Factor AB H-F Adjusted (Epsilon = ",hfeAB ,") for n =",n, "=", powerhfAB))}
      {print(paste("Power Factor AB G-G Adjusted (Epsilon = ", ggeAB,") for n =",n, "=", powerggAB))}
    }

  if (levels=="4"){
    if (!is.null(s)){
      s1.1<-s; s2.1<-s;s3.1<-s;s4.1<-s;s1.2<-s;s2.2<-s;s3.2<-s;s4.2<-s
      var1<-s^2; var2<-s^2;var3<-s^2;var4<-s^2;var5<-s^2;var6<-s^2;var7<-s^2;var8<-s^2}
    if (is.null(s)){var1<-s1.1^2; var2<-s2.1^2;var3<-s3.1^2;var4<-s4.1^2;var5<-s1.2^2;var6<-s2.2^2;var7<-s3.2^2;var8<-s4.2^2}
    r12<-r;r13<-r;r14<-r;r15<-r;r16<-r;r17<-r;r18<-r;r23<-r;r24<-r;r25<-r;r26<-r;r27<-r;r28<-r
    r34<-r;r35<-r;r36<-r;r37<-r;r38<-r;r45<-r;r46<-r;r47<-r;r48<-r;r56<-r;r57<-r;r58<-r
    r67<-r;r68<-r;r78<-r
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
    model<-ez::ezANOVA(data=out, dv=.(dv), wid=.(id), within = .(iv1,iv2), type=3, detailed=TRUE)
    dfA<-model$ANOVA$DFn[2]
    dfB<-model$ANOVA$DFn[3]
    dfAB<-model$ANOVA$DFn[4]
    dfWA<-model$ANOVA$DFd[2]
    dfWB<-model$ANOVA$DFd[3]
    dfWAB<-model$ANOVA$DFd[4]
    SSA<-model$ANOVA$SSn[2]
    SSB<-model$ANOVA$SSn[3]
    SSAB<-model$ANOVA$SSn[4]
    SSWA<-model$ANOVA$SSd[2]
    SSWB<-model$ANOVA$SSd[3]
    SSWAB<-model$ANOVA$SSd[4]
    eta2A<-SSA/(SSA+SSWA)
    eta2B<-SSB/(SSB+SSWB)
    eta2AB<-SSAB/(SSAB+SSWAB)
    f2A<-eta2A/(1-eta2A)
    f2B<-eta2B/(1-eta2B)
    f2AB<-eta2AB/(1-eta2AB)
    lambdaA<-f2A*dfWA
    lambdaB<-f2B*dfWB
    lambdaAB<-f2AB*dfWAB
    minusalpha<-1-alpha
    FtA<-qf(minusalpha, dfA, dfWA)
    FtB<-qf(minusalpha, dfB, dfWB)
    FtAB<-qf(minusalpha, dfAB, dfWAB)
    powerA<-round(1-pf(FtA, dfA,dfWA,lambdaA),3)
    powerB<-round(1-pf(FtB, dfB,dfWB,lambdaB),3)
    powerAB<-round(1-pf(FtAB, dfAB,dfWAB,lambdaAB),3)
    ggeA<-round(model$`Sphericity Corrections`$GGe[1],3)
    ggeAB<-round(model$`Sphericity Corrections`$GGe[2],3)
    hfeA<-round(model$`Sphericity Corrections`$HFe[1],3)
    hfeAB<-round(model$`Sphericity Corrections`$HFe[2],3)
    hfeA[hfeA>1]<-1
    hfeAB[hfeAB>1]<-1
    ggdfA<-ggeA*dfA
    ggdfAB<-ggeAB*dfAB
    ggdfWA<-ggeA*dfWA
    ggdfWAB<-ggeAB*dfWAB
    hfdfA<-hfeA*dfA
    hfdfAB<-hfeAB*dfAB
    hfdfWA<-hfeA*dfWA
    hfdfWAB<-hfeAB*dfWAB
    lambdaggA<-f2A*ggdfWA
    lambdaggAB<-f2AB*ggdfWAB
    lambdahfA<-f2A*hfdfWA
    lambdahfAB<-f2AB*hfdfWAB
    FtggA<-qf(minusalpha, ggdfA, ggdfWA)
    FtggAB<-qf(minusalpha, ggdfAB, ggdfWAB)
    FthfA<-qf(minusalpha, hfdfA, hfdfWA)
    FthfAB<-qf(minusalpha, hfdfAB, hfdfWAB)
    powerggA<-round(1-pf(FtggA, ggdfA,ggdfWA,lambdaggA),3)
    powerggAB<-round(1-pf(FtggAB, ggdfAB,ggdfWAB,lambdaggAB),3)
    powerhfA<-round(1-pf(FthfA, hfdfA,hfdfWA,lambdahfA),3)
    powerhfAB<-round(1-pf(FthfAB, hfdfAB,hfdfWAB,lambdahfAB),3)
    {print(paste("Power Factor A (Unadjusted) for n =",n,"=", powerA))}
    {print(paste("Power Factor A H-F Adjusted (Epsilon = ",hfeA ,") for n =",n, "=", powerhfA))}
    {print(paste("Power Factor A G-G Adjusted (Epsilon = ", ggeA,") for n =",n, "=", powerggA))}
    {print(paste("Power Factor B (Unadjusted) for n =",n,"=", powerB))}
    {print(paste("Power Factor B Adjusted - There is no adjustment when levels = 2"))}
    {print(paste("Power Factor AB (Unadjusted) for n =",n,"=", powerAB))}
    {print(paste("Power Factor AB H-F Adjusted (Epsilon = ",hfeAB ,") for n =",n, "=", powerhfAB))}
    {print(paste("Power Factor AB G-G Adjusted (Epsilon = ", ggeAB,") for n =",n, "=", powerggAB))}
  }}


