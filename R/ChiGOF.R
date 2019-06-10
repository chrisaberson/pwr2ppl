#'Compute power for an Chi Square Goodness of Fit
#'Takes proportions for up to six group. Alpha is .05 by default, alternative values may be entered by user
#'@param po1 Proportion observed Group 1
#'@param po2 Proportion observed Group 2
#'@param po3 Proportion observed Group 3
#'@param po4 Proportion observed Group 4
#'@param po5 Proportion observed Group 5
#'@param po6 Proportion observed Group 6
#'@param groups Number of groups
#'@param n Total sample size
#'@param alpha Type I error (default is .05)
#'@examples
#'ChiGOF(po1=.25, po2=.20, po3=.20, po4=.35, groups=4,n=100)
#'@return Power for Chi Square Goodness of Fit
#'@export
#'
#'

ChiGOF<-function(groups, po1, po2, po3=NULL, po4=NULL, po5=NULL, po6=NULL, n, alpha=.05)
{

  df<-groups-1
  if(groups==2){
    pe1<-1/groups
    pe2<-1/groups
    sum<-po1+po2
    lambda<-n*((((po1-pe1)^2)/pe1)+(((po2-pe2)^2)/pe2))
  }
  else if(groups==3)
  {
    pe1<-1/groups
    pe2<-1/groups
    pe3<-1/groups
    sum<-po1+po2+po3
    lambda<-n*((((po1-pe1)^2)/pe1)+(((po2-pe2)^2)/pe2)+(((po3-pe3)^2)/pe3))
  }
  else if(groups==4)
  {
    pe1<-1/groups
    pe2<-1/groups
    pe3<-1/groups
    pe4<-1/groups
    sum<-po1+po2+po3+po4
    lambda<-n*((((po1-pe1)^2)/pe1)+(((po2-pe2)^2)/pe2)+(((po3-pe3)^2)/pe3)+(((po4-pe4)^2)/pe4))
  }
  else if(groups==5)
  {
    pe1<-1/groups
    pe2<-1/groups
    pe3<-1/groups
    pe4<-1/groups
    pe5<-1/groups
    sum<-po1+po2+po3+po4+po5
    lambda<-n*((((po1-pe1)^2)/pe1)+(((po2-pe2)^2)/pe2)+(((po3-pe3)^2)/pe3)+(((po4-pe4)^2)/pe4)+(((po5-pe5)^2)/pe5))
  }
  else if(groups==6)
  {
    pe1<-1/groups
    pe2<-1/groups
    pe3<-1/groups
    pe4<-1/groups
    pe5<-1/groups
    pe6<-1/groups
    sum<-po1+po2+po3+po4+po5+po6
    lambda<-n*((((po1-pe1)^2)/pe1)+(((po2-pe2)^2)/pe2)+(((po3-pe3)^2)/pe3)+(((po4-pe4)^2)/pe4)+(((po5-pe5)^2)/pe5)+(((po6-pe6)^2)/pe6))
  }

  tabled<-stats::qchisq(1-alpha, df=df)
  power<-round(1-stats::pchisq(tabled, df=df, lambda),3)
  if(sum!=1.0){stop("Expected proportions must add to 1.0. Check input po values")
  }
  else message("Power for n of ", n, " = ", power)
  result <- data.frame(matrix(ncol = 2))
  colnames(result) <- c( "n","Power")
  result[, 1]<-n
  result[, 2]<-power
  output<-na.omit(result)
  rownames(output)<- c()
  invisible(output)

  }
