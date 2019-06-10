#'Compute Power for Logistic Regression with a Single Categorical Predictor
#'@param p0 Probability of a Desirable Outcome in the Control Condition
#'@param p1 Probability of a Desirable Outcome in the Treatment Condition
#'@param prop Proportion in the Treatment Condition
#'@param alpha Type I error (default is .05)
#'@param power Desired Power
#'@param R2 How Well Predictor of Interest is Explained by Other Predictors (default is 0)
#'@examples
#'LRcat(p0=.137,p1=.611,prop =.689,power=.95)
#'@return Power for Logistic Regression with a Single Categorical Predictor
#'@export
#'
#'



LRcat<-function(p0=NULL, p1=NULL, prop=.50, alpha=.05, power, R2=.00)
{
    R<-prop
    pbar<-((1-R)*p0)+(R*p1)
    zalpha<-stats::qnorm(1-alpha/2)
    zbeta<-stats::qnorm(power)
    num1<-zalpha*(((pbar*(1-pbar))/R)^.5)
    num2<-zbeta*(((p0*(1-p0))+((p1*(1-p1)*(1-R)))/R))^.5
    den<-((p0-p1)^2)*(1-R)
    n<-((num1+num2)^2/den)/(1-R2)
    nprint<-ceiling(n)
    OR<-round((p1/(1-p1))/(p0/(1-p0)),3)
    message("Sample Size = ", nprint, " for Odds Ratio = ", OR)
    result <- data.frame(matrix(ncol = 3))
    colnames(result) <- c("n", "OR","power")
    result[, 1]<-nprint
    result[, 2]<-OR
    result[, 3]<-power
    output<-na.omit(result)
    rownames(output)<- c()
    invisible(output)
    }

