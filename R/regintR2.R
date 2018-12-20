#'Compute Power for Regression Interaction (R2 Change Approach)
#'@param R2Mod Full Model R2
#'@param R2Ch Change in R2 Added by Interaction
#'@param mod_pred Full Model Number of Predictors
#'@param ch_pred Change Model Number of Predictors
#'@param nlow starting sample size
#'@param nhigh ending sample size
#'@param by incrimental increase in sample (e.g. nlow = 10, nhigh = 24, by = 2, produces estimates of 10, 12, and 14)
#'@param alpha Type I error (default is .05)
#'@return Power for Regression Interaction (R2 Change Approach)
#'@export
#'
#'


regintR2<-function(R2Mod, R2Ch, mod_pred, ch_pred,nlow, nhigh, by =1, alpha=.05)
{
for(n in seq(nlow,nhigh, by)){
df_denom <- (n - mod_pred)-1
f2 <- R2Ch/(1-R2Mod)
lambda = f2*df_denom
minusalpha<-1-alpha
Ft<-stats::qf(minusalpha, ch_pred, df_denom)
Power<-round(1-stats::pf(Ft, ch_pred, df_denom,lambda),4)
print(paste("Power with n = ", n, "= ", Power))}
}
