# First source mixedEffectFunctions.R found at https://osf.io/p74wk/
source('rep code/mixedEffectFunctions.R')
set.seed(126)
mydat <- myDataSim(fixedEffectCondition = 0.5, randomSlopeConditionSigma = 1, randomInterceptConditionSigma = 0.5,
                   nSub = 20, nCond = 2, nItem = 15)
# Take average y per subject, per condition
mydatAgg <- as.data.frame.table(tapply(mydat$y,list(mydat$id,mydat$cond),mean))
colnames(mydatAgg)=c("id","cond", "y")

library(rstanarm)
library(loo)

#### Full data

mod3_fulldata <- stan_glmer( y ~ (1|id), data = mydat, iter = 4000)
mod4_fulldata <- stan_glmer(y ~ cond + (1|id), data = mydat, iter = 4000)
mod5_fulldata <- stan_glmer(y ~ (1 + cond|id), data = mydat, iter = 4000)
mod6_fulldata <- stan_glmer(y ~ cond + (1 + cond|id), data = mydat, iter = 4000)
kfolds_grps <- kfold_split_grouped(K= 20,  x = mydat$id)

cv_fulldata_3 <- rstanarm::kfold(mod3_fulldata, K = 20, folds = kfolds_grps)
cv_fulldata_4 <- rstanarm::kfold(mod4_fulldata, K = 20, folds = kfolds_grps)
cv_fulldata_5 <- rstanarm::kfold(mod5_fulldata, K = 20, folds = kfolds_grps)
cv_fulldata_6 <- rstanarm::kfold(mod6_fulldata, K = 20, folds = kfolds_grps)

loo_compare(cv_fulldata_6,cv_fulldata_5,cv_fulldata_4,cv_fulldata_3)

##### Aggr Data

mod3_aggrdata <- stan_glmer( y ~ (1|id), data = mydatAgg, iter = 4000)
mod4_aggrdata <- stan_glmer(y ~ cond + (1|id), data = mydatAgg, iter = 4000)
mod5_aggrdata <- stan_glmer(y ~ (1 +cond|id), data = mydatAgg, iter = 4000)
mod6_aggrdata <- stan_glmer(y ~ cond + (1 + cond|id), data = mydatAgg, iter = 4000, adapt_delta = .99)
kfolds_grps <- kfold_split_grouped(K= 20,  x = mydatAgg$id)

cv_aggrdata_3 <- rstanarm::kfold(mod3_aggrdata, K = 20, folds = kfolds_grps)
cv_aggrdata_4 <- rstanarm::kfold(mod4_aggrdata, K = 20, folds = kfolds_grps)
cv_aggrdata_5 <- rstanarm::kfold(mod5_aggrdata, K = 20, folds = kfolds_grps)
cv_aggrdata_6 <- rstanarm::kfold(mod6_aggrdata, K = 20, folds = kfolds_grps)

loo_compare(cv_aggrdata_6,cv_aggrdata_5,cv_aggrdata_4,cv_aggrdata_3)

