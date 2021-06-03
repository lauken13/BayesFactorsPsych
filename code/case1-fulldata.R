

# First source mixedEffectFunctions.R found at https://osf.io/p74wk/
source(here::here('target-article-code/mixedEffectFunctions.R'))
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

saveRDS(list(mod3_fulldata, mod4_fulldata,mod5_fulldata, mod6_fulldata),file = "results/model-fits-fulldata.rds")

kfolds_grps <- kfold_split_grouped(K= 20,  x = mydat$id)

cv_fulldata_3 <- rstanarm::kfold(mod3_fulldata, K = 20, folds = kfolds_grps)
cv_fulldata_4 <- rstanarm::kfold(mod4_fulldata, K = 20, folds = kfolds_grps)
cv_fulldata_5 <- rstanarm::kfold(mod5_fulldata, K = 20, folds = kfolds_grps)
cv_fulldata_6 <- rstanarm::kfold(mod6_fulldata, K = 20, folds = kfolds_grps)

model_comp_grp <- loo_compare(cv_fulldata_6,cv_fulldata_5,cv_fulldata_4,cv_fulldata_3)

saveRDS(list(cv_fulldata_3,cv_fulldata_4,cv_fulldata_5,cv_fulldata_6, model_comp_grp),file = "results/kfold-grp.rds")

# Using the adjusted se formula from https://avehtari.github.io/modelselection/rats_kcv.html#53_Grouped_K-fold_for_leave-one-group-out

cvgfix <- function(cv, cvidx) {
  groupwise=numeric();
  K <- length(unique(cvidx))
  for (i in 1:K) { groupwise[i]=sum(cv$pointwise[cvidx==i,"elpd_kfold"])}
  cv$pointwise <- cbind(elpd_kfolds=groupwise)
  cv$se_elpd_kfold <- sd(groupwise)*sqrt(K)
  cv$estimates[2] <- cv$se_elpd_kfold
  cv
}


cvgg_3 <- cvgfix(cv_fulldata_3, kfolds_grps)
cvgg_4 <- cvgfix(cv_fulldata_4, kfolds_grps)
cvgg_5 <- cvgfix(cv_fulldata_5, kfolds_grps)
cvgg_6 <- cvgfix(cv_fulldata_6, kfolds_grps)


lcvgg_model_comp <- loo_compare(cvgg_6,cvgg_5,cvgg_4,cvgg_3)


saveRDS(list(cvgg_3, cvgg_4, cvgg_5, cvgg_6, lcvgg_model_comp),file = "results/kfold-grp-correctse.rds")

# What if we did normal loo?

cv_fulldata_3_ind <- rstanarm::loo(mod3_fulldata)
cv_fulldata_4_ind <- rstanarm::loo(mod4_fulldata)
cv_fulldata_5_ind <- rstanarm::loo(mod5_fulldata, cores = 4)
cv_fulldata_6_ind <- rstanarm::loo(mod6_fulldata, cores = 4)

model_comp_ind <- loo_compare(cv_fulldata_6_ind,cv_fulldata_5_ind,cv_fulldata_4_ind,cv_fulldata_3_ind)


saveRDS(list(cv_fulldata_3_ind,cv_fulldata_4_ind,cv_fulldata_5_ind,cv_fulldata_6_ind, model_comp_ind),file = "results/kfold-ind.rds")


