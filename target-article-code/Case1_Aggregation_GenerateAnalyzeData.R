# First source mixedEffectFunctions.R found at https://osf.io/p74wk/
source('mixedEffectFunctions.R')

set.seed(126)
mydat <- myDataSim(fixedEffectCondition = 0.5, randomSlopeConditionSigma = 1, randomInterceptConditionSigma = 0.5,
                   nSub = 20, nCond = 2, nItem = 15)
# Take average y per subject, per condition
mydatAgg <- as.data.frame.table(tapply(mydat$y,list(mydat$id,mydat$cond),mean))
colnames(mydatAgg)=c("id","cond", "y")

set.seed(123)
bfs <- generalTestBF(y ~ cond * id, mydat, whichRandom = c('id'), whichModels = 'all')@bayesFactor
bfs["cond + id", 1] - bfs["id", 1] # Model 4 vs Model 3

bfs["cond + id + cond:id", 1] - bfs["id + cond:id", 1] # Model 6 vs Model 5

bfs["cond + id + cond:id", 1] -  bfs["id", 1]  # Model 6 vs Model 3


bfs <- generalTestBF(y ~ cond * id, mydatAgg, whichRandom = c('id'), whichModels = 'all')@bayesFactor
bfs["cond + id", 1] - bfs["id", 1] # Model 4 vs Model 3

bfs["cond + id + cond:id", 1] - bfs["id + cond:id", 1] # Model 6 vs Model 5

bfs["cond + id + cond:id", 1] -  bfs["id", 1]  # Model 6 vs Model 3