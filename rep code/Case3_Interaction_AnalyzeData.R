# First source mixedEffectFunctions.R found at https://osf.io/p74wk/
source('rep code/mixedEffectFunctions.R')
library(xtable)
library(data.table)
set_sum_contrasts()

# Load data found at https://osf.io/p6cj2/
mydat <- read.table("data/Case3_Interaction_FullData.txt", header = TRUE)

mydat <- mydat[mydat$stim_type != "target", ]
mydat$stim_type <- droplevels(mydat$stim_type)

# apply exclusion criteria from Lukacs et al. (2020)
mydat <-  mydat[mydat$rt_start >= 150 & 
                  mydat$too_slow == 0 &
                  mydat$block_number %in% list(4, 5, 7, 8) &
                  mydat$incorrect == 0 &
                  mydat$rt_end <= 1100,] 
levels(mydat$handposition) <- c("Index", "Thumb")
levels(mydat$stim_type) <- c("Irrelevant", "Probe")
mydat <- data.frame(y = mydat$rt_end - mydat$rt_start,
                    id = factor(mydat$subject_id),
                    handpos = mydat$handpos,
                    type = droplevels(mydat$stim_type),
                    trial = factor(mydat$trial_number))
mydatAgg <- as.data.frame.table(tapply(mydat$y,list(mydat$id, mydat$handpos, mydat$type),mean))
colnames(mydatAgg)=c("id","handpos", "type", "y")

bfsObject <- generalTestBF(y  ~ handpos * type + id + type:id + handpos:id + handpos:id:type,
                           whichRandom = c('id', 'type:id', 'handpos:id'), data = mydat, whichModels = "all", multicore = TRUE)
bfsAggObject <- generalTestBF(y  ~ handpos * type + id + type:id + handpos:id+ handpos:id:type,
                              whichRandom = c('id', 'type:id', 'handpos:id'), data = mydatAgg, whichModels = "all", multicore = TRUE)
bfs <- bfsObject@bayesFactor
bfsAgg <- bfsAggObject@bayesFactor

# Minimal RM AOV BF Full Data
bfs["handpos + type + id + handpos:type", 1] - bfs["handpos + type + id", 1]

# RM AOV BF Full Data
bfs["handpos + type + id + handpos:type + type:id + handpos:id", 1] - 
  bfs["handpos + type + id + type:id + handpos:id", 1] 

# Oberauer BF Full Data
bfs["handpos + type + id + handpos:type + type:id + handpos:id + handpos:type:id", 1] - 
  bfs["handpos + type + id + type:id + handpos:id + handpos:type:id", 1]

# Rouder BF Full Data
bfs["handpos + type + id + handpos:type + type:id + handpos:id + handpos:type:id", 1] -  
  bfs["handpos + type + id + type:id + handpos:id", 1] 

# Minimal RM AOV BF Aggregate Data
bfsAgg["handpos + type + id + handpos:type", 1] - bfsAgg["handpos + type + id", 1]

# RM AOV BF Aggregate Data
bfsAgg["handpos + type + id + handpos:type + type:id + handpos:id", 1] - 
  bfsAgg["handpos + type + id + type:id + handpos:id", 1] 

# Oberauer BF Aggregate Data
bfsAgg["handpos + type + id + handpos:type + type:id + handpos:id + handpos:type:id", 1] - 
  bfsAgg["handpos + type + id + type:id + handpos:id + handpos:type:id", 1]

# Rouder BF Aggregate Data
bfsAgg["handpos + type + id + handpos:type + type:id + handpos:id + handpos:type:id", 1] -  
  bfsAgg["handpos + type + id + type:id + handpos:id", 1] 



myModels <- c(
  "handpos + type + id", # mod 1
  "handpos + type + id + handpos:type", #mod 2
  "handpos + type + id + type:id + handpos:id", # mod 3
  "handpos + type + id + handpos:type + type:id + handpos:id", # mod 4
  "handpos + type + id + type:id + handpos:id + handpos:type:id", # mod 5
  "handpos + type + id + handpos:type + type:id + handpos:id + handpos:type:id" # mod 6
)

# Make table with all model comparisons
bfsForTable <- bfs
bfsAggForTable <- bfsAgg
bfMatrix <- matrix(ncol = length(myModels), nrow = length(myModels), dimnames = list(myModels, myModels))

cc <- 0
for (modOne in myModels) {
  skipVec <- 1:cc
  for (modTwo in myModels[-skipVec]) {
    bfMatrix[modOne, modTwo] <- round(bfsForTable[modOne, 1] - bfsForTable[modTwo, 1], 2)
    bfMatrix[modTwo, modOne] <- -round(bfsAggForTable[modOne, 1] - bfsAggForTable[modTwo, 1], 2)
  }
  cc <- cc + 1
}
bfMatrix[1,1] <- 0
colnames(bfMatrix) <- paste0("(", (1:length(myModels)), ")" )

bfMatrix <- cbind(data.frame(Model = paste0("(", (1:(length(myModels))), ")" )),
                  bfMatrix)
bfMatrix[bfMatrix == 0] <- NA

print(xtable(bfMatrix), include.rownames=FALSE)
print(xtable(data.frame(Model = paste0("(", (1:(length(myModels))), ")" ),
                        Name = c( myModels))), include.rownames=FALSE)