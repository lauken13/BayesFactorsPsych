library(BayesFactor)
library(dplyr)

sim_case3_quantiles <- function(fixedEffect = 0, randomSlopeSigma = 0) {
  
  I <- 20
  J <- 2
  K <- 100
  batch <- round(K / round(exp(seq(0, log(K), length.out = 7))))
  
  subject <- rep(1:I, each = J)
  cond <- rep(c(-0.5, 0.5), I)
  
  alpha <- sample(qnorm(ppoints(I), sd = 0.5))
  theta <- sample(qnorm(ppoints(I), mean = fixedEffect, sd = randomSlopeSigma))
  y <- alpha[subject] + cond * theta[subject] + rnorm(I*J, mean = 0, sd = 1/sqrt(K))
  
  case2_data <- list()
  
  for(i in batch) {
    if(K/i == 1) {
      y_i <- y
    } else {
      y_i <- as.vector(sapply(y, function(x) x + scale(qnorm(ppoints(K/i), sd = 1)) * sqrt(1/i)))
    }
    
    case2_data[[paste0("trials_", K/i)]] <- data.frame(
      id = factor(rep(subject, each = K/i))
      , cond = factor(rep(cond, each = K/i), levels = c(-0.5, 0.5), labels = c("a", "b"))
      , y = y_i
    )
  }
  
  
  case2_subject_descriptives <- lapply(
    case2_data
    , function(x) {
      x %>%
        group_by(id, cond) %>%
        summarize(m = mean(y), sd = sd(y), n = length(y), se = sd/sqrt(n), .groups = "keep")
    }
  )
  
  # Same participant means and standard errors?
  same_p_means <- sapply(
    case2_subject_descriptives[-c(1:2)]
    , function(x) {
      all.equal(case2_subject_descriptives[[2]][, c("id", "cond", "m", "se")], x[, c("id","cond", "m", "se")])
    }
  )
  stopifnot(all(same_p_means))
  
  case2_condition_descriptives <- lapply(
    case2_data
    , function(x) {
      x %>%
        group_by(id, cond) %>%
        summarize(y = mean(y), .groups = "keep") %>%
        group_by(cond) %>%
        summarize(m = mean(y), sd = sd(y), n = length(y), se = sd/sqrt(n), .groups = "keep")
    }
  )
  
  # Same condition means and standard errors?
  same_c_means <- sapply(
    case2_condition_descriptives[-1]
    , function(x) {
      all.equal(case2_condition_descriptives[[1]][, c("cond", "m", "se")], x[, c("cond", "m", "se")])
    }
  )
  stopifnot(all(same_c_means))
  
  return(case2_data)
}


set.seed(1234)
nTimes <- 10
noiseSig <- 1
nItems <- 100
nSub <- 20
nTimes <- 10

counter <- 1

allFixed <- c(0, 0.3, 0.5, 1)
allSlope <- c(0, 0.3, 0.5)
allInter <- c(0.5)
# allInter <- allSlope <- allFixed <- 0

thisFixedEffect <- c(0.5)
thisRandomSlopeSigma <- c(0.5)
thisRandomInterceptSigma <- c(0.5)
cexAxis <- 1.3
cexLab <- 1.5

showAdjust <- FALSE
myCount <- 1

bfList <- list()

for (thisFixedEffect in c(0, 0.5)) {
  for (thisRandomSlopeSigma in c(0, 0.5)) {
    if (thisFixedEffect == 0 && thisRandomSlopeSigma == 0.5) next
    set.seed(1234)
    
    allDat <- sim_case3_quantiles(fixedEffect = thisFixedEffect,  
                                  randomSlopeSigma = thisRandomSlopeSigma)
    
    itemIndicator <-   names(allDat)
    allTitles <- rev(c("Full Data", "50 Items", "20 items", "10 Items", "5 Items", "2 Items",  "Fully Aggr."))
    
    allBF <- matrix(ncol = 3, nrow = length(itemIndicator))
    
    for (i in 1:length(itemIndicator)) {
      
      smallDat <- allDat[[i]]
      
      sigErr <-   1 / sqrt(100 / nItems[i])
      
      bfs <- generalTestBF(y ~ cond * id, smallDat, whichRandom = c('id'), whichModels = 'all', multicore = F,
                           rscaleFixed = 0.5, rscaleRandom = 1)@bayesFactor
      
      aovBF <- round(bfs["cond + id", 1] - bfs["id", 1], 3) # RM ANOVA
      klausBF <- round(bfs["cond + id + cond:id", 1] - bfs["id + cond:id", 1], 3) # Klaus
      rouderBF <- round(bfs["cond + id + cond:id", 1] -  bfs["id", 1], 3) # Rouder
      
      
      allBF[i, ] <- log(exp(c(aovBF, klausBF, rouderBF)))
    }
    
    colnames(allBF) <- c("rmAOV", "Oberauer", "Rouder")
    rownames(allBF) <- allTitles
    
    bfList[[myCount]] <- allBF
    myCount <- myCount + 1
  }
}
