library(MASS)
library(lme4)
library(BayesFactor)
library(afex)
library(lmerTest)

myDataSim <- function(fixedEffectItem = 0, randomSlopeItemSigma = 0, randomInterceptItemSigma = 0, randomRhoSlopeInterceptItem = 0,
                      fixedEffectCondition = 0, randomSlopeConditionSigma = 0, randomInterceptConditionSigma = 0, randomRhoSlopeInterceptCondition = 0,
                      errorSigma = 1, errorSigBetween = 0, nSub = 30, nItem = 2, nCond = 2) {
  
  I <- nSub # n per item per condition
  id <- 1:I
  intercept <- 0
  
  M <- nItem # items
  J <- nCond # conditions
  grandN <- I * J * M
  cond <- rep(1:J, each = M*I)
  item <- rep(1:M, J*I)
  id <- rep(rep(1:I, each = M), J)
  fixedEffectConditionVector <- seq(-1/nCond, 1/nCond, length.out = nCond) 
  fixedEffectItemVector <- seq(-1/nItem, 1/nItem, length.out = nItem) 
  
  # setup random effects  
  randomEffectMatItem <- diag(c(randomInterceptItemSigma, randomSlopeItemSigma))
  randomEffectMatItem[1, 2] <- randomEffectMatItem[2, 1] <- randomRhoSlopeInterceptItem * sqrt(randomInterceptItemSigma) * sqrt(randomSlopeItemSigma)
  
  randomEffectMatCondition <- diag(c(randomInterceptConditionSigma, randomSlopeConditionSigma))
  randomEffectMatCondition[1, 2] <- randomEffectMatCondition[2, 1] <- randomRhoSlopeInterceptCondition * sqrt(randomInterceptConditionSigma) * sqrt(randomSlopeConditionSigma)
  
  err <- rnorm(grandN, 0, errorSigma)

  randomEffectsItem <- MASS::mvrnorm(n = I, mu = c(0,0) , Sigma = randomEffectMatItem)
  randomEffectsCondition <- MASS::mvrnorm(n = I, mu = c(0,0) , Sigma = randomEffectMatCondition)
  colnames(randomEffectsItem) <- colnames(randomEffectsCondition) <- c("intercept", "slope")
  
  y <- intercept + randomEffectsItem[, "intercept"][id] + randomEffectsCondition[, "intercept"][id] +
    (fixedEffectCondition + randomEffectsCondition[, "slope"][id]) * fixedEffectConditionVector[cond] +
    (fixedEffectItem + randomEffectsItem[, "slope"][id]) * fixedEffectItemVector[item] + err
  
  mydat <- data.frame(y = y,
                      cond = as.factor(cond),
                      item = as.factor(item),
                      id = as.factor(id))
  return(mydat)
}


plotBoxplotLines <- function(mydat, xlab = "Condition", ylab = "y", myCex = 1.3, myCexLab  = 1.3,
                             mylwd = 2, myAlpha = 0.4) {
  par(mfrow = c(1, 2), cex = myCex, cex.lab  = myCexLab, cex.main = myCexLab)
  mydatAgg <- as.data.frame.table(tapply(mydat$y,list(mydat$id,mydat$cond),mean))
  colnames(mydatAgg)=c("id","cond", "y")
  
  nCond <- nlevels(mydat$cond)
  nItem <- nlevels(mydat$item)
  nSub <- nlevels(mydat$id)
  palette(rainbow(nSub))
  allCols <- adjustcolor(palette(), alpha.f = myAlpha)
  plot(mydat$cond, mydat$y, bty = "n", las = 1, type ="n", xlab = xlab, ylab = ylab, main = "Full Data")
  
  points(as.numeric(mydat$cond), mydat$y)
  for (i in 1:nSub) {
    ss <- subset(mydat, id == levels(mydat$id)[i])
    for (j in 1:nItem) {
      sss <- subset(ss, item == levels(mydat$item)[j])
      lines(1:nCond, sss$y, col = allCols[i], lwd = mylwd)
    }
  }
  
  plot(mydat$y ~ mydat$cond, bty = "n", las = 1, type ="n", 
          xlab = xlab, ylab = ylab, main = "Aggregated Data")
  boxplot(mydatAgg$y ~ mydatAgg$cond, add = TRUE, las = 1)
  points(as.numeric(mydatAgg$cond), mydatAgg$y)
  for (i in 1:nSub) {
    ss <- subset(mydatAgg, id == levels(mydat$id)[i])
    lines(1:nCond, tapply(ss$y, ss$cond, mean), col = allCols[i], lwd = mylwd)
  }
}



analyzeSamples <- function(nSub, nItem, nCond, fixEffect, errorSigma, randInterSigma, randSlopeSigma) {
  afex::set_sum_contrasts()
  allBfAgg <- allBf <- AIClmer <- AIClmerAgg <- data.frame(rInter = NA, rInterCond = NA, rSlopInter = NA, full = NA)

  mydat <- myDataSim(nSub = nSub, nItem = nItem, nCond = nCond, fixedEffectCondition = fixEffect, 
                     errorSigma = errorSigma, randomInterceptConditionSigma = randInterSigma,
                     randomSlopeConditionSigma = randSlopeSigma)
  mydatAgg <- as.data.frame.table(tapply(mydat$y,list(mydat$id,mydat$cond),mean)) # same as mydatUnivar
  colnames(mydatAgg)=c("id","cond", "y")
  
  bfs <- BayesFactor::generalTestBF(y ~ cond*id, mydat, whichRandom = c('id'), whichModels = 'all')
  idbf <- bfs@bayesFactor$bf[rownames(bfs@bayesFactor) == "id"]
  condbf <- bfs@bayesFactor$bf[rownames(bfs@bayesFactor) == "cond + id"]
  ranSlopeBf <- bfs@bayesFactor$bf[rownames(bfs@bayesFactor) == "id + cond:id"]
  fullBf <- bfs@bayesFactor$bf[rownames(bfs@bayesFactor) == "cond + id + cond:id"]
  allBf[1, ] <- c(idbf, condbf, ranSlopeBf, fullBf) 
  
  bfs <- BayesFactor::generalTestBF(y ~ cond*id, mydatAgg, whichRandom = c('id'), whichModels = 'all')
  idbf <- bfs@bayesFactor$bf[rownames(bfs@bayesFactor) == "id"]
  condbf <- bfs@bayesFactor$bf[rownames(bfs@bayesFactor) == "cond + id"]
  ranSlopeBf <- bfs@bayesFactor$bf[rownames(bfs@bayesFactor) == "id + cond:id"]
  fullBf <- bfs@bayesFactor$bf[rownames(bfs@bayesFactor) == "cond + id + cond:id"]
  allBfAgg[1, ] <- c(idbf, condbf, ranSlopeBf, fullBf) 
  
  fullMod <- lmerTest::lmer(y ~ cond + (cond | id), data = mydat)
  ranSlope <- lmerTest::lmer(y ~ (cond | id), data = mydat)
  condId <- lmerTest::lmer(y ~ cond + (1 | id), data = mydat)
  idOnly <- lmerTest::lmer(y ~ (1 | id), data = mydat)
  AIClmer[1, ] <- AIC(idOnly, condId, ranSlope, fullMod)[,2]
  
  baseModPval <- summary(condId)$coefficients["cond2", "Pr(>|t|)"]
  fullModPval <- summary(fullMod)$coefficients["cond2", "Pr(>|t|)"]
  
  condId <- lmerTest::lmer(y ~ cond + (1 | id), data = mydatAgg)
  idOnly <- lmerTest::lmer(y ~ (1 | id), data = mydatAgg)
  AIClmerAgg[1, ] <- c(AIC(idOnly, condId)[,2], NA, NA)
  
  baseModPvalAgg <- summary(condId)$coefficients["cond2", "Pr(>|t|)"]
  
  return(list(allBfAgg = allBfAgg, allBf = allBf, AIClmer = AIClmer, AIClmerAgg = AIClmerAgg, 
              fullModPval = fullModPval, baseModPval = baseModPval, baseModPvalAgg = baseModPvalAgg))
}

