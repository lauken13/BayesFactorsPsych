
setwd('/mnt/lustre/projects/Mona0070/lkennedy/BFCogSci')

library(BayesFactor)
library(dplyr)
library(rstanarm)
library(loo)

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

gpcvlist <- list()
indcvlist <- list()

for (thisFixedEffect in c(0, 0.5)) {
  for (thisRandomSlopeSigma in c(0, 0.5)) {
    if (thisFixedEffect == 0 && thisRandomSlopeSigma == 0.5) next
    set.seed(1234)
    print(myCount)
    allDat <- sim_case3_quantiles(fixedEffect = thisFixedEffect,  
                                  randomSlopeSigma = thisRandomSlopeSigma)
    
    itemIndicator <-   names(allDat)
    allTitles <- rev(c("Full Data", "50 Items", "20 items", "10 Items", "5 Items", "2 Items",  "Fully Aggr."))
    
    allcv <- data.frame(matrix(ncol = 8,nrow=0)) # 4 candidate models
    
    colnames(allcv) <- c('elpd_diff^mod4','elpd_diff^mod5','elpd_diff^mod6','elpd_diff^mod3','se_diff^mod4','se_diff^mod6','se_diff^mod5','se_diff^mod3')
    
    indcv <- data.frame(matrix(ncol = 8,nrow=0)) # 4 candidate models
    
    colnames(indcv) <- c('elpd_diff^mod4','elpd_diff^mod5','elpd_diff^mod6','elpd_diff^mod3','se_diff^mod4','se_diff^mod6','se_diff^mod5','se_diff^mod3')
    
    
    for (i in 1:length(itemIndicator)) {
      print(i)
      smallDat <- allDat[[i]]
      
      sigErr <-   1 / sqrt(100 / nItems[i])
      
      mod3 <- stan_glmer( y ~ (1|id), data = smallDat, iter = 4000, cores = 4)
      mod4 <- stan_glmer(y ~ cond + (1|id), data = smallDat, iter = 4000, cores = 4)
      mod5 <- stan_glmer(y ~ (1 + cond|id), data = smallDat, iter = 4000, cores = 4)
      mod6 <- stan_glmer(y ~ cond + (1 + cond|id), data = smallDat, iter = 4000, cores = 4)
      
      saveRDS(list(mod3, mod4,mod5, mod6),file = paste0("results/model-fits-",allTitles[i],
                                                        "-FE-",thisFixedEffect,"-SigmaSlope-",thisFixedEffect,".rds"))
      
      kfolds_grps <- kfold_split_grouped(K= 20,  x = smallDat$id)
      
      cvgfix <- function(cv, cvidx) {
        groupwise=numeric();
        K <- length(unique(cvidx))
        for (i in 1:K) { groupwise[i]=sum(cv$pointwise[cvidx==i,"elpd_kfold"])}
        cv$pointwise <- cbind(elpd_kfolds=groupwise)
        cv$se_elpd_kfold <- sd(groupwise)*sqrt(K)
        cv$estimates[2] <- cv$se_elpd_kfold
        cv
      }
      
      cv_3_gp <- rstanarm::kfold(mod3,  K=20, folds = kfolds_grps, cores = 20)
      cv_4_gp <- rstanarm::kfold(mod4,  K=20, folds = kfolds_grps, cores = 20)
      cv_5_gp <- rstanarm::kfold(mod5,  K=20, folds = kfolds_grps, cores = 20)
      cv_6_gp <- rstanarm::kfold(mod6,  K=20, folds = kfolds_grps, cores = 20)
      
      cv_3_gp_fix <- cvgfix(cv_3_gp,kfolds_grps)
      cv_4_gp_fix <- cvgfix(cv_4_gp,kfolds_grps)
      cv_5_gp_fix <- cvgfix(cv_5_gp,kfolds_grps)
      cv_6_gp_fix <- cvgfix(cv_6_gp,kfolds_grps)
      
      model_comp <- data.frame(elpd_diff = c(sum(cv_6_gp_fix$pointwise-cv_6_gp_fix$pointwise),
                                             sum(cv_5_gp_fix$pointwise-cv_6_gp_fix$pointwise),
                                             sum(cv_4_gp_fix$pointwise-cv_6_gp_fix$pointwise),
                                             sum(cv_3_gp_fix$pointwise-cv_6_gp_fix$pointwise)),
                               se_diff =  c(sqrt(20)*sd(cv_6_gp_fix$pointwise-cv_6_gp_fix$pointwise),
                                            sqrt(20)*sd(cv_5_gp_fix$pointwise-cv_6_gp_fix$pointwise),
                                            sqrt(20)*sd(cv_4_gp_fix$pointwise-cv_6_gp_fix$pointwise),
                                            sqrt(20)*sd(cv_3_gp_fix$pointwise-cv_6_gp_fix$pointwise)))
      rownames(model_comp) <- c("mod6","mod5","mod4","mod3")                         
                               
      model_comp <- model_comp %>% data.frame()%>% 
        select("elpd_diff","se_diff")%>%
        mutate(model = rownames(.))%>%
        tidyr::pivot_wider(c("elpd_diff","se_diff"),
                           names_from = model,
                           names_sep = "^",
                           values_from = c("elpd_diff","se_diff"))
      allcv <- rbind(allcv[,c(5:8,1:4)],model_comp)
      
        cv_3_ind <- rstanarm::loo(mod3, cores = 20)
        cv_4_ind <- rstanarm::loo(mod4, cores = 20)
        cv_5_ind <- rstanarm::loo(mod5, cores = 20)
        cv_6_ind <- rstanarm::loo(mod6, cores = 20)
        
        model_comp_ind <- data.frame(elpd_diff = c(sum(cv_6_ind$pointwise[,1]-cv_6_ind$pointwise[,1]),
                                               sum(cv_5_ind$pointwise[,1]-cv_6_ind$pointwise[,1]),
                                               sum(cv_4_ind$pointwise[,1]-cv_6_ind$pointwise[,1]),
                                               sum(cv_3_ind$pointwise[,1]-cv_6_ind$pointwise[,1])),
                                 se_diff =  c(sqrt(nrow(smallDat))*sd(cv_6_ind$pointwise[,1]-cv_6_ind$pointwise[,1]),
                                              sqrt(nrow(smallDat))*sd(cv_5_ind$pointwise[,1]-cv_6_ind$pointwise[,1]),
                                              sqrt(nrow(smallDat))*sd(cv_4_ind$pointwise[,1]-cv_6_ind$pointwise[,1]),
                                              sqrt(nrow(smallDat))*sd(cv_3_ind$pointwise[,1]-cv_6_ind$pointwise[,1])))
        
        rownames(model_comp_ind) <- c("mod6","mod5","mod4","mod3")  
        
        model_comp_ind <- model_comp_ind%>% data.frame()%>% 
          select("elpd_diff","se_diff")%>%
          mutate(model = rownames(.))%>%
          tidyr::pivot_wider(c("elpd_diff","se_diff"),
                             names_from = model,
                             names_sep = "^",
                             values_from = c("elpd_diff","se_diff"))
        
        
        indcv <- rbind(indcv[,c(5:8,1:4)],model_comp_ind)
        
    }
    
    allcv$data <- allTitles
    indcv$data <- allTitles
    
    gpcvlist[[myCount]] <- allcv
    indcvlist[[myCount]] <- indcv
    myCount <- myCount + 1
  }
}

saveRDS(list(gpcvlist,indcvlist), file = "results/final_case2_results.rds")

