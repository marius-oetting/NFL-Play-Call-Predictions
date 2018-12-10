# This R-Code fits a hidden Markov model (HMM) to play-call data of a chosen NFL team and predicts play-calls of the
# 2017 season.
# Author: Marius Oetting
# Both the training and test data are stored in .RData files.
# Running the code includes the following 3 steps:
# 1. choose a team and subset the corresponding data (training data: data.list; test data: data.list.17)
# 2. optimize the HMM likelihood by using the mle.2states() function
# 3. forecast play-calls of the 2017 season by using the hmm.forecast() function

library(data.table)
# load training and test data
load("TrainingDataList.RData")
load("TestDataList.RData")
# list of all available teams:
names(data.list.all)


# choose one of the teams above here, e.g. the Arizona Cardinals:
chosen.team <- "ARI"
# training data
data.list <- data.list.all[[chosen.team]]
# test data 
data.list.17 <- data.list.17.all[[chosen.team]]
data.17.df <- as.data.frame(rbindlist(data.list.17))


# HMM likelihood ----------------------------------------------------------
L.hmm <- function(theta.star, x, N){
  Gamma1 <- diag(N)
  Gamma1[!Gamma1] <- exp(theta.star[1:2])
  Gamma1 <- Gamma1/rowSums(Gamma1)
  
  Gamma2 <- diag(N)
  Gamma2[!Gamma2] <- exp(theta.star[3:4])
  Gamma2 <- Gamma2/rowSums(Gamma2)
  
  pis.beta <- theta.star[5:6]  
  beta.shotgun <- theta.star[7]
  beta.down1 <- theta.star[8]
  beta.down2 <- theta.star[9]
  beta.down3 <- theta.star[10]
  beta.ydstogo <- theta.star[11]
  beta.nohuddle <- theta.star[12]
  beta.scorediff <- theta.star[13]
  beta.nrdrive <- theta.star[14]
  beta.playtype <- theta.star[15]
  beta.lasttwo <- theta.star[16]
  beta.goaltogo <- theta.star[17]
  beta.shotgun_ydstogo <- theta.star[18]
  beta.shotgun_down1 <- theta.star[19]
  beta.shotgun_down2 <- theta.star[20]
  beta.shotgun_down3 <- theta.star[21]
  beta.scorediff_lasttwo <- theta.star[22]
  beta.down1_ydstogo <- theta.star[23]
  beta.down2_ydstogo <- theta.star[24]
  beta.down3_ydstogo <- theta.star[25]
  delta <- c(plogis(theta.star[26]), 1 - plogis(theta.star[26]))
  
  l.all <- 0
  for(k in 1:length(x)){ 
    allprobs <- matrix(1, nrow(x[[k]]), N)
    ind <- which(!is.na(x[[k]]$PlayType.new)) 
    
    # covariates
    shotgun <- x[[k]]$shotgun
    goaltogo <- x[[k]]$GoalToGo
    nohuddle <- x[[k]]$nohuddle
    down1 <- as.numeric(x[[k]]$down==1)
    down2 <- as.numeric(x[[k]]$down==2)
    down3 <- as.numeric(x[[k]]$down==3)
    ydstogo <- x[[k]]$ydstogo
    scorediff <- x[[k]]$ScoreDiff
    nrdrive <- x[[k]]$NrDrive
    playmin1 <- c(0, x[[k]]$PlayType.new[1:(nrow(x[[k]]) - 1)]) 
    lasttwo <- x[[k]]$LastTwoMin
    
    for(j in 1:N){
      allprobs[ind, j] <- dbinom(x[[k]]$PlayType.new[ind], size = 1, plogis(pis.beta[j] + beta.shotgun * shotgun[ind] + 
                                                                              beta.nohuddle * nohuddle[ind] + 
                                                                              beta.down1 * down1[ind] + 
                                                                              beta.down2 * down2[ind] +
                                                                              beta.down3 * down3[ind] + 
                                                                              beta.ydstogo * ydstogo[ind] + 
                                                                              beta.scorediff * scorediff[ind] + 
                                                                              beta.nrdrive * nrdrive[ind] + 
                                                                              beta.playtype * playmin1[ind] + 
                                                                              beta.lasttwo * lasttwo[ind] + 
                                                                              beta.goaltogo * goaltogo[ind] +
                                                                              beta.shotgun_ydstogo * shotgun[ind] * ydstogo[ind] +
                                                                              beta.shotgun_down1 * shotgun[ind] * down1[ind] + 
                                                                              beta.shotgun_down2 * shotgun[ind] * down2[ind] + 
                                                                              beta.shotgun_down3 * shotgun[ind] * down3[ind] + 
                                                                              beta.scorediff_lasttwo * lasttwo[ind] * scorediff[ind] + 
                                                                              beta.down1_ydstogo * down1[ind] * ydstogo[ind] + 
                                                                              beta.down2_ydstogo * down2[ind] * ydstogo[ind] +
                                                                              beta.down3_ydstogo * down3[ind] * ydstogo[ind]))
    }
    foo <- delta %*% diag(allprobs[1,])
    l <- log(sum(foo))
    phi <- foo/sum(foo)
    for(t in 2:nrow(x[[k]])){
      if(x[[k]]$Drive[t] == x[[k]]$Drive[t-1]){
        foo <- phi %*% Gamma1 %*% diag(allprobs[t,])
        l <- l + log(sum(foo))
        phi <- foo/sum(foo)
      }
      else{
        foo <- phi %*% Gamma2 %*% diag(allprobs[t,])
        l <- l + log(sum(foo))
        phi <- foo/sum(foo)
      }
    }
    l.all <- l.all + l
  }
  return(-l.all)
}

# starting values
theta.star <- c(rep(-2, 2), rep(-2.3, 2), qlogis(0.4), qlogis(0.3), qlogis(0.55), qlogis(0.53), qlogis(0.45), qlogis(0.55),
                qlogis(0.6), qlogis(0.4), qlogis(0.5), qlogis(0.7), qlogis(0.4), qlogis(0.5), qlogis(0.4), qlogis(0.5), qlogis(0.4), 
                qlogis(0.55), qlogis(0.52), qlogis(0.5), qlogis(0.5), qlogis(0.5), qlogis(0.5),
                qlogis(0.7))
L.hmm(theta.star = theta.star, data.list, 2)

# function for optimizing the likelihood
mle.2states <- function(theta, N){
  theta.star <- theta
  mod <- nlm(L.hmm, theta.star, x = data.list, N = N, print.level = 2, iterlim = 1000)
  hess.mat <- mod$hessian
  result <- c(mod$estimate[1:26])
  list(pi.betas=result[5:6], beta.shotgun=result[7], beta.down1=result[8], 
       beta.down2=result[9], beta.down3=result[10], beta.ydstogo=result[11], beta.nohuddle=result[12],
       beta.scorediff=result[13], beta.nrdrive=result[14], beta.playtype=result[15],
       beta.lasttwo=result[16], beta.goaltogo=result[17], beta.shotgun_ydstogo=result[18],
       beta.shotgun_down1=result[19], beta.shotgun_down2=result[20], beta.shotgun_down3=result[21],
       beta.scorediff_lasttwo=result[22], beta.down1_ydstogo=result[23], beta.down2_ydstogo=result[24],
       beta.down3_ydstogo=result[25], 
       delta = c(plogis(result[26]), 1 - plogis(result[26])), gamma1=result[1:2], gamma2=result[3:4],
       wp = mod$estimate, AIC=2*(mod$minimum+length(mod$estimate)), llk=mod$minimum, hessian=hess.mat)
}

# optimize likelihood
mod <- mle.2states(theta.star, 2)


# function for forecast ---------------------------------------------------
hmm.forecast <- function(xf, x.new, mod){
  n <- length(x.new)
  nr.plays <- sum(unlist(lapply(x.new, nrow)))
  nxf <- length(xf)
  dxf.fin <- matrix(ncol = 2)
  
  for(k in 1:n){
    dxf <- matrix(0, nrow = nrow(x.new[[k]]), ncol = nxf)
    
    allprobs <- matrix(1, nrow(x.new[[k]]), 2)
    ind <- which(!is.na(x.new[[k]]$PlayType.new))
    
    # covariates
    shotgun <- x.new[[k]]$shotgun
    goaltogo <- x.new[[k]]$GoalToGo
    down1 <- as.numeric(x.new[[k]]$down==1)
    down2 <- as.numeric(x.new[[k]]$down==2)
    down3 <- as.numeric(x.new[[k]]$down==3)
    ydstogo <- x.new[[k]]$ydstogo
    nohuddle <- x.new[[k]]$nohuddle
    scorediff <- x.new[[k]]$ScoreDiff
    nrdrive <- x.new[[k]]$NrDrive
    playmin1 <- c(0, x.new[[k]]$PlayType.new[1:(nrow(x.new[[k]]) - 1)])
    lasttwo <- x.new[[k]]$LastTwoMin
    
    # t.p.m
    Gamma1 <- diag(2)
    Gamma1[!Gamma1] <- exp(mod$gamma1[1:2])
    Gamma1 <- Gamma1/rowSums(Gamma1)
    
    Gamma2 <- diag(2)
    Gamma2[!Gamma2] <- exp(mod$gamma2[1:2])
    Gamma2 <- Gamma2/rowSums(Gamma2)
    
    
    for(j in 1:2){
      allprobs[ind, j] <- dbinom(x.new[[k]]$PlayType.new[ind], size = 1, prob = plogis(mod$pi.betas[j] + mod$beta.shotgun * shotgun[ind] + 
                                                                                         mod$beta.goaltogo * goaltogo[ind] +
                                                                                         mod$beta.nohuddle * nohuddle[ind] + 
                                                                                         mod$beta.down1 * down1[ind] + 
                                                                                         mod$beta.down2 * down2[ind] +
                                                                                         mod$beta.down3 * down3[ind] + 
                                                                                         mod$beta.ydstogo * ydstogo[ind] + 
                                                                                         mod$beta.scorediff * scorediff[ind] + 
                                                                                         mod$beta.nrdrive * nrdrive[ind] + 
                                                                                         mod$beta.playtype * playmin1[ind] + 
                                                                                         mod$beta.lasttwo * lasttwo[ind] + 
                                                                                         mod$beta.shotgun_ydstogo * shotgun[ind] * ydstogo[ind] +
                                                                                         mod$beta.shotgun_down1 * shotgun[ind] * down1[ind] + 
                                                                                         mod$beta.shotgun_down2 * shotgun[ind] * down2[ind] + 
                                                                                         mod$beta.shotgun_down3 * shotgun[ind] * down3[ind] + 
                                                                                         mod$beta.scorediff_lasttwo * scorediff[ind] * lasttwo[ind] + 
                                                                                         mod$beta.down1_ydstogo * down1[ind] * ydstogo[ind] + 
                                                                                         mod$beta.down2_ydstogo * down2[ind] * ydstogo[ind] + 
                                                                                         mod$beta.down3_ydstogo * down3[ind] * ydstogo[ind]))
    }
    
    foo <- mod$delta %*% diag(allprobs[1,])
    if(x.new[[k]]$PlayType.new[ind][1] == 1) dxf[1,] <- c(1 - sum(foo), sum(foo))
    else dxf[1,] <- c(sum(foo), 1 - sum(foo))
    sumfoo <- sum(foo)
    lscale <- log(sumfoo)
    foo <- foo/sumfoo
    
    
    for(i in 2:nrow(x.new[[k]])){
      if(x.new[[k]]$Drive[i] == x.new[[k]]$Drive[i-1]){
        foo <- foo %*% Gamma1
      }
      else foo <- foo %*% Gamma2
      
      for(j in 1:2){
        dxf[i,] <- dxf[i,] + foo[j] * dbinom(xf, size = 1, prob = plogis(mod$pi.betas[j] + 
                                                                           mod$beta.shotgun * x.new[[k]]$shotgun[i] + 
                                                                           mod$beta.goaltogo * x.new[[k]]$GoalToGo[i] +
                                                                           mod$beta.nohuddle * x.new[[k]]$nohuddle[i] + 
                                                                           mod$beta.down1 * as.numeric(x.new[[k]]$down[i]==1) + 
                                                                           mod$beta.down2 * as.numeric(x.new[[k]]$down[i]==2) +
                                                                           mod$beta.down3 * as.numeric(x.new[[k]]$down[i]==3) + 
                                                                           mod$beta.ydstogo * x.new[[k]]$ydstogo[i] + 
                                                                           mod$beta.scorediff * x.new[[k]]$ScoreDiff[i] + 
                                                                           mod$beta.nrdrive * x.new[[k]]$NrDrive[i] + 
                                                                           mod$beta.playtype * x.new[[k]]$playmin1[i] +
                                                                           mod$beta.lasttwo * x.new[[k]]$LastTwo[i] + 
                                                                           mod$beta.shotgun_ydstogo * x.new[[k]]$shotgun[i] * x.new[[k]]$ydstogo[i] +
                                                                           mod$beta.shotgun_down1 * x.new[[k]]$shotgun[i] * as.numeric(x.new[[k]]$down[i]==1) + 
                                                                           mod$beta.shotgun_down2 * x.new[[k]]$shotgun[i] * as.numeric(x.new[[k]]$down[i]==2) + 
                                                                           mod$beta.shotgun_down3 * x.new[[k]]$shotgun[i] * as.numeric(x.new[[k]]$down[i]==3) + 
                                                                           mod$beta.scorediff_lasttwo * x.new[[k]]$ScoreDiff[i] * x.new[[k]]$LastTwo[i] + 
                                                                           mod$beta.down1_ydstogo * as.numeric(x.new[[k]]$down[i]==1) * x.new[[k]]$ydstogo[i] + 
                                                                           mod$beta.down2_ydstogo * as.numeric(x.new[[k]]$down[i]==2) * x.new[[k]]$ydstogo[i] + 
                                                                           mod$beta.down3_ydstogo * as.numeric(x.new[[k]]$down[i]==3) * x.new[[k]]$ydstogo[i]))
      }
    }
    dxf.fin <- rbind(dxf.fin, dxf)
  }
  dxf.fin <- dxf.fin[-1,]
  return(dxf.fin)
}


forecast.probs <- hmm.forecast(c(0:1), data.list.17, mod = mod)
forecast.plays <- apply(forecast.probs, 1, which.max)

forecast.plays[forecast.plays==1] <- 0
forecast.plays[forecast.plays==2] <- 1

# true plays in test data
true.play <- data.17.df$PlayType.new

# PCP
length(which(forecast.plays == true.play)) / length(true.play)
# TPR
sum(forecast.plays == 1 & true.play == 1) / sum(true.play == 1)
# FPR
sum(forecast.plays == 1 & true.play == 0) / sum(true.play == 0)
# PPV
sum(forecast.plays == 1 & true.play == 1) / (sum(forecast.plays == 1))
# NPV
sum(forecast.plays == 0 & true.play == 0) / (sum(forecast.plays == 0))
