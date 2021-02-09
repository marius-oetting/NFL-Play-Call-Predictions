### R code for the paper
### "Predicting play-calls in the National Football League using hidden Markov models
### author: Marius Ötting


## load packages
library(dplyr)
library(ggplot2)
library(lemon)


# import data
all.teams.df <- read.csv("nfl_data.csv")


# figures -----------------------------------------------------------------
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## Figure 1

plot.shotgun <- all.teams.df %>% ggplot(aes(x = down, y = PlayType.new, fill = factor(shotgun))) + 
  geom_bar(position = "dodge", stat = "summary", fun.y = "mean") + ylab("pass proportion") +
  scale_fill_manual(name = "Shotgun formation", 
                    labels = c("no", 
                               "yes"), 
                    values = cbbPalette[6:7]) + theme_minimal() +
  theme(text = element_text(size = 22))
plot.shotgun


## Figure 2
all.teams.df <- all.teams.df %>% mutate(ScoreDiffCat = cut(ScoreDiff, breaks = c(-Inf, -8, 0, 8, Inf))) %>% 
  mutate(ScoreDiff1 = if_else(ScoreDiffCat == "(-Inf,-8]", 1, 0)) %>% 
  mutate(ScoreDiff2 = if_else(ScoreDiffCat == "(-8,0]", 1, 0)) %>% 
  mutate(ScoreDiff3 = if_else(ScoreDiffCat == "(0,8]", 1, 0)) %>% 
  mutate(ScoreDiff4 = if_else(ScoreDiffCat == "(8, Inf]", 1, 0))

plot.scorediff <- all.teams.df %>% group_by(ydstogo, ScoreDiffCat) %>% summarize(meanpass = mean(PlayType.new), count = n()) %>% 
  filter(ydstogo <= 25) %>% 
  ggplot(aes(x = ydstogo, y = meanpass, color = ScoreDiffCat)) + 
  geom_pointline() + ylab("pass proportion") + xlab("yards to go") +
  scale_color_manual(name = "Categorized\nscore difference", 
                     labels = c(expression(phantom(x) < "-7"), 
                                expression(paste(" ", phantom(x) >= "-7", " &", phantom(x) <= "0")), 
                                expression(paste(" ", phantom(x) > "0", " &", phantom(x) <= "7")), 
                                expression(phantom(x) > "7")), 
                     values = cbbPalette) + theme_minimal() +
  theme(text = element_text(size = 22))
plot.scorediff



# model fitting -----------------------------------------------------------

# first, different functions are defined:
# 1) HMM likelihood
# 2) function to fit the model using nlm()
# 3) two functions for forecasting
# 4) a function which performs data preprocessing, model fitting, AIC forward selection,
#    and forecasting for a selected team

# TODO!!
# Using these functions, models are fitted to each team separately (see below) and
# summary statistics on the forecasts are computed

L.choosecov <- function(theta.star, x, N, covariates){
  nr.covariates <- length(covariates) + 1 # 10 + 1
  beta <- matrix(theta.star[1:((N - 1) * N * nr.covariates)], nrow = N * (N - 1), ncol = nr.covariates, byrow = FALSE)
  pis.beta <- theta.star[((N - 1) * N * nr.covariates + 1):(((N - 1) * N * nr.covariates) + N)]
  delta <- exp(c(theta.star[(((N - 1) * N * nr.covariates) + N + 1):(((N - 1) * N * nr.covariates) + N + N - 1)], 0)) / 
    sum(exp(c(theta.star[(((N - 1) * N * nr.covariates) + N + 1):(((N - 1) * N * nr.covariates) + N + N - 1)], 0)))
  
  l.all <- 0
  for(k in 1:length(x)){
    
    idx.covariates <- which(colnames(x[[k]]) %in% covariates)
    idx.response <- which(colnames(x[[k]]) == "PlayType.new")
    zwerg <- lapply(x, "[", c(idx.covariates, idx.response)) # response is included
    
    covariate.mat <- zwerg[[k]][!names(zwerg[[k]]) %in% c("PlayType.new")]
    covariate.mat <- covariate.mat[, order(match(colnames(covariate.mat), covariates))] # sort columns
    covariate.mat <- cbind(c(rep(1, nrow(zwerg[[k]]))), covariate.mat) # response is not included
    covariate.mat <- as.matrix(covariate.mat)
    
    ind <- which(!is.na(zwerg[[k]]$PlayType.new))
    allprobs <- matrix(1, nrow(zwerg[[k]]), N)
    
    
    for(j in 1:N){
      allprobs[ind, j] <- dbinom(zwerg[[k]]$PlayType.new[ind], size = 1, plogis(pis.beta[j])) 
    }
    foo <- delta %*% diag(allprobs[1,])
    l <- log(sum(foo))
    phi <- foo/sum(foo)
    for(t in 2:nrow(zwerg[[k]])){
      eta <- as.vector(beta %*% covariate.mat[t,])
      Gamma <- diag(N)
      Gamma[!Gamma] <- exp(eta)
      Gamma <- Gamma/rowSums(Gamma)
      foo <- phi %*% Gamma %*% diag(allprobs[t,])
      
      l <- l + log(sum(foo))
      phi <- foo/sum(foo)
    }
    l.all <- l.all + l
  }
  return(-l.all)
}


mle <- function(theta, N, data, covariates){
  theta.star <- theta
  mod <- nlm(L.choosecov, theta.star, x = data, N = N, covariates = covariates, print.level = 2, iterlim = 10000, hessian = FALSE)
  hess.mat <- mod$hessian
  result <- mod$estimate
  nr.covariates <- length(covariates) + 1
  
  list(beta = matrix(result[1:((N - 1) * N * nr.covariates)], nrow = N * (N - 1), ncol = nr.covariates, byrow = FALSE), 
       pis = plogis(result[((N - 1) * N * nr.covariates + 1):(((N - 1) * N * nr.covariates) + N)]), 
       delta = c(exp(c(result[(((N - 1) * N * nr.covariates) + N + 1):(((N - 1) * N * nr.covariates) + N + N - 1)], 0)) / 
                   sum(exp(c(result[(((N - 1) * N * nr.covariates) + N + 1):(((N - 1) * N * nr.covariates) + N + N - 1)], 0)))),
       wp = mod$estimate, AIC = 2 * (mod$minimum + length(mod$estimate)), llk = mod$minimum, hessian = hess.mat)
}



hmm.forecast <- function(xf, x.new, mod, covariates = selected.covariates){
  n <- length(x.new)
  nr.plays <- sum(unlist(lapply(x.new, nrow)))
  nxf <- length(xf)
  dxf.fin <- matrix(ncol = 3)
  
  for(k in 1:n){
    dxf <- matrix(0, nrow = nrow(x.new[[k]]), ncol = nxf)
    
    allprobs <- matrix(1, nrow(x.new[[k]]), N)
    ind <- which(!is.na(x.new[[k]]$PlayType.new))
    
    # covariates
    idx.covariates <- which(colnames(x.new[[k]]) %in% covariates)
    idx.response <- which(colnames(x.new[[k]]) == "PlayType.new")
    zwerg <- lapply(x.new, "[", c(idx.covariates, idx.response)) # response is included
    
    covariate.mat <- zwerg[[k]][!names(zwerg[[k]]) %in% c("PlayType.new")]
    covariate.mat <- covariate.mat[, order(match(colnames(covariate.mat), covariates))] # sort columns
    covariate.mat <- cbind(c(rep(1, nrow(zwerg[[k]]))), covariate.mat) # response is not included
    covariate.mat <- as.matrix(covariate.mat)
    
    # covariate effects
    beta <- mod$beta
    
    foo <- mod$delta 
    sumfoo <- sum(foo)
    lscale <- log(sumfoo)
    foo <- foo/sumfoo
    for(j in 1:N){
      dxf[1,] <- dxf[1,] + foo[j] * dbinom(xf, size = 1, prob = mod$pis[j])
    }    
    
    for(i in 2:nrow(x.new[[k]])){
      eta <- as.vector(beta %*% covariate.mat[i,])
      Gamma1 <- diag(N)
      Gamma1[!Gamma1] <- exp(eta)
      Gamma1 <- Gamma1/rowSums(Gamma1)
      foo <- foo %*% Gamma1
      
      for(j in 1:N){
        dxf[i,] <- dxf[i,] + foo[j] * dbinom(xf, size = 1, prob = mod$pis[j])
      }
    }
    dxf <- cbind(dxf, k)
    dxf.fin <- rbind(dxf.fin, dxf)
  }
  dxf.fin <- dxf.fin[-1, ]
  return(dxf.fin)
}


hmm.forecast.step <- function(xf, x.new, mod, H = 1, n, N, covariates = selected.covariates){
  # input: one time series, i.e. one match
  nxf <- length(xf)
  
  allprobs <- matrix(1, n, nxf)
  
  # covariates
  idx.covariates <- which(colnames(x.new) %in% covariates)
  idx.response <- which(colnames(x.new) == "PlayType.new")
  zwerg <- x.new[, c(idx.covariates, idx.response)] # response is included
  
  covariate.mat <- zwerg[!names(zwerg) %in% c("PlayType.new")]
  covariate.mat <- covariate.mat[, order(match(colnames(covariate.mat), covariates))] # sort columns
  covariate.mat <- cbind(c(rep(1, nrow(zwerg))), covariate.mat) # response is not included
  covariate.mat <- as.matrix(covariate.mat)
  
  # covariate effects
  beta <- mod$beta
  
  for(j in 1:N){
    allprobs[, j] <- dbinom(x.new$PlayType.new[1:n], size = 1, mod$pis[j])
  }
  foo <- mod$delta %*% diag(allprobs[1, ])
  
  sumfoo <- sum(foo)
  lscale <- log(sumfoo)
  foo <- foo / sumfoo
  
  for(i in 2:n){
    eta <- as.vector(beta %*% covariate.mat[i,])
    Gamma1 <- diag(N)
    Gamma1[!Gamma1] <- exp(eta)
    Gamma1 <- Gamma1 / rowSums(Gamma1)
    
    foo <- foo %*% Gamma1 * allprobs[i, ]
    sumfoo <- sum(foo)
    lscale <- lscale + log(sumfoo)
    foo <- foo / sumfoo
  }
  xi <- matrix(NA, nrow = N, ncol = H)
  allprobs.new <- matrix(0, nrow = H, ncol = nxf)
  
  for(i in (n + 1):(n + H)){
    eta <- as.vector(beta %*% covariate.mat[i,])
    Gamma1 <- diag(N)
    Gamma1[!Gamma1] <- exp(eta)
    Gamma1 <- Gamma1 / rowSums(Gamma1)
    
    foo <- foo %*% Gamma1
    xi[, (i - n)] <- foo
    
    for(j in 1:N){
      allprobs.new[(i - n), ] <- allprobs.new[(i - n),] + foo[j] * dbinom(xf, size = 1, mod$pis[j])
    }
  }
  return(allprobs.new)
}




fit.team <- function(team){
  akt.team <- team
  gameid.49 <- unique(all.teams.df[which(all.teams.df$posteam == akt.team | all.teams.df$defteam == akt.team), ]$game_id)
  data.team <- all.teams.df %>% filter(game_id %in% gameid.49)
  data.team <- data.team %>% filter(posteam==akt.team)
  data.team$GameIDDrive <- paste0(data.team$game_id, data.team$drive)
  
  # create yards per run / yards per pass column
  gameids <- unique(data.team$game_id)
  data.team$prev_run <- NA
  data.team$prev_pass <- NA
  zwerg.list <- list()
  
  for(i in 1:length(gameids)){ 
    cur.match <- data.team %>% filter(game_id == gameids[i])
    cur.match$prev_pass[1] <- 0
    cur.match$prev_run[1] <- 0
    for(j in 2:nrow(cur.match)){
      zwerg_prev_run <- cur.match[1:(j - 1),] %>% filter(play_type == "run") %>% .[["yards_gained"]]
      if(length(zwerg_prev_run) == 0) cur.match$prev_run[j] <- 0
      else cur.match$prev_run[j] <- mean(zwerg_prev_run)
      
      zwerg_prev_pass <- cur.match[1:(j - 1),] %>% filter(play_type == "pass") %>% .[["yards_gained"]]
      if(length(zwerg_prev_pass) == 0) cur.match$prev_pass[j] <- 0
      else cur.match$prev_pass[j] <- mean(zwerg_prev_pass)
    }
    zwerg.list[[i]] <- cur.match
    #print(i)
  }
  
  data.team <- bind_rows(zwerg.list) %>% as.data.frame()
  
  # one time series = one match
  data.list <- split(data.team, f = data.team$game_id)
  
  # remove time series with less than 3 observations
  idx.delete <- sapply(data.list, nrow) < 3
  idx.delete <- names(which(idx.delete == TRUE))
  data.team <- data.team[!(data.team$GameIDDrive %in% idx.delete),]
  
  # playtype as numeric
  data.team$PlayType.new <- ifelse(data.team$play_type == "pass", 1, 0)
  
  # generate scaled covariates
  data.team$ScoreDiffscale <- as.vector(scale(data.team$ScoreDiff))
  data.team$ydstogoscale <- as.vector(scale(data.team$ydstogo))
  data.team$downscale <- as.vector(scale(data.team$down))
  data.team$ydstogo_down <- data.team$ydstogo * data.team$down
  data.team$ydstogo_down_scale <- as.vector(scale(data.team$ydstogo_down))
  data.team$prev_pass_scale <- as.vector(scale(data.team$prev_pass)) 
  data.team$prev_run_scale <- as.vector(scale(data.team$prev_run)) 
  
  data.team$down1 <- ifelse(data.team$down == 1, 1, 0)
  data.team$down2 <- ifelse(data.team$down == 2, 1, 0)
  data.team$down3 <- ifelse(data.team$down == 3, 1, 0)
  
  # generate interaction columns
  data.team$interact_down1_yds <- data.team$down1 * data.team$ydstogoscale
  data.team$interact_down2_yds <- data.team$down2 * data.team$ydstogoscale
  data.team$interact_down3_yds <- data.team$down3 * data.team$ydstogoscale
  
  data.team$interact_shotgun_ydstogo <- data.team$shotgun * data.team$ydstogoscale
  data.team$interact_nohuddle_ScoreDiff <- data.team$nohuddle * data.team$ScoreDiffscale
  data.team$interact_nohuddle_shotgun <- data.team$nohuddle * data.team$shotgun
  
  
  # test data = data for season 2018/19
  data.team.18.short <- data.team[data.team$season == "2018/19",]
  # training data = data until season 2017/18
  data.team <- data.team[data.team$season != "2018/19",]
  
  # training data as list
  data.list <- split(data.team, f = data.team$game_id)
  # test data as list
  data.list.18 <- split(data.team.18.short, f = data.team.18.short$game_id)
  
  
  ### forward selection
  considered.covariates <- c("ScoreDiffscale", "ydstogoscale", "shotgun", "down1", "down2", "down3", "Home", "yardline90",
                             "interact_down1_yds", "interact_down2_yds", "interact_down3_yds", "interact_shotgun_ydstogo",
                             "nohuddle", "interact_nohuddle_ScoreDiff", "interact_nohuddle_shotgun")
  aic.vec <- rep(NA, length(considered.covariates) + 1)
  N <- 2
  
  # step 0: fit model without any covariates
  selected.covariates <- c()
  nr.covariates <- length(selected.covariates) + 1
  theta.star <- c(runif(N * (N - 1) * nr.covariates, -0.1, 0.1), qlogis(runif(N, 0, 1)), qlogis(runif(N - 1, 0, 1)))
  mod <- mle(theta = theta.star, N = N, data = data.list, covariates = selected.covariates)
  aic.vec[1] <- mod$AIC
  
  # step 1: fit model with covariate 1
  # step 2: add covariate 2, if AIC increases, keep covariate 2, ...
  # ...
  
  for(i in 1:length(considered.covariates)){
    selected.covariates <- c(selected.covariates, considered.covariates[i])
    # choose initial values
    nr.covariates <- length(selected.covariates) + 1
    theta.star <- c(runif(N * (N - 1) * nr.covariates, -0.1, 0.1), qlogis(runif(N, 0, 1)), qlogis(runif(N - 1, 0, 1)))
    # fit model
    mod <- tryCatch(mle(theta = theta.star, N = N, data = data.list, covariates = selected.covariates),
                    error = function(e) NA)
    
    # if nlm failed, run again with different starting values
    if(is.na(mod[1])){
      for(a in 1:5){
        theta.star <- c(runif(N * (N - 1) * nr.covariates, -0.1, 0.1), qlogis(runif(N, 0, 1)), qlogis(runif(N - 1, 0, 1)))
        mod <- tryCatch(mle(theta = theta.star, N = N, data = data.list, covariates = selected.covariates),
                        error = function(e) NA)
        if(!is.na(mod[1])) break;
      }
    }
    
    # store and compare AIC
    aic.vec[i + 1] <- mod$AIC
    if(mod$AIC < aic.vec[i]) selected.covariates <- selected.covariates
    else selected.covariates <- selected.covariates[selected.covariates != considered.covariates[i]]
  }
  
  ### fit model with selected.covariates
  
  N <- 2
  nr.covariates <- length(selected.covariates) + 1
  theta.star <- c(runif(N * (N - 1) * nr.covariates, -0.1, 0.1), qlogis(runif(N, 0, 1)), qlogis(runif(N - 1, 0, 1)))
  
  mod <- mle(theta = theta.star, N = N, data = data.list, covariates = selected.covariates)
  
  # forecast
  nr.matches <- length(data.list.18)
  forecast.all <- matrix(NA, ncol = 3)
  for(i in 1:nr.matches){
    cur.match <- data.list.18[[i]]
    forecast.probs.match <- matrix(NA, nrow = nrow(cur.match), ncol = 2)
    for(j in 2:(nrow(cur.match) - 1)){
      forecast.probs.step <- hmm.forecast.step(c(0:1), cur.match, mod = mod, H = 1, n = j, N = N, covariates = selected.covariates)
      forecast.probs.match[j + 1, ] <- forecast.probs.step
    }
    forecast.probs.firsttwo <- hmm.forecast(c(0:1), data.list.18, mod = mod, covariates = selected.covariates) %>% 
      as.data.frame() %>% filter(k == i) %>% head(2) %>% select(V1, V2)
    forecast.probs.match[1:2, ] <- forecast.probs.firsttwo %>% as.matrix()
    forecast.probs.match <- cbind(forecast.probs.match, i)
    forecast.all <- rbind(forecast.all, forecast.probs.match)
    # dynamically update forecasts throughout season
    data.list[[length(data.list) + 1]] <- cur.match
    mod <- mle(theta = theta.star, N = N, data = data.list, covariates = selected.covariates)
  }
  forecast.all <- forecast.all[-1, ] 
  
  forecast.plays <- apply(forecast.all[, 1:2], 1, which.max)
  forecast.plays[forecast.plays == 1] <- 0
  forecast.plays[forecast.plays == 2] <- 1
  
  true.play <- data.team.18.short$PlayType.new
  # prediction accuracy
  pred.accuracy.step <- length(which(forecast.plays == true.play)) / length(na.omit(true.play))
  # precision pass
  prec.pass <- sum(forecast.plays == 1 & true.play == 1) / (sum(forecast.plays == 1 & true.play == 1) + sum(forecast.plays == 1 & true.play == 0))
  # precision run
  prec.run <- sum(forecast.plays == 0 & true.play == 0) / (sum(forecast.plays == 0 & true.play == 0) + sum(forecast.plays == 0 & true.play == 1))
  # recall pass
  reca.pass <- sum(forecast.plays == 1 & true.play == 1) / (sum(forecast.plays == 1 & true.play == 1) + sum(forecast.plays == 0 & true.play == 1))
  # recall run
  reca.run <- sum(forecast.plays == 0 & true.play == 0) / (sum(forecast.plays == 0 & true.play == 0) + sum(forecast.plays == 1 & true.play == 0))
  
  # nr observations
  nr.obs <- length(forecast.plays)
  
  zwerg.team <- data.frame(akt.team, nr.obs, pred.accuracy.step, prec.pass, prec.run, reca.pass, reca.run)
  return(zwerg.team)
}


### loop over all teams
all.teams <- all.teams.df$home_team %>% as.character() %>% unique

res <- data.frame(team = NA, count = NA, pred.acc = NA, precision.pass = NA, precision.run = NA,
                  recall.pass = NA, recall.run = NA)

set.seed(123)

for(i in 1:length(all.teams)){  
  N <- 2
  cur.team <- all.teams[i]
  if(cur.team == "BUF" | cur.team == "DET" | cur.team == "CHI") N <- 3
  zwerg <- fit.team(cur.team)

  colnames(zwerg) <- c("team", "count", "pred.acc", "precision.pass", "precision.run", "recall.pass", "recall.run")
  res <- rbind(res, zwerg) 
  print(paste("Finished iteration", i))
}

### "res" includes summary statistics on the prediction accuracy etc. for all teams
View(res)
