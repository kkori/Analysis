install.packages('GA')
install.packages("pROC")
install.packages("ROCR")

library(ROCR)
library(pROC)
library(GA)
library(dummies)
library(moments)
library(dplyr)
library(ggplot2)
library(GGally)
library(corrplot)
library(psych)
library(glmnet)
library(car)

df <- read.csv("Life Expectancy Data.csv")
head(df)
table(is.na(df))
# 결측치를 평균값으로 대체
for (i in 4:22){
  df[i][is.na(df[i])] <- mean(df[,i], na.rm=TRUE)
}

table(is.na(df))
# 명목형 변수 변환
# Country제외
nCou <- nrow(df)
nVar <- ncol(df)

no_idx <- c(1)
category_id <- 3

dummy_s <- rep(0, nCou)
dummy_d <- rep(0, nCou)

s_idx <- which(df$Status == 'Developed')
d_idx <- which(df$Status == 'Developing')

dummy_s[s_idx] <- 1
dummy_d[d_idx] <- 1

sta <- data.frame(dummy_s, dummy_d)
names(sta) <- c("Developed", "Developing")

df <- cbind(df[,-c(no_idx,category_id)], sta)
colnames(df)

# 종속변수를 맨 마지막 컬럼으로
arrange.vars <- function(data, vars){
  ##stop if not a data.frame (but should work for matrices as well)
  stopifnot(is.data.frame(data))
  
  ##sort out inputs
  data.nms <- names(data)
  var.nr <- length(data.nms)
  var.nms <- names(vars)
  var.pos <- vars
  ##sanity checks
  stopifnot( !any(duplicated(var.nms)), 
             !any(duplicated(var.pos)) )
  stopifnot( is.character(var.nms), 
             is.numeric(var.pos) )
  stopifnot( all(var.nms %in% data.nms) )
  stopifnot( all(var.pos > 0), 
             all(var.pos <= var.nr) )
  
  ##prepare output
  out.vec <- character(var.nr)
  out.vec[var.pos] <- var.nms
  out.vec[-var.pos] <- data.nms[ !(data.nms %in% var.nms) ]
  stopifnot( length(out.vec)==var.nr )
  
  ##re-arrange vars by position
  data <- data[ , out.vec]
  return(data)
}
df <- arrange.vars(df, c("Life.expectancy"=22))
colnames(df)

#데이터셋 분할
set.seed(365)
ndf <- nrow(df)
ndf
df_trn_idx <- sample(1:ndf, round(0.7*ndf))
df_trn <- df[df_trn_idx,]
df_val <- df[-df_trn_idx,]

# RMSE, MAE, MAPE 계산
perf_reg <- function(tg_y, pr_y){
  #RMSE
  rmse <- sqrt(mean((tg_y - pr_y)^2))
  #MAE
  mae <- mean(abs(tg_y - pr_y))
  #MAPE
  mape <- 100*mean(abs((tg_y - pr_y)/tg_y))
  return(c(rmse, mae, mape))
}

# matrix
Perf_Table = matrix(0, nrow = 4, ncol = 3)

rownames(Perf_Table) = c("Forward","Bakcward","Stepwise","GA")
colnames(Perf_Table) = c("RMSE", "MAE", "MAPE")
Perf_Table

# [Q1]
# Forward
tmp_x <- paste(colnames(df_trn)[-22], collapse=' + ')
tmp_xy <- paste("Life.expectancy ~ ", tmp_x, collapse = "")
as.formula(tmp_xy)

start_time <- proc.time()
forward_model <- step(lm(Life.expectancy ~ 1, data = df_trn),
                      scope = list(upper = as.formula(tmp_xy),
                                   lower = Life.expectancy ~ 1), direction = 'forward', trace =1)
end_time <- proc.time()
end_time - start_time

summary(forward_model)
fore_coeff <- as.matrix(forward_model$coefficients,22,1)
fore_coeff

forward_model_prey <- predict(forward_model, newdata=df_val)
Perf_Table[1,] <- perf_reg(df_val$Life.expectancy,forward_model_prey)
Perf_Table

#Backward
start_time <- proc.time()
backward_model <- step(lm(Life.expectancy ~ ., data = df_trn),
                       scope = list(upper=as.formula(tmp_xy),
                                    lower=Life.expectancy ~ 1), direction = 'backward', trace = 1)
end_time <- proc.time()
end_time - start_time

summary(backward_model)
back_coeff <- as.matrix(backward_model$coefficients,22,1)
back_coeff

backward_model_prey <- predict(backward_model, newdata=df_val)
Perf_Table[2,] <- perf_reg(df_val$Life.expectancy, backward_model_prey)
Perf_Table

#Stepwise
start_time <- proc.time()
stepwise_model <- step(lm(Life.expectancy ~ 1, data = df_trn),
                       scope = list(upper=as.formula(tmp_xy),
                                    lower=Life.expectancy ~ 1), direction="both", trace=1)
end_time <- proc.time()
end_time - start_time
summary(stepwise_model)
step_coeff <- as.matrix(stepwise_model$coefficients,22,1)
step_coeff

step_model_prey <- predict(stepwise_model, newdata=df_val)
Perf_Table[3,] <- perf_reg(df_val$Life.expectancy,step_model_prey)
Perf_Table

#[Q2]
#GA
#fitness function
fit_adj2 <- function(string){
  sel_var_idx <- which(string == 1)
  sel_x <- x[, sel_var_idx]
  xy <- data.frame(sel_x, y)
  
  ga_ml <- lm(y ~ ., data=xy)
  adj2 <- summary(ga_ml)$adj.r.squared
  
  return(adj2)
}

x <- as.matrix(df_trn[,-22])
y <- df_trn[,22]

start_time <- proc.time()
ga_adj2 <- ga(type="binary", fitness = fit_adj2, nBits = ncol(x),
              names = colnames(x),popSize = 50, pcrossover = 0.5, pmutation = 0.01,
              maxiter = 100, elitism = 2, seed = 2021)
end_time <- proc.time()
end_time - start_time

best_var_idx <- which(ga_adj2@solution == 1)
best_var_idx
GA_trn <- df_trn[,c(best_var_idx,22)]
GA_tst <- df_val[,c(best_var_idx,22)]

GA_model <- lm(Life.expectancy ~ ., GA_trn)
summary(GA_model)

GA_coeff <- as.matrix(GA_model$coefficients,22,1)
GA_coeff

GA_model_prey <- predict(GA_model, newdata = GA_tst)
Perf_Table[4,] <- perf_reg(GA_tst$Life.expectancy, GA_model_prey)
Perf_Table

# [Q3]
# 비교를 위한 matrix 생성
Perf_Table_ga <- matrix(0, nrow = 27, ncol = 3)
Perf_Table_ga_adj <- matrix(0, nrow=27, ncol = 1)
Perf_Time <- matrix(0, nrow=27, ncol =1)
rownames(Perf_Table_ga) <- c(1:27)
colnames(Perf_Table_ga) <- c("RMSE", "MAE", "MAPE")
rownames(Perf_Table_ga_adj) <- c(1:27)
colnames(Perf_Table_ga_adj) <- c("Adj-R2")
rownames(Perf_Time) <- c(1:27)
colnames(Perf_Time) <- c("Time")


# Popsize는 (10,100,500), pcrossover는 (0.1,0,4,0.7), pmutation은 (0.05,0.1,0.3) 이렇게 3가지씩으로 총 27개의 조합을 만든다.
# 10, 0.1, 0.05
start_time <- proc.time()
ga_adj2 <- ga(type="binary", fitness = fit_adj2, nBits = ncol(x),
              names = colnames(x),popSize = 10, pcrossover = 0.1, pmutation = 0.05,
              maxiter = 100, elitism = 2, seed = 2021, monitor = plot)
end_time <- proc.time()
tt <- end_time - start_time
Perf_Time[1,] <- tt[3]

best_var_idx <- which(ga_adj2@solution == 1)
GA_trn <- df_trn[,c(best_var_idx,22)]
GA_tst <- df_val[,c(best_var_idx,22)]

GA_model <- lm(Life.expectancy ~ ., GA_trn)
summary(GA_model)

GA_coeff <- as.matrix(GA_model$coefficients,22,1)
GA_coeff
GA_model_prey <- predict(GA_model, newdata = GA_tst)
Perf_Table_ga[1,] <- perf_reg(GA_tst$Life.expectancy, GA_model_prey)
Perf_Table_ga_adj[1,] <- summary(GA_model)$adj.r.squared

# 10, 0.1, 0.1
start_time <- proc.time()
ga_adj2 <- ga(type="binary", fitness = fit_adj2, nBits = ncol(x),
              names = colnames(x),popSize = 10, pcrossover = 0.1, pmutation = 0.1,
              maxiter = 100, elitism = 2, seed = 2021, monitor = plot)
end_time <- proc.time()
tt <- end_time - start_time
Perf_Time[2,] <- tt[3]

best_var_idx <- which(ga_adj2@solution == 1)
GA_trn <- df_trn[,c(best_var_idx,22)]
GA_tst <- df_val[,c(best_var_idx,22)]

GA_model <- lm(Life.expectancy ~ ., GA_trn)
summary(GA_model)

GA_coeff <- as.matrix(GA_model$coefficients,22,1)
GA_coeff
GA_model_prey <- predict(GA_model, newdata = GA_tst)
Perf_Table_ga[2,] <- perf_reg(GA_tst$Life.expectancy, GA_model_prey)
Perf_Table_ga_adj[2,] <- summary(GA_model)$adj.r.squared

# 10, 0.1, 0.3
start_time <- proc.time()
ga_adj2 <- ga(type="binary", fitness = fit_adj2, nBits = ncol(x),
              names = colnames(x),popSize = 10, pcrossover = 0.1, pmutation = 0.3,
              maxiter = 100, elitism = 2, seed = 2021, monitor = plot)
end_time <- proc.time()
tt <- end_time - start_time
Perf_Time[3,] <- tt[3]


best_var_idx <- which(ga_adj2@solution == 1)
GA_trn <- df_trn[,c(best_var_idx,22)]
GA_tst <- df_val[,c(best_var_idx,22)]

GA_model <- lm(Life.expectancy ~ ., GA_trn)
summary(GA_model)

GA_coeff <- as.matrix(GA_model$coefficients,22,1)
GA_coeff
GA_model_prey <- predict(GA_model, newdata = GA_tst)
Perf_Table_ga[3,] <- perf_reg(GA_tst$Life.expectancy, GA_model_prey)
Perf_Table_ga_adj[3,] <- summary(GA_model)$adj.r.squared

# 10, 0.4, 0.05
start_time <- proc.time()
ga_adj2 <- ga(type="binary", fitness = fit_adj2, nBits = ncol(x),
              names = colnames(x),popSize = 10, pcrossover = 0.4, pmutation = 0.05,
              maxiter = 100, elitism = 2, seed = 2021, monitor = plot)
end_time <- proc.time()
tt <- end_time - start_time
Perf_Time[4,] <- tt[3]


best_var_idx <- which(ga_adj2@solution == 1)
GA_trn <- df_trn[,c(best_var_idx,22)]
GA_tst <- df_val[,c(best_var_idx,22)]

GA_model <- lm(Life.expectancy ~ ., GA_trn)
summary(GA_model)

GA_coeff <- as.matrix(GA_model$coefficients,22,1)
GA_coeff
GA_model_prey <- predict(GA_model, newdata = GA_tst)
Perf_Table_ga[4,] <- perf_reg(GA_tst$Life.expectancy, GA_model_prey)
Perf_Table_ga_adj[4,] <- summary(GA_model)$adj.r.squared

# 10, 0.4, 0.1
start_time <- proc.time()
ga_adj2 <- ga(type="binary", fitness = fit_adj2, nBits = ncol(x),
              names = colnames(x),popSize = 10, pcrossover = 0.4, pmutation = 0.1,
              maxiter = 100, elitism = 2, seed = 2021, monitor = plot)
end_time <- proc.time()
tt <- end_time - start_time
Perf_Time[5,] <- tt[3]


best_var_idx <- which(ga_adj2@solution == 1)
GA_trn <- df_trn[,c(best_var_idx,22)]
GA_tst <- df_val[,c(best_var_idx,22)]

GA_model <- lm(Life.expectancy ~ ., GA_trn)
summary(GA_model)

GA_coeff <- as.matrix(GA_model$coefficients,22,1)
GA_coeff
GA_model_prey <- predict(GA_model, newdata = GA_tst)
Perf_Table_ga[5,] <- perf_reg(GA_tst$Life.expectancy, GA_model_prey)
Perf_Table_ga_adj[5,] <- summary(GA_model)$adj.r.squared

# 10, 0.4, 0.3
start_time <- proc.time()
ga_adj2 <- ga(type="binary", fitness = fit_adj2, nBits = ncol(x),
              names = colnames(x),popSize = 10, pcrossover = 0.4, pmutation = 0.3,
              maxiter = 100, elitism = 2, seed = 2021, monitor = plot)
end_time <- proc.time()
tt <- end_time - start_time
Perf_Time[6,] <- tt[3]


best_var_idx <- which(ga_adj2@solution == 1)
GA_trn <- df_trn[,c(best_var_idx,22)]
GA_tst <- df_val[,c(best_var_idx,22)]

GA_model <- lm(Life.expectancy ~ ., GA_trn)
summary(GA_model)

GA_coeff <- as.matrix(GA_model$coefficients,22,1)
GA_coeff
GA_model_prey <- predict(GA_model, newdata = GA_tst)
Perf_Table_ga[6,] <- perf_reg(GA_tst$Life.expectancy, GA_model_prey)
Perf_Table_ga_adj[6,] <- summary(GA_model)$adj.r.squared

# 10, 0.7, 0.05
start_time <- proc.time()
ga_adj2 <- ga(type="binary", fitness = fit_adj2, nBits = ncol(x),
              names = colnames(x),popSize = 10, pcrossover = 0.7, pmutation = 0.05,
              maxiter = 100, elitism = 2, seed = 2021, monitor = plot)
end_time <- proc.time()
tt <- end_time - start_time
Perf_Time[7,] <- tt[3]


best_var_idx <- which(ga_adj2@solution == 1)
GA_trn <- df_trn[,c(best_var_idx,22)]
GA_tst <- df_val[,c(best_var_idx,22)]

GA_model <- lm(Life.expectancy ~ ., GA_trn)
summary(GA_model)

GA_coeff <- as.matrix(GA_model$coefficients,22,1)
GA_coeff
GA_model_prey <- predict(GA_model, newdata = GA_tst)
Perf_Table_ga[7,] <- perf_reg(GA_tst$Life.expectancy, GA_model_prey)
Perf_Table_ga_adj[7,] <- summary(GA_model)$adj.r.squared

# 10, 0.7, 0.1
start_time <- proc.time()
ga_adj2 <- ga(type="binary", fitness = fit_adj2, nBits = ncol(x),
              names = colnames(x),popSize = 10, pcrossover = 0.7, pmutation = 0.1,
              maxiter = 100, elitism = 2, seed = 2021, monitor = plot)
end_time <- proc.time()
tt <- end_time - start_time
Perf_Time[8,] <- tt[3]

best_var_idx <- which(ga_adj2@solution == 1)
GA_trn <- df_trn[,c(best_var_idx,22)]
GA_tst <- df_val[,c(best_var_idx,22)]

GA_model <- lm(Life.expectancy ~ ., GA_trn)
summary(GA_model)

GA_coeff <- as.matrix(GA_model$coefficients,22,1)
GA_coeff
GA_model_prey <- predict(GA_model, newdata = GA_tst)
Perf_Table_ga[8,] <- perf_reg(GA_tst$Life.expectancy, GA_model_prey)
Perf_Table_ga_adj[8,] <- summary(GA_model)$adj.r.squared

# 10, 0.7, 0.3
start_time <- proc.time()
ga_adj2 <- ga(type="binary", fitness = fit_adj2, nBits = ncol(x),
              names = colnames(x),popSize = 10, pcrossover = 0.7, pmutation = 0.3,
              maxiter = 100, elitism = 2, seed = 2021, monitor = plot)
end_time <- proc.time()
tt <- end_time - start_time
Perf_Time[9,] <- tt[3]


best_var_idx <- which(ga_adj2@solution == 1)
GA_trn <- df_trn[,c(best_var_idx,22)]
GA_tst <- df_val[,c(best_var_idx,22)]

GA_model <- lm(Life.expectancy ~ ., GA_trn)
summary(GA_model)

GA_coeff <- as.matrix(GA_model$coefficients,22,1)
GA_coeff
GA_model_prey <- predict(GA_model, newdata = GA_tst)
Perf_Table_ga[9,] <- perf_reg(GA_tst$Life.expectancy, GA_model_prey)
Perf_Table_ga_adj[9,] <- summary(GA_model)$adj.r.squared

# 100, 0.1, 0.05
start_time <- proc.time()
ga_adj2 <- ga(type="binary", fitness = fit_adj2, nBits = ncol(x),
              names = colnames(x),popSize = 100, pcrossover = 0.1, pmutation = 0.05,
              maxiter = 100, elitism = 2, seed = 2021, monitor = plot)
end_time <- proc.time()
tt <- end_time - start_time
Perf_Time[10,] <- tt[3]


best_var_idx <- which(ga_adj2@solution == 1)
GA_trn <- df_trn[,c(best_var_idx,22)]
GA_tst <- df_val[,c(best_var_idx,22)]

GA_model <- lm(Life.expectancy ~ ., GA_trn)
summary(GA_model)

GA_coeff <- as.matrix(GA_model$coefficients,22,1)
GA_coeff
GA_model_prey <- predict(GA_model, newdata = GA_tst)
Perf_Table_ga[10,] <- perf_reg(GA_tst$Life.expectancy, GA_model_prey)
Perf_Table_ga_adj[10,] <- summary(GA_model)$adj.r.squared

# 100, 0.1, 0.1
start_time <- proc.time()
ga_adj2 <- ga(type="binary", fitness = fit_adj2, nBits = ncol(x),
              names = colnames(x),popSize = 100, pcrossover = 0.1, pmutation = 0.1,
              maxiter = 100, elitism = 2, seed = 2021, monitor = plot)
end_time <- proc.time()
tt <- end_time - start_time
Perf_Time[11,] <- tt[3]


best_var_idx <- which(ga_adj2@solution == 1)
GA_trn <- df_trn[,c(best_var_idx,22)]
GA_tst <- df_val[,c(best_var_idx,22)]

GA_model <- lm(Life.expectancy ~ ., GA_trn)
summary(GA_model)

GA_coeff <- as.matrix(GA_model$coefficients,22,1)
GA_coeff
GA_model_prey <- predict(GA_model, newdata = GA_tst)
Perf_Table_ga[11,] <- perf_reg(GA_tst$Life.expectancy, GA_model_prey)
Perf_Table_ga_adj[11,] <- summary(GA_model)$adj.r.squared

# 100,0.1, 0.3
start_time <- proc.time()
ga_adj2 <- ga(type="binary", fitness = fit_adj2, nBits = ncol(x),
              names = colnames(x),popSize = 100, pcrossover = 0.1, pmutation = 0.3,
              maxiter = 100, elitism = 2, seed = 2021, monitor = plot)
end_time <- proc.time()
tt <- end_time - start_time
Perf_Time[12,] <- tt[3]


best_var_idx <- which(ga_adj2@solution == 1)
GA_trn <- df_trn[,c(best_var_idx,22)]
GA_tst <- df_val[,c(best_var_idx,22)]

GA_model <- lm(Life.expectancy ~ ., GA_trn)
summary(GA_model)

GA_coeff <- as.matrix(GA_model$coefficients,22,1)
GA_coeff
GA_model_prey <- predict(GA_model, newdata = GA_tst)
Perf_Table_ga[12,] <- perf_reg(GA_tst$Life.expectancy, GA_model_prey)
Perf_Table_ga_adj[12,] <- summary(GA_model)$adj.r.squared

# 100, 0.4, 0.05
start_time <- proc.time()
ga_adj2 <- ga(type="binary", fitness = fit_adj2, nBits = ncol(x),
              names = colnames(x),popSize = 100, pcrossover = 0.4, pmutation = 0.05,
              maxiter = 100, elitism = 2, seed = 2021, monitor = plot)
end_time <- proc.time()
tt <- end_time - start_time
Perf_Time[13,] <- tt[3]


best_var_idx <- which(ga_adj2@solution == 1)
GA_trn <- df_trn[,c(best_var_idx,22)]
GA_tst <- df_val[,c(best_var_idx,22)]

GA_model <- lm(Life.expectancy ~ ., GA_trn)
summary(GA_model)

GA_coeff <- as.matrix(GA_model$coefficients,22,1)
GA_coeff
GA_model_prey <- predict(GA_model, newdata = GA_tst)
Perf_Table_ga[13,] <- perf_reg(GA_tst$Life.expectancy, GA_model_prey)
Perf_Table_ga_adj[13,] <- summary(GA_model)$adj.r.squared

# 100, 0.4, 0.1
start_time <- proc.time()
ga_adj2 <- ga(type="binary", fitness = fit_adj2, nBits = ncol(x),
              names = colnames(x),popSize = 100, pcrossover = 0.4, pmutation = 0.1,
              maxiter = 100, elitism = 2, seed = 2021, monitor = plot)
end_time <- proc.time()
tt <- end_time - start_time
Perf_Time[14,] <- tt[3]


best_var_idx <- which(ga_adj2@solution == 1)
GA_trn <- df_trn[,c(best_var_idx,22)]
GA_tst <- df_val[,c(best_var_idx,22)]

GA_model <- lm(Life.expectancy ~ ., GA_trn)
summary(GA_model)

GA_coeff <- as.matrix(GA_model$coefficients,22,1)
GA_coeff
GA_model_prey <- predict(GA_model, newdata = GA_tst)
Perf_Table_ga[14,] <- perf_reg(GA_tst$Life.expectancy, GA_model_prey)
Perf_Table_ga_adj[14,] <- summary(GA_model)$adj.r.squared

# 100, 0.4, 0.3
start_time <- proc.time()
ga_adj2 <- ga(type="binary", fitness = fit_adj2, nBits = ncol(x),
              names = colnames(x),popSize = 100, pcrossover = 0.4, pmutation = 0.3,
              maxiter = 100, elitism = 2, seed = 2021, monitor = plot)
end_time <- proc.time()
tt <- end_time - start_time
Perf_Time[15,] <- tt[3]


best_var_idx <- which(ga_adj2@solution == 1)
GA_trn <- df_trn[,c(best_var_idx,22)]
GA_tst <- df_val[,c(best_var_idx,22)]

GA_model <- lm(Life.expectancy ~ ., GA_trn)
summary(GA_model)

GA_coeff <- as.matrix(GA_model$coefficients,22,1)
GA_coeff
GA_model_prey <- predict(GA_model, newdata = GA_tst)
Perf_Table_ga[15,] <- perf_reg(GA_tst$Life.expectancy, GA_model_prey)
Perf_Table_ga_adj[15,] <- summary(GA_model)$adj.r.squared

# 100, 0.7, 0.05
start_time <- proc.time()
ga_adj2 <- ga(type="binary", fitness = fit_adj2, nBits = ncol(x),
              names = colnames(x),popSize = 100, pcrossover = 0.7, pmutation = 0.05,
              maxiter = 100, elitism = 2, seed = 2021, monitor = plot)
end_time <- proc.time()
tt <- end_time - start_time
Perf_Time[16,] <- tt[3]


best_var_idx <- which(ga_adj2@solution == 1)
GA_trn <- df_trn[,c(best_var_idx,22)]
GA_tst <- df_val[,c(best_var_idx,22)]

GA_model <- lm(Life.expectancy ~ ., GA_trn)
summary(GA_model)

GA_coeff <- as.matrix(GA_model$coefficients,22,1)
GA_coeff
GA_model_prey <- predict(GA_model, newdata = GA_tst)
Perf_Table_ga[16,] <- perf_reg(GA_tst$Life.expectancy, GA_model_prey)
Perf_Table_ga_adj[16,] <- summary(GA_model)$adj.r.squared

# 100, 0.7, 0.1
start_time <- proc.time()
ga_adj2 <- ga(type="binary", fitness = fit_adj2, nBits = ncol(x),
              names = colnames(x),popSize = 100, pcrossover = 0.7, pmutation = 0.1,
              maxiter = 100, elitism = 2, seed = 2021, monitor = plot)
end_time <- proc.time()
tt <- end_time - start_time
Perf_Time[17,] <- tt[3]


best_var_idx <- which(ga_adj2@solution == 1)
GA_trn <- df_trn[,c(best_var_idx,22)]
GA_tst <- df_val[,c(best_var_idx,22)]

GA_model <- lm(Life.expectancy ~ ., GA_trn)
summary(GA_model)

GA_coeff <- as.matrix(GA_model$coefficients,22,1)
GA_coeff
GA_model_prey <- predict(GA_model, newdata = GA_tst)
Perf_Table_ga[17,] <- perf_reg(GA_tst$Life.expectancy, GA_model_prey)
Perf_Table_ga_adj[17,] <- summary(GA_model)$adj.r.squared

# 100, 0.7, 0.3
start_time <- proc.time()
ga_adj2 <- ga(type="binary", fitness = fit_adj2, nBits = ncol(x),
              names = colnames(x),popSize = 100, pcrossover = 0.7, pmutation = 0.3,
              maxiter = 100, elitism = 2, seed = 2021, monitor = plot)
end_time <- proc.time()
tt <- end_time - start_time
Perf_Time[18,] <- tt[3]


best_var_idx <- which(ga_adj2@solution == 1)
GA_trn <- df_trn[,c(best_var_idx,22)]
GA_tst <- df_val[,c(best_var_idx,22)]

GA_model <- lm(Life.expectancy ~ ., GA_trn)
summary(GA_model)

GA_coeff <- as.matrix(GA_model$coefficients,22,1)
GA_coeff
GA_model_prey <- predict(GA_model, newdata = GA_tst)
Perf_Table_ga[18,] <- perf_reg(GA_tst$Life.expectancy, GA_model_prey)
Perf_Table_ga_adj[18,] <- summary(GA_model)$adj.r.squared

# 500, 0.1, 0.05
start_time <- proc.time()
ga_adj2 <- ga(type="binary", fitness = fit_adj2, nBits = ncol(x),
              names = colnames(x),popSize = 500, pcrossover = 0.1, pmutation = 0.05,
              maxiter = 100, elitism = 2, seed = 2021, monitor = plot)
end_time <- proc.time()
tt <- end_time - start_time
Perf_Time[19,] <- tt[3]


best_var_idx <- which(ga_adj2@solution == 1)
GA_trn <- df_trn[,c(best_var_idx,22)]
GA_tst <- df_val[,c(best_var_idx,22)]

GA_model <- lm(Life.expectancy ~ ., GA_trn)
summary(GA_model)

GA_coeff <- as.matrix(GA_model$coefficients,22,1)
GA_coeff
GA_model_prey <- predict(GA_model, newdata = GA_tst)
Perf_Table_ga[19,] <- perf_reg(GA_tst$Life.expectancy, GA_model_prey)
Perf_Table_ga_adj[19,] <- summary(GA_model)$adj.r.squared

# 500, 0.1, 0.1
start_time <- proc.time()
ga_adj2 <- ga(type="binary", fitness = fit_adj2, nBits = ncol(x),
              names = colnames(x),popSize = 500, pcrossover = 0.1, pmutation = 0.1,
              maxiter = 100, elitism = 2, seed = 2021, monitor = plot)
end_time <- proc.time()
tt <- end_time - start_time
Perf_Time[20,] <- tt[3]


best_var_idx <- which(ga_adj2@solution == 1)
GA_trn <- df_trn[,c(best_var_idx,22)]
GA_tst <- df_val[,c(best_var_idx,22)]

GA_model <- lm(Life.expectancy ~ ., GA_trn)
summary(GA_model)

GA_coeff <- as.matrix(GA_model$coefficients,22,1)
GA_coeff
GA_model_prey <- predict(GA_model, newdata = GA_tst)
Perf_Table_ga[20,] <- perf_reg(GA_tst$Life.expectancy, GA_model_prey)
Perf_Table_ga_adj[20,] <- summary(GA_model)$adj.r.squared

# 500, 0.1, 0.3
start_time <- proc.time()
ga_adj2 <- ga(type="binary", fitness = fit_adj2, nBits = ncol(x),
              names = colnames(x),popSize = 500, pcrossover = 0.1, pmutation = 0.3,
              maxiter = 100, elitism = 2, seed = 2021, monitor = plot)
end_time <- proc.time()
tt <- end_time - start_time
Perf_Time[21,] <- tt[3]


best_var_idx <- which(ga_adj2@solution == 1)
GA_trn <- df_trn[,c(best_var_idx,22)]
GA_tst <- df_val[,c(best_var_idx,22)]

GA_model <- lm(Life.expectancy ~ ., GA_trn)
summary(GA_model)

GA_coeff <- as.matrix(GA_model$coefficients,22,1)
GA_coeff
GA_model_prey <- predict(GA_model, newdata = GA_tst)
Perf_Table_ga[21,] <- perf_reg(GA_tst$Life.expectancy, GA_model_prey)
Perf_Table_ga_adj[21,] <- summary(GA_model)$adj.r.squared

# 500, 0.4, 0.05
start_time <- proc.time()
ga_adj2 <- ga(type="binary", fitness = fit_adj2, nBits = ncol(x),
              names = colnames(x),popSize = 500, pcrossover = 0.4, pmutation = 0.05,
              maxiter = 100, elitism = 2, seed = 2021, monitor = plot)
end_time <- proc.time()
tt <- end_time - start_time
Perf_Time[22,] <- tt[3]


best_var_idx <- which(ga_adj2@solution == 1)
GA_trn <- df_trn[,c(best_var_idx,22)]
GA_tst <- df_val[,c(best_var_idx,22)]

GA_model <- lm(Life.expectancy ~ ., GA_trn)
summary(GA_model)

GA_coeff <- as.matrix(GA_model$coefficients,22,1)
GA_coeff
GA_model_prey <- predict(GA_model, newdata = GA_tst)
Perf_Table_ga[22,] <- perf_reg(GA_tst$Life.expectancy, GA_model_prey)
Perf_Table_ga_adj[22,] <- summary(GA_model)$adj.r.squared

# 500, 0.4, 0.1
start_time <- proc.time()
ga_adj2 <- ga(type="binary", fitness = fit_adj2, nBits = ncol(x),
              names = colnames(x),popSize = 500, pcrossover = 0.4, pmutation = 0.1,
              maxiter = 100, elitism = 2, seed = 2021, monitor = plot)
end_time <- proc.time()
tt <- end_time - start_time
Perf_Time[23,] <- tt[3]


best_var_idx <- which(ga_adj2@solution == 1)
GA_trn <- df_trn[,c(best_var_idx,22)]
GA_tst <- df_val[,c(best_var_idx,22)]

GA_model <- lm(Life.expectancy ~ ., GA_trn)
summary(GA_model)

GA_coeff <- as.matrix(GA_model$coefficients,22,1)
GA_coeff
GA_model_prey <- predict(GA_model, newdata = GA_tst)
Perf_Table_ga[23,] <- perf_reg(GA_tst$Life.expectancy, GA_model_prey)
Perf_Table_ga_adj[23,] <- summary(GA_model)$adj.r.squared

# 500, 0.4, 0.3
start_time <- proc.time()
ga_adj2 <- ga(type="binary", fitness = fit_adj2, nBits = ncol(x),
              names = colnames(x),popSize = 500, pcrossover = 0.4, pmutation = 0.3,
              maxiter = 100, elitism = 2, seed = 2021, monitor = plot)
end_time <- proc.time()
tt <- end_time - start_time
Perf_Time[24,] <- tt[3]


best_var_idx <- which(ga_adj2@solution == 1)
GA_trn <- df_trn[,c(best_var_idx,22)]
GA_tst <- df_val[,c(best_var_idx,22)]

GA_model <- lm(Life.expectancy ~ ., GA_trn)
summary(GA_model)

GA_coeff <- as.matrix(GA_model$coefficients,22,1)
GA_coeff
GA_model_prey <- predict(GA_model, newdata = GA_tst)
Perf_Table_ga[24,] <- perf_reg(GA_tst$Life.expectancy, GA_model_prey)
Perf_Table_ga_adj[24,] <- summary(GA_model)$adj.r.squared

# 500, 0.7, 0.05
start_time <- proc.time()
ga_adj2 <- ga(type="binary", fitness = fit_adj2, nBits = ncol(x),
              names = colnames(x),popSize = 500, pcrossover = 0.7, pmutation = 0.05,
              maxiter = 100, elitism = 2, seed = 2021, monitor = plot)
end_time <- proc.time()
tt <- end_time - start_time
Perf_Time[25,] <- tt[3]


best_var_idx <- which(ga_adj2@solution == 1)
GA_trn <- df_trn[,c(best_var_idx,22)]
GA_tst <- df_val[,c(best_var_idx,22)]

GA_model <- lm(Life.expectancy ~ ., GA_trn)
summary(GA_model)

GA_coeff <- as.matrix(GA_model$coefficients,22,1)
GA_coeff
GA_model_prey <- predict(GA_model, newdata = GA_tst)
Perf_Table_ga[25,] <- perf_reg(GA_tst$Life.expectancy, GA_model_prey)
Perf_Table_ga_adj[25,] <- summary(GA_model)$adj.r.squared

# 500, 0.7, 0.1
start_time <- proc.time()
ga_adj2 <- ga(type="binary", fitness = fit_adj2, nBits = ncol(x),
              names = colnames(x),popSize = 500, pcrossover = 0.7, pmutation = 0.1,
              maxiter = 100, elitism = 2, seed = 2021, monitor = plot)
end_time <- proc.time()
tt <- end_time - start_time
Perf_Time[26,] <- tt[3]


best_var_idx <- which(ga_adj2@solution == 1)
GA_trn <- df_trn[,c(best_var_idx,22)]
GA_tst <- df_val[,c(best_var_idx,22)]

GA_model <- lm(Life.expectancy ~ ., GA_trn)
summary(GA_model)

GA_coeff <- as.matrix(GA_model$coefficients,22,1)
GA_coeff
GA_model_prey <- predict(GA_model, newdata = GA_tst)
Perf_Table_ga[26,] <- perf_reg(GA_tst$Life.expectancy, GA_model_prey)
Perf_Table_ga_adj[26,] <- summary(GA_model)$adj.r.squared

# 500, 0.7, 0.3
start_time <- proc.time()
ga_adj2 <- ga(type="binary", fitness = fit_adj2, nBits = ncol(x),
              names = colnames(x),popSize = 500, pcrossover = 0.7, pmutation = 0.3,
              maxiter = 100, elitism = 2, seed = 2021, monitor = plot)
end_time <- proc.time()
tt <- end_time - start_time
Perf_Time[27,] <- tt[3]


best_var_idx <- which(ga_adj2@solution == 1)
GA_trn <- df_trn[,c(best_var_idx,22)]
GA_tst <- df_val[,c(best_var_idx,22)]

GA_model <- lm(Life.expectancy ~ ., GA_trn)
summary(GA_model)

GA_coeff <- as.matrix(GA_model$coefficients,22,1)
GA_coeff
GA_model_prey <- predict(GA_model, newdata = GA_tst)
Perf_Table_ga[27,] <- perf_reg(GA_tst$Life.expectancy, GA_model_prey)
Perf_Table_ga_adj[27,] <- summary(GA_model)$adj.r.squared

Perf_Table_ga
Perf_Table_ga_adj
Perf_Time

# Logistic Regression
#[Q4]
data <- read.csv("heart_failure_clinical_records_dataset.csv")
head(data)
table(is.na(data))

perf_eval <- function(cm){
  TPR <- cm[2,2]/sum(cm[2,])
  PRE <- cm[2,2]/sum(cm[,2])
  TNR <- cm[1,1]/sum(cm[1,])
  ACC <- (cm[1,1]+cm[2,2])/sum(cm)
  BCR <- sqrt(TPR*TNR)
  F1 <- 2*TPR*PRE/(TPR+PRE)
  
  return(c(ACC, BCR, F1))
}

#AUROC
AUROC <- function(lr_response, lr_target){
  res_sor <- sort(lr_response, decreasing=TRUE)
  AUROC = 0                                                  
  x <- c(0)                                                  
  y <- c(0)
  
  for (i in c(2:length(res_sor)+1)){                            
    lr_target_1 <- lr_target
    lr_pred_1 <- rep(0, length(lr_target_1))                    
    lr_pred_1[which(lr_response > res_sor[i-1])] <- 1
    
    cm <- table(lr_target_1, lr_pred_1)
    TPR_1 <- cm[2,2]/sum(cm[2,])
    FPR_1 <- cm[1,2]/sum(cm[1,])
    
    if (i==length(res_sor)+1){                                  
      FPR_2 <- 1
      TPR_2 <- 1   
    } else{                                               
      lr_pred_2 <- rep(0, length(lr_target_1))
      lr_pred_2[which(lr_response > res_sor[i])] <- 1
      
      cm_2 <- table(lr_target_1, lr_pred_2)
      TPR_2 <- cm_2[2,2]/sum(cm_2[2,])
      FPR_2 <- cm_2[1,2]/sum(cm_2[1,])
    }
    
    if(FPR_2 == FPR_1){
      AUROC = AUROC
    }else{
      AUROC = AUROC + ((FPR_2-FPR_1)*(TPR_2))
    }
    x <- c(x,FPR_2)                                          
    y <- c(y,TPR_2)
  }
  plot(x, y, type='l',xlab="False Positive Rate(FPR)", ylab="True Positive Rate(TPR)")
  
  return(AUROC)
}

# 비교를 위한 Matrix생성
Perf_Table_lr <- matrix(0, nrow=4, ncol=3)
Perf_tr_auroc <- matrix(0, nrow=4, ncol=1)
Perf_ts_auroc <- matrix(0, nrow=4, ncol=1)
rownames(Perf_Table_lr) <- c("Forward", "Backward","Stepwise","GA")
colnames(Perf_Table_lr) <- c("ACC","BCR","F1")
rownames(Perf_tr_auroc) <- c("Forward", "Backward","Stepwise","GA")
colnames(Perf_tr_auroc) <- c("Train AUROC")
rownames(Perf_ts_auroc) <- c("Forward", "Backward","Stepwise","GA")
colnames(Perf_ts_auroc) <- c("Test AUROC")

# 수치형 변수에 대해서 nomarlization 실행
input_num_idx <- c(1,3,5,7,8,9,12)
input_cat_idx <- c(2,4,6,10,11)

data_input <- data[,input_num_idx]
data_cat_input <- data[,input_cat_idx]

data_input <- scale(data_input, center=TRUE, scale = TRUE)

data_input <- data.frame(data_input, data_cat_input)
data_target <- data[,13]

heart_df <- data.frame(data_input, data_target)
heart_df <- rename(heart_df,"DEATH_EVENT"="data_target")

# 분할
set.seed(506)
trn_idx <- sample(1:nrow(heart_df), round(0.7*nrow(heart_df)))
data_trn <- heart_df[trn_idx,]
data_val <- heart_df[-trn_idx,]


# Forward
tmp_x_lr <- paste(colnames(data_trn)[-13], collapse="+")
tmp_xy_lr <- paste("DEATH_EVENT ~ ", tmp_x_lr, collapse="")
as.formula(tmp_xy_lr)
start_time <- proc.time()
forward_model_lr <- step(glm(DEATH_EVENT ~ 1, data = data_trn),
                         scope=list(upper=as.formula(tmp_xy_lr),
                                    lower=DEATH_EVENT ~ 1), direction="forward", trace=1)
end_time <- proc.time()
end_time - start_time
summary(forward_model_lr)
forward_lr_coef  <- as.matrix(forward_model_lr$coefficients,13,1)
forward_lr_coef

# Train AUROC
forward_lr_prob_t <- predict(forward_model_lr, type="response", newdata=data_trn)
forward_lr_prey_t <- rep(0,nrow(data_trn))
forward_lr_prey_t[which(forward_lr_prob_t >= 0.5)] <- 1
fs_a <- AUROC(forward_lr_prob_t, data_trn$DEATH_EVENT)
fs_a
Perf_tr_auroc[1,] <- fs_a

# predict
forward_lr_prob <- predict(forward_model_lr, type="response", newdata=data_val)
forward_lr_prey <- rep(0,nrow(data_val))
forward_lr_prey[which(forward_lr_prob >= 0.5)] <- 1
fs_a_ts <- AUROC(forward_lr_prob,data_val$DEATH_EVENT)
Perf_ts_auroc[1,] <- fs_a_ts
forward_lr_cm <- table(data_val$DEATH_EVENT, forward_lr_prey)
forward_lr_cm
Perf_Table_lr[1,] <- perf_eval(forward_lr_cm)

# Backward
start_time <- proc.time()
backward_model_lr <- step(glm(DEATH_EVENT ~ ., family = binomial, data_trn),
                          scope=list(upper = as.formula(tmp_xy_lr),
                                     lower=DEATH_EVENT ~ 1), direction="backward", trace = 1)
end_time <- proc.time()
end_time - start_time
summary(backward_model_lr)
backward_lr_coef <- as.matrix(backward_model_lr$coefficients,13,1)
backward_lr_coef

backward_lr_prob_t <- predict(backward_model_lr, type="response", newdata=data_trn)
be_a <- AUROC(backward_lr_prob_t,data_trn$DEATH_EVENT)
be_a
Perf_tr_auroc[2,] <- be_a

backward_lr_prob <- predict(backward_model_lr, type="response",newdata=data_val)
backward_lr_prey <- rep(0, nrow(data_val))
backward_lr_prey[which(backward_lr_prob >= 0.5)] <- 1
be_a_ts <- AUROC(backward_lr_prob,data_val$DEATH_EVENT)
Perf_ts_auroc[2,] <- be_a_ts

backward_lr_cm <- table(data_val$DEATH_EVENT, backward_lr_prey)
backward_lr_cm
Perf_Table_lr[2,] <- perf_eval(backward_lr_cm)

# Stepwise
start_time <- proc.time()
stepwise_model_lr <- step(glm(DEATH_EVENT ~ 1, data = data_trn),
                          scope=list(upper=as.formula(tmp_xy_lr),
                                     lower=DEATH_EVENT ~ 1), direction = "both", trace=1)
end_time <- proc.time()
end_time - start_time
summary(stepwise_model_lr)
stepwise_lr_coef <- as.matrix(stepwise_model_lr$coefficients,13,1)
stepwise_lr_coef

stepwise_lr_prob_t <- predict(stepwise_model_lr, type="response", newdata=data_trn)
ss_a <- AUROC(stepwise_lr_prob_t,data_trn$DEATH_EVENT)
ss_a
Perf_tr_auroc[3,] <- ss_a

stepwise_lr_prob <- predict(stepwise_model_lr, type="response",newdata=data_val)
stepwise_lr_prey <- rep(0, nrow(data_val))
stepwise_lr_prey[which(stepwise_lr_prob >= 0.5)] <- 1
ss_a_ts <- AUROC(stepwise_lr_prob,data_val$DEATH_EVENT)
Perf_ts_auroc[3,] <- ss_a_ts

stepwise_lr_cm <- table(data_val$DEATH_EVENT, stepwise_lr_prey)
stepwise_lr_cm
Perf_Table_lr[3,] <- perf_eval(stepwise_lr_cm)
Perf_Table_lr
Perf_tr_auroc
Perf_ts_auroc

#[Q5]
# 그래프 없는 AUROC
AUROC_2 <- function(lr_response, lr_target){
  res_sor <- sort(lr_response, decreasing=TRUE)
  AUROC = 0
  x <- c(0)
  y <- c(0)

  for (i in c(2:length(res_sor)+1)){
    lr_target_1 <- lr_target
    lr_pred_1 <- rep(0, length(lr_target_1))
    lr_pred_1[which(lr_response > res_sor[i-1])] <- 1

    cm <- table(lr_target_1, lr_pred_1)
    TPR_1 <- cm[2,2]/sum(cm[2,])
    FPR_1 <- cm[1,2]/sum(cm[1,])

    if (i==length(res_sor)+1){
      FPR_2 <- 1
      TPR_2 <- 1
    } else{
      lr_pred_2 <- rep(0, length(lr_target_1))
      lr_pred_2[which(lr_response > res_sor[i])] <- 1

      cm_2 <- table(lr_target_1, lr_pred_2)
      TPR_2 <- cm_2[2,2]/sum(cm_2[2,])
      FPR_2 <- cm_2[1,2]/sum(cm_2[1,])
    }

    if(FPR_2 == FPR_1){
      AUROC = AUROC
    }else{
      AUROC = AUROC + ((FPR_2-FPR_1)*(TPR_2))
    }
    x <- c(x,FPR_2)
    y <- c(y,TPR_2)
  }

  return(AUROC)
}


fit_AUROC <- function(string){
  sel_var_idx <- which(string == 1)
  sel_x <- x[,sel_var_idx]
  xy <- data.frame(sel_x,y)
  GA_lr <- glm(y ~ ., family = binomial, data=xy)
  GA_lr_prob <- predict(GA_lr, type="response",newdata=xy)
  GA_auroc <- AUROC_2(GA_lr_prob,y)
  return(GA_auroc)
}


x <- as.matrix(data_trn[,-13])
y <- data_trn[,13]


start_time <- proc.time()
GA_model_lr <- ga(type="binary", fitness = fit_AUROC, nBits = ncol(x),
                  names = colnames(x), popSize = 50, pcrossover = 0.5, pmutation = 0.01,
                  maxiter = 100, elitism = 2, seed = 2021)
end_time <- proc.time()
end_time - start_time

best_var_idx <- which(GA_model_lr@solution == 1)
GA_trn <- data_trn[,c(best_var_idx,13)]
GA_tst <- data_val[,c(best_var_idx,13)]

GA_model <- glm(DEATH_EVENT ~ ., family = binomial, GA_trn)
summary(GA_model)

GA_coeff <- as.matrix(GA_model$coefficients,13,1)
GA_coeff

#Make prediction
GA_model_prob <- predict(GA_model, type="response", newdata = GA_tst)
GA_model_prey <- rep(0, nrow(data_val))
GA_model_prey[which(GA_model_prob >= 0.5)] <- 1

GA_a_ts <- AUROC_2(GA_model_prob, data_val$DEATH_EVENT)
GA_a_ts
Perf_ts_auroc[4,] <- GA_a_ts

GA_lr_cm <- table(data_val$DEATH_EVENT, GA_model_prey)
GA_lr_cm
Perf_Table_lr[4,] <- perf_eval(GA_lr_cm)
Perf_Table_lr
Perf_ts_auroc


#[Q6]
# 비교를 위한 matrix 생성
Perf_Table_ga <- matrix(0, nrow = 27, ncol = 3)
Perf_lr_time <- matrix(0, nrow=27, ncol = 1)
Perf_lr_auroc <- matrix(0, nrow=27, ncol = 1)
rownames(Perf_Table_ga) <- c(1:27)
colnames(Perf_Table_ga) <- c("Accuracy", "BCR", "F1-Measure")
rownames(Perf_lr_time) <- c(1:27)
colnames(Perf_lr_time) <- c("Time")
rownames(Perf_lr_auroc) <- c(1:27)
colnames(Perf_lr_auroc) <- c("AUROC")

# maxiter(50,100,300), pcrossover는 (0.1,0,4,0.7), pmutation은 (0.05,0.1,0.3) 이렇게 3가지씩으로 총 27개의 조합을 만든다.
# 50, 0.1, 0.05
start_time <- proc.time()
GA_model_lr <- ga(type="binary", fitness = fit_AUROC, nBits = ncol(x),
              names = colnames(x),popSize = 50, pcrossover = 0.1, pmutation = 0.05,
              maxiter = 50, elitism = 2, seed = 2021, monitor = plot)
end_time <- proc.time()
tt <- end_time - start_time

best_var_idx <- which(GA_model_lr@solution == 1)
GA_trn <- data_trn[,c(best_var_idx,13)]
GA_tst <- data_val[,c(best_var_idx,13)]

GA_model <- glm(DEATH_EVENT ~ ., family = binomial, GA_trn)
summary(GA_model)

GA_coeff <- as.matrix(GA_model$coefficients,13,1)
GA_coeff

GA_model_prob <- predict(GA_model, type="response", newdata = GA_tst)
GA_model_prey <- rep(0, nrow(data_val))
GA_model_prey[which(GA_model_prob >= 0.5)] <- 1

GA_a_ts <- AUROC_2(GA_model_prob, data_val$DEATH_EVENT)
GA_lr_cm <- table(data_val$DEATH_EVENT, GA_model_prey)

Perf_lr_time[1,] <- tt[3]
Perf_lr_auroc[1,] <- GA_a_ts
Perf_Table_ga[1,] <- perf_eval(GA_lr_cm)

# 50, 0.1, 0.1
start_time <- proc.time()
GA_model_lr <- ga(type="binary", fitness = fit_AUROC, nBits = ncol(x),
                  names = colnames(x),popSize = 50, pcrossover = 0.1, pmutation = 0.1,
                  maxiter = 50, elitism = 2, seed = 2021, monitor = plot)
end_time <- proc.time()
tt <- end_time - start_time

best_var_idx <- which(GA_model_lr@solution == 1)
GA_trn <- data_trn[,c(best_var_idx,13)]
GA_tst <- data_val[,c(best_var_idx,13)]

GA_model <- glm(DEATH_EVENT ~ ., family = binomial, GA_trn)
summary(GA_model)

GA_coeff <- as.matrix(GA_model$coefficients,13,1)
GA_coeff

GA_model_prob <- predict(GA_model, type="response", newdata = GA_tst)
GA_model_prey <- rep(0, nrow(data_val))
GA_model_prey[which(GA_model_prob >= 0.5)] <- 1

GA_a_ts <- AUROC_2(GA_model_prob, data_val$DEATH_EVENT)
GA_lr_cm <- table(data_val$DEATH_EVENT, GA_model_prey)

Perf_lr_time[2,] <- tt[3]
Perf_lr_auroc[2,] <- GA_a_ts
Perf_Table_ga[2,] <- perf_eval(GA_lr_cm)

# 50, 0.1, 0.3
start_time <- proc.time()
GA_model_lr <- ga(type="binary", fitness = fit_AUROC, nBits = ncol(x),
                  names = colnames(x),popSize = 50, pcrossover = 0.1, pmutation = 0.3,
                  maxiter = 50, elitism = 2, seed = 2021, monitor = plot)
end_time <- proc.time()
tt <- end_time - start_time

best_var_idx <- which(GA_model_lr@solution == 1)
GA_trn <- data_trn[,c(best_var_idx,13)]
GA_tst <- data_val[,c(best_var_idx,13)]

GA_model <- glm(DEATH_EVENT ~ ., family = binomial, GA_trn)
summary(GA_model)

GA_coeff <- as.matrix(GA_model$coefficients,13,1)
GA_coeff

GA_model_prob <- predict(GA_model, type="response", newdata = GA_tst)
GA_model_prey <- rep(0, nrow(data_val))
GA_model_prey[which(GA_model_prob >= 0.5)] <- 1

GA_a_ts <- AUROC_2(GA_model_prob, data_val$DEATH_EVENT)
GA_lr_cm <- table(data_val$DEATH_EVENT, GA_model_prey)

Perf_lr_time[3,] <- tt[3]
Perf_lr_auroc[3,] <- GA_a_ts
Perf_Table_ga[3,] <- perf_eval(GA_lr_cm)

# 50, 0.4, 0.05
start_time <- proc.time()
GA_model_lr <- ga(type="binary", fitness = fit_AUROC, nBits = ncol(x),
                  names = colnames(x),popSize = 50, pcrossover = 0.4, pmutation = 0.05,
                  maxiter = 50, elitism = 2, seed = 2021, monitor = plot)
end_time <- proc.time()
tt <- end_time - start_time

best_var_idx <- which(GA_model_lr@solution == 1)
GA_trn <- data_trn[,c(best_var_idx,13)]
GA_tst <- data_val[,c(best_var_idx,13)]

GA_model <- glm(DEATH_EVENT ~ ., family = binomial, GA_trn)
summary(GA_model)

GA_coeff <- as.matrix(GA_model$coefficients,13,1)
GA_coeff

GA_model_prob <- predict(GA_model, type="response", newdata = GA_tst)
GA_model_prey <- rep(0, nrow(data_val))
GA_model_prey[which(GA_model_prob >= 0.5)] <- 1

GA_a_ts <- AUROC_2(GA_model_prob, data_val$DEATH_EVENT)
GA_lr_cm <- table(data_val$DEATH_EVENT, GA_model_prey)

Perf_lr_time[4,] <- tt[3]
Perf_lr_auroc[4,] <- GA_a_ts
Perf_Table_ga[4,] <- perf_eval(GA_lr_cm)

# 50, 0.4, 0.1
start_time <- proc.time()
GA_model_lr <- ga(type="binary", fitness = fit_AUROC, nBits = ncol(x),
                  names = colnames(x),popSize = 50, pcrossover = 0.4, pmutation = 0.1,
                  maxiter = 50, elitism = 2, seed = 2021, monitor = plot)
end_time <- proc.time()
tt <- end_time - start_time

best_var_idx <- which(GA_model_lr@solution == 1)
GA_trn <- data_trn[,c(best_var_idx,13)]
GA_tst <- data_val[,c(best_var_idx,13)]

GA_model <- glm(DEATH_EVENT ~ ., family = binomial, GA_trn)
summary(GA_model)

GA_coeff <- as.matrix(GA_model$coefficients,13,1)
GA_coeff

GA_model_prob <- predict(GA_model, type="response", newdata = GA_tst)
GA_model_prey <- rep(0, nrow(data_val))
GA_model_prey[which(GA_model_prob >= 0.5)] <- 1

GA_a_ts <- AUROC_2(GA_model_prob, data_val$DEATH_EVENT)
GA_lr_cm <- table(data_val$DEATH_EVENT, GA_model_prey)

Perf_lr_time[5,] <- tt[3]
Perf_lr_auroc[5,] <- GA_a_ts
Perf_Table_ga[5,] <- perf_eval(GA_lr_cm)

# 50, 0.4, 0.3
start_time <- proc.time()
GA_model_lr <- ga(type="binary", fitness = fit_AUROC, nBits = ncol(x),
                  names = colnames(x),popSize = 50, pcrossover = 0.4, pmutation = 0.3,
                  maxiter = 50, elitism = 2, seed = 2021, monitor = plot)
end_time <- proc.time()
tt <- end_time - start_time

best_var_idx <- which(GA_model_lr@solution == 1)
GA_trn <- data_trn[,c(best_var_idx,13)]
GA_tst <- data_val[,c(best_var_idx,13)]

GA_model <- glm(DEATH_EVENT ~ ., family = binomial, GA_trn)
summary(GA_model)

GA_coeff <- as.matrix(GA_model$coefficients,13,1)
GA_coeff

GA_model_prob <- predict(GA_model, type="response", newdata = GA_tst)
GA_model_prey <- rep(0, nrow(data_val))
GA_model_prey[which(GA_model_prob >= 0.5)] <- 1

GA_a_ts <- AUROC_2(GA_model_prob, data_val$DEATH_EVENT)
GA_lr_cm <- table(data_val$DEATH_EVENT, GA_model_prey)

Perf_lr_time[6,] <- tt[3]
Perf_lr_auroc[6,] <- GA_a_ts
Perf_Table_ga[6,] <- perf_eval(GA_lr_cm)

# 50, 0.7, 0.05
start_time <- proc.time()
GA_model_lr <- ga(type="binary", fitness = fit_AUROC, nBits = ncol(x),
                  names = colnames(x),popSize = 50, pcrossover = 0.7, pmutation = 0.05,
                  maxiter = 50, elitism = 2, seed = 2021, monitor = plot)
end_time <- proc.time()
tt <- end_time - start_time

best_var_idx <- which(GA_model_lr@solution == 1)
GA_trn <- data_trn[,c(best_var_idx,13)]
GA_tst <- data_val[,c(best_var_idx,13)]

GA_model <- glm(DEATH_EVENT ~ ., family = binomial, GA_trn)
summary(GA_model)

GA_coeff <- as.matrix(GA_model$coefficients,13,1)
GA_coeff

GA_model_prob <- predict(GA_model, type="response", newdata = GA_tst)
GA_model_prey <- rep(0, nrow(data_val))
GA_model_prey[which(GA_model_prob >= 0.5)] <- 1

GA_a_ts <- AUROC_2(GA_model_prob, data_val$DEATH_EVENT)
GA_lr_cm <- table(data_val$DEATH_EVENT, GA_model_prey)

Perf_lr_time[7,] <- tt[3]
Perf_lr_auroc[7,] <- GA_a_ts
Perf_Table_ga[7,] <- perf_eval(GA_lr_cm)

# 50, 0.7, 0.1
start_time <- proc.time()
GA_model_lr <- ga(type="binary", fitness = fit_AUROC, nBits = ncol(x),
                  names = colnames(x),popSize = 50, pcrossover = 0.7, pmutation = 0.1,
                  maxiter = 50, elitism = 2, seed = 2021, monitor = plot)
end_time <- proc.time()
tt <- end_time - start_time


best_var_idx <- which(GA_model_lr@solution == 1)
GA_trn <- data_trn[,c(best_var_idx,13)]
GA_tst <- data_val[,c(best_var_idx,13)]

GA_model <- glm(DEATH_EVENT ~ ., family = binomial, GA_trn)
summary(GA_model)

GA_coeff <- as.matrix(GA_model$coefficients,13,1)
GA_coeff

GA_model_prob <- predict(GA_model, type="response", newdata = GA_tst)
GA_model_prey <- rep(0, nrow(data_val))
GA_model_prey[which(GA_model_prob >= 0.5)] <- 1

GA_a_ts <- AUROC_2(GA_model_prob, data_val$DEATH_EVENT)
GA_lr_cm <- table(data_val$DEATH_EVENT, GA_model_prey)

Perf_lr_time[8,] <- tt[3]
Perf_lr_auroc[8,] <- GA_a_ts
Perf_Table_ga[8,] <- perf_eval(GA_lr_cm)

# 50, 0.7, 0.3
start_time <- proc.time()
GA_model_lr <- ga(type="binary", fitness = fit_AUROC, nBits = ncol(x),
                  names = colnames(x),popSize = 50, pcrossover = 0.7, pmutation = 0.3,
                  maxiter = 50, elitism = 2, seed = 2021, monitor = plot)
end_time <- proc.time()
tt <- end_time - start_time

best_var_idx <- which(GA_model_lr@solution == 1)

GA_trn <- data_trn[,c(best_var_idx,13)]
GA_tst <- data_val[,c(best_var_idx,13)]

GA_model <- glm(DEATH_EVENT ~ ., family = binomial, GA_trn)
summary(GA_model)

GA_coeff <- as.matrix(GA_model$coefficients,13,1)
GA_coeff

GA_model_prob <- predict(GA_model, type="response", newdata = GA_tst)
GA_model_prey <- rep(0, nrow(data_val))
GA_model_prey[which(GA_model_prob >= 0.5)] <- 1

GA_a_ts <- AUROC_2(GA_model_prob, data_val$DEATH_EVENT)
GA_lr_cm <- table(data_val$DEATH_EVENT, GA_model_prey)

Perf_lr_time[9,] <- tt[3]
Perf_lr_auroc[9,] <- GA_a_ts
Perf_Table_ga[9,] <- perf_eval(GA_lr_cm)
# 100, 0.1, 0.05
start_time <- proc.time()
GA_model_lr <- ga(type="binary", fitness = fit_AUROC, nBits = ncol(x),
                  names = colnames(x),popSize = 50, pcrossover = 0.1, pmutation = 0.05,
                  maxiter = 100, elitism = 2, seed = 2021, monitor = plot)
end_time <- proc.time()
tt <- end_time - start_time

best_var_idx <- which(GA_model_lr@solution == 1)

GA_trn <- data_trn[,c(best_var_idx,13)]
GA_tst <- data_val[,c(best_var_idx,13)]

GA_model <- glm(DEATH_EVENT ~ ., family = binomial, GA_trn)
summary(GA_model)

GA_coeff <- as.matrix(GA_model$coefficients,13,1)
GA_coeff

GA_model_prob <- predict(GA_model, type="response", newdata = GA_tst)
GA_model_prey <- rep(0, nrow(data_val))
GA_model_prey[which(GA_model_prob >= 0.5)] <- 1

GA_a_ts <- AUROC_2(GA_model_prob, data_val$DEATH_EVENT)
GA_lr_cm <- table(data_val$DEATH_EVENT, GA_model_prey)

Perf_lr_time[10,] <- tt[3]
Perf_lr_auroc[10,] <- GA_a_ts
Perf_Table_ga[10,] <- perf_eval(GA_lr_cm)

# 100, 0.1, 0.1
start_time <- proc.time()
GA_model_lr <- ga(type="binary", fitness = fit_AUROC, nBits = ncol(x),
                  names = colnames(x),popSize = 50, pcrossover = 0.1, pmutation = 0.1,
                  maxiter = 100, elitism = 2, seed = 2021, monitor = plot)
end_time <- proc.time()
tt <- end_time - start_time


best_var_idx <- which(GA_model_lr@solution == 1)
GA_trn <- data_trn[,c(best_var_idx,13)]
GA_tst <- data_val[,c(best_var_idx,13)]

GA_model <- glm(DEATH_EVENT ~ ., family = binomial, GA_trn)
summary(GA_model)

GA_coeff <- as.matrix(GA_model$coefficients,13,1)
GA_coeff

GA_model_prob <- predict(GA_model, type="response", newdata = GA_tst)
GA_model_prey <- rep(0, nrow(data_val))
GA_model_prey[which(GA_model_prob >= 0.5)] <- 1

GA_a_ts <- AUROC_2(GA_model_prob, data_val$DEATH_EVENT)
GA_lr_cm <- table(data_val$DEATH_EVENT, GA_model_prey)

Perf_lr_time[11,] <- tt[3]
Perf_lr_auroc[11,] <- GA_a_ts
Perf_Table_ga[11,] <- perf_eval(GA_lr_cm)

# 100,0.1, 0.3
start_time <- proc.time()
GA_model_lr <- ga(type="binary", fitness = fit_AUROC, nBits = ncol(x),
                  names = colnames(x),popSize = 50, pcrossover = 0.1, pmutation = 0.3,
                  maxiter = 100, elitism = 2, seed = 2021, monitor = plot)
end_time <- proc.time()
tt <- end_time - start_time


best_var_idx <- which(GA_model_lr@solution == 1)
GA_trn <- data_trn[,c(best_var_idx,13)]
GA_tst <- data_val[,c(best_var_idx,13)]

GA_model <- glm(DEATH_EVENT ~ ., family = binomial, GA_trn)
summary(GA_model)

GA_coeff <- as.matrix(GA_model$coefficients,13,1)
GA_coeff

GA_model_prob <- predict(GA_model, type="response", newdata = GA_tst)
GA_model_prey <- rep(0, nrow(data_val))
GA_model_prey[which(GA_model_prob >= 0.5)] <- 1

GA_a_ts <- AUROC_2(GA_model_prob, data_val$DEATH_EVENT)
GA_lr_cm <- table(data_val$DEATH_EVENT, GA_model_prey)

Perf_lr_time[12,] <- tt[3]
Perf_lr_auroc[12,] <- GA_a_ts
Perf_Table_ga[12,] <- perf_eval(GA_lr_cm)

# 100, 0.4, 0.05
start_time <- proc.time()
GA_model_lr <- ga(type="binary", fitness = fit_AUROC, nBits = ncol(x),
                  names = colnames(x),popSize = 50, pcrossover = 0.4, pmutation = 0.05,
                  maxiter = 100, elitism = 2, seed = 2021, monitor = plot)
end_time <- proc.time()
tt <- end_time - start_time


best_var_idx <- which(GA_model_lr@solution == 1)
GA_trn <- data_trn[,c(best_var_idx,13)]
GA_tst <- data_val[,c(best_var_idx,13)]

GA_model <- glm(DEATH_EVENT ~ ., family = binomial, GA_trn)
summary(GA_model)

GA_coeff <- as.matrix(GA_model$coefficients,13,1)
GA_coeff

GA_model_prob <- predict(GA_model, type="response", newdata = GA_tst)
GA_model_prey <- rep(0, nrow(data_val))
GA_model_prey[which(GA_model_prob >= 0.5)] <- 1

GA_a_ts <- AUROC_2(GA_model_prob, data_val$DEATH_EVENT)
GA_lr_cm <- table(data_val$DEATH_EVENT, GA_model_prey)

Perf_lr_time[13,] <- tt[3]
Perf_lr_auroc[13,] <- GA_a_ts
Perf_Table_ga[13,] <- perf_eval(GA_lr_cm)

# 100, 0.4, 0.1
start_time <- proc.time()
GA_model_lr <- ga(type="binary", fitness = fit_AUROC, nBits = ncol(x),
                  names = colnames(x),popSize = 50, pcrossover = 0.4, pmutation = 0.1,
                  maxiter = 100, elitism = 2, seed = 2021, monitor = plot)
end_time <- proc.time()
tt <- end_time - start_time


best_var_idx <- which(GA_model_lr@solution == 1)
GA_trn <- data_trn[,c(best_var_idx,13)]
GA_tst <- data_val[,c(best_var_idx,13)]

GA_model <- glm(DEATH_EVENT ~ ., family = binomial, GA_trn)
summary(GA_model)

GA_coeff <- as.matrix(GA_model$coefficients,13,1)
GA_coeff

GA_model_prob <- predict(GA_model, type="response", newdata = GA_tst)
GA_model_prey <- rep(0, nrow(data_val))
GA_model_prey[which(GA_model_prob >= 0.5)] <- 1

GA_a_ts <- AUROC_2(GA_model_prob, data_val$DEATH_EVENT)
GA_lr_cm <- table(data_val$DEATH_EVENT, GA_model_prey)

Perf_lr_time[14,] <- tt[3]
Perf_lr_auroc[14,] <- GA_a_ts
Perf_Table_ga[14,] <- perf_eval(GA_lr_cm)

# 100, 0.4, 0.3
start_time <- proc.time()
GA_model_lr <- ga(type="binary", fitness = fit_AUROC, nBits = ncol(x),
                  names = colnames(x),popSize = 50, pcrossover = 0.4, pmutation = 0.3,
                  maxiter = 100, elitism = 2, seed = 2021, monitor = plot)
end_time <- proc.time()
tt <- end_time - start_time


best_var_idx <- which(GA_model_lr@solution == 1)
GA_trn <- data_trn[,c(best_var_idx,13)]
GA_tst <- data_val[,c(best_var_idx,13)]

GA_model <- glm(DEATH_EVENT ~ ., family = binomial, GA_trn)
summary(GA_model)

GA_coeff <- as.matrix(GA_model$coefficients,13,1)
GA_coeff

GA_model_prob <- predict(GA_model, type="response", newdata = GA_tst)
GA_model_prey <- rep(0, nrow(data_val))
GA_model_prey[which(GA_model_prob >= 0.5)] <- 1

GA_a_ts <- AUROC_2(GA_model_prob, data_val$DEATH_EVENT)
GA_lr_cm <- table(data_val$DEATH_EVENT, GA_model_prey)

Perf_lr_time[15,] <- tt[3]
Perf_lr_auroc[15,] <- GA_a_ts
Perf_Table_ga[15,] <- perf_eval(GA_lr_cm)

# 100, 0.7, 0.05
start_time <- proc.time()
GA_model_lr <- ga(type="binary", fitness = fit_AUROC, nBits = ncol(x),
                  names = colnames(x),popSize = 50, pcrossover = 0.7, pmutation = 0.05,
                  maxiter = 100, elitism = 2, seed = 2021, monitor = plot)
end_time <- proc.time()
tt <- end_time - start_time


best_var_idx <- which(GA_model_lr@solution == 1)
GA_trn <- data_trn[,c(best_var_idx,13)]
GA_tst <- data_val[,c(best_var_idx,13)]

GA_model <- glm(DEATH_EVENT ~ ., family = binomial, GA_trn)
summary(GA_model)

GA_coeff <- as.matrix(GA_model$coefficients,13,1)
GA_coeff

GA_model_prob <- predict(GA_model, type="response", newdata = GA_tst)
GA_model_prey <- rep(0, nrow(data_val))
GA_model_prey[which(GA_model_prob >= 0.5)] <- 1

GA_a_ts <- AUROC_2(GA_model_prob, data_val$DEATH_EVENT)
GA_lr_cm <- table(data_val$DEATH_EVENT, GA_model_prey)

Perf_lr_time[16,] <- tt[3]
Perf_lr_auroc[16,] <- GA_a_ts
Perf_Table_ga[16,] <- perf_eval(GA_lr_cm)

# 100, 0.7, 0.1
start_time <- proc.time()
GA_model_lr <- ga(type="binary", fitness = fit_AUROC, nBits = ncol(x),
                  names = colnames(x),popSize = 50, pcrossover = 0.7, pmutation = 0.1,
                  maxiter = 100, elitism = 2, seed = 2021, monitor = plot)
end_time <- proc.time()
tt <- end_time - start_time


best_var_idx <- which(GA_model_lr@solution == 1)
GA_trn <- data_trn[,c(best_var_idx,13)]
GA_tst <- data_val[,c(best_var_idx,13)]

GA_model <- glm(DEATH_EVENT ~ ., family = binomial, GA_trn)
summary(GA_model)

GA_coeff <- as.matrix(GA_model$coefficients,13,1)
GA_coeff

GA_model_prob <- predict(GA_model, type="response", newdata = GA_tst)
GA_model_prey <- rep(0, nrow(data_val))
GA_model_prey[which(GA_model_prob >= 0.5)] <- 1

GA_a_ts <- AUROC_2(GA_model_prob, data_val$DEATH_EVENT)
GA_lr_cm <- table(data_val$DEATH_EVENT, GA_model_prey)

Perf_lr_time[17,] <- tt[3]
Perf_lr_auroc[17,] <- GA_a_ts
Perf_Table_ga[17,] <- perf_eval(GA_lr_cm)

# 100, 0.7, 0.3
start_time <- proc.time()
GA_model_lr <- ga(type="binary", fitness = fit_AUROC, nBits = ncol(x),
                  names = colnames(x),popSize = 50, pcrossover = 0.7, pmutation = 0.3,
                  maxiter = 100, elitism = 2, seed = 2021, monitor = plot)
end_time <- proc.time()
tt <- end_time - start_time


best_var_idx <- which(GA_model_lr@solution == 1)
GA_trn <- data_trn[,c(best_var_idx,13)]
GA_tst <- data_val[,c(best_var_idx,13)]

GA_model <- glm(DEATH_EVENT ~ ., family = binomial, GA_trn)
summary(GA_model)

GA_coeff <- as.matrix(GA_model$coefficients,13,1)
GA_coeff

GA_model_prob <- predict(GA_model, type="response", newdata = GA_tst)
GA_model_prey <- rep(0, nrow(data_val))
GA_model_prey[which(GA_model_prob >= 0.5)] <- 1

GA_a_ts <- AUROC_2(GA_model_prob, data_val$DEATH_EVENT)
GA_lr_cm <- table(data_val$DEATH_EVENT, GA_model_prey)

Perf_lr_time[18,] <- tt[3]
Perf_lr_auroc[18,] <- GA_a_ts
Perf_Table_ga[18,] <- perf_eval(GA_lr_cm)

# 200, 0.1, 0.05
start_time <- proc.time()
GA_model_lr <- ga(type="binary", fitness = fit_AUROC, nBits = ncol(x),
                  names = colnames(x),popSize = 50, pcrossover = 0.7, pmutation = 0.3,
                  maxiter = 200, elitism = 2, seed = 2021, monitor = plot)
end_time <- proc.time()
tt <- end_time - start_time


best_var_idx <- which(GA_model_lr@solution == 1)
GA_trn <- data_trn[,c(best_var_idx,13)]
GA_tst <- data_val[,c(best_var_idx,13)]

GA_model <- glm(DEATH_EVENT ~ ., family = binomial, GA_trn)
summary(GA_model)

GA_coeff <- as.matrix(GA_model$coefficients,13,1)
GA_coeff

GA_model_prob <- predict(GA_model, type="response", newdata = GA_tst)
GA_model_prey <- rep(0, nrow(data_val))
GA_model_prey[which(GA_model_prob >= 0.5)] <- 1

GA_a_ts <- AUROC_2(GA_model_prob, data_val$DEATH_EVENT)
GA_lr_cm <- table(data_val$DEATH_EVENT, GA_model_prey)

Perf_lr_time[19,] <- tt[3]
Perf_lr_auroc[19,] <- GA_a_ts
Perf_Table_ga[19,] <- perf_eval(GA_lr_cm)

# 200, 0.1, 0.1
start_time <- proc.time()
GA_model_lr <- ga(type="binary", fitness = fit_AUROC, nBits = ncol(x),
                  names = colnames(x),popSize = 50, pcrossover = 0.1, pmutation = 0.1,
                  maxiter = 200, elitism = 2, seed = 2021, monitor = plot)
end_time <- proc.time()
tt <- end_time - start_time


best_var_idx <- which(GA_model_lr@solution == 1)
GA_trn <- data_trn[,c(best_var_idx,13)]
GA_tst <- data_val[,c(best_var_idx,13)]

GA_model <- glm(DEATH_EVENT ~ ., family = binomial, GA_trn)
summary(GA_model)

GA_coeff <- as.matrix(GA_model$coefficients,13,1)
GA_coeff

GA_model_prob <- predict(GA_model, type="response", newdata = GA_tst)
GA_model_prey <- rep(0, nrow(data_val))
GA_model_prey[which(GA_model_prob >= 0.5)] <- 1

GA_a_ts <- AUROC_2(GA_model_prob, data_val$DEATH_EVENT)
GA_lr_cm <- table(data_val$DEATH_EVENT, GA_model_prey)

Perf_lr_time[20,] <- tt[3]
Perf_lr_auroc[20,] <- GA_a_ts
Perf_Table_ga[20,] <- perf_eval(GA_lr_cm)

# 200, 0.1, 0.3
start_time <- proc.time()
GA_model_lr <- ga(type="binary", fitness = fit_AUROC, nBits = ncol(x),
                  names = colnames(x),popSize = 50, pcrossover = 0.1, pmutation = 0.3,
                  maxiter = 200, elitism = 2, seed = 2021, monitor = plot)
end_time <- proc.time()
tt <- end_time - start_time


best_var_idx <- which(GA_model_lr@solution == 1)
GA_trn <- data_trn[,c(best_var_idx,13)]
GA_tst <- data_val[,c(best_var_idx,13)]

GA_model <- glm(DEATH_EVENT ~ ., family = binomial, GA_trn)
summary(GA_model)

GA_coeff <- as.matrix(GA_model$coefficients,13,1)
GA_coeff

GA_model_prob <- predict(GA_model, type="response", newdata = GA_tst)
GA_model_prey <- rep(0, nrow(data_val))
GA_model_prey[which(GA_model_prob >= 0.5)] <- 1

GA_a_ts <- AUROC_2(GA_model_prob, data_val$DEATH_EVENT)
GA_lr_cm <- table(data_val$DEATH_EVENT, GA_model_prey)

Perf_lr_time[21,] <- tt[3]
Perf_lr_auroc[21,] <- GA_a_ts
Perf_Table_ga[21,] <- perf_eval(GA_lr_cm)

# 200, 0.4, 0.05
start_time <- proc.time()
GA_model_lr <- ga(type="binary", fitness = fit_AUROC, nBits = ncol(x),
                  names = colnames(x),popSize = 50, pcrossover = 0.4, pmutation = 0.05,
                  maxiter = 200, elitism = 2, seed = 2021, monitor = plot)
end_time <- proc.time()
tt <- end_time - start_time


best_var_idx <- which(GA_model_lr@solution == 1)
GA_trn <- data_trn[,c(best_var_idx,13)]
GA_tst <- data_val[,c(best_var_idx,13)]

GA_model <- glm(DEATH_EVENT ~ ., family = binomial, GA_trn)
summary(GA_model)

GA_coeff <- as.matrix(GA_model$coefficients,13,1)
GA_coeff

GA_model_prob <- predict(GA_model, type="response", newdata = GA_tst)
GA_model_prey <- rep(0, nrow(data_val))
GA_model_prey[which(GA_model_prob >= 0.5)] <- 1

GA_a_ts <- AUROC_2(GA_model_prob, data_val$DEATH_EVENT)
GA_lr_cm <- table(data_val$DEATH_EVENT, GA_model_prey)

Perf_lr_time[22,] <- tt[3]
Perf_lr_auroc[22,] <- GA_a_ts
Perf_Table_ga[22,] <- perf_eval(GA_lr_cm)

# 200, 0.4, 0.1
start_time <- proc.time()
GA_model_lr <- ga(type="binary", fitness = fit_AUROC, nBits = ncol(x),
                  names = colnames(x),popSize = 50, pcrossover = 0.4, pmutation = 0.1,
                  maxiter = 200, elitism = 2, seed = 2021, monitor = plot)
end_time <- proc.time()
tt <- end_time - start_time


best_var_idx <- which(GA_model_lr@solution == 1)
GA_trn <- data_trn[,c(best_var_idx,13)]
GA_tst <- data_val[,c(best_var_idx,13)]

GA_model <- glm(DEATH_EVENT ~ ., family = binomial, GA_trn)
summary(GA_model)

GA_coeff <- as.matrix(GA_model$coefficients,13,1)
GA_coeff

GA_model_prob <- predict(GA_model, type="response", newdata = GA_tst)
GA_model_prey <- rep(0, nrow(data_val))
GA_model_prey[which(GA_model_prob >= 0.5)] <- 1

GA_a_ts <- AUROC_2(GA_model_prob, data_val$DEATH_EVENT)
GA_lr_cm <- table(data_val$DEATH_EVENT, GA_model_prey)

Perf_lr_time[23,] <- tt[3]
Perf_lr_auroc[23,] <- GA_a_ts
Perf_Table_ga[23,] <- perf_eval(GA_lr_cm)

# 200, 0.4, 0.3
start_time <- proc.time()
GA_model_lr <- ga(type="binary", fitness = fit_AUROC, nBits = ncol(x),
                  names = colnames(x),popSize = 50, pcrossover = 0.4, pmutation = 0.3,
                  maxiter = 200, elitism = 2, seed = 2021, monitor = plot)
end_time <- proc.time()
tt <- end_time - start_time


best_var_idx <- which(GA_model_lr@solution == 1)
GA_trn <- data_trn[,c(best_var_idx,13)]
GA_tst <- data_val[,c(best_var_idx,13)]

GA_model <- glm(DEATH_EVENT ~ ., family = binomial, GA_trn)
summary(GA_model)

GA_coeff <- as.matrix(GA_model$coefficients,13,1)
GA_coeff

GA_model_prob <- predict(GA_model, type="response", newdata = GA_tst)
GA_model_prey <- rep(0, nrow(data_val))
GA_model_prey[which(GA_model_prob >= 0.5)] <- 1

GA_a_ts <- AUROC_2(GA_model_prob, data_val$DEATH_EVENT)
GA_lr_cm <- table(data_val$DEATH_EVENT, GA_model_prey)

Perf_lr_time[24,] <- tt[3]
Perf_lr_auroc[24,] <- GA_a_ts
Perf_Table_ga[24,] <- perf_eval(GA_lr_cm)

# 200, 0.7, 0.05
start_time <- proc.time()
GA_model_lr <- ga(type="binary", fitness = fit_AUROC, nBits = ncol(x),
                  names = colnames(x),popSize = 50, pcrossover = 0.7, pmutation = 0.05,
                  maxiter = 200, elitism = 2, seed = 2021, monitor = plot)
end_time <- proc.time()
tt <- end_time - start_time


best_var_idx <- which(GA_model_lr@solution == 1)
GA_trn <- data_trn[,c(best_var_idx,13)]
GA_tst <- data_val[,c(best_var_idx,13)]

GA_model <- glm(DEATH_EVENT ~ ., family = binomial, GA_trn)
summary(GA_model)

GA_coeff <- as.matrix(GA_model$coefficients,13,1)
GA_coeff

GA_model_prob <- predict(GA_model, type="response", newdata = GA_tst)
GA_model_prey <- rep(0, nrow(data_val))
GA_model_prey[which(GA_model_prob >= 0.5)] <- 1

GA_a_ts <- AUROC_2(GA_model_prob, data_val$DEATH_EVENT)
GA_lr_cm <- table(data_val$DEATH_EVENT, GA_model_prey)

Perf_lr_time[25,] <- tt[3]
Perf_lr_auroc[25,] <- GA_a_ts
Perf_Table_ga[25,] <- perf_eval(GA_lr_cm)

# 200, 0.7, 0.1
start_time <- proc.time()
GA_model_lr <- ga(type="binary", fitness = fit_AUROC, nBits = ncol(x),
                  names = colnames(x),popSize = 50, pcrossover = 0.7, pmutation = 0.1,
                  maxiter = 200, elitism = 2, seed = 2021, monitor = plot)
end_time <- proc.time()
tt <- end_time - start_time


best_var_idx <- which(GA_model_lr@solution == 1)
GA_trn <- data_trn[,c(best_var_idx,13)]
GA_tst <- data_val[,c(best_var_idx,13)]

GA_model <- glm(DEATH_EVENT ~ ., family = binomial, GA_trn)
summary(GA_model)

GA_coeff <- as.matrix(GA_model$coefficients,13,1)
GA_coeff

GA_model_prob <- predict(GA_model, type="response", newdata = GA_tst)
GA_model_prey <- rep(0, nrow(data_val))
GA_model_prey[which(GA_model_prob >= 0.5)] <- 1

GA_a_ts <- AUROC_2(GA_model_prob, data_val$DEATH_EVENT)
GA_lr_cm <- table(data_val$DEATH_EVENT, GA_model_prey)

Perf_lr_time[26,] <- tt[3]
Perf_lr_auroc[26,] <- GA_a_ts
Perf_Table_ga[26,] <- perf_eval(GA_lr_cm)

# 200, 0.7, 0.3
start_time <- proc.time()
GA_model_lr <- ga(type="binary", fitness = fit_AUROC, nBits = ncol(x),
                  names = colnames(x),popSize = 50, pcrossover = 0.7, pmutation = 0.3,
                  maxiter = 200, elitism = 2, seed = 2021, monitor = plot)
end_time <- proc.time()
tt <- end_time - start_time


best_var_idx <- which(GA_model_lr@solution == 1)
GA_trn <- data_trn[,c(best_var_idx,13)]
GA_tst <- data_val[,c(best_var_idx,13)]

GA_model <- glm(DEATH_EVENT ~ ., family = binomial, GA_trn)
summary(GA_model)

GA_coeff <- as.matrix(GA_model$coefficients,13,1)
GA_coeff

GA_model_prob <- predict(GA_model, type="response", newdata = GA_tst)
GA_model_prey <- rep(0, nrow(data_val))
GA_model_prey[which(GA_model_prob >= 0.5)] <- 1

GA_a_ts <- AUROC_2(GA_model_prob, data_val$DEATH_EVENT)
GA_lr_cm <- table(data_val$DEATH_EVENT, GA_model_prey)

Perf_lr_time[27,] <- tt[3]
Perf_lr_auroc[27,] <- GA_a_ts
Perf_Table_ga[27,] <- perf_eval(GA_lr_cm)
Perf_lr_time
Perf_lr_auroc
Perf_Table_ga
