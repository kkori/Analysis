
library(pROC)
library(ROCR)
library(GA)
library(tree)
library(party)
library(nnet) 
library(ggplot2)

data <- read.csv("earthquake.csv")
str(data)
data <- data[,-1]
table(is.na(data))
dim(data)

# [Q1], numeric확인 / Bar chart
no_num <- c()
n_instance <- dim(data)[1]
n_var <- dim(data)[2]

for (i in c(1:n_var)){
  if (is.numeric(data[,i]) == FALSE){
    no_num <- c(no_num, i)
  }
}

no_num

# bar chart
par(mfrow = c(2,2))
for (i in no_num){
  no_num_data <- table(data[,i])
  barplot(no_num_data, main = colnames(data[i]), col = "#1b98e0", border = "#1b98e0")
  print(no_num_data)
}


# [Q2], 1-of-C coding
data_no_num <- data[,no_num]

for (i in c(1:ncol(data_no_num))){
  for(x in unique(data_no_num[,i])){
    data_no_num[paste(colnames(data_no_num[i]),x, sep = "_")] <- ifelse(data_no_num[,i] == x, 1, 0)
  }
}
data_no_num <- data_no_num[,-c(1:8)]

# 데이터 정규화 및 정리
input_num <- data[,-n_var]
num_data <- input_num[,-no_num]
data_input <- scale(num_data, center =TRUE, scale = TRUE)
data_target <- class.ind(data[,n_var])

data_normalized <- data.frame(data_input, data_no_num, Class = data_target)

# 데이터 분할
set.seed(2021)
full_idx <- c(1:n_instance)
trn_idx <- sample(1:n_instance, size = 200000)
val_idx <- sample(trn_idx, size = 50000)
data_tst <- data_normalized[-trn_idx,]
data_val <- data_normalized[val_idx,]
data_trn <- data_normalized[trn_idx[!(trn_idx%in%val_idx)],]

trn_input <- data_trn[,-(69:71)]
trn_target <- data_trn[,69:71]
val_input <- data_val[,-(69:71)]
val_target <- data_val[,69:71]

# Performance evaluation function for multi-class classification 
perf_eval_multi <- function(cm){ 
  # Simple Accuracy
  ACC = sum(diag(cm))/sum(cm)
  # Balanced Correction Rate
  BCR = 1
  for (i in 1:dim(cm)[1]){
    BCR = BCR*(cm[i,i]/sum(cm[i,])) 
  }
  BCR = BCR^(1/dim(cm)[1])
  return(c(ACC, BCR))
}

perf_eval_multi_2 <- function(tar, pred){
  tar_pred <- cbind(tar, pred)
  colnames(tar_pred) <- c("target", "predict")
  
  tp_df <- as.data.frame(tar_pred)
  ACC <- length(which(tp_df$target==tp_df$predict)) / nrow(tp_df)
  BCR = 1
  for (i in unique(tp_df$target)){
    BCR <- BCR*(length(which(tp_df$target==i & tp_df$predict==i)) / length(which(tp_df$target==i)))
  }
  
  BCR = BCR^(1/length(unique(tp_df$target)))
  
  return(c(ACC, BCR))
}

# [Q3]
# 하이퍼파라미터 선정
nh <- seq(from=7, to=15, by=1)
mi <- c(50,100,200)
val_perf <- matrix(0, nrow=length(nh)*length(mi), ncol=4)
val_perf_temp <- matrix(0, nrow=5, ncol=4)

num <- 1
start_time <- proc.time()

for (i in 1:length(nh)) {
  cat("Training ANN: the number of hidden nodes:", nh[i], "\n")

  for (j in 1:length(mi)) {
    hp_nnet <- nnet(trn_input,trn_target, size = nh[i], decay = 5e-4, maxit = mi[j], MaxNWts = 10000)
    val_perf[num, 1:2] <- c(nh[i], mi[j])
    val_perf[num, 3:4] <- perf_eval_multi_2(max.col(val_target), max.col(predict(hp_nnet, val_input)))    
    num <- num + 1
  }
}

proc.time() - start_time

# ACC 기준 정렬
ordered_perf_acc <- val_perf[order(val_perf[,3], decreasing = TRUE),]    
colnames(ordered_perf_acc) <- c("nh", "max_it", "ACC", "BCR")
ordered_perf_acc

# BCR 기준 정렬
ordered_perf_bcr <- val_perf[order(val_perf[,4], decreasing = TRUE),]    
colnames(ordered_perf_bcr) <- c("nh", "max_it", "ACC", "BCR")
ordered_perf_bcr


# [Q4]
rang <- seq(0.1, 0.9, 0.2)
perf_rang <- matrix(0,length(rang),3)
colnames(perf_rang) <- c("rang", "ACC", "BCR")


start_time <- proc.time()
num <- 1

for (i in 1:length(rang)){
  rang_nnet <- nnet(trn_input, trn_target, size = 8  , rang = rang[i], decay = 5e-4, maxit = 50, MaxNWts = 10000)
  rang_prey <- predict(rang_nnet, val_input)
  perf_rang[num, 1] <- rang[i]
  perf_rang[num, 2:3] <- perf_eval_multi_2(max.col(val_target), max.col(rang_prey))
  num <- num + 1
  
}

proc.time() - start_time
perf_rang
ordered_perf_rang <- perf_rang[order(perf_rang[,3], decreasing = TRUE),]    
colnames(ordered_perf_rang) <- c("rang", "ACC", "BCR")
ordered_perf_rang


# [Q5], best ANN 선정, 데이터셋 결합 후 학습
# nh = 8, max_iter = 50, rang = 0.9
data_trn_com <- rbind(data_trn, data_val)

trn_input_com <- data_trn_com[,-(69:71)]
trn_target_com <- data_trn_com[,69:71]
tst_input <- data_tst[,-(69:71)]
tst_target <- data_tst[,69:71]

perf_best <- matrix(0, 10, 2)
colnames(perf_best) <- c("ACC", "BCR")

num <- 1
start_time <- proc.time()

for (i in 1:10){
  best_nnet <- nnet(trn_input_com, trn_target_com, size = 8 , rang = 0.9, decay = 5e-4, maxit = 50, MaxNWts = 10000 )
  best_prey <- predict(best_nnet, tst_input)
  perf_best[num, 1:2] <- perf_eval_multi_2(max.col(tst_target), max.col(best_prey))
  num <- num + 1 
}
proc.time() - start_time
perf_best



# [Q6]

# fitness function
fit_BCR <- function(string){
  sel_var_idx <- which(string == 1)
  sel_x <- x[,sel_var_idx]
  GA_nnet <- nnet(sel_x, y, size = 8, rang = 0.9, decay = 5e-4, maxit = 50, MaxNWts = 10000)
  GA_nnet_prey <- predict(GA_nnet, ga_tst_input)
  tst_cm <- table(factor(max.col(ga_tst_target), levels = c("1","2","3")),factor(max.col(GA_nnet_prey), levels = c("1","2","3")))
  GA_bcr <- perf_eval_multi(tst_cm)[2]
  
  return(GA_bcr)
}

ndf <- nrow(data_trn_com)
x_idx <- sample(1:ndf, round(0.5*ndf))
ga_trn <- data_trn_com[x_idx,]
x <- as.matrix(ga_trn[,-(69:71)])
y <- ga_trn[,69:71]

ntdf <- nrow(data_tst)
tst_idx <- sample(1:ntdf, round(0.5*ntdf))
ga_tst <- data_tst[tst_idx,]
ga_tst_input <- ga_tst[,-(69:71)]
ga_tst_target <- ga_tst[,69:71]


# max_iter 50, cross 0.1, mutation 0.05, popsize = 25
start_time <- proc.time()
GA_model_nnet_1 <- ga(type="binary", fitness = fit_BCR, nBits = ncol(x),
                  names = colnames(x), popSize = 25, pcrossover = 0.1, pmutation = 0.05,
                  maxiter = 50, elitism = 2, seed = 2021, monitor=plot)

proc.time() - start_time
best_var_idx_1 <- which(GA_model_nnet_1@solution[1,] == 1)
best_var_idx_1

# max_iter 50, cross 0.1, mutation 0.05, popsize = 25
start_time <- proc.time()
GA_model_nnet_2 <- ga(type="binary", fitness = fit_BCR, nBits = ncol(x),
                  names = colnames(x), popSize = 25, pcrossover = 0.1, pmutation = 0.05,
                  maxiter = 50, elitism = 2, seed = 202, monitor=plot)

proc.time() - start_time
best_var_idx_2 <- which(GA_model_nnet_2@solution[1,] == 1)
best_var_idx_2

# max_iter 50, cross 0.1, mutation 0.05, popsize = 25
start_time <- proc.time()
GA_model_nnet_3 <- ga(type="binary", fitness = fit_BCR, nBits = ncol(x),
                    names = colnames(x), popSize = 25, pcrossover = 0.1, pmutation = 0.05,
                    maxiter = 50, elitism = 2, seed = 20, monitor=plot)

proc.time() - start_time
best_var_idx_3 <- which(GA_model_nnet_3@solution[1,] == 1)
best_var_idx_3


# [Q7], 선택된 변수 사용
select_idx <- c(1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 14, 15, 16, 17, 19, 20, 21, 23, 24, 25, 27, 29, 30, 31, 33, 34, 35, 36, 37, 38, 39, 40, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 56, 57, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68)
trn_input_com <- data_trn_com[,-(69:71)]
trn_input_com <- trn_input_com[,select_idx]
trn_target_com <- data_trn_com[,69:71]
tst_input <- data_tst[,-(69:71)]
tst_input <- tst_input[,select_idx]
tst_target <- data_tst[,69:71]

perf_select <- matrix(0, 1, 2)
colnames(perf_select) <- c("ACC", "BCR")

select_nnet <- nnet(trn_input_com, trn_target_com, size = 8 , rang = 0.9, decay = 5e-4, maxit = 50, MaxNWts = 10000 )
select_prey <- predict(select_nnet, tst_input)
perf_select[1, 1:2] <- perf_eval_multi_2(max.col(tst_target), max.col(select_prey))
perf_select 


# [Q8], Decision Tree
# eval함수 지정
perf_eval <- function(cm){
  
  # True positive rate: TPR (Recall)
  TPR <- cm[2,2]/sum(cm[2,])
  # Precision
  PRE <- cm[2,2]/sum(cm[,2])
  # True negative rate: TNR
  TNR <- cm[1,1]/sum(cm[1,])
  # Simple Accuracy
  ACC <- (cm[1,1]+cm[2,2])/sum(cm)
  # Balanced Correction Rate
  BCR <- sqrt(TPR*TNR)
  # F1-Measure
  F1 <- 2*TPR*PRE/(TPR+PRE)
  
  return(c(TPR, PRE, TNR, ACC, BCR, F1))
}


df <- read.csv("earthquake.csv")
# id 제거
df <- df[,-1]

df_no_num <- df[,no_num]

for (i in c(1:ncol(df_no_num))){
  for(x in unique(df_no_num[,i])){
    df_no_num[paste(colnames(df_no_num[i]),x, sep = "_")] <- ifelse(df_no_num[,i] == x, 1, 0)
  }
}

df_no_num <- df_no_num[,-c(1:8)]
df_input <- df[,-(39)]
df_input <- df_input[,-no_num]
df_input <- data.frame(df_input, df_no_num)
df_target <- as.factor(df[,39])


CT_trn <- data.frame(df_input[trn_idx[!(trn_idx%in%val_idx)],], dfYN = df_target[trn_idx[!(trn_idx%in%val_idx)]])
CT_val <- data.frame(df_input[val_idx,], dfYN = df_target[val_idx])
CT_tst <- data.frame(df_input[-trn_idx,], dfYN = df_target[-trn_idx])

# 하이퍼파라미터 후보값
min_criterion <- c(0.7, 0.9, 0.95)
min_split <- c(5, 10, 50, 100)
max_depth <- c(5, 10, 15, 20, 0)


CT_pre_search_result <- matrix(0, length(min_criterion)*length(min_split)*length(max_depth),11)
colnames(CT_pre_search_result) <- c("min_criterion", "min_split", "max_depth", 
                                     "TPR", "Precision", "TNR", "ACC", "BCR", "F1", "AUROC", "N_leaves")

iter_cnt = 1

for (i in 1:length(min_criterion)){
  for ( j in 1:length(min_split)){
    for ( k in 1:length(max_depth)){
      
      cat("CART Min criterion:", min_criterion[i], ", Min split:", min_split[j], ", Max depth:", max_depth[k], "\n")
      tmp_control = ctree_control(mincriterion = min_criterion[i], minsplit = min_split[j], maxdepth = max_depth[k])
      tmp_tree <- ctree(dfYN ~ ., data = CT_trn, controls = tmp_control)
      tmp_tree_val_prediction <- predict(tmp_tree, newdata = CT_val)
      tmp_tree_val_response <- treeresponse(tmp_tree, newdata = CT_val)
      tmp_tree_val_prob <- 1-unlist(tmp_tree_val_response, use.names=F)[seq(1,nrow(CT_val)*2,2)]
      tmp_tree_val_rocr <- multiclass.roc(CT_val$dfYN,tmp_tree_val_prob)
      # Confusion matrix for the validation dataset
      tmp_tree_val_cm <- table(CT_val$dfYN, tmp_tree_val_prediction)
      
      # parameters
      CT_pre_search_result[iter_cnt,1] = min_criterion[i]
      CT_pre_search_result[iter_cnt,2] = min_split[j]
      CT_pre_search_result[iter_cnt,3] = max_depth[k]
      # Performances from the confusion matrix
      CT_pre_search_result[iter_cnt,4:9] = perf_eval(tmp_tree_val_cm)
      # AUROC
      CT_pre_search_result[iter_cnt,10] = multiclass.roc(CT_val$dfYN,tmp_tree_val_prob)$auc
      # Number of leaf nodes
      CT_pre_search_result[iter_cnt,11] = length(nodes(tmp_tree, unique(where(tmp_tree))))
      iter_cnt = iter_cnt + 1
    }
  }
}

# best 하이퍼파라미터 찾기
CT_pre_search_result <- CT_pre_search_result[order(CT_pre_search_result[,10], decreasing = T),]
CT_pre_search_result
best_criterion <- CT_pre_search_result[1,1]
best_split <- CT_pre_search_result[1,2]
best_depth <- CT_pre_search_result[1,3]

tree_control = ctree_control(mincriterion = best_criterion, minsplit = best_split, maxdepth = best_depth)

# training과 validation 결합된 데이터셋 사용하여 학습
CT_trn <- rbind(CT_trn, CT_val)
CT_pre <- ctree(dfYN ~ ., data = CT_trn, controls = tree_control)
CT_pre_prediction <- predict(CT_pre, newdata = CT_tst)
CT_pre_response <- treeresponse(CT_pre, newdata = CT_tst)
CT_pre_cm <- table(CT_tst$dfYN, CT_pre_prediction)

# Plot, *node가 너무 많아 plot이 보기 좋지 않게 나와서 보고서에서는 생략하였습니다.
plot(CT_pre)
plot(CT_pre, type="simple")

perf_tree <- matrix(0, nrow = 1, ncol = 2)
rownames(perf_tree) <- c("Decision Tree")
colnames(perf_tree) <- c("Acc", "BCR")

perf_tree[1,] <- perf_eval_multi(CT_pre_cm)
perf_tree

# [Q9], Multinomial logistic regression
input_num <- data[,-n_var]
num_data <- input_num[,-no_num]
data_input <- scale(num_data, center =TRUE, scale = TRUE)
data_target <- as.factor(data[,39])
data_normalized <- data.frame(data_input, data_no_num, Class = data_target)

perf_mlr <- matrix(0, nrow=1, ncol=2)
colnames(perf_mlr) <- c("ACC", "BCR")
rownames(perf_mlr) <- c("Multi_Logit")

mlr_trn <- data_normalized[trn_idx,]
mlr_tst <- data_normalized[-trn_idx,]

ml_logit <- multinom(Class ~ ., data = mlr_trn)

summary(ml_logit)
a <- t(summary(ml_logit)$coefficients)

# Predict
ml_logit_prey <- predict(ml_logit, newdata = mlr_tst)
cfmatrix <- table(mlr_tst$Class, ml_logit_prey)
cfmatrix

perf_mlr[1,] <- perf_eval_multi(cfmatrix)
perf_mlr


# [추가 분석]
# hidden node 후보 값 변경
nh <- seq(from=16, to=25, by=3)
mi <- c(50,100,200)
val_perf <- matrix(0, nrow=length(nh)*length(mi), ncol=4)
val_perf_temp <- matrix(0, nrow=5, ncol=4)

num <- 1

for (i in 1:length(nh)) {
  cat("Training ANN: the number of hidden nodes:", nh[i], "\n")
  
  for (j in 1:length(mi)) {
    hp_nnet <- nnet(trn_input,trn_target, size = nh[i], decay = 5e-4, maxit = mi[j], MaxNWts = 10000)
    val_perf[num, 1:2] <- c(nh[i], mi[j])
    val_perf[num, 3:4] <- perf_eval_multi_2(max.col(val_target), max.col(predict(hp_nnet, val_input)))    
    num <- num + 1
  }
}

# BCR 기준 정렬
ordered_perf_bcr <- val_perf[order(val_perf[,4], decreasing = TRUE),]    
colnames(ordered_perf_bcr) <- c("nh", "max_it", "ACC", "BCR")
ordered_perf_bcr

# rang
rang <- seq(0.1, 0.9, 0.2)
perf_rang <- matrix(0,length(rang),3)
colnames(perf_rang) <- c("rang", "ACC", "BCR")


start_time <- proc.time()
num <- 1

for (i in 1:length(rang)){
  rang_nnet <- nnet(trn_input, trn_target, size = 22  , rang = rang[i], decay = 5e-4, maxit = 200, MaxNWts = 10000)
  rang_prey <- predict(rang_nnet, val_input)
  perf_rang[num, 1] <- rang[i]
  perf_rang[num, 2:3] <- perf_eval_multi_2(max.col(val_target), max.col(rang_prey))
  num <- num + 1
  
}

proc.time() - start_time
perf_rang
ordered_perf_rang <- perf_rang[order(perf_rang[,3], decreasing = TRUE),]    
colnames(ordered_perf_rang) <- c("rang", "ACC", "BCR")
ordered_perf_rang

# nh = 22, max_iter = 200, rang = 0.1, 10회 반복수행
data_trn_com <- rbind(data_trn, data_val)

trn_input_com <- data_trn_com[,-(69:71)]
trn_target_com <- data_trn_com[,69:71]
tst_input <- data_tst[,-(69:71)]
tst_target <- data_tst[,69:71]

perf_best <- matrix(0, 10, 2)
colnames(perf_best) <- c("ACC", "BCR")

num <- 1
start_time <- proc.time()

for (i in 1:10){
  best_nnet <- nnet(trn_input_com, trn_target_com, size = 22 , rang = 0.1, decay = 5e-4, maxit = 200, MaxNWts = 10000 )
  best_prey <- predict(best_nnet, tst_input)
  perf_best[num, 1:2] <- perf_eval_multi_2(max.col(tst_target), max.col(best_prey))
  num <- num + 1 
}
proc.time() - start_time
perf_best
ordered_perf_best <- perf_best[order(perf_best[,2], decreasing = TRUE),]    
colnames(ordered_perf_best) <- c("ACC", "BCR")
ordered_perf_best
