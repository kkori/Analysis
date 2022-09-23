
install.packages("tree")
install.packages("party")
install.packages("ROCR")
install.packages("dummies")
install.packages("MASS")
install.packages("dplyr")

library(ROCR)
library(tree)
library(party)
library(dummies)
library(MASS)
library(dplyr)

# heart_failure
# 데이터 불러오기
df <- read.csv("heart_failure_clinical_records_dataset.csv")
str(df)
# 결측치 확인
table(is.na(df))

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

# Performance table
Perf_Table <- matrix(0, nrow = 4, ncol = 6)
rownames(Perf_Table) <- c("Without-Pruning","Post-Pruning", "Pre-Pruning","Logistic Regression")
colnames(Perf_Table) <- c("TPR", "Precision", "TNR", "Accuracy", "BCR", "F1-Measure")
Perf_Table

# [Q2]
# factor 변환
df$anaemia <- as.factor(df$anaemia)
df$diabetes <- as.factor(df$diabetes)
df$high_blood_pressure <- as.factor(df$high_blood_pressure)
df$sex <- as.factor(df$sex)
df$smoking <- as.factor(df$smoking)

input_idx <- c(1:12)
target_idx <- 13

df_input <- df[,input_idx]
df_target <- as.factor(df[,target_idx])

# train/test => 80:20
set.seed(2021)
trn_idx <- sample(1:nrow(df), round(0.8*nrow(df)))

CT_trn <- data.frame(df_input[trn_idx,], dfYN = df_target[trn_idx])
CT_tst <- data.frame(df_input[-trn_idx,], dfYN = df_target[-trn_idx])

# Training
CT_ori <- tree(dfYN ~., CT_trn)
summary(CT_ori)

# Plot
plot(CT_ori)
text(CT_ori, pretty=1)

# prediction
CT_ori_prey <- predict(CT_ori, CT_tst, type="class")
CT_ori_cm <- table(CT_tst$dfYN, CT_ori_prey)
CT_ori_cm
Perf_Table[1,] <- perf_eval(CT_ori_cm)
Perf_Table

# AUROC
CT_ori_prey <- predict(CT_ori, CT_tst, type="vector")[,2]
CT_ori_rocr <- prediction(CT_ori_prey, CT_tst$dfYN)
CT_ori_perf <- performance(CT_ori_rocr,"tpr",'fpr')
plot(CT_ori_perf, col=5,lwd=3)
performance(CT_ori_rocr, "auc")@y.values[[1]]


# [Q3], Post-Pruning
# Find the best tree
set.seed(2021)
CT_post_cv <- cv.tree(CT_ori, FUN = prune.misclass)

# Plot
plot(CT_post_cv$size, CT_post_cv$dev, type="b")
CT_post_cv

# final model, best = 7
CT_post_pruned <- prune.misclass(CT_ori, best=7)
plot(CT_post_pruned)
text(CT_post_pruned,pretty=1)

# Prediction
CT_post_prey <- predict(CT_post_pruned, CT_tst, type = "class")
CT_post_cm <- table(CT_tst$dfYN, CT_post_prey)
CT_post_cm

Perf_Table[2,] <- perf_eval(CT_post_cm)
Perf_Table

# AUROC
CT_post_prey <- predict(CT_post_pruned, CT_tst, type="vector")[,2]
CT_post_rocr <- prediction(CT_post_prey, CT_tst$dfYN)
CT_post_perf <- performance(CT_post_rocr,"tpr",'fpr')
plot(CT_post_perf, col=5,lwd=3)
performance(CT_post_rocr, "auc")@y.values[[1]]


# [Q4], 최적의 하이퍼파라미터 찾기
# train:val:test = 6:2:2
set.seed(2021)
trn_idx <- sample(1:nrow(df), round(0.6*nrow(df)))
val_idx <- sample(1:nrow(df[-trn_idx,]), round(0.2*nrow(df)))
ful <- c(1:nrow(df))
ful <- ful[-trn_idx]
tst_idx <- ful[-val_idx]

CT_trn <- data.frame(df_input[trn_idx,], dfYN = df_target[trn_idx])
CT_val <- data.frame(df_input[val_idx,], dfYN = df_target[val_idx])
CT_tst <- data.frame(df_input[tst_idx,], dfYN = df_target[tst_idx])

min_criterion <- c(0.6, 0.8, 0.9, 0.95, 0.99)
min_split <- c(2,5,10,50,100,150)
max_depth <- c(3,5,10,15,20,0)

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
      tmp_tree_val_rocr <- prediction(tmp_tree_val_prob, CT_val$dfYN)
      # Confusion matrix for the validation dataset
      tmp_tree_val_cm <- table(CT_val$dfYN, tmp_tree_val_prediction)
      
      # parameters
      CT_pre_search_result[iter_cnt,1] = min_criterion[i]
      CT_pre_search_result[iter_cnt,2] = min_split[j]
      CT_pre_search_result[iter_cnt,3] = max_depth[k]
      # Performances from the confusion matrix
      CT_pre_search_result[iter_cnt,4:9] = perf_eval(tmp_tree_val_cm)
      # AUROC
      CT_pre_search_result[iter_cnt,10] = unlist(performance(tmp_tree_val_rocr, "auc")@y.values)
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

# [Q5], 최적의 모델 구축
tree_control = ctree_control(mincriterion = best_criterion, minsplit = best_split, maxdepth = best_depth)
CT_best <- ctree(dfYN ~., data = CT_trn, controls = tree_control)
summary(CT_best)

# Plot
plot(CT_best)
plot(CT_best, type="simple")

# [Q6], training과 validation 결합
CT_trn <- rbind(CT_trn, CT_val)
CT_pre <- ctree(dfYN ~ ., data = CT_trn, controls = tree_control)
CT_pre_prediction <- predict(CT_pre, newdata = CT_tst)
CT_pre_response <- treeresponse(CT_pre, newdata = CT_tst)

CT_pre_cm <- table(CT_tst$dfYN, CT_pre_prediction)
CT_pre_cm
Perf_Table[3,] <- perf_eval(CT_pre_cm)
Perf_Table

# roc
CT_pre_prob <- 1-unlist(CT_pre_response,use.names = F)[seq(1,nrow(CT_tst)*2,2)]
CT_pre_rocr <- prediction(CT_pre_prob,CT_tst$dfYN)
CT_pre_perf <- performance(CT_pre_rocr,"tpr","fpr")
plot(CT_pre_perf, col =5 ,lwd=3)
unlist(performance(CT_pre_rocr, "auc")@y.values)

# plot tree
plot(CT_pre)
plot(CT_pre,type="simple")

# -----------2번째 데이터 사용
df <- read.csv("diabetes.csv")
str(df)
# 결측치 확인
table(is.na(df))

Perf_Table <- matrix(0, nrow = 4, ncol = 6)
rownames(Perf_Table) <- c("Without-Pruning","Post-Pruning", "Pre-Pruning","Logistic Regression")
colnames(Perf_Table) <- c("TPR", "Precision", "TNR", "Accuracy", "BCR", "F1-Measure")
Perf_Table

# [Q2]
input_idx <- c(1:8)
target_idx <- 9

df_input <- df[,input_idx]
df_target <- as.factor(df[,target_idx])

# train/test => 80:20
set.seed(2021)
trn_idx <- sample(1:nrow(df), round(0.8*nrow(df)))

CT_trn <- data.frame(df_input[trn_idx,], dfYN = df_target[trn_idx])
CT_tst <- data.frame(df_input[-trn_idx,], dfYN = df_target[-trn_idx])

# Training
CT_ori <- tree(dfYN ~., CT_trn)
summary(CT_ori)

# Plot
plot(CT_ori)
text(CT_ori, pretty=1)

# prediction
CT_ori_prey <- predict(CT_ori, CT_tst, type="class")
CT_ori_cm <- table(CT_tst$dfYN, CT_ori_prey)
CT_ori_cm
Perf_Table[1,] <- perf_eval(CT_ori_cm)
Perf_Table

# AUROC
CT_ori_prey <- predict(CT_ori, CT_tst, type="vector")[,2]
CT_ori_rocr <- prediction(CT_ori_prey, CT_tst$dfYN)
CT_ori_perf <- performance(CT_ori_rocr,"tpr",'fpr')
plot(CT_ori_perf, col=5,lwd=3)
performance(CT_ori_rocr, "auc")@y.values[[1]]


# [Q3], post-pruning
# Find the best tree
set.seed(2021)
CT_post_cv <- cv.tree(CT_ori, FUN = prune.misclass)

# Plot
plot(CT_post_cv$size, CT_post_cv$dev, type="b")
CT_post_cv

# final model, best = 6
CT_post_pruned <- prune.misclass(CT_ori, best=6)
plot(CT_post_pruned)
text(CT_post_pruned,pretty=1)

# Prediction
CT_post_prey <- predict(CT_post_pruned, CT_tst, type = "class")
CT_post_cm <- table(CT_tst$dfYN, CT_post_prey)
CT_post_cm

Perf_Table[2,] <- perf_eval(CT_post_cm)
Perf_Table

# AUROC
CT_post_prey <- predict(CT_post_pruned, CT_tst, type="vector")[,2]
CT_post_rocr <- prediction(CT_post_prey, CT_tst$dfYN)
CT_post_perf <- performance(CT_post_rocr,"tpr",'fpr')
plot(CT_post_perf, col=5,lwd=3)
performance(CT_post_rocr, "auc")@y.values[[1]]

# [Q4], Pre-pruning, 최적의 하이퍼파라미터 찾기
# train:val:test = 6:2:2
set.seed(2021)
trn_idx <- sample(1:nrow(df), round(0.6*nrow(df)))
val_idx <- sample(1:nrow(df[-trn_idx,]), round(0.2*nrow(df)))
ful <- c(1:nrow(df))
ful <- ful[-trn_idx]
tst_idx <- ful[-val_idx]

CT_trn <- data.frame(df_input[trn_idx,], dfYN = df_target[trn_idx])
CT_val <- data.frame(df_input[val_idx,], dfYN = df_target[val_idx])
CT_tst <- data.frame(df_input[tst_idx,], dfYN = df_target[tst_idx])

min_criterion <- c(0.6, 0.8, 0.9, 0.95,0.99)
min_split <- c(2,5,10,50,100,150)
max_depth <- c(3,5,10,15,20,0)

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
      tmp_tree_val_rocr <- prediction(tmp_tree_val_prob, CT_val$dfYN)
      # Confusion matrix for the validation dataset
      tmp_tree_val_cm <- table(CT_val$dfYN, tmp_tree_val_prediction)
      
      # parameters
      CT_pre_search_result[iter_cnt,1] = min_criterion[i]
      CT_pre_search_result[iter_cnt,2] = min_split[j]
      CT_pre_search_result[iter_cnt,3] = max_depth[k]
      # Performances from the confusion matrix
      CT_pre_search_result[iter_cnt,4:9] = perf_eval(tmp_tree_val_cm)
      # AUROC
      CT_pre_search_result[iter_cnt,10] = unlist(performance(tmp_tree_val_rocr, "auc")@y.values)
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

# [Q5], 최적의 모델 구축
tree_control = ctree_control(mincriterion = best_criterion, minsplit = best_split, maxdepth = best_depth)
CT_best <- ctree(dfYN ~., data = CT_trn, controls = tree_control)
summary(CT_best)

# Plot
plot(CT_best)
plot(CT_best, type="simple")

# [Q6], training과 validation 결합하여 학습, 위에서 구한 하이퍼파라미터 값 이용
CT_trn <- rbind(CT_trn, CT_val)
CT_pre <- ctree(dfYN ~ ., data = CT_trn, controls = tree_control)
CT_pre_prediction <- predict(CT_pre, newdata = CT_tst)
CT_pre_response <- treeresponse(CT_pre, newdata = CT_tst)

CT_pre_cm <- table(CT_tst$dfYN, CT_pre_prediction)
CT_pre_cm
Perf_Table[3,] <- perf_eval(CT_pre_cm)
Perf_Table

# roc
CT_pre_prob <- 1-unlist(CT_pre_response,use.names = F)[seq(1,nrow(CT_tst)*2,2)]
CT_pre_rocr <- prediction(CT_pre_prob,CT_tst$dfYN)
CT_pre_perf <- performance(CT_pre_rocr,"tpr","fpr")
plot(CT_pre_perf, col =5 ,lwd=3)
unlist(performance(CT_pre_rocr, "auc")@y.values)

# plot tree
plot(CT_pre)
plot(CT_pre,type="simple")


#[Q7] Logistic Regression
# AUROC
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

# Dataset1 불러오기
df <- read.csv("heart_failure_clinical_records_dataset.csv")
str(df)
# 결측치 확인
table(is.na(df))

# factor로 변환
df$anaemia <- as.factor(df$anaemia)
df$diabetes <- as.factor(df$diabetes)
df$high_blood_pressure <- as.factor(df$high_blood_pressure)
df$sex <- as.factor(df$sex)
df$smoking <- as.factor(df$smoking)
df$DEATH_EVENT <- as.factor(df$DEATH_EVENT)

# 수치형 변수에 대해서 nomarlization 실행
input_num_idx <- c(1,3,5,7,8,9,12)
input_cat_idx <- c(2,4,6,10,11)

df_input <- df[,input_num_idx]
df_cat_input <- df[,input_cat_idx]

df_input <- scale(df_input, center=TRUE, scale = TRUE)

df_input <- data.frame(df_input, df_cat_input)
df_target <- df[,13]

heart_df <- data.frame(df_input, df_target)
names(heart_df)[13] <- c("DEATH_EVENT")

# 시드 설정 후, 80:20 비율로 train data와 test data 분할
set.seed(365)
trn_idx <- sample(1:nrow(heart_df), round(0.8*nrow(heart_df)))
heart_trn <- heart_df[trn_idx,]
heart_tst <- heart_df[-trn_idx,]

full_lr <- glm(DEATH_EVENT ~ ., family = binomial, heart_trn)
summary(full_lr)

tst_response <- predict(full_lr, type="response", newdata=heart_tst)
tst_target <- heart_tst$DEATH_EVENT
tst_predicted <- rep(0,length(tst_target))
tst_predicted[which(tst_response >= 0.5)] <- 1
cm_tst_full <- table(tst_target, tst_predicted)
cm_tst_full
Perf_Table[4,] <- perf_eval(cm_tst_full)
Perf_Table

AUROC(tst_response, tst_target)

# Dataset2, Logistic Regression
df <- read.csv("diabetes.csv")
str(df)
# 결측치 확인
table(is.na(df))

input_idx <- c(1:8)
target_idx <- 9

df_input <- df[,input_idx]
df_input <- scale(df_input, center=TRUE, scale = TRUE)

df_target <- as.factor(df[,target_idx])

df <- data.frame(df_input, df_target)
names(df)[9] <- c("Outcome")

# 시드 설정 후, 80:20 비율로 train data와 test data 분할
set.seed(365)
trn_idx <- sample(1:nrow(df), round(0.8*nrow(df)))
df_trn <- df[trn_idx,]
df_tst <- df[-trn_idx,]

full_lr <- glm(Outcome ~ ., family = binomial, df_trn)
summary(full_lr)

tst_response <- predict(full_lr, type="response", newdata=df_tst)
tst_target <- df_tst$Outcome
tst_predicted <- rep(0,length(tst_target))
tst_predicted[which(tst_response >= 0.5)] <- 1
cm_tst_full <- table(tst_target, tst_predicted)
cm_tst_full
Perf_Table[4,] <- perf_eval(cm_tst_full)
Perf_Table

AUROC(tst_response, tst_target)
