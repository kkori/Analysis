
.libPaths("C:/Users/cjoi0/OneDrive/바탕 화면/library")
.libPaths()

install.packages('caret', dependencies = TRUE)
install.packages("doParallel")
install.packages("randomForest")
install.packages("ada")
install.packages("gbm")
install.packages("ipred")
install.packages("tree")
install.packages("party")
install.packages("nnet")
install.packages("ggplot2")


library(ipred)
library(caret)
library(doParallel)
library(ada)
library(gbm)
library(randomForest)
library(tree)
library(party)
library(nnet) 
library(ggplot2)

data <- read.csv("Earthquake_Damage.csv")
str(data)
data <- data[,-1]
table(is.na(data))
dim(data)

# 사전 작업, 1-hot encoding
no_num <- c()
n_instance <- dim(data)[1]
n_var <- dim(data)[2]

for (i in c(1:n_var)){
  if (is.numeric(data[,i]) == FALSE){
    no_num <- c(no_num, i)
  }
}
no_num

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

# 사전작업, 데이터 분할
set.seed(365)
full_idx <- c(1:n_instance)
trn_idx <- sample(1:n_instance, size = 200000)
val_idx <- sample(trn_idx, size = 50000)
data_tst <- data_normalized[-trn_idx,]
data_val <- data_normalized[val_idx,]
data_trn <- data_normalized[trn_idx[!(trn_idx%in%val_idx)],]


# 평가함수
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

# [Q1]
perf_single <- matrix(0, nrow=3,ncol=2)
rownames(perf_single) <- c("Multinomial","CART","ANN")
colnames(perf_single) <- c("ACC", "BCR")


# Multinomial Logistic Regression
input_num <- data[,-n_var]
num_data <- input_num[,-no_num]
data_input <- scale(num_data, center =TRUE, scale = TRUE)
data_target <- as.factor(data[,39])
data_normalized <- data.frame(data_input, data_no_num, Class = data_target)

mlr_trn <- data_normalized[trn_idx,]
mlr_tst <- data_normalized[-trn_idx,]

ml_logit <- multinom(Class ~ ., data = mlr_trn)

summary(ml_logit)
t(summary(ml_logit)$coefficients)

# Predict
ml_logit_prey <- predict(ml_logit, newdata = mlr_tst)
mlr_cm <- table(mlr_tst$Class, ml_logit_prey)
mlr_cm

perf_single[1,] <- perf_eval_multi(cfmatrix)
perf_single

# CART
df <- read.csv("Earthquake_Damage.csv")
# id 제거
df <- df[,-1]
# damage_grade가 1인 개수, 25124개만 추출
df_target_1 <- df[df$damage_grade == 1,]
df_target_2 <- df[df$damage_grade == 2,]
df_target_2 <- df_target_2[sample(nrow(df_target_2), 25124),]
df_target_3 <- df[df$damage_grade == 3,]
df_target_3 <- df_target_3[sample(nrow(df_target_3), 25124),]
df <- rbind(df_target_1,df_target_2,df_target_3)

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

# train, validatoin, test : 3:1:1.2 정도의 비율로 다시 설정
n_instance <- dim(df)[1]
n_var <- dim(df)[2]
set.seed(365)
full_idx <- c(1:n_instance)
trn_idx <- sample(1:n_instance, size = 57977)
val_idx <- sample(trn_idx, size = 14494)
data_tst <- data_normalized[-trn_idx,]
data_val <- data_normalized[val_idx,]
data_trn <- data_normalized[trn_idx[!(trn_idx%in%val_idx)],]

CT_trn <- data.frame(df_input[trn_idx[!(trn_idx%in%val_idx)],], dfYN = df_target[trn_idx[!(trn_idx%in%val_idx)]])
CT_val <- data.frame(df_input[val_idx,], dfYN = df_target[val_idx])
CT_tst <- data.frame(df_input[-trn_idx,], dfYN = df_target[-trn_idx])

min_split <- c(5,10,15,25,30)
max_depth <- c(5,15,25,35,0)
cp <- seq(0.01,0.03,0.01)

perf_CT <- matrix(0,length(min_split)*length(max_depth)*length(cp), 5)
colnames(perf_CT) <- c("min_split", "max_depth", "cp", "ACC", "BCR")

start_time <- proc.time()
num <- 1
for (i in 1:length(min_split)) {
  cat("Training CART: the number of min_split:", min_split[i], "\n")
  for (j in 1:length(max_depth)) {
    cat("Training CART: the number of max_depth: ", max_depth[j], "\n")
    for (k in 1:length(cp)){
      tree_control <- rpart.control(min_split[i],max_depth[j],cp[k])
      CT_model <- rpart(dfYN ~ ., CT_trn, control = tree_control)
      prey <- predict(CT_model, CT_val, type = "class")
      perf_CT[num, 1:3] <- c(min_split[i],max_depth[j],cp[k])
      CT_cm <- table(CT_val$dfYN, prey)
      perf_CT[num, 4:5] <- perf_eval_multi_2(CT_val$dfYN,prey)    
      num <- num + 1
    }
  }
}
proc.time() - start_time

CT_hp_result <- perf_CT[order(perf_CT[,5], decreasing = TRUE),]    
colnames(CT_hp_result) <- c("min_split", "max_depth","cp", "ACC", "BCR")
best_min_split <- CT_hp_result[1,1]
best_max_depth <- CT_hp_result[1,2]
best_cp <- CT_hp_result[1,3]

# predict
CT_trn <- rbind(CT_trn, CT_val)

tree_control <- rpart.control(minsplit=best_min_split, maxdepth=best_max_depth, cp=best_cp)
CT_best <- rpart(dfYN ~ ., CT_trn, control=tree_control)
CT_best_prey <- predict(CT_best, newdata= CT_tst, type = "class")
CT_best_cm <- table(CT_tst$dfYN, CT_best_prey)
CT_best_cm
perf_single[2,] <- perf_eval_multi_2(CT_tst$dfYN, CT_best_prey)
perf_single


# ANN
set.seed(365)
full_idx <- c(1:n_instance)
trn_idx <- sample(1:n_instance, size = 200000)
val_idx <- sample(trn_idx, size = 50000)
data_target <- class.ind(data[,n_var])
data_normalized <- data.frame(data_input, data_no_num, Class = data_target)

data_tst <- data_normalized[-trn_idx,]
data_val <- data_normalized[val_idx,]
data_trn <- data_normalized[trn_idx[!(trn_idx%in%val_idx)],]

trn_input <- data_trn[,-(69:71)]
trn_target <- data_trn[,69:71]
val_input <- data_val[,-(69:71)]
val_target <- data_val[,69:71]

# 하이퍼파라미터 선정
nh <- seq(from=7, to=20, by=2)
mi <- c(50,100,200)
rang <- c(0.1, 0.5, 0.9)
val_perf <- matrix(0, nrow=length(nh)*length(mi)*length(rang), ncol=5)
colnames(val_perf) <- c("nh","max_iter","rang","ACC","BCR")

num <- 1
start_time <- proc.time()

for (i in 1:length(nh)) {
  cat("Training ANN: the number of hidden nodes:", nh[i], "\n")
  
  for (j in 1:length(mi)) {
    for (k in 1:length(rang)) {
      hp_nnet <- nnet(trn_input, trn_target, size = nh[i], decay = 5e-4, maxit = mi[j],rang=rang[k], MaxNWts = 10000)
      val_perf[num, 1:3] <- c(nh[i], mi[j],rang[k])
      val_perf[num, 4:5] <- perf_eval_multi_2(max.col(val_target), max.col(predict(hp_nnet, val_input)))    
      num <- num + 1
    }
  }
}

proc.time() - start_time


# ACC 기준 정렬
ordered_perf_acc <- val_perf[order(val_perf[,4], decreasing = TRUE),]    
colnames(ordered_perf_acc) <- c("nh", "max_it","rang", "ACC", "BCR")
ordered_perf_acc

# BCR 기준 정렬
ordered_perf_bcr <- val_perf[order(val_perf[,5], decreasing = TRUE),]    
colnames(ordered_perf_bcr) <- c("nh", "max_it", "rang","ACC", "BCR")
ordered_perf_bcr

# predict
data_trn_com <- rbind(data_trn, data_val)

trn_input_com <- data_trn_com[,-(69:71)]
trn_target_com <- data_trn_com[,69:71]
tst_input <- data_tst[,-(69:71)]
tst_target <- data_tst[,69:71]

best_nnet <- nnet(trn_input_com, trn_target_com, size = 19  , rang = 0.1, decay = 5e-4, maxit = 200, MaxNWts = 10000 )
best_prey <- predict(best_nnet, tst_input)
perf_single[3,] <- perf_eval_multi_2(max.col(tst_target), max.col(best_prey))

perf_single

# Assign core
cl <- makeCluster(4)
registerDoParallel(cl)


# [Q2], CART Bagging
CT_trn <- data.frame(df_input[trn_idx[!(trn_idx%in%val_idx)],], dfYN = df_target[trn_idx[!(trn_idx%in%val_idx)]])
CT_val <- data.frame(df_input[val_idx,], dfYN = df_target[val_idx])
CT_tst <- data.frame(df_input[-trn_idx,], dfYN = df_target[-trn_idx])

perf_CT_bag <- matrix(0, nrow=10, ncol=3)
colnames(perf_CT_bag) <- c("Bootstrap","ACC","BCR")
boot <- seq(30,300,30)

start_time <- proc.time()
num <- 1
for (i in 1:10) {
  cat("Bootstrap : ", boot[i], "\n")
  tree_control <- rpart.control(minsplit = 5, maxdepth = 5, cp = 0.01)
  CT_bagging <- bagging(dfYN ~ ., CT_trn, control = tree_control, nbagg = boot[i], coob=TRUE)
  prey <- predict(CT_bagging, CT_val, type = "class")
  perf_CT_bag[num,1] <- boot[i]
  bag_cm <- table(CT_val$dfYN, prey)
  perf_CT_bag[num, 2:3] <- perf_eval_multi_2(CT_val$dfYN, prey)
  num <- num + 1
}
proc.time() - start_time

perf_CT_bag

# predict
CT_trn <- rbind(CT_trn, CT_val)

CT_bagging_best <- bagging(dfYN ~ ., CT_trn, control = tree_control, nbagg = 210, coob=TRUE)
prey <- predict(CT_bagging, CT_tst, type = "class")
CT_bag_best_cm <- table(CT_tst$dfYN, prey)
CT_bag_best_cm
perf_eval_multi_2(CT_tst$dfYN, prey)



# [Q3], Random Forest
CT_trn <- data.frame(df_input[trn_idx[!(trn_idx%in%val_idx)],], dfYN = df_target[trn_idx[!(trn_idx%in%val_idx)]])
CT_val <- data.frame(df_input[val_idx,], dfYN = df_target[val_idx])
CT_tst <- data.frame(df_input[-trn_idx,], dfYN = df_target[-trn_idx])

RF.trn <- CT_trn
RF.tst <- CT_tst
RF.val <- CT_val

ntr <- seq(30,300,30)

perf_rf <-  matrix(0,nrow=10,ncol=3)
colnames(perf_rf) <- c("Tree의 수", "ACC", "BCR")


start_time <- proc.time()
num <- 1
for (i in 1:length(ntr)){
  RF.model <- randomForest(dfYN ~ ., data = RF.trn, ntree = ntr[i], importance = TRUE, do.trace = TRUE)
  
  
  # Variable importance
  Var.imp <- importance(RF.model)
  barplot(Var.imp[order(Var.imp[,4], decreasing = TRUE),4])
  
  # Prediction
  RF.prey <- predict(RF.model, newdata = RF.val, type = "class")
  RF.cfm <- table(RF.prey, RF.val$dfYN)
  
  perf_rf[num,1] <- ntr[i]
  perf_rf[num,2:3] <- perf_eval_multi(RF.cfm)    
  
  num <- num + 1
} 
proc.time() - start_time

perf_rf

# predict
RF.trn <- rbind(RF.trn, RF.val)

RF.model <- randomForest(dfYN ~ ., data = RF.trn, ntree = 270 , importance = TRUE, do.trace = TRUE)

print(RF.model)
plot(RF.model)  

Var.imp <- importance(RF.model)
barplot(Var.imp[order(Var.imp[,4], decreasing = TRUE),4])

RF.prey <- predict(RF.model, newdata = RF.tst, type = "class")
RF.cm <- table(RF.tst$dfYN,RF.prey)
RF.cm
perf_eval_multi_2(RF.tst$dfYN,RF.prey)

# plot
boots <- seq(30,300,30)
CT_bag_acc <- c(0.6061129,0.6021112,0.6021112,0.6021112,0.6021112,0.6021112,0.6081827,0.6021112,0.6021112,0.6021112)
CT_bag_bcr <- c(0.5474177,0.5361400,0.5361400,0.5361400,0.5361400,0.5361400,0.5573837,0.5361400,0.5361400,0.5361400)
RF_acc <- c(0.70728,0.70958,0.71202,0.71362,0.71276,0.71178,0.71354,0.71318,0.71380,0.71390)
RF_bcr <- c(0.7063273,0.7083739,0.7125012,0.7158831,0.7141888,0.7125026,0.7155706,0.7144077,0.7169550,0.7152452)

plot(CT_bag_acc, axes=F, main = "ACC",xlab="bootstrap", ylab="ACC", type='o', col='red', ylim = c(0.5,0.8))
axis(1, at = 1:10, lab=boots)
axis(2, ylim=c(0.5,0.8))
lines(RF_acc, type='o', col='blue')

plot(CT_bag_bcr, axes=F, main = "BCR", xlab="bootstrap", ylab="BCR", type='o', col='red', ylim = c(0.4,0.8))
axis(1, at = 1:10, lab=boots)
axis(2, ylim=c(0.4,0.8))
lines(RF_bcr, type='o', col='blue')

# [Q4], ANN 반복수행
perf_ann_rep <- matrix(0,30,2)
colnames(perf_ann_rep) <- c("ACC", "BCR")

start_time <- proc.time()
num <- 1
for (i in 1:30){
  best_nnet <- nnet(trn_input_com, trn_target_com, size = 19 , rang =0.1 , decay = 5e-4, maxit =200 , MaxNWts = 10000 )
  best_prey <- predict(best_nnet, tst_input)
  perf_ann_rep[num,] <- perf_eval_multi_2(max.col(tst_target), max.col(best_prey))
  num <- num + 1
}
proc.time() - start_time
perf_ann_rep

me_var <- matrix(0,2,2)
rownames(me_var) <- c("ACC","BCR")
colnames(me_var) <- c("Mean","Variance")
me_var[1,] <- c(mean(perf_ann_rep[,1]), var(perf_ann_rep[,1]))
me_var[2,] <- c(mean(perf_ann_rep[,2]), var(perf_ann_rep[,2]))
me_var

# [Q5], ANN Bagging
# 데이터 축소
set.seed(365)
data <- read.csv("Earthquake_Damage.csv")
data <- data[,-1]

no_num <- c()
n_instance <- dim(data)[1]
n_var <- dim(data)[2]

for (i in c(1:n_var)){
  if (is.numeric(data[,i]) == FALSE){
    no_num <- c(no_num, i)
  }
}
no_num

data_no_num <- data[,no_num]

for (i in c(1:ncol(data_no_num))){
  for(x in unique(data_no_num[,i])){
    data_no_num[paste(colnames(data_no_num[i]),x, sep = "_")] <- ifelse(data_no_num[,i] == x, 1, 0)
  }
}
data_no_num <- data_no_num[,-c(1:8)]

df_input <- data[,-(39)]
df_input <- df_input[,-no_num]
data_target <- data[,n_var]
data <- data.frame(df_input, data_no_num, Class = data_target)

df_target_1 <- data[data$Class == 1,]
df_target_1 <- df_target_1[sample(nrow(df_target_1), 1000),]
df_target_2 <- data[data$Class == 2,]
df_target_2 <- df_target_2[sample(nrow(df_target_2), 6000),]
df_target_3 <- data[data$Class == 3,]
df_target_3 <- df_target_3[sample(nrow(df_target_3), 3500),]

data <- rbind(df_target_1,df_target_2,df_target_3)


set.seed(365)
trn_idx_s <- sample(1:nrow(data), size=8076)
val_idx_s <- sample(trn_idx_s, size=2019)
data_tst_s <- data[-trn_idx_s,]
data_val_s <- data[val_idx_s,]
data_trn_s <- data[trn_idx_s[!(trn_idx_s%in%val_idx_s)],]

trn_input_s <- data_trn_s[,-69]
trn_target_s <- class.ind(data_trn_s[,69])
val_input_s <- data_val_s[,-69]
val_target_s <- class.ind(data_val_s[,69])
tst_input_s <- data_tst_s[,-69]
tst_target_s <- class.ind(data_tst_s[,69])

boots <- seq(30,300,30)
perf_ann_bag <- matrix(0, nrow=300,ncol=3)
colnames(perf_ann_bag) <- c("Bootstrap", "ACC", "BCR")

num <- 1
start_time <- proc.time()
for (i in 1:10) {
  for (j in (1:30)){
    cat("ANN Bagging Bootstrap : ", i,"_", j, "\n")
    Bagging.ANN.model <- avNNet(trn_input_s, trn_target_s, size = 19, decay = 5e-4,rang=0.1, maxit = 200,
                                repeats = boots[i], bag = TRUE, allowParallel = TRUE, trace = TRUE,MaxNWts = 10000)
    
    Bagging.ANN.prey <- predict(Bagging.ANN.model, newdata = val_input_s)
    Bagging.ANN.cfm <- table(factor(max.col(val_target_s),levels=c("1","2","3")), factor(max.col(Bagging.ANN.prey),levels=c("1","2","3")))
    
    
    perf_ann_bag[num, 1] <- boots[i]
    perf_ann_bag[num, 2:3] <- perf_eval_multi(Bagging.ANN.cfm)    
    num <- num + 1
    
  } 
}
proc.time() - start_time
perf_ann_bag
perf_ann_boot_acc <- matrix(0, nrow= 10, ncol=3)
colnames(perf_ann_boot_acc) <- c("bootstrap","mean","sd")

j <- 0
for (i in 1:10) {
  perf_ann_boot_acc[i,1] <- boots[i]
  perf_ann_boot_acc[i,2] <- mean(perf_ann_bag[(i+j):(i+j+29),2])
  perf_ann_boot_acc[i,3] <- sd(perf_ann_bag[(i+j):(i+j+29),2])
  j <- j + 29
}
perf_ann_boot_acc

perf_ann_boot_bcr <- matrix(0, nrow= 10, ncol=3)
colnames(perf_ann_boot_bcr) <- c("bootstrap","mean","sd")

j <- 0
for (i in 1:10) {
  perf_ann_boot_bcr[i,1] <- boots[i]
  perf_ann_boot_bcr[i,2] <- mean(perf_ann_bag[(i+j):(i+j+29),3])
  perf_ann_boot_bcr[i,3] <- sd(perf_ann_bag[(i+j):(i+j+29),3])
  j <- j + 29
}
perf_ann_boot_bcr

# predict
trn_input_s <- rbind(trn_input_s, val_input_s)
trn_target_s <- rbind(trn_target_s, val_target_s)
Bagging.ANN.model <- avNNet(trn_input_s, trn_target_s, size = 19, decay = 5e-4,rang=0.1, maxit = 200,
                            repeats = 270, bag = TRUE, allowParallel = TRUE, trace = TRUE,MaxNWts = 10000)

Bagging.ANN.prey <- predict(Bagging.ANN.model, newdata = tst_input_s)
Bagging.ANN.cfm <- table(factor(max.col(tst_target_s),levels=c("1","2","3")), factor(max.col(Bagging.ANN.prey),levels=c("1","2","3")))
Bagging.ANN.cfm

perf_ann_bag_best <- matrix(0,1,3)
rownames(perf_ann_bag_best) <- c("ANN.bagging")
colnames(perf_ann_bag_best) <- c("Bootstrap","ACC","BCR")

perf_ann_bag_best[1, 1] <- 270
perf_ann_bag_best[1, 2:3] <- perf_eval_multi(Bagging.ANN.cfm)
perf_ann_bag_best


# [Q6], AdaBoost
iter_num <- c(50, 100, 150, 200)
frac <- c(0.3,0.5,0.7)
perf_ada <- matrix(0, )

CT_trn <- data.frame(df_input[trn_idx[!(trn_idx%in%val_idx)],], dfYN = df_target[trn_idx[!(trn_idx%in%val_idx)]])
CT_val <- data.frame(df_input[val_idx,], dfYN = df_target[val_idx])
CT_tst <- data.frame(df_input[-trn_idx,], dfYN = df_target[-trn_idx])

Ada.trn <- CT_trn
Ada.val <- CT_val
Ada.tst <- CT_tst


for (i in 1:length(iter_num)) {
  cat("Training AdaBoost: iter: ", iter[i], "\n")
  for (j in 1:length(frac)) {
    AdaBoost.model <- ada(Ada.trn[,-69],Ada.trn[,69], loss = "exponential",
                          iter=iter_num[i], bag.frac = frac[j], verbose = TRUE)
    AdaBoost.prey <- predict(AdaBoost.model, Ada.tst[,-69])
    
    
  }
  
}
AdaBoost.model <- ada(trn_input,trn_target, loss = "exponential",
)