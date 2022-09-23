
# 라이브러리 불러오기
library(dummies)
library(moments)
library(dplyr)
library(ggplot2)
library(GGally)
library(corrplot)
library(psych)
library(glmnet)
library(car)

# 데이터 불러오기
df <- read.csv("heart_failure_clinical_records_dataset.csv")
str(df)
# 결측치 확인
table(is.na(df))

# [Q3] 단변량 통계량, 정규성 검정
info <- function(x){
  Mean <- mean(x)
  Standard_deviation <- sd(x)
  Skewness <- skewness(x)
  Kurtosis <- kurtosis(x)
  result <- data.frame(Mean, Standard_deviation, Skewness, Kurtosis)
  
  return(result)
}

# age
info(df$age)
boxplot(df$age, col="blue", main='Boxplot_age')
# 이상치 존재x
qqnorm(df$age, main="Q-Q plot_age")
qqline(df$age, col="red")
hist(df$age, freq=FALSE,main="Normal curve over histogram")
curve(dnorm(x, mean=info(df$age)$Mean, sd=info(df$age)$Standard_deviation), col="darkblue",lwd=2, add=TRUE,yaxt="n")
shapiro.test(df$age)

# CPK
info(df$creatinine_phosphokinase)
boxplot(df$creatinine_phosphokinase, col="blue", main='Boxplot_CPK')
# 이상치 존재
qqnorm(df$creatinine_phosphokinase, main="Q-Q plot_CPK")
qqline(df$creatinine_phosphokinase, col="red")
hist(df$creatinine_phosphokinase, freq=FALSE,main="Normal curve over histogram")
curve(dnorm(x, mean=info(df$creatinine_phosphokinase)$Mean, sd=info(df$creatinine_phosphokinase)$Standard_deviation), col="darkblue",lwd=2, add=TRUE,yaxt="n")
shapiro.test(df$creatinine_phosphokinase)

# ejection fraction
info(df$ejection_fraction)
boxplot(df$ejection_fraction, col="blue", main='Boxplot_ejection_fraction')
# 이상치 존재
qqnorm(df$ejection_fraction, main="Q-Q plot_ejection_fraction")
qqline(df$ejection_fraction, col="red")
hist(df$ejection_fraction, freq=FALSE,main="Normal curve over histogram")
curve(dnorm(x, mean=info(df$ejection_fraction)$Mean, sd=info(df$ejection_fraction)$Standard_deviation), col="darkblue",lwd=2, add=TRUE,yaxt="n")
shapiro.test(df$ejection_fraction)

# platelets
info(df$platelets)
boxplot(df$platelets, col="blue", main='Boxplot_platelets')
# 이상치 존재
qqnorm(df$platelets, main="Q-Q plot_platelets")
qqline(df$platelets, col="red")
hist(df$platelets, freq=FALSE,main="Normal curve over histogram")
curve(dnorm(x, mean=info(df$platelets)$Mean, sd=info(df$platelets)$Standard_deviation), col="darkblue",lwd=2, add=TRUE,yaxt="n")
shapiro.test(df$platelets)

# serum creatinine
info(df$serum_creatinine)
boxplot(df$serum_creatinine, col="blue", main='Boxplot_serum_creatinine')
# 이상치 존재
qqnorm(df$serum_creatinine, main="Q-Q plot_serum_creatinine")
qqline(df$serum_creatinine, col="red")
hist(df$serum_creatinine, freq=FALSE,main="Normal curve over histogram")
curve(dnorm(x, mean=info(df$serum_creatinine)$Mean, sd=info(df$serum_creatinine)$Standard_deviation), col="darkblue",lwd=2, add=TRUE,yaxt="n")
shapiro.test(df$serum_creatinine)

# serum sodium
info(df$serum_sodium)
boxplot(df$serum_sodium, col="blue", main='Boxplot_serum_sodium')
# 이상치 존재
qqnorm(df$serum_sodium, main="Q-Q plot_serum_sodium")
qqline(df$serum_sodium, col="red")
hist(df$serum_sodium, freq=FALSE,main="Normal curve over histogram")
curve(dnorm(x, mean=info(df$serum_sodium)$Mean, sd=info(df$serum_sodium)$Standard_deviation), col="darkblue",lwd=2, add=TRUE,yaxt="n")
shapiro.test(df$serum_sodium)

# time
info(df$time)
boxplot(df$time, col="blue", main='Boxplot_time')

qqnorm(df$time, main="Q-Q plot_time")
qqline(df$time, col="red")
hist(df$time, freq=FALSE,main="Normal curve over histogram")
curve(dnorm(x, mean=info(df$time)$Mean, sd=info(df$time)$Standard_deviation), col="darkblue",lwd=2, add=TRUE,yaxt="n")
shapiro.test(df$time)

# [Q4], 위 boxplot을 토대로 이상치가 존재한 변수들에 대해서 이상치 제거
# CPK
df <- df[-which(df$creatinine_phosphokinase>(summary(df$creatinine_phosphokinase)[5] + 1.5*IQR(df$creatinine_phosphokinase))),]
boxplot(df$creatinine_phosphokinase,main="Boxplot_CPK", col="blue")

# ejection fraction
df <- df[-which(df$ejection_fraction>(summary(df$ejection_fraction)[5] + 1.5*IQR(df$ejection_fraction))),]
boxplot(df$ejection_fraction,main="Boxplot_ejection_fraction", col="blue")

# platelets
df <- df[-which(df$platelets>(summary(df$platelets)[5] + 1.5*IQR(df$platelets))),]
df <- df[-which(df$platelets<(summary(df$platelets)[2] - 1.5*IQR(df$platelets))),]
boxplot(df$platelets,main="Boxplot_platelets", col="blue")

# serum creatinine
df <- df[-which(df$serum_creatinine>(summary(df$serum_creatinine)[5] + 1.5*IQR(df$serum_creatinine))),]
boxplot(df$serum_creatinine,main="Boxplot_serum_creatinine", col="blue")

# serum sodium
df <- df[-which(df$serum_sodium<(summary(df$serum_sodium)[2] - 1.5*IQR(df$serum_sodium))),]
boxplot(df$serum_sodium,main="Boxplot_serum_sodium", col="blue")

# 기존에 이상치가 없던 변수들의 box plot
boxplot(df$age, col="blue", main='Boxplot_age')
boxplot(df$time, col="blue", main='Boxplot_time')

# [Q5], Scatter plot / Correlation plot
col_names <- colnames(df)
col_names

# age
par(mfrow = c(3,4))
for (i in c(2:12)){
  plot(df[,1], df[,i], xlab=col_names[1], ylab=col_names[i], col="blue")
}

# anaemia
par(mfrow = c(3,4))
for (i in c(1,3:12)){
  plot(df[,2], df[,i], xlab=col_names[2], ylab=col_names[i], col="blue")
}

# CPK
par(mfrow = c(3,4))
for (i in c(1:2,4:12)){
  plot(df[,3], df[,i], xlab=col_names[3], ylab=col_names[i], col="blue")
}

# diabetes
par(mfrow = c(3,4))
for (i in c(1:3,5:12)){
  plot(df[,4], df[,i], xlab=col_names[4], ylab=col_names[i], col="blue")
}

# ejection fraction
par(mfrow = c(3,4))
for (i in c(1:4,6:12)){
  plot(df[,5], df[,i], xlab=col_names[5], ylab=col_names[i], col="blue")
}

# high blood pressure
par(mfrow = c(3,4))
for (i in c(1:5,7:12)){
  plot(df[,6], df[,i], xlab=col_names[6], ylab=col_names[i], col="blue")
}

# platelets
par(mfrow = c(3,4))
for (i in c(1:6,8:12)){
  plot(df[,7], df[,i], xlab=col_names[7], ylab=col_names[i], col="blue")
}

# serum creatinine
par(mfrow = c(3,4))
for (i in c(1:7,9:12)){
  plot(df[,8], df[,i], xlab=col_names[8], ylab=col_names[i], col="blue")
}

# serum sodium
par(mfrow = c(3,4))
for (i in c(1:8,10:12)){
  plot(df[,9], df[,i], xlab=col_names[9], ylab=col_names[i], col="blue")
}

# sex
par(mfrow = c(3,4))
for (i in c(1:9,11:12)){
  plot(df[,10], df[,i], xlab=col_names[10], ylab=col_names[i], col="blue")
}

# smoking
par(mfrow = c(3,4))
for (i in c(1:10,12)){
  plot(df[,11], df[,i], xlab=col_names[11], ylab=col_names[i], col="blue")
}

# time
par(mfrow = c(3,4))
for (i in c(1:11)){
  plot(df[,12], df[,i], xlab=col_names[12], ylab=col_names[i], col="blue")
}

# correlation plot
par(mfrow = c(1,1))
col <- colorRampPalette(c("#BB4444","#EE9988","#FFFFFF","#77AADD","#4477AA"))
r <- cor(df, method = "pearson")
corrplot(round(r,2),
         method = "color",
         col = col(200),
         type="lower",
         order="hclust",
         number.cex=.7,
         addCoef.col="black",
         tl.col="black",
         tl.srt=15,
         sig.level=0.01,
         insig="blank",
         diag=FALSE)
# [Q6]
perf_eval <- function(cm){
  TPR <- cm[2,2]/sum(cm[2,])
  PRE <- cm[2,2]/sum(cm[,2])
  TNR <- cm[1,1]/sum(cm[1,])
  ACC <- (cm[1,1]+cm[2,2])/sum(cm)
  BCR <- sqrt(TPR*TNR)
  F1 <- 2*TPR*PRE/(TPR+PRE)
  
  return(c(ACC, BCR, F1))
}

perf_mat <- matrix(0,1,3)
colnames(perf_mat) <- c("ACC","BCR","F1")
rownames(perf_mat) <- "Logistic Regression"

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
heart_df <- rename(heart_df,"DEATH_EVENT"="df_target")

# 시드 설정 후, 70:30 비율로 train data와 test data 분할
set.seed(365)
trn_idx <- sample(1:nrow(heart_df), round(0.7*nrow(heart_df)))
heart_trn <- heart_df[trn_idx,]
heart_tst <- heart_df[-trn_idx,]

# model 학습
full_lr <- glm(DEATH_EVENT ~ ., family = binomial, heart_trn)
summary(full_lr)
# vif계수 확인
vif(full_lr)

# train data, confusion matrix
tr_response <- predict(full_lr, type="response", newdata=heart_trn)
tr_target <- heart_trn$DEATH_EVENT
tr_predicted <- rep(0,length(tr_target))
tr_predicted[which(tr_response >= 0.5)] <- 1
cm_trn_full <- table(tr_target, tr_predicted)
cm_trn_full

perf_mat[1,] <- perf_eval(cm_trn_full)
perf_mat

# test data, confusion matrix
tst_response <- predict(full_lr, type="response", newdata=heart_tst)
tst_target <- heart_tst$DEATH_EVENT
tst_predicted <- rep(0,length(tst_target))
tst_predicted[which(tst_response >= 0.5)] <- 1
cm_tst_full <- table(tst_target, tst_predicted)
cm_tst_full

perf_mat[1,] <- perf_eval(cm_tst_full)
perf_mat

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

AUROC(tr_response, tr_target)
AUROC(tst_response, tst_target)

# [Q7]
heart_df_ns <- heart_df[,-11]
set.seed(365)
trn_idx_ns <- sample(1:nrow(heart_df_ns), round(0.7*nrow(heart_df_ns)))
heart_trn_ns <- heart_df_ns[trn_idx_ns,]
heart_tst_ns <- heart_df_ns[-trn_idx_ns,]

# model 학습
full_lr_ns <- glm(DEATH_EVENT ~ ., family = binomial, heart_trn_ns)
summary(full_lr_ns)

# train data, confusion matrix
tr_response <- predict(full_lr_ns, type="response", newdata=heart_trn_ns)
tr_target <- heart_trn_ns$DEATH_EVENT
tr_predicted <- rep(0,length(tr_target))
tr_predicted[which(tr_response >= 0.5)] <- 1
cm_trn_ns <- table(tr_target, tr_predicted)
cm_trn_ns

perf_mat[1,] <- perf_eval(cm_trn_ns)
perf_mat

# test data, confusion matrix
tst_response <- predict(full_lr_ns, type="response", newdata=heart_tst_ns)
tst_target <- heart_tst_ns$DEATH_EVENT
tst_predicted <- rep(0,length(tst_target))
tst_predicted[which(tst_response >= 0.5)] <- 1
cm_tst_ns <- table(tst_target, tst_predicted)
cm_tst_ns

perf_mat[1,] <- perf_eval(cm_tst_ns)
perf_mat

# AUROC
AUROC(tr_response, tr_target)
AUROC(tst_response, tst_target)

#[Q8], smoking변수 제거
heart_df_nsm <- heart_df[,-12]
set.seed(365)
trn_idx_nsm <- sample(1:nrow(heart_df_nsm), round(0.7*nrow(heart_df_nsm)))
heart_trn_nsm <- heart_df_nsm[trn_idx_nsm,]
heart_tst_nsm <- heart_df_nsm[-trn_idx_nsm,]

# model 학습
full_lr_nsm <- glm(DEATH_EVENT ~ ., family = binomial, heart_trn_nsm)
summary(full_lr_nsm)

# train data, confusion matrix
tr_response <- predict(full_lr_nsm, type="response", newdata=heart_trn_nsm)
tr_target <- heart_trn_nsm$DEATH_EVENT
tr_predicted <- rep(0,length(tr_target))
tr_predicted[which(tr_response >= 0.5)] <- 1
cm_trn_nsm <- table(tr_target, tr_predicted)
cm_trn_nsm

perf_mat[1,] <- perf_eval(cm_trn_nsm)
perf_mat

# test data, confusion matrix
tst_response <- predict(full_lr_nsm, type="response", newdata=heart_tst_nsm)
tst_target <- heart_tst_nsm$DEATH_EVENT
tst_predicted <- rep(0,length(tst_target))
tst_predicted[which(tst_response >= 0.5)] <- 1
cm_tst_nsm <- table(tst_target, tst_predicted)
cm_tst_nsm

perf_mat[1,] <- perf_eval(cm_tst_nsm)
perf_mat

# AUROC
AUROC(tr_response, tr_target)
AUROC(tst_response, tst_target)

# [Q8], seed에 따른 변화
# 100
set.seed(100)
trn_idx <- sample(1:nrow(heart_df), round(0.7*nrow(heart_df)))
heart_trn <- heart_df[trn_idx,]
heart_tst <- heart_df[-trn_idx,]

# model 학습
full_lr <- glm(DEATH_EVENT ~ ., family = binomial, heart_trn)
summary(full_lr)

# test data, confusion matrix
tst_response <- predict(full_lr, type="response", newdata=heart_tst)
tst_target <- heart_tst$DEATH_EVENT
tst_predicted <- rep(0,length(tst_target))
tst_predicted[which(tst_response >= 0.5)] <- 1
cm_tst_full <- table(tst_target, tst_predicted)
cm_tst_full

perf_mat[1,] <- perf_eval(cm_tst_full)
perf_mat

AUROC(tst_response, tst_target)

# 200
set.seed(200)
trn_idx <- sample(1:nrow(heart_df), round(0.7*nrow(heart_df)))
heart_trn <- heart_df[trn_idx,]
heart_tst <- heart_df[-trn_idx,]

# model 학습
full_lr <- glm(DEATH_EVENT ~ ., family = binomial, heart_trn)
summary(full_lr)

# test data, confusion matrix
tst_response <- predict(full_lr, type="response", newdata=heart_tst)
tst_target <- heart_tst$DEATH_EVENT
tst_predicted <- rep(0,length(tst_target))
tst_predicted[which(tst_response >= 0.5)] <- 1
cm_tst_full <- table(tst_target, tst_predicted)
cm_tst_full

perf_mat[1,] <- perf_eval(cm_tst_full)
perf_mat

AUROC(tst_response, tst_target)

# 300
set.seed(300)
trn_idx <- sample(1:nrow(heart_df), round(0.7*nrow(heart_df)))
heart_trn <- heart_df[trn_idx,]
heart_tst <- heart_df[-trn_idx,]

# model 학습
full_lr <- glm(DEATH_EVENT ~ ., family = binomial, heart_trn)
summary(full_lr)

# test data, confusion matrix
tst_response <- predict(full_lr, type="response", newdata=heart_tst)
tst_target <- heart_tst$DEATH_EVENT
tst_predicted <- rep(0,length(tst_target))
tst_predicted[which(tst_response >= 0.5)] <- 1
cm_tst_full <- table(tst_target, tst_predicted)
cm_tst_full

perf_mat[1,] <- perf_eval(cm_tst_full)
perf_mat
AUROC(tst_response, tst_target)

# 400
set.seed(400)
trn_idx <- sample(1:nrow(heart_df), round(0.7*nrow(heart_df)))
heart_trn <- heart_df[trn_idx,]
heart_tst <- heart_df[-trn_idx,]

# model 학습
full_lr <- glm(DEATH_EVENT ~ ., family = binomial, heart_trn)
summary(full_lr)

# test data, confusion matrix
tst_response <- predict(full_lr, type="response", newdata=heart_tst)
tst_target <- heart_tst$DEATH_EVENT
tst_predicted <- rep(0,length(tst_target))
tst_predicted[which(tst_response >= 0.5)] <- 1
cm_tst_full <- table(tst_target, tst_predicted)
cm_tst_full

perf_mat[1,] <- perf_eval(cm_tst_full)
perf_mat

AUROC(tst_response, tst_target)

# 500
set.seed(500)
trn_idx <- sample(1:nrow(heart_df), round(0.7*nrow(heart_df)))
heart_trn <- heart_df[trn_idx,]
heart_tst <- heart_df[-trn_idx,]

# model 학습
full_lr <- glm(DEATH_EVENT ~ ., family = binomial, heart_trn)
summary(full_lr)

# test data, confusion matrix
tst_response <- predict(full_lr, type="response", newdata=heart_tst)
tst_target <- heart_tst$DEATH_EVENT
tst_predicted <- rep(0,length(tst_target))
tst_predicted[which(tst_response >= 0.5)] <- 1
cm_tst_full <- table(tst_target, tst_predicted)
cm_tst_full

perf_mat[1,] <- perf_eval(cm_tst_full)
perf_mat
AUROC(tst_response, tst_target)

