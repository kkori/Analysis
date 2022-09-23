
install.packages("dummies")
install.packages("moments")
install.packages("dplyr")
install.packages("ggplot2")
install.pa
library(dummies)
library(moments)
library(dplyr)
library(ggplot2)
library(GGally)
library(corrplot)
library(psych)
library(glmnet)
library(car)ckages("GGally")
install.packages("corrplot")
install.packages("psych")
install.packages("glmnet")
install.packages("car")


#[Q1]
df <- read.csv("Life Expectancy Data.csv")

# 결측치 확인
table(is.na(df))

# 결측치를 평균값으로 대체
for (i in 4:22){
  df[i][is.na(df[i])] <- mean(df[,i], na.rm=TRUE)
}

# 다시 결측치 확인
table(is.na(df))

# 더미형 변수로 변환
nCou <- nrow(df)
nVar <- ncol(df)

no_idx <- c(1,2,6,22)
category_id <- 3

dummy_s <- rep(0, nCou)
dummy_d <- rep(0, nCou)

s_idx <- which(df$Status == 'Developed')
d_idx <- which(df$Status == 'Developing')

dummy_s[s_idx] <- 1
dummy_d[d_idx] <- 1

sta <- data.frame(dummy_s, dummy_d)
names(sta) <- c("Developed", "Developing")

df <- cbind(df[,-c(no_idx, category_id)], sta)
colnames(df)

#[Q3]
# 변수별 단변량 통계 및 정규성 검정
info <- function(x){
  Mean <- mean(x)
  Standard_deviation <- sd(x)
  Skewness <- skewness(x)
  Kurtosis <- kurtosis(x)
  result <- data.frame(Mean, Standard_deviation, Skewness, Kurtosis)
  
  return(result)
}

# Adult Mortality
info(df$Adult.Mortality)
boxplot(df$Adult.Mortality, col="blue", main='Boxplot_Adult.Mortality')
# 이상치 존재
qqnorm(df$Adult.Mortality, main="Q-Q plot_Adult.Mortality")
qqline(df$Adult.Mortality, col="red")
hist(df$Adult.Mortality, freq=FALSE,main="Normal curve over histogram")
curve(dnorm(x, mean=info(df$Adult.Mortality)$Mean, sd=info(df$Adult.Mortality)$Standard_deviation), col="darkblue",lwd=2, add=TRUE,yaxt="n")
shapiro.test(df$Adult.Mortality)
# P-value < 0.05 => 기각 => 정규성x

# Alcohol
info(df$Alcohol)
boxplot(df$Alcohol, col="blue", main='Boxplot_Alcohol')
#이상치 존재
qqnorm(df$Alcohol, main="Q-Q plot_Alcohol")
qqline(df$Alcohol, col="red")
hist(df$Alcohol, freq=FALSE,main="Normal curve over histogram")
curve(dnorm(x, mean=info(df$Alcohol)$Mean, sd=info(df$Alcohol)$Standard_deviation), col="darkblue",lwd=2, add=TRUE,yaxt="n")
shapiro.test(df$Alcohol)

# percentage expenditure
info(df$percentage.expenditure)
boxplot(df$percentage.expenditure, col="blue", main='Boxplot_percentage.expenditure')
#이상치 존재
qqnorm(df$percentage.expenditure, main="Q-Q plot_percentage.expenditure")
qqline(df$percentage.expenditure, col="red")
hist(df$percentage.expenditure, freq=FALSE,main="Normal curve over histogram")
curve(dnorm(x, mean=info(df$percentage.expenditure)$Mean,
            sd=info(df$percentage.expenditure)$Standard_deviation),col="darkblue",lwd=2,
      add=TRUE,yaxt="n")
shapiro.test(df$percentage.expenditure)

# Hepatitis B
info(df$Hepatitis.B)
boxplot(df$Hepatitis.B, col="blue", main='Boxplot_Hepatitis.B')
#이상치 존재
qqnorm(df$Hepatitis.B, main="Q-Q plot_Hepatitis.B")
qqline(df$Hepatitis.B, col="red")
hist(df$Hepatitis.B, freq=FALSE,main="Normal curve over histogram")
curve(dnorm(x, mean=info(df$Hepatitis.B)$Mean,
            sd=info(df$Hepatitis.B)$Standard_deviation),col="darkblue",lwd=2,
      add=TRUE,yaxt="n")
shapiro.test(df$Hepatitis.B)

# Measles
info(df$Measles)
boxplot(df$Measles, col="blue", main='Boxplot_Measles')
#이상치 존재
qqnorm(df$Measles, main="Q-Q plot_Measles")
qqline(df$Measles, col="red")
hist(df$Measles, freq=FALSE,main="Normal curve over histogram")
curve(dnorm(x, mean=info(df$Measles)$Mean,
            sd=info(df$Measles)$Standard_deviation),col="darkblue",lwd=2,
      add=TRUE,yaxt="n")
shapiro.test(df$Measles)

# BMI
info(df$BMI)
boxplot(df$BMI, col="blue", main='Boxplot_BMI')
#이상치 x
qqnorm(df$BMI, main="Q-Q plot_BMI")
qqline(df$BMI, col="red")
hist(df$BMI, freq=FALSE,main="Normal curve over histogram")
curve(dnorm(x, mean=info(df$BMI)$Mean,
            sd=info(df$BMI)$Standard_deviation),col="darkblue",lwd=2,
      add=TRUE,yaxt="n")
shapiro.test(df$BMI)

# Polio
info(df$Polio)
boxplot(df$Polio, col="blue", main='Boxplot_Polio')
#이상치 존재
qqnorm(df$Polio, main="Q-Q plot_Polio")
qqline(df$Polio, col="red")
hist(df$Polio, freq=FALSE,main="Normal curve over histogram")
curve(dnorm(x, mean=info(df$Polio)$Mean,
            sd=info(df$Polio)$Standard_deviation),col="darkblue",lwd=2,
      add=TRUE,yaxt="n")
shapiro.test(df$Polio)

# Total expenditure
info(df$Total.expenditure)
boxplot(df$Total.expenditure, col="blue", main='Boxplot_Total.expenditure')
#이상치 존재
qqnorm(df$Total.expenditure, main="Q-Q plot_Total.expenditure")
qqline(df$Total.expenditure, col="red")
hist(df$Total.expenditure, freq=FALSE,main="Normal curve over histogram")
curve(dnorm(x, mean=info(df$Total.expenditure)$Mean,
            sd=info(df$Total.expenditure)$Standard_deviation),col="darkblue",lwd=2,
      add=TRUE,yaxt="n")
shapiro.test(df$Total.expenditure)

# Diphtheria
info(df$Diphtheria)
boxplot(df$Diphtheria, col="blue", main='Boxplot_Diphtheria')
#이상치 존재
qqnorm(df$Diphtheria, main="Q-Q plot_Diphtheria")
qqline(df$Diphtheria, col="red")
hist(df$Diphtheria, freq=FALSE,main="Normal curve over histogram")
curve(dnorm(x, mean=info(df$Diphtheria)$Mean,
            sd=info(df$Diphtheria)$Standard_deviation),col="darkblue",lwd=2,
      add=TRUE,yaxt="n")
shapiro.test(df$Diphtheria)

# HIV/AIDS
info(df$HIV.AIDS)
boxplot(df$HIV.AIDS, col="blue", main='Boxplot_HIV.AIDS')
#이상치 존재
qqnorm(df$HIV.AIDS, main="Q-Q plot_HIV.AIDS")
qqline(df$HIV.AIDS, col="red")
hist(df$HIV.AIDS, freq=FALSE,main="Normal curve over histogram")
curve(dnorm(x, mean=info(df$HIV.AIDS)$Mean,
            sd=info(df$HIV.AIDS)$Standard_deviation),col="darkblue",lwd=2,
      add=TRUE,yaxt="n")
shapiro.test(df$HIV.AIDS)

# GDP
info(df$GDP)
boxplot(df$GDP, col="blue", main='Boxplot_GDP')
#이상치 존재
qqnorm(df$GDP, main="Q-Q plot_GDP")
qqline(df$GDP, col="red")
hist(df$GDP, freq=FALSE,main="Normal curve over histogram")
curve(dnorm(x, mean=info(df$GDP)$Mean,
            sd=info(df$GDP)$Standard_deviation),col="darkblue",lwd=2,
      add=TRUE,yaxt="n")
shapiro.test(df$GDP)

# Population
info(df$Population)
boxplot(df$Population, col="blue", main='Boxplot_Population')
#이상치 존재
qqnorm(df$Population, main="Q-Q plot_Population")
qqline(df$Population, col="red")
hist(df$Population, freq=FALSE,main="Normal curve over histogram")
curve(dnorm(x, mean=info(df$Population)$Mean,
            sd=info(df$Population)$Standard_deviation),col="darkblue",lwd=2,
      add=TRUE,yaxt="n")
shapiro.test(df$Population)

# thinness 1-19 years
info(df$thinness..1.19.years)
boxplot(df$thinness..1.19.years, col="blue", main='Boxplot_thinness..1.19.years')
#이상치 존재
qqnorm(df$thinness..1.19.years, main="Q-Q plot_thinness..1.19.years")
qqline(df$thinness..1.19.years, col="red")
hist(df$thinness..1.19.years, freq=FALSE,main="Normal curve over histogram")
curve(dnorm(x, mean=info(df$thinness..1.19.years)$Mean,
            sd=info(df$thinness..1.19.years)$Standard_deviation),col="darkblue",lwd=2,
      add=TRUE,yaxt="n")
shapiro.test(df$thinness..1.19.years)

# thinness 5-9 years
info(df$thinness.5.9.years)
boxplot(df$thinness.5.9.years, col="blue", main='Boxplot_thinness.5.9.years')
#이상치 존재
qqnorm(df$thinness.5.9.years, main="Q-Q plot_thinness.5.9.years")
qqline(df$thinness.5.9.years, col="red")
hist(df$thinness.5.9.years, freq=FALSE,main="Normal curve over histogram")
curve(dnorm(x, mean=info(df$thinness.5.9.years)$Mean,
            sd=info(df$thinness.5.9.years)$Standard_deviation),col="darkblue",lwd=2,
      add=TRUE,yaxt="n")
shapiro.test(df$thinness.5.9.years)

# Income composition of resource
info(df$Income.composition.of.resources)
boxplot(df$Income.composition.of.resources, col="blue", main='Boxplot_Income.composition.of.resources')
#이상치 존재
qqnorm(df$Income.composition.of.resources, main="Q-Q plot_Income.composition.of.resources")
qqline(df$Income.composition.of.resources, col="red")
hist(df$Income.composition.of.resources, freq=FALSE,main="Normal curve over histogram")
curve(dnorm(x, mean=info(df$Income.composition.of.resources)$Mean,
            sd=info(df$Income.composition.of.resources)$Standard_deviation),col="darkblue",lwd=2,
      add=TRUE,yaxt="n")
shapiro.test(df$Income.composition.of.resources)

# Developed
info(df$Developed)
boxplot(df$Developed, col="blue", main='Boxplot_Developed')
qqnorm(df$Developed, main="Q-Q plot_Developed")
qqline(df$Developed, col="red")
hist(df$Developed, freq=FALSE,main="Normal curve over histogram")
curve(dnorm(x, mean=info(df$Developed)$Mean,
            sd=info(df$Developed)$Standard_deviation),col="darkblue",lwd=2,
      add=TRUE,yaxt="n")
shapiro.test(df$Developed)

# Developing
info(df$Developing)
boxplot(df$Developing, col="blue", main='Boxplot_Developing')
qqnorm(df$Developing, main="Q-Q plot_Developing")
qqline(df$Developing, col="red")
hist(df$Developing, freq=FALSE,main="Normal curve over histogram")
curve(dnorm(x, mean=info(df$Developing)$Mean,
            sd=info(df$Developing)$Standard_deviation),col="darkblue",lwd=2,
      add=TRUE,yaxt="n")
shapiro.test(df$Developing)

# Under.five.deaths
info(df$under.five.deaths)
boxplot(df$under.five.deaths, col="blue", main='Boxplot_under.five.deaths')
qqnorm(df$under.five.deaths, main="Q-Q plot_under.five.deaths")
qqline(df$under.five.deaths, col="red")
hist(df$under.five.deaths, freq=FALSE,main="Normal curve over histogram")
curve(dnorm(x, mean=info(df$under.five.deaths)$Mean,
            sd=info(df$under.five.deaths)$Standard_deviation),col="darkblue",lwd=2,
      add=TRUE,yaxt="n")
shapiro.test(df$under.five.deaths)

# [Q4]
# 이상치 제거, 기존 Box plot을 확인 후에 위,아래 혹은 전체 부분의 이상치를 제거한다.
# A.M
df <- df[-which(df$Adult.Mortality>(summary(df$Adult.Mortality)[5] + 1.5*IQR(df$Adult.Mortality))),]
boxplot(df$Adult.Mortality, main="Boxplot_Adult.Mortality")

# Alcohol
df <- df[-which(df$Alcohol>summary(df$Alcohol)[5] + 1.5*IQR(df$Alcohol)),]
boxplot(df$Alcohol,main="Boxplot_Alcohol")

# P.E, 많이 진행하면 데이터 손실
df <- df[-which(df$percentage.expenditure>(summary(df$percentage.expenditure)[5] + 1.5*IQR(df$percentage.expenditure))),]
boxplot(df$percentage.expenditure, main="Boxplot_Percentage.expenditure")

# H.B
df <- df[-which(df$Hepatitis.B<summary(df$Hepatitis.B)[2] - 1.5*IQR(df$Hepatitis.B)),]
boxplot(df$Hepatitis.B,main="Boxplot_Hepatitis.B")

# Measles, 많이 진행하면 데이터 손실
df <- df[-which(df$Measles>summary(df$Measles)[5] + 1.5*IQR(df$Measles)),]
boxplot(df$Measles,main="Boxplot_Measles")

# Polio
df <- df[-which(df$Polio<summary(df$Polio)[2] - 1.5*IQR(df$Polio)),]
boxplot(df$Polio,main="Boxplot_Polio")

# T.E
df <- df[-which(df$Total.expenditure<summary(df$Total.expenditure)[2] - 1.5*IQR(df$Total.expenditure)),]
boxplot(df$Total.expenditure,main="Boxplot_Total.expenditure")

df <- df[-which(df$Total.expenditure>summary(df$Total.expenditure)[5] + 1.5*IQR(df$Total.expenditure)),]
boxplot(df$Total.expenditure,main="Boxplot_Total.expenditure")

# Diphtheria
df <- df[-which(df$Diphtheria<summary(df$Diphtheria)[2] - 1.5*IQR(df$Diphtheria)),]
boxplot(df$Diphtheria,main="Boxplot_Diphtheria")

# HIV/AIDS
df <- df[-which(df$HIV.AIDS>summary(df$HIV.AIDS)[5] + 1.5*IQR(df$HIV.AIDS)),]
boxplot(df$HIV.AIDS,main="HIV/AIDS")

# GDP
df <- df[-which(df$GDP>summary(df$GDP)[5] + 1.5*IQR(df$GDP)),]
boxplot(df$GDP,main="Boxplot_GDP")

# Population
df <- df[-which(df$Population>summary(df$Population)[5] + 1.5*IQR(df$Population)),]
boxplot(df$Population,main="Boxplot_Population")

# Thinness 1-19
df <- df[-which(df$thinness..1.19.years>summary(df$thinness..1.19.years)[5] + 1.5*IQR(df$thinness..1.19.years)),]
boxplot(df$thinness..1.19.years,main="Boxplot_thinness.1-19")

# Thinness 5-9
df <- df[-which(df$thinness.5.9.years>summary(df$thinness.5.9.years)[5] + 1.5*IQR(df$thinness.5.9.years)),]
boxplot(df$thinness.5.9.years,main="Boxplot_thinness.5-9")

# I.C.R
df <- df[-which(df$Income.composition.of.resources<summary(df$Income.composition.of.resources)[2] - 1.5*IQR(df$Income.composition.of.resources)),]
boxplot(df$Income.composition.of.resources, main="Boxplot_Income.composition.of.resources")

# U.F.D
df <- df[-which(df$under.five.deaths>summary(df$under.five.deaths)[5] + 1.5*IQR(df$under.five.deaths)),]
boxplot(df$under.five.deaths,main="Boxplot_Under.five.deaths")

# [Q5]
# Scatter Plot

# A.M
col_name <- colnames(df)
col_name
par(mfrow = c(4,5))
for (i in c(3:19)){
  plot(df[,2], df[,i], xlab=col_name[2], ylab=col_name[i], col="blue")
}

# Alcohol
par(mfrow = c(4,5))
for (i in c(2,4:19)){
  plot(df[,3], df[,i], xlab=col_name[3], ylab=col_name[i], col="blue")
}

# P.E
par(mfrow = c(4,5))
for (i in c(2:3,5:19)){
  plot(df[,4], df[,i], xlab=col_name[4], ylab=col_name[i], col="blue")
}

# H.B
par(mfrow = c(4,5))
for (i in c(2:4,6:19)){
  plot(df[,5], df[,i], xlab=col_name[5], ylab=col_name[i], col="blue")
}

# Measles
par(mfrow = c(4,5))
for (i in c(2:5, 7:19)){
  plot(df[,6], df[,i], xlab=col_name[6], ylab=col_name[i], col="blue")
}

# BMI
par(mfrow = c(4,5))
for (i in c(2:6, 8:19)){
  plot(df[,7], df[,i], xlab=col_name[7], ylab=col_name[i], col="blue")
}

# U.F
par(mfrow = c(4,5))
for (i in c(2:7, 9:19)){
  plot(df[,8], df[,i], xlab=col_name[8], ylab=col_name[i], col="blue")
}

# Polio
par(mfrow = c(4,5))
for (i in c(2:8, 10:19)){
  plot(df[,9], df[,i], xlab=col_name[9], ylab=col_name[i], col="blue")
}

# T.E
par(mfrow = c(4,5))
for (i in c(2:9, 11:19)){
  plot(df[,10], df[,i], xlab=col_name[10], ylab=col_name[i], col="blue")
}

# Diphtheria
par(mfrow = c(4,5))
for (i in c(2:10, 12:19)){
  plot(df[,11], df[,i], xlab=col_name[11], ylab=col_name[i], col="blue")
}

# HIV.AIDS
par(mfrow = c(4,5))
for (i in c(2:11,13:19)){
  plot(df[,12], df[,i], xlab=col_name[12], ylab=col_name[i], col="blue")
}

# GDP
par(mfrow = c(4,5))
for (i in c(2:12, 14:19)){
  plot(df[,13], df[,i], xlab=col_name[13], ylab=col_name[i], col="blue")
}

# Population
par(mfrow = c(4,5))
for (i in c(2:13, 15:19)){
  plot(df[,14], df[,i], xlab=col_name[14], ylab=col_name[i], col="blue")
}

# Thinness 1-19
par(mfrow = c(4,5))
for (i in c(2:14,16:19)){
  plot(df[,15], df[,i], xlab=col_name[15], ylab=col_name[i], col="blue")
}

# Thinness 5-9
par(mfrow = c(4,5))
for (i in c(2:15,17:19)){
  plot(df[,16], df[,i], xlab=col_name[16], ylab=col_name[i], col="blue")
}

# I.C.R
par(mfrow = c(4,5))
for (i in c(2:16,18:19)){
  plot(df[,17], df[,i], xlab=col_name[17], ylab=col_name[i], col="blue")
}

# Developed
par(mfrow = c(4,5))
for (i in c(2:17,19)){
  plot(df[,18], df[,i], xlab=col_name[18], ylab=col_name[i], col="blue")
}

#Developing
par(mfrow = c(4,5))
for (i in c(2:18)){
  plot(df[,19], df[,i], xlab=col_name[19], ylab=col_name[i], col="blue")
}

# Correlation Plot
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
# 데이터셋 분리
set.seed(365)
ndf <- nrow(df)
ndf
df_trn_idx <- sample(1:ndf, round(0.7*ndf))
df_trn_data <- df[df_trn_idx,]
df_val_data <- df[-df_trn_idx,]

# MLR_Train
mlr_df <- lm(Life.expectancy ~ ., data = df_trn_data)
mlr_df
summary(mlr_df)
plot(mlr_df)

#resid check
df_resid <- resid(mlr_df)
m <- mean(df_resid)
std <- sqrt(var(df_resid))

hist(df_resid, density=20, breaks=50, prob=TRUE,xlab="x_variable", main="Normal curve over histogram")
curve(dnorm(x, mean=m, sd=std), col="darkblue", lwd=2, add=TRUE, yaxt="n")

skewness(df_resid)
kurtosis(df_resid)

#[Q8]
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

perf_mat <- matrix(0, nrow=1, ncol=3)
rownames(perf_mat) <- c("Life.Expectancy")
colnames(perf_mat) <- c("RMSE", "MAE", "MAPE")

mlr_df_haty <- predict(mlr_df, newdata = df_val_data)
perf_mat[1,] <- perf_reg(df_val_data$Life.expectancy, mlr_df_haty)
perf_mat

# [Q10]
# 일정변수만 선택 후, 모델링
mlr_df_new <- lm(Life.expectancy~ Income.composition.of.resources+
                   Adult.Mortality + thinness.5.9.years+Total.expenditure+
                   Population+Developed+BMI, data=df_trn_data)
summary(mlr_df_new)
plot(mlr_df_new)

df_resid_1 <- resid(mlr_df_new)
m_1 <- mean(df_resid_1)
std_1 <- sqrt(var(df_resid_1))
hist(df_resid_1, density=20, breaks=50, prob=TRUE,xlab="x_variable", main="Normal curve over histogram")
curve(dnorm(x, mean=m_1, sd=std_1), col="darkblue", lwd=2, add=TRUE, yaxt="n")

skewness(df_resid_1)
kurtosis(df_resid_1)

perf_mat_new <- matrix(0, nrow=1, ncol=3)
rownames(perf_mat_new) <- c("Life.Expectancy")
colnames(perf_mat_new) <- c("RMSE", "MAE", "MAPE")

mlr_df_haty_new <- predict(mlr_df_new, newdata = df_val_data)
perf_mat_new[1,] <- perf_reg(df_val_data$Life.expectancy, mlr_df_haty_new)
perf_mat_new

# [Extra Question]
# Ridge regression
lambdas <- seq(0.1, 100, by = .5)
df_x <- as.matrix(df[,-1])
df_y <- as.matrix(df[,1])
rd_fit <- cv.glmnet(df_x[df_trn_idx,],df_y[df_trn_idx],alpha=0, lambda=lambdas)
summary(rd_fit)
plot(rd_fit)

# 최적의 lambda결정
opt_lambda <- rd_fit$lambda.min
opt_lambda
fin <- glmnet(df_x[df_trn_idx,],df_y[df_trn_idx],alpha=0,lambda = opt_lambda)
coef(fin)

rd_pr <- predict(fin, s=opt_lambda, newx=df_x[-df_trn_idx,])
rd_haty <- predict(fin, newx=df_x[-df_trn_idx,])

perf_mat_rd <- matrix(0, nrow=1, ncol=3)
rownames(perf_mat_rd) <- c("Life.Expectancy")
colnames(perf_mat_rd) <- c("RMSE", "MAE", "MAPE")
perf_mat_rd[1,] <- perf_reg(df_y[-df_trn_idx,1], rd_haty)
perf_mat_rd

# Lasso regression
ls_fit <- cv.glmnet(df_x[df_trn_idx,],df_y[df_trn_idx], alpha=1,lambda=lambdas)
opt_lam_2 <- ls_fit$lambda.min
opt_lam_2
fin_ls <- glmnet(df_x[df_trn_idx,],df_y[df_trn_idx],alpha=1,lambda = opt_lam_2)
coef(fin_ls)

ls_pr <- predict(fin_ls, s=opt_lam_2, newx=df_x[-df_trn_idx,])
ls_haty <- predict(fin_ls, newx=df_x[-df_trn_idx,])

perf_mat_ls <- matrix(0, nrow=1, ncol=3)
rownames(perf_mat_ls) <- c("Life.Expectancy")
colnames(perf_mat_ls) <- c("RMSE", "MAE", "MAPE")
perf_mat_ls[1,] <- perf_reg(df_y[-df_trn_idx,1], ls_haty)
perf_mat_ls

# vif 확인
mlr_df_vi <- lm(Life.expectancy ~ .-Developing, data = df_trn_data)
summary(mlr_df_vi)
vif(mlr_df_vi)

# 다중공선성 제거 모델
mlr_df_fn <- lm(Life.expectancy ~ .-(Developing+thinness..1.19.years),data=df_trn_data)
summary(mlr_df_fn)
plot(mlr_df_fn)

# resid check
df_resid_fn <- resid(mlr_df_fn)
m <- mean(df_resid_fn)
std <- sqrt(var(df_resid_fn))

hist(df_resid_fn, density=20, breaks=50, prob=TRUE,xlab="x_variable", main="Normal curve over histogram")
curve(dnorm(x, mean=m, sd=std), col="darkblue", lwd=2, add=TRUE, yaxt="n")

skewness(df_resid_fn)
kurtosis(df_resid_fn)

perf_mat_fn <- matrix(0, nrow=1, ncol=3)
rownames(perf_mat_fn) <- c("Life.Expectancy")
colnames(perf_mat_fn) <- c("RMSE", "MAE", "MAPE")

mlr_df_hat_fn <- predict(mlr_df_fn, newdata = df_val_data)
perf_mat_fn[1,] <- perf_reg(df_val_data$Life.expectancy, mlr_df_hat_fn)
perf_mat_fn

