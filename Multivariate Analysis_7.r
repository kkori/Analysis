
.libPaths("C:/Users/cjoi0/OneDrive/바탕 화면")
.libPaths()

install.packages("arules")
install.packages("Rgraphviz")
install.packages("arulesViz")
install.packages("wordcloud")

install.packages("BiocManager")
BiocManager::install("Rgraphviz")

library(wordcloud)
library(dplyr)
library(arules)
library(arulesViz)
library(Rgraphviz)

mooc_dataset <- read.csv("big_student_clear_third_version.csv")
str(mooc_dataset)
head(mooc_dataset)

Institute <- mooc_dataset["institute"]
Course <- mooc_dataset["course_id"]
Region <- mooc_dataset["final_cc_cname_DI"]
Degree <- mooc_dataset["LoE_DI"]

str(Region$final_cc_cname_DI)
Region<- gsub('\\s', '', Region$final_cc_cname_DI)
head(Region)
str(Region)

RawTransactions <- paste(Institute$institute, Course$course_id, Region, Degree$LoE_DI,sep="_")
str(RawTransactions)
head(RawTransactions)

MOOC_transactions <- paste(mooc_dataset$userid_DI, RawTransactions, sep = " ")
str(MOOC_transactions)
head(MOOC_transactions)

write.table(MOOC_transactions, "MOOC_User_Course.csv", row.names = FALSE, col.names = FALSE, sep=',', quote=FALSE)


# [Q2]
single_format <- read.transactions("MOOC_User_Course.csv", format='single', cols=c(1,2), rm.duplicates=TRUE)
inspect(single_format)
summary(single_format)


# 워드클라우드
item_name <- itemLabels(single_format)
item_count <- itemFrequency(single_format)*nrow(single_format)
col <- brewer.pal(8, "Pastel1")
wordcloud(words = item_name, freq = item_count, min.freq = 500, scale = c(3,0.1), col = col, random.order = FALSE)

itemFrequencyPlot(single_format, support = 0.01, cex.names = 0.6,col = "#50bcdf", main = "Support 0.01")

itemFrequencyPlot(single_format, support=0.01, cex.names=0.8, col="#50bcdf", topN=5, main="Top5")

# [Q3]
#rule 생성, 최소 10개 이상
rule1 <- apriori(single_format, parameter = list(support = 0.001, confidence = 0.001))
rule2 <- apriori(single_format, parameter = list(support = 0.001, confidence = 0.005))
rule3 <- apriori(single_format, parameter = list(support = 0.001, confidence = 0.01))
rule4 <- apriori(single_format, parameter = list(support = 0.001, confidence = 0.02))
rule5 <- apriori(single_format, parameter = list(support = 0.005, confidence = 0.001))
rule6 <- apriori(single_format, parameter = list(support = 0.005, confidence = 0.005))
rule7 <- apriori(single_format, parameter = list(support = 0.005, confidence = 0.01))
rule8 <- apriori(single_format, parameter = list(support = 0.005, confidence = 0.02))
rule9 <- apriori(single_format, parameter = list(support = 0.01, confidence = 0.001))
rule10 <- apriori(single_format, parameter = list(support = 0.01, confidence = 0.005))
rule11 <- apriori(single_format, parameter = list(support = 0.01, confidence = 0.01))
rule12 <- apriori(single_format, parameter = list(support = 0.01, confidence = 0.02))
rule13 <- apriori(single_format, parameter = list(support = 0.02, confidence = 0.001))
rule14 <- apriori(single_format, parameter = list(support = 0.02, confidence = 0.005))
rule15 <- apriori(single_format, parameter = list(support = 0.02, confidence = 0.01))
rule16 <- apriori(single_format, parameter = list(support = 0.02, confidence = 0.02))

# support = 0.001, confidence = 0.05
rules <- apriori(single_format, parameter = list(support = 0.001, confidence = 0.05))

summary(rules)
inspect(rules)
inspect(sort(rule, by = "support"))
inspect(sort(rules, by="confidence"))
inspect(sort(rules, by="lift"))

rules_df <- as.data.frame(inspect(rule3))
scl <- rules_df$support * rules_df*confidence
scl <- scl * rules_df$lift

rules_n <- data.frame(rules_df, scl)
head(rules_n[rev(order(rules_n$measure)),])

plot(rules, method="graph")

