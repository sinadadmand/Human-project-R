library(readr)
library(dplyr)
library(ROCR)
library(caret)

df <- read_csv("Disease probability.csv")
df <- df[df$Status=='reviewed',]
df$`Involvement in disease (source: UniProt)`[!is.na(df$`Involvement in disease (source: UniProt)`)] <- 1
df$`Involvement in disease (source: UniProt)`[is.na(df$`Involvement in disease (source: UniProt)`)] <- 0
df$`Involvement in disease (source: UniProt)` <- factor(df$`Involvement in disease (source: UniProt)`)

set.seed(100)
trainDataIndex <- createDataPartition(df$`Involvement in disease (source: UniProt)`, p=0.8, list = F)  # 70% training data
trainData <- df[trainDataIndex, ]
testData <- df[-trainDataIndex, ]
down_train <- downSample(x = df$`Mutual information estimation (in pits)`,
                         y = df$`Involvement in disease (source: UniProt)`)

model <- glm(`Involvement in disease (source: UniProt)`~ `Mutual information estimation (in pits)` + `Number of interactions`
             + `Network length` ,family = 'binomial', data= df)

summary(model)
confint(model)

pred <- predict(model, type = "response")
y_pred_num <- ifelse(pred > 0.5, 1, 0)
y_pred <- factor(y_pred_num, levels=c(0, 1))
y_act <- df$`Involvement in disease (source: UniProt)`
mean(y_pred == y_act)

table_mat <- table(df$`Involvement in disease (source: UniProt)`, pred > 0.5)
table_mat


ROCRpred <- prediction(pred, df$`Involvement in disease (source: UniProt)`)
ROCRperf <- performance(ROCRpred, 'tpr', 'fpr')
plot(ROCRperf, colorize = TRUE, text.adj = c(-0.2, 1.7))
