source("~/GIT/CPGARCH-HAR/R/scp.R")

# HAR

library(tidyverse)
library(caret)
library(xts)
library(highfrequency)
library(AdaptiveConformal)
library(rugarch)

set.seed(123)
RVSPY <- as.xts(SPYRM$RV5, order.by = SPYRM$DT)

id <- 1:nrow(RVSPY)
treino_id <- 1:ceiling(nrow(RVSPY)*0.6)
val_id <- c((min(id[-treino_id])-22):(min(id[-treino_id])-1), id[-treino_id])

treino <- RVSPY[treino_id, ]
val <- RVSPY[val_id, ]

val_matriz <- data.frame(y = val) %>% 
  mutate(`(Intercept)` = 1,
         RV1 = lag(y, 1),
         RV5 = rollmean(RV1, k = 5, align = "right", fill = NA),
         RV22 = rollmean(RV1, k = 22, align = "right", fill = NA)) %>% 
  na.omit() %>% 
  as.matrix()

fit <- HARmodel(data = treino, periods = c(1,5,22), RVest = c("rCov"),
                type = "HAR", h = 1, transform = NULL, inputType = "RM")

treino_matriz <- fit[["model"]] %>% 
  mutate(`(Intercept)` = 1) %>% 
  select(`(Intercept)`, RV1, RV5, RV22) %>% 
  as.matrix()

betas <- coef(fit)

treino_pred <- treino_matriz %*% betas
val_pred <- as.numeric(val_matriz[,-1] %*% betas)
val_gab <- as.numeric(RVSPY[-treino_id,])

cal_id <- 1:ceiling(length(val_pred)*0.5)

cal_gab <- val_gab[cal_id]
test_gab <- val_gab[-cal_id]

cal_pred <- val_pred[cal_id]
test_pred <- val_pred[-cal_id]

test_scp <- scp(cal = cal_pred, y = cal_gab, test = test_pred)

ic_test <- test_scp$ic_test
ic_test$gab <- test_gab

ic_test <- ic_test %>% 
  mutate(check = ifelse(gab <= up & gab >= down, T, F))

mean(ic_test$check)
mean(ic_test$amp)

# REG

n <- 10000

e <- rnorm(n, sd = 2)

x <- rnorm(n, 10, 4)

y <- 10 + 5 * x + e

data <- data.frame(y, x)
treino <- data[1:7000, ]
cal <- data[7001:8500, ]
test <- data[8501:10000, ]

fit <- lm(y ~ x, treino)
summary(fit)

cal_pred <- predict(fit, data.frame(x = cal$x))
test_pred <- predict(fit, data.frame(x = test$x))

cal_gab <- cal$y
test_gab <- test$y

test_scp <- scp(cal = cal_pred, y = cal_gab, test = test_pred)

ic_test <- test_scp$ic_test
ic_test$gab <- test_gab

ic_test <- ic_test %>% 
  mutate(check = ifelse(gab <= up & gab >= down, T, F))

mean(ic_test$check)
mean(ic_test$amp)

# GARCH

data(dmbp)

treino <- dmbp[1:1382, 1]
cal <- dmbp[1383:1678, 1]
test <- dmbp[1679:1974, 1]

spec <- ugarchspec()
fit <- ugarchfit(data = dmbp[,1], spec = spec, out.sample = 592)

val_fore <- ugarchforecast(fit, n.ahead = 1, n.roll = 591)

val_pred <- val_fore@forecast$sigmaFor %>% t() %>% as.numeric()

gab <- dmbp$V1^2

cal_pred <- val_pred[1:296]
test_pred <- val_pred[297:592]

cal_gab <- gab[1383:1678]
test_gab <- gab[1679:1974]

test_scp <- scp(cal = cal_pred, y = cal_gab, test = test_pred)

ic_test <- test_scp$ic_test
ic_test$gab <- test_gab

ic_test <- ic_test %>% 
  mutate(check = ifelse(gab <= up & gab >= down, T, F))

mean(ic_test$check)
mean(ic_test$amp)
