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

scores <- abs(cal_gab - cal_pred)/sqrt(cal_pred)

n_cal <- length(scores)

q_cal <- quantile(scores, probs = ceiling((1-0.05) * (n_cal+1))/n_cal)

test_scp <- data.frame(test = test_pred, gab = test_gab, 
                       down = test_pred - q_cal * sqrt(test_pred), up = test_pred + q_cal * sqrt(test_pred))

test_scp$amp <- test_scp$up - test_scp$down

ic_test <- test_scp
ic_test$gab <- test_gab

ic_test <- test_scp %>% 
  mutate(check = ifelse(gab <= up & gab >= down, T, F))

ic_test %>%
  ggplot(aes(x = seq(1, nrow(ic_test)))) +
  geom_line(aes(y = gab)) +
  geom_line(aes(y = down), color = "red") +
  geom_line(aes(y = up), color = "blue") +
  labs(title = "HAR: RVSPY data", x = "Index", y = "Realized volatility") +
  theme_minimal()

mean(ic_test$check)
mean(ic_test$amp)

test_scp <- scp(cal = cal_pred, y = cal_gab, test = test_pred)

ic_test <- test_scp$ic_test
ic_test$gab <- test_gab

ic_test <- ic_test %>% 
  mutate(check = ifelse(gab <= up & gab >= down, T, F))

ic_test %>%
  ggplot(aes(x = seq(1, nrow(ic_test)))) +
  geom_line(aes(y = gab)) +
  geom_line(aes(y = down), color = "red") +
  geom_line(aes(y = up), color = "blue") +
  labs(title = "HAR: RVSPY data", x = "Index", y = "Realized volatility") +
  theme_minimal()

mean(ic_test$check)
mean(ic_test$amp)

# REG

n <- 10000

e <- rnorm(n, sd = 2)

x <- rnorm(n, 10, 4)

y <- 10 + 5 * x + e

data <- data.frame(y, x)
treino <- data[1:round(n * 0.7), ]
cal <- data[(round(n * 0.7)+1):(round(n * 0.7) + round(n * 0.3/2)), ]
test <- data[(round(n * 0.7) + round(n * 0.3/2)+1):n, ]

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

ic_test %>% 
  mutate(x = x[8501:10000]) %>% 
  ggplot(aes(x = x)) +
  geom_line(aes(y = gab)) +
  geom_line(aes(y = down), color = "red") +
  geom_line(aes(y = up), color = "blue") +
  labs(title = "Regression: simulated data", x = "x", y = "y") +
  theme_minimal()

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
val_pred <- val_pred^2

gab <- dmbp$V1^2

cal_pred <- val_pred[1:296]
test_pred <- val_pred[297:592]

cal_gab <- gab[1383:1678]
test_gab <- gab[1679:1974]

scores <- abs(cal_gab - cal_pred)/sqrt(cal_pred)

n_cal <- length(scores)

q_cal <- quantile(scores, probs = ceiling((1-0.05) * (n_cal+1))/n_cal)

test_scp <- data.frame(test = test_pred, gab = test_gab, 
                       down = test_pred - q_cal * sqrt(test_pred), up = test_pred + q_cal * sqrt(test_pred))

test_scp$amp <- test_scp$up - test_scp$down

ic_test <- test_scp
ic_test$gab <- test_gab

ic_test <- test_scp %>% 
  mutate(check = ifelse(gab <= up & gab >= down, T, F))

ic_test %>%
  ggplot(aes(x = seq(1, nrow(ic_test)))) +
  geom_line(aes(y = gab)) +
  geom_line(aes(y = down), color = "red") +
  geom_line(aes(y = up), color = "blue") +
  labs(title = "GARCH: dmbp data", x = "Index", y = "Squared returns") +
  theme_minimal()

mean(ic_test$check)
mean(ic_test$amp)

# GARCH SIM

library(fGarch)

spec <- garchSpec(model = list())
sim <- garchSim(spec, n = 10000)
