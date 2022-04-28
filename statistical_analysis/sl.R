
# Data preparation --------------------------------------------------------

packages = c("devtools",
             "usethis",
             "here",
             "readr",
             "readxl",
             "tidyverse",
             "tidylog",
             "lubridate",
             "ggplot2",
             "tidylog",
             "ggplotgui",
             "ggthemes",
             "arsenal",
             "SuperLearner",
             "glmnet",
             "xgboost",
             "earth",
             "drtmle")
package.check <- lapply(packages, FUN = function(x){
  if (!require(x, character.only = TRUE)){
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})

# install.packages("remotes")
# remotes::install_github("ecpolley/SuperLearner")
# devtools::install_github("benkeser/drtmle")

getwd()
rm(list=ls())

set.seed(1116)
df <- read_csv("input/data.csv")

df %>% glimpse()
df %>% colnames()

df <- df %>% 
  mutate_all(.funs = ~ as.numeric(.)) %>% 
  mutate(eosi_a = wbc * eosi_p/100) %>% 
  select(eosi_a, bun, rr, ams, hr, hot, wheeze, adl, death, steroid) 

df_c <- df %>%
  drop_na()

df_x <- df_c %>% 
  select(-death)

df_x0 <- df_x %>% 
  mutate(steroid = 0)
df_x1 <- df_x %>% 
  mutate(steroid = 1)

df_xx <- df_x %>% 
  select(-steroid)

# G-computation & IPTW ----------------------------------------------------

## G-computation
fit_or <- glm(df_c$death ~ .,
              data = df_x,
              family = binomial())
fit_or

Qbar0W <- predict(fit_or,
                  newdata = df_x0,
                  type = "response")
Qbar0W <- if_else(Qbar0W > .5, 1, 0)

Qbar1W <- predict(fit_or,
                  newdata = df_x1,
                  type = "response")

psi_nG0 <- mean(Qbar0W)
psi_nG1 <- mean(Qbar1W)

gamma_nG <- psi_nG1 - psi_nG0
gamma_nG

M <- 500
gammaVec_nG <- replicate(M, {
  ind <- sample (1:nrow(df_x), replace = TRUE )
  fit_or <- glm(df_c$death[ind] ~ ., data = df_x[ind ,])
  Qbar0W <- predict(fit_or, newdata = df_x0[ind ,])
  Qbar1W <- predict(fit_or, newdata = df_x1[ind ,])
  psi_nG1 <- mean(Qbar1W)
  psi_nG0 <- mean(Qbar0W)
  gamma_nG <- psi_nG1 - psi_nG0
  return (gamma_nG)
})

quantile(gammaVec_nG, c(0.025, 0.975))

## IPTW
fit_ps <- glm(df_x$steroid ~ .,
              data = df_xx,
              family = binomial())
fit_ps

g1W <- predict(fit_ps,
               newdata = df_xx,
               type = "response")

psi_nIPTW0 <- mean(as.numeric(df_x0$steroid) / (1 - g1W) * df_x$steroid)
psi_nIPTW1 <- mean(as.numeric(df_x1$steroid) / g1W * df_x$steroid)

gamma_nIPTW <- psi_nIPTW1 - psi_nIPTW0
gamma_nIPTW

M <- 500
gammaVec_ps <- replicate(M, {
  ind <- sample (1:nrow(df_x), replace = TRUE )
  fit_ps <- glm(df_x$steroid ~ .,
                data = df_xx[ind ,],
                family = binomial())
  g1W <- predict(fit_ps,
                 newdata = df_xx[ind ,],
                 type = "response")
  psi_nIPTW0 <- mean(as.numeric(df_x0$steroid) / (1 - g1W) * df_x$steroid)
  psi_nIPTW1 <- mean(as.numeric(df_x1$steroid) / g1W * df_x$steroid)
  gamma_nIPTW <- psi_nIPTW1 - psi_nIPTW0
  return (gamma_nIPTW)
})

quantile(gammaVec_ps, c(0.025, 0.975))

# SL ----------------------------------------------------------------------

## Basic SuperLearner
sl_fit1 <- SuperLearner(Y = df_c$death,
                        X = df_x,
                        SL.library = c("SL.glmnet", "SL.xgboost"),
                        family = binomial(),
                        method = "method.CC_nloglik",
                        cvControl = list(V=4),
                        verbose = FALSE)
sl_fit1

pred0 <- predict(sl_fit1, df_x0)
pred1 <- predict(sl_fit1, df_x1)

pred0$pred
pred1$pred

sl_cv <- CV.SuperLearner(Y = df_c$death,
                        X = df_x,
                        SL.library = c("SL.glmnet", "SL.xgboost", "SL.glm"),
                        family = binomial(),
                        method = "method.CC_nloglik",
                        cvControl = list(V = 2),
                        innerCvControl = list(list(V = 2),
                                            list(V = 2)))
plot(sl_cv)

## AIPTW
fit_or <- SuperLearner(Y = df_c$death,
                       X = df_x,
                       SL.library = c("SL.glmnet", "SL.xgboost"),
                       family = binomial(),
                       method = "method.CC_nloglik",
                       cvControl = list(V=4),
                       verbose = FALSE)
fit_ps <- SuperLearner(Y = df_c$steroid,
                       X = df_xx,
                       SL.library = c("SL.glmnet", "SL.xgboost"),
                       family = binomial(),
                       method = "method.CC_nloglik",
                       cvControl = list(V=4),
                       verbose = FALSE)

Qbar0 <- predict(fit_or,
                 newdata = df_x0)
Qbar1 <- predict(fit_or,
                 newdata = df_x1)
g1W <- fit_ps$SL.predict

psi_nAIPTW0 <- mean(Qbar0$pred) + mean(as.numeric(df_x0$steroid) / (1-g1W) * (df_c$death - Qbar0$pred))
psi_nAIPTW1 <- mean(Qbar1$pred) + mean(as.numeric(df_x1$steroid) / (g1W) * (df_c$death - Qbar1$pred))

gamma_nAIPTW <- psi_nAIPTW1 - psi_nAIPTW0
gamma_nAIPTW

tau2_n0 <- mean((as.numeric(df_x0$steroid) / (1-g1W) * (df_c$death - Qbar0$pred) + Qbar0$pred - psi_nAIPTW0)^2)
ci0 <- c(psi_nAIPTW0 - 1.96 * sqrt (tau2_n0/nrow(df_c)),
         psi_nAIPTW0 + 1.96 * sqrt (tau2_n0/nrow(df_c)))
tau2_n1 <- mean((as.numeric(df_x1$steroid) / g1W * (df_c$death - Qbar1$pred) + Qbar1$pred - psi_nAIPTW1)^2)
ci1 <- c(psi_nAIPTW1 - 1.96 * sqrt (tau2_n1/nrow(df_c)),
         psi_nAIPTW1 + 1.96 * sqrt (tau2_n1/nrow(df_c)))

## TMLE

tmle <- drtmle(W = df_xx,
               A = df_x$steroid,
               Y = df_c$death,
               a_0 = c(0,1),
               family = binomial(),
               SL_g = c("SL.glm", "SL.earth"),
               SL_Q = c("SL.glm", "SL.earth"),
               SL_gr = c("SL.glm", "SL.earth"),
               SL_Qr = c("SL.glm", "SL.earth"),
               stratify = FALSE)
tmle
ci(tmle)
ci(tmle, contrast = c(-1, 1))

### manual approach

QbarA <- fit_or$SL.predict
Z1 <- df_x$steroid/g1W; Z0 <- (1-df_x$steroid)/(1 - g1W)
l <- min(df_c$death); u <- max(df_c$death)
Ystar <- (df_c$death - l) / (u-l)
logistic_fit <- glm(Ystar ~ -1 + offset(qlogis(QbarA)) + Z0 + Z1 ,
                       family = binomial ())
alpha <- coef(logistic_fit)

Qbarstar0 <- (u-l)* plogis(qlogis(Qbar0$pred) + alpha [1]/(1 - g1W)) + l
Qbarstar1 <- (u-l)* plogis(qlogis(Qbar1$pred) + alpha [2]/ g1W) + l

tau2_n0 <- mean((as.numeric(df_x0$steroid) / (1-g1W) * (df_c$death - Qbarstar0) + Qbarstar0 - psi_nAIPTW0)^2)
ci0 <- c(psi_nAIPTW0 - 1.96 * sqrt (tau2_n0/nrow(df_c)),
         psi_nAIPTW0 + 1.96 * sqrt (tau2_n0/nrow(df_c)))
tau2_n1 <- mean((as.numeric(df_x1$steroid) / g1W * (df_c$death - Qbarstar1) + Qbarstar1 - psi_nAIPTW1)^2)
ci1 <- c(psi_nAIPTW1 - 1.96 * sqrt (tau2_n1/nrow(df_c)),
         psi_nAIPTW1 + 1.96 * sqrt (tau2_n1/nrow(df_c)))

ci1 - ci0
