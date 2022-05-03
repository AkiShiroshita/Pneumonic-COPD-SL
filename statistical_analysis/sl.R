
### Caution!!!###
## Just for demonstration ##
## I wrote a "rough" code ##

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
             "survey",
             "broom",
             "cobalt",
             "WeightIt",
             "MatchIt",
             "SuperLearner",
             "glmnet",
             "xgboost",
             "earth",
             "drtmle",
             "mice")
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

# pick up the first dataset just for my simulations
df_m0 <- df %>% 
  mutate(ams = as.logical(ams),
         hot = as.logical(hot),
         wheeze = as.logical(wheeze),
         adl = as.logical(adl),
         death = as.logical(death),
         steroid = as.logical(steroid))
m <- 100 
imp <- mice(df_m0, m = m ,seed=12345 ,maxit=50 ,printFlag=FALSE)
df_m <- complete(imp,action="long") %>% 
  as_tibble()

df_m1 <- df_m %>% 
  filter(.imp == 1) %>% 
  select(-.imp, -.id) %>% 
  mutate_all(.funs = ~ as.numeric(.))

df_x <- df_m1 %>% 
  select(-death)

df_x0 <- df_x %>% 
  mutate(steroid = 0)
df_x1 <- df_x %>% 
  mutate(steroid = 1)

df_xx <- df_x %>% 
  select(-steroid)

## for complete case analysis

#df_c <- df %>%
#  drop_na()
#
#df_x <- df_c %>% 
#  select(-death)
#
#df_x0 <- df_x %>% 
#  mutate(steroid = 0)
#df_x1 <- df_x %>% 
#  mutate(steroid = 1)
#
#df_xx <- df_x %>% 
#  select(-steroid)

# G-computation & IPTW ----------------------------------------------------

## G-computation
fit_or <- glm(df_m1$death ~ .,
              data = df_x,
              family = binomial())
fit_or
plot(fit_or)

Qbar0W <- predict(fit_or,
                  newdata = df_x0,
                  type = "response")
Qbar0W <- if_else(Qbar0W > .5, 1, 0)

Qbar1W <- predict(fit_or,
                  newdata = df_x1,
                  type = "response")
Qbar1W <- if_else(Qbar1W > .5, 1, 0)

psi_nG0 <- mean(Qbar0W)
psi_nG1 <- mean(Qbar1W)

gamma_nG <- psi_nG1 - psi_nG0

M <- 500
gammaVec_nG <- replicate(M, {
  ind <- sample (1:nrow(df_x), replace = TRUE )
  fit_or <- glm(df_m1$death[ind] ~ .,
                data = df_x[ind ,],
                family = binomial())
  Qbar0W <- predict(fit_or, newdata = df_x0[ind ,])
  Qbar1W <- predict(fit_or, newdata = df_x1[ind ,])
  Qbar0W <- if_else(Qbar0W > .5, 1, 0)
  Qbar1W <- if_else(Qbar1W > .5, 1, 0)
  psi_nG0 <- mean(Qbar0W)
  psi_nG1 <- mean(Qbar1W)
  gamma_nG <- psi_nG1 - psi_nG0
  return (gamma_nG)
})

exp(quantile(gammaVec_nG, c(0.025, 0.975)))

## IPTW
fit_ps <- glm(df_x$steroid ~ .,
              data = df_xx,
              family = binomial())
fit_ps

g1W <- predict(fit_ps,
               newdata = df_xx,
               type = "response")

psi_nIPTW0 <- mean(as.numeric(df_x0$steroid) / (1 - g1W) * df_m1$death)
psi_nIPTW1 <- mean(as.numeric(df_x1$steroid) / g1W * df_m1$death)

gamma_nIPTW <- psi_nIPTW1 - psi_nIPTW0
gamma_nIPTW/100

M <- 500
gammaVec_ps <- replicate(M, {
  ind <- sample (1:nrow(df_x), replace = TRUE )
  fit_ps <- glm(df_x$steroid ~ .,
                data = df_xx[ind ,],
                family = binomial())
  g1W <- predict(fit_ps,
                 newdata = df_xx[ind ,],
                 type = "response")
  psi_nIPTW0 <- mean(as.numeric(df_x0$steroid) / (1 - g1W) * df_m1$death)
  psi_nIPTW1 <- mean(as.numeric(df_x1$steroid) / g1W * df_m1$death)
  gamma_nIPTW <- psi_nIPTW1 - psi_nIPTW0
  return (gamma_nIPTW)
})

quantile(gammaVec_ps, c(0.025, 0.975))/100

# Extra matching and weighting ----------------------------------------------------------------

match_var <- colnames(select(df_m1, -steroid, -death))
out_var <-  colnames(select(df_m1, -steroid, -death))
trt_form <- f.build("steroid", match_var)
out_form <- f.build("death", out_var)

df_m1 <- df_m1 %>% 
  mutate(steroid = case_when(steroid == 0 ~ "Non-steroid",
                             steroid == 1 ~ "Steroid")) 
set.cobalt.options(binary = "std",
                   un = TRUE,
                   disp.v.ratio = TRUE,
                   disp.ks = TRUE)

t.test(death ~ steroid, data = df_m1) %>% tidy

e.match <- matchit(trt_form,
                   data = df_m1,
                   method = "exact")
summary(e.match)
e.data <- match.data(e.match)

glm(death ~ steroid,
    data = e.data,
    family = binomial(),
    weights = weights) %>% 
  tidy() 

## entropy balancing

ebal.out <- weightit(trt_form,
                     method = "ebal",
                     moments = 3,
                     data = df_m1,
                     estimand = "ATT")
summary(ebal.out)
love.plot(ebal.out,
          stas= c("mean", "var", "ks"),
          var.order = "unadjusted",
          thresholds = c(.1, 2, .05),
          line = TRUE)

df_m1 <- df_m1 %>% 
  mutate(ebal_wt = get.w(ebal.out))

glm(death ~ steroid,
    data = df_m1,
    family = binomial(),
    weights = ebal_wt) %>% 
  tidy() 

# SL ----------------------------------------------------------------------

## Basic SuperLearner
sl_fit1 <- SuperLearner(Y = df_m1$death,
                        X = df_x,
                        SL.library = c("SL.glmnet", "SL.xgboost"),
                        family = binomial(),
                        method = "method.CC_nloglik",
                        cvControl = list(V=4),
                        verbose = FALSE)
sl_fit1

pred0 <- predict(sl_fit1, df_x0)
pred0 <- if_else(pred0$pred > .5, 1, 0)

pred1 <- predict(sl_fit1, df_x1)
pred1 <- if_else(pred1$pred > .5, 1, 0)

mpred0 <- mean(pred0)
mpred1 <- mean(pred1)

gamma_nG <- mpred1 - mpred0

sl_cv <- CV.SuperLearner(Y = df_m1$death,
                         X = df_x,
                         SL.library = c("SL.glmnet", "SL.xgboost", "SL.glm"),
                         family = binomial(),
                         method = "method.CC_nloglik",
                         cvControl = list(V = 2),
                         innerCvControl = list(list(V = 2),
                                            list(V = 2)))
plot(sl_cv)

## AIPTW

df_m1 <- df_m1 %>% 
  mutate_all(steroid = as.numeric(as.factor(steroid)))

fit_or <- SuperLearner(Y = df_m1$death,
                       X = df_x,
                       SL.library = c("SL.glmnet", "SL.xgboost"),
                       family = binomial(),
                       method = "method.CC_nloglik",
                       cvControl = list(V=4),
                       verbose = FALSE)
fit_ps <- SuperLearner(Y = df_m1$steroid,
                       X = df_xx,
                       SL.library = c("SL.glmnet", "SL.xgboost"),
                       family = binomial(),
                       method = "method.CC_nloglik",
                       cvControl = list(V=4),
                       verbose = FALSE)

Qbar0 <- predict(fit_or,
                 newdata = df_x0)
Qbar0 <- if_else(Qbar0$pred > .5, 1, 0)

Qbar1 <- predict(fit_or,
                 newdata = df_x1)
Qbar1 <- if_else(Qbar1$pred > .5, 1, 0)

g1W <- fit_ps$SL.predict

psi_nAIPTW0 <- mean(Qbar0) + mean(as.numeric(df_x0$steroid) / (1-g1W) * (df_m1$death - Qbar0))
psi_nAIPTW1 <- mean(Qbar1) + mean(as.numeric(df_x1$steroid) / (g1W) * (df_m1$death - Qbar1))

gamma_nAIPTW <- psi_nAIPTW1 - psi_nAIPTW0
gamma_nAIPTW

tau2_n0 <- mean((as.numeric(df_x0$steroid) / (1-g1W) * (df_m1$death - Qbar0) + Qbar0 - psi_nAIPTW0)^2)
ci0 <- c(psi_nAIPTW0 - 1.96 * sqrt (tau2_n0/nrow(df_m1)),
         psi_nAIPTW0 + 1.96 * sqrt (tau2_n0/nrow(df_m1)))
ci0

tau2_n1 <- mean((as.numeric(df_x1$steroid) / g1W * (df_m1$death - Qbar1) + Qbar1 - psi_nAIPTW1)^2)
ci1 <- c(psi_nAIPTW1 - 1.96 * sqrt (tau2_n1/nrow(df_m1)),
         psi_nAIPTW1 + 1.96 * sqrt (tau2_n1/nrow(df_m1)))
ci1

ci1 - ci0

## TMLE
## doubly robust inference

tmle <- drtmle(W = df_xx,
               A = df_x$steroid,
               Y = df_m1$death,
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
l <- min(df_m1$death); u <- max(df_m1$death)
Ystar <- (df_m1$death - l) / (u-l)
logistic_fit <- glm(Ystar ~ -1 + offset(qlogis(QbarA)) + Z0 + Z1 ,
                       family = binomial ())
alpha <- coef(logistic_fit)

Qbarstar0 <- (u-l)* plogis(qlogis(Qbar0) + alpha [1]/(1 - g1W)) + l
Qbarstar1 <- (u-l)* plogis(qlogis(Qbar1) + alpha [2]/ g1W) + l

tau2_n0 <- mean((as.numeric(df_x0$steroid) / (1-g1W) * (df_m1$death - Qbarstar0) + Qbarstar0 - psi_nAIPTW0)^2)
ci0 <- c(psi_nAIPTW0 - 1.96 * sqrt (tau2_n0/nrow(df_m1)),
         psi_nAIPTW0 + 1.96 * sqrt (tau2_n0/nrow(df_m1)))
tau2_n1 <- mean((as.numeric(df_x1$steroid) / g1W * (df_m1$death - Qbarstar1) + Qbarstar1 - psi_nAIPTW1)^2)
ci1 <- c(psi_nAIPTW1 - 1.96 * sqrt (tau2_n1/nrow(df_m1)),
         psi_nAIPTW1 + 1.96 * sqrt (tau2_n1/nrow(df_m1)))

ci1 - ci0
