
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
remotes::install_github("ecpolley/SuperLearner")
devtools::install_github("benkeser/drtmle")

getwd()
rm(list=ls())

# SL ----------------------------------------------------------------------

set.seed(1116)
df <- read_csv("input/data.csv")

df %>% glimpse()
df %>% colnames()

df_c <- df %>%
  drop_na()

df_x <- df_c %>% 
  select(2:4, 9:17, 25)

## Basic SuperLearner
sl_fit1 <- SuperLearner(Y = df_c$death,
                        X = df_x,
                        SL.library = c("SL.glmnet", "SL.xgboost"),
                        family = binomial(),
                        method = "method.CC_nloglik",
                        cvControl = list(V=4),
                        verbose = FALSE)
sl_fit1

## TMLE
df_xx <- df_x %>% 
  select(-steroid)

tmle <- drtmle(W = df_xx,
               A = df_x$steroid,
               Y = df_c$death,
               a_0 = c(0,1),
               family = binomial(),
               SL_g = c("SL.glmnet", "SL.xgboost"),
               SL_Q = c("SL.glmnet", "SL.xgboost"),
               SL_gr = c("SL.glmnet", "SL.xgboost"),
               SL_Qr = c("SL.glmnet", "SL.xgboost"),
               stratify = FALSE)
tmle
ci(tmle)
ci(tmle, contrast = c(-1, 1))
