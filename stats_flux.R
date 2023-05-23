#Load packages
library(tidyverse)

#Loading dataset
df.cml.last <- read.csv("cml_respired_lastcyc.csv", row.names = 1)
df.cml.last$treatment <- as.factor(df.cml.last$treatment)
df.cml.last$temperature <- as.factor(df.cml.last$temperature)

df.c.trt <- c("water-extractable PyOM-C","non-water-extractable PyOM-C", "SOC")
df.c.prim <- c("SOC", "SOC(no PyOM)")

#Ungrouping the two temp datasets
df.c.350 <- df.cml.last %>%
  filter(treatment %in% df.c.trt) %>%
  filter(temperature == 350)
df.c.550 <- df.cml.last %>%
  filter(treatment %in% df.c.trt) %>%
  filter(temperature == 550)

#Ungrouping temp and treatments
df.c.350.wx <- df.cml.last %>%
  filter(temperature == 350 & treatment == "water-extractable PyOM-C")
df.c.350.nwx <- df.cml.last %>%
  filter(temperature == 350 & treatment == "non-water-extractable PyOM-C")
df.c.350.s <- df.cml.last %>%
  filter(temperature == 350 & treatment == "SOC")

#Testing for normality
shapiro.test(df.c.350.s$C_mg_norm_cml)
# If p-value is greater than 0.05; null hypothesis accepted; variable is normally distributed
# Normalized C values within each temp-treatment are normally distrubuted 

#Testing for homogeneity of variance
#Bartlett's test
res <- bartlett.test(C_mg_norm_cml ~ treatment, data = df.c.550)
res 
#If p-value is greater than 0.05; null hypothesis is accepted; data are normally distributed
#Fails the homogeneity assumption for both 350 and 550 data sets

#Testing for homogeneity of variance using Levene's test 
library(car)
leveneTest(C_mg_cml ~ treatment, data = df.c.350)
#If p-value is greater than 0.05; null hypothesis is accepted; variances of the two groups are equal
#Fails the homogeneity assumption for both 350 and 550 data sets

#Performing Welch's ANOVa for groups with unequal variances
model <- oneway.test(C_mg_norm_cml ~ treatment, data = df.c.350, var.equal = FALSE)
model
#If p-value is less than 0.05, null hypothesis rejected, cumulative C different between the groups

#Performing Games-Howell post-hoc test; follow up test for Welch's ANOVA
library("rstatix")
s = games_howell_test(data = df.c.550, C_mg_norm_cml ~ treatment, conf.level = 0.95, detailed  = FALSE)
s
