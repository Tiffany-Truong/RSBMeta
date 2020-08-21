# Clean environment
rm(list = ls())
# Clear Console
cat("\014")

# Load library
library(metafor)
library(plotly)
library(robumeta)

#setwd("~/git/RSBMeta")
filePath <- "~/git/RSBMeta/Analyses/20-08-2020/output.txt"

# Redirect console to a log file
outputFile <- file(filePath, "w")
sink(outputFile, append=FALSE, type= c("output", "message"), split = TRUE)

#### Description ####
# rsbtype
# 1 = unprotected sex
# 2 = multiple life partners 
# 3 = hazardous sex
# 4 = sexual initiation 
# id = Study Number
# article = Author(s) and year
# ss = total study sample size 
# age.m = age mean 
# age.sd = age standard deviation 
# female = percentage of sample female 
# caucasian = percentage of sample caucasian 
# r = correlation coefficient 
# var = variance
# n = sample size utilized 
# Z = Fisher's Z transformation 
# Z.var = Fisher's Z variance transformation 
# back.transform = Fisher's Z back to r 


######################## 1 - UNPROTECTED SEX ######################## 

#### Data Prep ####

# pull data 
fullset = read.csv("/Users/tiffanytruong/git/RSBMeta/Analyses/20-08-2020/dataset.csv", header = T, sep = ",")


# set working directoy to subfolder - holds plots 
setwd("/Users/tiffanytruong/git/RSBMeta/Analyses/20-08-2020/1_UnprotectedSex")

# subset RSB
ds = subset(fullset, rsbtype == 1)

# Fisher's Z transformation
ds$Z <- .5*log((1 + ds$r) / (1 - ds$r))

# Variance of Fisher's Z
ds$Z.var <- 1 / (ds$n - 3)

####  Forest Plot ####

# Create and print default forest plot
# yi: effect size
# vi: effect size variance
# data: data set 

forestplot_ds <- rma.uni(yi = Z, vi = Z.var, data = ds)
forest(forestplot_ds, slab = ds$article, order = order(ds$Z), showweights = TRUE)


##################### Multivariate/Dependency Analyses ##################### 

#### RE model ####
# RE model in robumeta with correlational weights
model_rve_c <- robu(Z ~ 1, 
                    data = ds,
                    studynum = id,
                    var.eff.size = Z.var, 
                    modelweights = "CORR",
                    rho = .80)
model_rve_c


# Checking model sensitivity of correlational weights (i.e., choice of rho)
sensitivity(model_rve_c)


# forest plot
jpeg('Forest_Plot.jpg', width = 610, height = 880)
forest.robu(model_rve_c,
            es.lab = "rsbtype",
            study.lab = "id")
dev.off()

#### multivariate moderator ####

# Conditional RE model in robumeta with one moderator using correlational weights 
# gender

model_rve_c_female <- robu(Z ~ female, 
                           data = ds,
                           studynum = id,
                           var.eff.size = Z.var, 
                           modelweights = "CORR",
                           rho = .80)
model_rve_c_female

# Conditional RE model in robumeta with one moderator using correlational weights 
# race

model_rve_c_caucasian <- robu(Z ~ caucasian, 
                              data = ds,
                              studynum = id,
                              var.eff.size = Z.var, 
                              modelweights = "CORR",
                              rho = .80)
model_rve_c_caucasian

# Conditional RE model in robumeta with one moderator using correlational weights 
# age

model_rve_c_age <- robu(Z ~ age.m, 
                        data = ds,
                        studynum = id,
                        var.eff.size = Z.var, 
                        modelweights = "CORR",
                        rho = .80)
model_rve_c_age


##################### INDEPENDENT RE MODEL ##################### 

#### RE Model Full Model #### 
RE.model <- rma.uni(yi = Z, vi = Z.var, data = ds, method = "REML")
RE.model

#### back transformation ####
RE_r <- (exp(2*RE.model$beta) - 1) / (exp(2*RE.model$beta) + 1)
RE_r

#### Sensitivity Analysis ####
# Print various influence point measures
options(max.print=1000000)
influence(RE.model)

# Plots of influence measures
plot(influence(RE.model))

leave1out(RE.model)

#### Publication Bias ####

#Funnel Plots 

par(mfrow = c(1,1))

# Funnel plot w/ SE as vertical axis
funnel(RE.model)

# Egger's Regression Test
regtest(x = ds$Z, vi = ds$Z.var, model = "lm")

# Basic trim-and-fill code
trimfill(RE.model, side = "right", verbose = T)
trimfill(RE.model, side = "left", verbose = T)


# Plot both versions of trim-and-fill (optional)
par(mfrow = c(2,1))
funnel(trimfill(RE.model, side = "right"), main = "Right-Side Imputation")
funnel(trimfill(RE.model, side = "left"), main = "Left-Side Imputation")

# Fail safe N 
fsn(yi = Z, vi = Z.var, data = ds)

# #### Gender Moderator ####
# # % female 
# 
# # Scatterplot without weights
# plot_ly(ds, x = ~female, y = ~Z, text = ~id, type = 'scatter')
# 
# # Mixed-effects meta-regression
# gender_MM_RE <- rma(yi = Z, vi = Z.var, mods = ~ female, data = ds)
# gender_MM_RE
# 
# #### Caucasian Moderator ####
# # % Caucasian 
# 
# # Scatterplot without weights
# plot_ly(ds, x = ~caucasian, y = ~Z, text = ~id, type = 'scatter')
# 
# # Mixed-effects meta-regression
# caucasian_MM_RE <- rma(yi = Z, vi = Z.var, mods = ~ caucasian, data = ds)
# caucasian_MM_RE
# 
# #### Age Moderator ####
# 
# # Scatterplot without weights
# plot_ly(ds, x = ~age.m, y = ~Z, text = ~id, type = 'scatter')
# 
# # Mixed-effects meta-regression
# age_MM_RE <- rma(yi = Z, vi = Z.var, mods = ~ age.m, data = ds)
# age_MM_RE

#### EXTRA: Back Transform ####
ds$back.transform <- (exp(2*ds$Z) - 1) / (exp(2*ds$Z) + 1)



#### Clean evironment #### 
rm(list = ls())



######################## 2 - MULTIPLE LIFETIME PARTNERS ######################## 

#### Data Prep ####

# pull data 
fullset = read.csv("/Users/tiffanytruong/git/RSBMeta/Analyses/20-08-2020/dataset.csv", header = T, sep = ",")


# set working directoy to subfolder - holds plots 
setwd("/Users/tiffanytruong/git/RSBMeta/Analyses/20-08-2020/2_MultipleLTPartners")

# subset RSB
ds = subset(fullset, rsbtype == 2)

# Fisher's Z transformation
ds$Z <- .5*log((1 + ds$r) / (1 - ds$r))

# Variance of Fisher's Z
ds$Z.var <- 1 / (ds$n - 3)

####  Forest Plot ####

# Create and print default forest plot
# yi: effect size
# vi: effect size variance
# data: data set 

forestplot_ds <- rma.uni(yi = Z, vi = Z.var, data = ds)
forest(forestplot_ds, slab = ds$article, order = order(ds$Z), showweights = TRUE)

##################### Multivariate/Dependency Analyses ##################### 

#### RE model ####
# RE model in robumeta with correlational weights
model_rve_c <- robu(Z ~ 1, 
                    data = ds,
                    studynum = id,
                    var.eff.size = Z.var, 
                    modelweights = "CORR",
                    rho = .80)
model_rve_c


# Checking model sensitivity of correlational weights (i.e., choice of rho)
sensitivity(model_rve_c)


# forest plot
jpeg('Forest_Plot.jpg', width = 610, height = 880)
forest.robu(model_rve_c,
            es.lab = "rsbtype",
            study.lab = "id")
dev.off()

#### multivariate moderator ####

# Conditional RE model in robumeta with one moderator using correlational weights 
# gender

model_rve_c_female <- robu(Z ~ female, 
                           data = ds,
                           studynum = id,
                           var.eff.size = Z.var, 
                           modelweights = "CORR",
                           rho = .80)
model_rve_c_female

# Conditional RE model in robumeta with one moderator using correlational weights 
# race

model_rve_c_caucasian <- robu(Z ~ caucasian, 
                              data = ds,
                              studynum = id,
                              var.eff.size = Z.var, 
                              modelweights = "CORR",
                              rho = .80)
model_rve_c_caucasian

# Conditional RE model in robumeta with one moderator using correlational weights 
# age

model_rve_c_age <- robu(Z ~ age.m, 
                        data = ds,
                        studynum = id,
                        var.eff.size = Z.var, 
                        modelweights = "CORR",
                        rho = .80)
model_rve_c_age


##################### INDEPENDENT RE MODEL ##################### 

#### RE Model Full Model #### 
RE.model <- rma.uni(yi = Z, vi = Z.var, data = ds, method = "REML")
RE.model

#### back transformation ####
RE_r <- (exp(2*RE.model$beta) - 1) / (exp(2*RE.model$beta) + 1)
RE_r

#### Sensitivity Analysis ####
# Print various influence point measures
options(max.print=1000000)
influence(RE.model)

# Plots of influence measures
plot(influence(RE.model))

leave1out(RE.model)

#### Publication Bias ####

#Funnel Plots 

par(mfrow = c(1,1))

# Funnel plot w/ SE as vertical axis
funnel(RE.model)

# Egger's Regression Test
regtest(x = ds$Z, vi = ds$Z.var, model = "lm")

# Basic trim-and-fill code
trimfill(RE.model, side = "right", verbose = T)
trimfill(RE.model, side = "left", verbose = T)


# Plot both versions of trim-and-fill (optional)
par(mfrow = c(2,1))
funnel(trimfill(RE.model, side = "right"), main = "Right-Side Imputation")
funnel(trimfill(RE.model, side = "left"), main = "Left-Side Imputation")

# Fail safe N 
fsn(yi = Z, vi = Z.var, data = ds)

# #### Gender Moderator ####
# # % female 
# 
# # Scatterplot without weights
# plot_ly(ds, x = ~female, y = ~Z, text = ~id, type = 'scatter')
# 
# # Mixed-effects meta-regression
# gender_MM_RE <- rma(yi = Z, vi = Z.var, mods = ~ female, data = ds)
# gender_MM_RE
# 
# #### Caucasian Moderator ####
# # % Caucasian 
# 
# # Scatterplot without weights
# plot_ly(ds, x = ~caucasian, y = ~Z, text = ~id, type = 'scatter')
# 
# # Mixed-effects meta-regression
# caucasian_MM_RE <- rma(yi = Z, vi = Z.var, mods = ~ caucasian, data = ds)
# caucasian_MM_RE
# 
# #### Age Moderator ####
# 
# # Scatterplot without weights
# plot_ly(ds, x = ~age.m, y = ~Z, text = ~id, type = 'scatter')
# 
# # Mixed-effects meta-regression
# age_MM_RE <- rma(yi = Z, vi = Z.var, mods = ~ age.m, data = ds)
# age_MM_RE

#### EXTRA: Back Transform ####
ds$back.transform <- (exp(2*ds$Z) - 1) / (exp(2*ds$Z) + 1)



#### Clean evironment #### 
rm(list = ls())



######################## 3 - HAZARDOUS SEX ######################## 

#### Data Prep ####

# pull data 
fullset = read.csv("/Users/tiffanytruong/git/RSBMeta/Analyses/20-08-2020/dataset.csv", header = T, sep = ",")


# set working directoy to subfolder - holds plots 
setwd("/Users/tiffanytruong/git/RSBMeta/Analyses/20-08-2020/3_HazardousSex")

# subset RSB
ds = subset(fullset, rsbtype == 3)

# Fisher's Z transformation
ds$Z <- .5*log((1 + ds$r) / (1 - ds$r))

# Variance of Fisher's Z
ds$Z.var <- 1 / (ds$n - 3)

####  Forest Plot ####

# Create and print default forest plot
# yi: effect size
# vi: effect size variance
# data: data set 

forestplot_ds <- rma.uni(yi = Z, vi = Z.var, data = ds)
forest(forestplot_ds, slab = ds$article, order = order(ds$Z), showweights = TRUE)

##################### Multivariate/Dependency Analyses ##################### 

#### RE model ####
# RE model in robumeta with correlational weights
model_rve_c <- robu(Z ~ 1, 
                    data = ds,
                    studynum = id,
                    var.eff.size = Z.var, 
                    modelweights = "CORR",
                    rho = .80)
model_rve_c


# Checking model sensitivity of correlational weights (i.e., choice of rho)
sensitivity(model_rve_c)


# forest plot
jpeg('Forest_Plot.jpg', width = 610, height = 880)
forest.robu(model_rve_c,
            es.lab = "rsbtype",
            study.lab = "id")
dev.off()

#### multivariate moderator ####

# Conditional RE model in robumeta with one moderator using correlational weights 
# gender

model_rve_c_female <- robu(Z ~ female, 
                           data = ds,
                           studynum = id,
                           var.eff.size = Z.var, 
                           modelweights = "CORR",
                           rho = .80)
model_rve_c_female

# Conditional RE model in robumeta with one moderator using correlational weights 
# race

model_rve_c_caucasian <- robu(Z ~ caucasian, 
                              data = ds,
                              studynum = id,
                              var.eff.size = Z.var, 
                              modelweights = "CORR",
                              rho = .80)
model_rve_c_caucasian

# Conditional RE model in robumeta with one moderator using correlational weights 
# age

model_rve_c_age <- robu(Z ~ age.m, 
                        data = ds,
                        studynum = id,
                        var.eff.size = Z.var, 
                        modelweights = "CORR",
                        rho = .80)
model_rve_c_age


##################### INDEPENDENT RE MODEL ##################### 


#### RE Model Full Model #### 
RE.model <- rma.uni(yi = Z, vi = Z.var, data = ds, method = "REML")
RE.model

#### back transformation ####
RE_r <- (exp(2*RE.model$beta) - 1) / (exp(2*RE.model$beta) + 1)
RE_r

#### Sensitivity Analysis ####
# Print various influence point measures
options(max.print=1000000)
influence(RE.model)

# Plots of influence measures
plot(influence(RE.model))

leave1out(RE.model)

#### Publication Bias ####

#Funnel Plots 

par(mfrow = c(1,1))

# Funnel plot w/ SE as vertical axis
funnel(RE.model)

# Egger's Regression Test
regtest(x = ds$Z, vi = ds$Z.var, model = "lm")

# Basic trim-and-fill code
trimfill(RE.model, side = "right", verbose = T)
trimfill(RE.model, side = "left", verbose = T)


# Plot both versions of trim-and-fill (optional)
par(mfrow = c(2,1))
funnel(trimfill(RE.model, side = "right"), main = "Right-Side Imputation")
funnel(trimfill(RE.model, side = "left"), main = "Left-Side Imputation")

# Fail safe N 
fsn(yi = Z, vi = Z.var, data = ds)

# #### Gender Moderator ####
# # % female 
# 
# # Scatterplot without weights
# plot_ly(ds, x = ~female, y = ~Z, text = ~id, type = 'scatter')
# 
# # Mixed-effects meta-regression
# gender_MM_RE <- rma(yi = Z, vi = Z.var, mods = ~ female, data = ds)
# gender_MM_RE
# 
# #### Caucasian Moderator ####
# # % Caucasian 
# 
# # Scatterplot without weights
# plot_ly(ds, x = ~caucasian, y = ~Z, text = ~id, type = 'scatter')
# 
# # Mixed-effects meta-regression
# caucasian_MM_RE <- rma(yi = Z, vi = Z.var, mods = ~ caucasian, data = ds)
# caucasian_MM_RE
# 
# #### Age Moderator ####
# 
# # Scatterplot without weights
# plot_ly(ds, x = ~age.m, y = ~Z, text = ~id, type = 'scatter')
# 
# # Mixed-effects meta-regression
# age_MM_RE <- rma(yi = Z, vi = Z.var, mods = ~ age.m, data = ds)
# age_MM_RE

#### EXTRA: Back Transform ####
ds$back.transform <- (exp(2*ds$Z) - 1) / (exp(2*ds$Z) + 1)



#### Clean evironment #### 
rm(list = ls())



######################## 4 - SEX INITIATION ######################## 

#### Data Prep ####

# pull data 
fullset = read.csv("/Users/tiffanytruong/git/RSBMeta/Analyses/20-08-2020/dataset.csv", header = T, sep = ",")


# set working directoy to subfolder - holds plots 
setwd("/Users/tiffanytruong/git/RSBMeta/Analyses/20-08-2020/4_SexInitiation")

# subset RSB
ds = subset(fullset, rsbtype == 4)

# Fisher's Z transformation
ds$Z <- .5*log((1 + ds$r) / (1 - ds$r))

# Variance of Fisher's Z
ds$Z.var <- 1 / (ds$n - 3)

####  Forest Plot ####

# Create and print default forest plot
# yi: effect size
# vi: effect size variance
# data: data set 

forestplot_ds <- rma.uni(yi = Z, vi = Z.var, data = ds)
forest(forestplot_ds, slab = ds$article, order = order(ds$Z), showweights = TRUE)

##################### Multivariate/Dependency Analyses ##################### 

#### RE model ####
# RE model in robumeta with correlational weights
model_rve_c <- robu(Z ~ 1, 
                    data = ds,
                    studynum = id,
                    var.eff.size = Z.var, 
                    modelweights = "CORR",
                    rho = .80)
model_rve_c


# Checking model sensitivity of correlational weights (i.e., choice of rho)
sensitivity(model_rve_c)


# forest plot
jpeg('Forest_Plot.jpg', width = 610, height = 880)
forest.robu(model_rve_c,
            es.lab = "rsbtype",
            study.lab = "id")
dev.off()

#### multivariate moderator ####

# Conditional RE model in robumeta with one moderator using correlational weights 
# gender

model_rve_c_female <- robu(Z ~ female, 
                           data = ds,
                           studynum = id,
                           var.eff.size = Z.var, 
                           modelweights = "CORR",
                           rho = .80)
model_rve_c_female

# Conditional RE model in robumeta with one moderator using correlational weights 
# race

model_rve_c_caucasian <- robu(Z ~ caucasian, 
                              data = ds,
                              studynum = id,
                              var.eff.size = Z.var, 
                              modelweights = "CORR",
                              rho = .80)
model_rve_c_caucasian

# Conditional RE model in robumeta with one moderator using correlational weights 
# age

model_rve_c_age <- robu(Z ~ age.m, 
                        data = ds,
                        studynum = id,
                        var.eff.size = Z.var, 
                        modelweights = "CORR",
                        rho = .80)
model_rve_c_age


##################### INDEPENDENT RE MODEL ##################### 



#### RE Model Full Model #### 
RE.model <- rma.uni(yi = Z, vi = Z.var, data = ds, method = "REML")
RE.model

#### back transformation ####
RE_r <- (exp(2*RE.model$beta) - 1) / (exp(2*RE.model$beta) + 1)
RE_r

#### Sensitivity Analysis ####
# Print various influence point measures
options(max.print=1000000)
influence(RE.model)

# Plots of influence measures
plot(influence(RE.model))

leave1out(RE.model)

#### Publication Bias ####

#Funnel Plots 

par(mfrow = c(1,1))

# Funnel plot w/ SE as vertical axis
funnel(RE.model)

# Egger's Regression Test
regtest(x = ds$Z, vi = ds$Z.var, model = "lm")

# Basic trim-and-fill code
trimfill(RE.model, side = "right", verbose = T)
trimfill(RE.model, side = "left", verbose = T)


# Plot both versions of trim-and-fill (optional)
par(mfrow = c(2,1))
funnel(trimfill(RE.model, side = "right"), main = "Right-Side Imputation")
funnel(trimfill(RE.model, side = "left"), main = "Left-Side Imputation")

# Fail safe N 
fsn(yi = Z, vi = Z.var, data = ds)

# #### Gender Moderator ####
# # % female 
# 
# # Scatterplot without weights
# plot_ly(ds, x = ~female, y = ~Z, text = ~id, type = 'scatter')
# 
# # Mixed-effects meta-regression
# gender_MM_RE <- rma(yi = Z, vi = Z.var, mods = ~ female, data = ds)
# gender_MM_RE
# 
# #### Caucasian Moderator ####
# # % Caucasian 
# 
# # Scatterplot without weights
# plot_ly(ds, x = ~caucasian, y = ~Z, text = ~id, type = 'scatter')
# 
# # Mixed-effects meta-regression
# caucasian_MM_RE <- rma(yi = Z, vi = Z.var, mods = ~ caucasian, data = ds)
# caucasian_MM_RE
# 
# #### Age Moderator ####
# 
# # Scatterplot without weights
# plot_ly(ds, x = ~age.m, y = ~Z, text = ~id, type = 'scatter')
# 
# # Mixed-effects meta-regression
# age_MM_RE <- rma(yi = Z, vi = Z.var, mods = ~ age.m, data = ds)
# age_MM_RE

#### EXTRA: Back Transform ####
ds$back.transform <- (exp(2*ds$Z) - 1) / (exp(2*ds$Z) + 1)



#### Clean evironment #### 
rm(list = ls())



######################## 5 - VIRGINITY STATUS ######################## 

#### Data Prep ####

# pull data 
fullset = read.csv("/Users/tiffanytruong/git/RSBMeta/Analyses/20-08-2020/dataset.csv", header = T, sep = ",")


# set working directoy to subfolder - holds plots 
setwd("/Users/tiffanytruong/git/RSBMeta/Analyses/20-08-2020/5_VirginityStatus")

# subset RSB
ds = subset(fullset, rsbtype == 5)

# Fisher's Z transformation
ds$Z <- .5*log((1 + ds$r) / (1 - ds$r))

# Variance of Fisher's Z
ds$Z.var <- 1 / (ds$n - 3)

####  Forest Plot ####

# Create and print default forest plot
# yi: effect size
# vi: effect size variance
# data: data set 

forestplot_ds <- rma.uni(yi = Z, vi = Z.var, data = ds)
forest(forestplot_ds, slab = ds$article, order = order(ds$Z), showweights = TRUE)

##################### Multivariate/Dependency Analyses ##################### 

#### RE model ####
# RE model in robumeta with correlational weights
model_rve_c <- robu(Z ~ 1, 
                    data = ds,
                    studynum = id,
                    var.eff.size = Z.var, 
                    modelweights = "CORR",
                    rho = .80)
model_rve_c


# Checking model sensitivity of correlational weights (i.e., choice of rho)
sensitivity(model_rve_c)


# forest plot
jpeg('Forest_Plot.jpg', width = 610, height = 880)
forest.robu(model_rve_c,
            es.lab = "rsbtype",
            study.lab = "id")
dev.off()

#### multivariate moderator ####

# Conditional RE model in robumeta with one moderator using correlational weights 
# gender

model_rve_c_female <- robu(Z ~ female, 
                           data = ds,
                           studynum = id,
                           var.eff.size = Z.var, 
                           modelweights = "CORR",
                           rho = .80)
model_rve_c_female

# Conditional RE model in robumeta with one moderator using correlational weights 
# race

model_rve_c_caucasian <- robu(Z ~ caucasian, 
                              data = ds,
                              studynum = id,
                              var.eff.size = Z.var, 
                              modelweights = "CORR",
                              rho = .80)
model_rve_c_caucasian

# Conditional RE model in robumeta with one moderator using correlational weights 
# age

model_rve_c_age <- robu(Z ~ age.m, 
                        data = ds,
                        studynum = id,
                        var.eff.size = Z.var, 
                        modelweights = "CORR",
                        rho = .80)
model_rve_c_age


##################### INDEPENDENT RE MODEL ##################### 




#### RE Model Full Model #### 
RE.model <- rma.uni(yi = Z, vi = Z.var, data = ds, method = "REML")
RE.model

#### back transformation ####
RE_r <- (exp(2*RE.model$beta) - 1) / (exp(2*RE.model$beta) + 1)
RE_r

#### Sensitivity Analysis ####
# Print various influence point measures
options(max.print=1000000)
influence(RE.model)

# Plots of influence measures
plot(influence(RE.model))

leave1out(RE.model)

#### Publication Bias ####

#Funnel Plots 

par(mfrow = c(1,1))

# Funnel plot w/ SE as vertical axis
funnel(RE.model)

# Egger's Regression Test
regtest(x = ds$Z, vi = ds$Z.var, model = "lm")

# Basic trim-and-fill code
trimfill(RE.model, side = "right", verbose = T)
trimfill(RE.model, side = "left", verbose = T)


# Plot both versions of trim-and-fill (optional)
par(mfrow = c(2,1))
funnel(trimfill(RE.model, side = "right"), main = "Right-Side Imputation")
funnel(trimfill(RE.model, side = "left"), main = "Left-Side Imputation")

# Fail safe N 
fsn(yi = Z, vi = Z.var, data = ds)

# #### Gender Moderator ####
# # % female 
# 
# # Scatterplot without weights
# plot_ly(ds, x = ~female, y = ~Z, text = ~id, type = 'scatter')
# 
# # Mixed-effects meta-regression
# gender_MM_RE <- rma(yi = Z, vi = Z.var, mods = ~ female, data = ds)
# gender_MM_RE
# 
# #### Caucasian Moderator ####
# # % Caucasian 
# 
# # Scatterplot without weights
# plot_ly(ds, x = ~caucasian, y = ~Z, text = ~id, type = 'scatter')
# 
# # Mixed-effects meta-regression
# caucasian_MM_RE <- rma(yi = Z, vi = Z.var, mods = ~ caucasian, data = ds)
# caucasian_MM_RE
# 
# #### Age Moderator ####
# 
# # Scatterplot without weights
# plot_ly(ds, x = ~age.m, y = ~Z, text = ~id, type = 'scatter')
# 
# # Mixed-effects meta-regression
# age_MM_RE <- rma(yi = Z, vi = Z.var, mods = ~ age.m, data = ds)
# age_MM_RE

#### EXTRA: Back Transform ####
ds$back.transform <- (exp(2*ds$Z) - 1) / (exp(2*ds$Z) + 1)



#### Clean evironment #### 
rm(list = ls())



######################## 6 - STD HISTORY ######################## 

#### Data Prep ####

# pull data 
fullset = read.csv("/Users/tiffanytruong/git/RSBMeta/Analyses/20-08-2020/dataset.csv", header = T, sep = ",")


# set working directoy to subfolder - holds plots 
setwd("/Users/tiffanytruong/git/RSBMeta/Analyses/20-08-2020/6_STDHistory")

# subset RSB
ds = subset(fullset, rsbtype == 6)

# Fisher's Z transformation
ds$Z <- .5*log((1 + ds$r) / (1 - ds$r))

# Variance of Fisher's Z
ds$Z.var <- 1 / (ds$n - 3)

####  Forest Plot ####

# Create and print default forest plot
# yi: effect size
# vi: effect size variance
# data: data set 

forestplot_ds <- rma.uni(yi = Z, vi = Z.var, data = ds)
forest(forestplot_ds, slab = ds$article, order = order(ds$Z), showweights = TRUE)


##################### INDEPENDENT RE MODEL ##################### 





#### RE Model Full Model #### 
RE.model <- rma.uni(yi = Z, vi = Z.var, data = ds, method = "REML")
RE.model

#### back transformation ####
RE_r <- (exp(2*RE.model$beta) - 1) / (exp(2*RE.model$beta) + 1)
RE_r

#### Sensitivity Analysis ####
# Print various influence point measures
options(max.print=1000000)
influence(RE.model)

# Plots of influence measures
plot(influence(RE.model))

leave1out(RE.model)

#### Publication Bias ####

#Funnel Plots 

par(mfrow = c(1,1))

# Funnel plot w/ SE as vertical axis
funnel(RE.model)

# Egger's Regression Test
regtest(x = ds$Z, vi = ds$Z.var, model = "lm")

# Basic trim-and-fill code
trimfill(RE.model, side = "right", verbose = T)
trimfill(RE.model, side = "left", verbose = T)


# Plot both versions of trim-and-fill (optional)
par(mfrow = c(2,1))
funnel(trimfill(RE.model, side = "right"), main = "Right-Side Imputation")
funnel(trimfill(RE.model, side = "left"), main = "Left-Side Imputation")

# Fail safe N 
fsn(yi = Z, vi = Z.var, data = ds)

# #### Gender Moderator ####
# # % female 
# 
# # Scatterplot without weights
# plot_ly(ds, x = ~female, y = ~Z, text = ~id, type = 'scatter')
# 
# # Mixed-effects meta-regression
# gender_MM_RE <- rma(yi = Z, vi = Z.var, mods = ~ female, data = ds)
# gender_MM_RE
# 
# #### Caucasian Moderator ####
# # % Caucasian 
# 
# # Scatterplot without weights
# plot_ly(ds, x = ~caucasian, y = ~Z, text = ~id, type = 'scatter')
# 
# # Mixed-effects meta-regression
# caucasian_MM_RE <- rma(yi = Z, vi = Z.var, mods = ~ caucasian, data = ds)
# caucasian_MM_RE
# 
# #### Age Moderator ####
# 
# # Scatterplot without weights
# plot_ly(ds, x = ~age.m, y = ~Z, text = ~id, type = 'scatter')
# 
# # Mixed-effects meta-regression
# age_MM_RE <- rma(yi = Z, vi = Z.var, mods = ~ age.m, data = ds)
# age_MM_RE

#### EXTRA: Back Transform ####
ds$back.transform <- (exp(2*ds$Z) - 1) / (exp(2*ds$Z) + 1)



#### Clean evironment #### 
rm(list = ls())

# Normal output
sink(file = NULL, type = c("output", "message"))
