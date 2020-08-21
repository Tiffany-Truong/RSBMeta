#### Run script using command - paste to console
# source("RSBMetascript_overall.R", echo=TRUE, max.deparse.length=10000)


# Clean environment
rm(list = ls())
# Clear Console
cat("\014")

# Load library
library(metafor)
library(plotly)
library(robumeta)

filePath <- "~/git/RSBMeta/Analyses/21-08-2020/Overalloutput.txt"

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


#### Data Prep #### 

ds = read.csv("/Users/tiffanytruong/git/RSBMeta/Analyses/20-08-2020/dataset.csv", header = T, sep = ",")
setwd("~/git/RSBMeta/Analyses/21-08-2020/")


# rsb type subsetting
# fullset = read.csv("dataset.csv", header = T, sep = ",")
# ds =
# subset(fullset, rsbtype == 1)
# subset(fullset, rsbtype == 2)
# subset(fullset, rsbtype == 3) 
# subset(fullset, rsbtype == 4)
# subset(fullset, rsbtype == 5)
# subset(fullset, rsbtype == 6)


# Fisher's Z transformation
ds$Z <- .5*log((1 + ds$r) / (1 - ds$r))

# Variance of Fisher's Z
ds$Z.var <- 1 / (ds$n - 3)

####  Forest Plot ####

# Create and print default forest plot
# yi: effect size
# vi: effect size variance
# data: data set 

jpeg('Forest_overall.jpg', width = 610, height = 880)
forestplot_ds <- rma.uni(yi = Z, vi = Z.var, data = ds)
forest(forestplot_ds, slab = ds$article, order = order(ds$Z), showweights = TRUE)
dev.off()

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
jpeg('Forest_Plot_multivariate.jpg', width = 610, height = 880)
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
jpeg('influence_overall.jpg', width = 610, height = 880)
plot(influence(RE.model))
dev.off()


leave1out(RE.model)

#### Publication Bias ####

#Funnel Plots 

par(mfrow = c(1,1))

# Funnel plot w/ SE as vertical axis
jpeg('funnel_overall.jpg', width = 610, height = 880)
funnel(RE.model)
dev.off()


# Egger's Regression Test
regtest(x = ds$Z, vi = ds$Z.var, model = "lm")

# Basic trim-and-fill code
trimfill(RE.model, side = "right", verbose = T)
trimfill(RE.model, side = "left", verbose = T)


# Plot both versions of trim-and-fill (optional)
par(mfrow = c(2,1), jpeg('tnf_overall.jpg'))

funnel(trimfill(RE.model, side = "right"), main = "Right-Side Imputation")
funnel(trimfill(RE.model, side = "left"), main = "Left-Side Imputation")
dev.off()

# Fail safe N 
fsn(yi = Z, vi = Z.var, data = ds)

#### EXTRA: Back Transform ####
ds$back.transform <- (exp(2*ds$Z) - 1) / (exp(2*ds$Z) + 1)






#### Clean evironment #### 
rm(list = ls())

# Normal output
sink(file = NULL, type = c("output", "message"))

# Re-set WD for some reason 
setwd("~/git/RSBMeta")

