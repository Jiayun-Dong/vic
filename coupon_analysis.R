library(dplyr)
library(plotly)
library(plyr)
library(dummies)

setwd("/Users/jiayun/Dropbox/VID_local/coupon")

rm(list = ls())

myData = read.csv("coupon.csv", header = TRUE)[,-1]

# focus on coffee house only and filter out missing data and those without a car
df <- dplyr::select(myData, destination, passanger, weather, temperature, time, coupon,
                    expiration, gender, age, education,
                    occupation, income, car, CoffeeHouse, toCoupon_GEQ5min, toCoupon_GEQ15min,
                    toCoupon_GEQ25min, direction_same, Y) %>% 
  filter(coupon == "Coffee House") %>%
  filter(CoffeeHouse != "") %>%
  filter(car == 1)

df = df[ , -which(names(df) %in% c("coupon","car","toCoupon_GEQ5min","toCoupon_GEQ15min","toCoupon_GEQ25min"))]
df = df[ , -which(names(df) %in% c("temperature", "occupation"))]

df = cbind(df, dummy(df$destination, sep = "_"))
df = cbind(df, dummy(df$passanger, sep = "_"))
df = cbind(df, dummy(df$weather, sep = "_"))
df = cbind(df, dummy(df$time, sep = "_"))
df = cbind(df, dummy(df$expiration, sep = "_"))
df = cbind(df, dummy(df$gender, sep = "_"))
df = cbind(df, dummy(df$age, sep = "_"))
df = cbind(df, dummy(df$education, sep = "_"))
df = cbind(df, dummy(df$income, sep = "_"))
df = cbind(df, dummy(df$CoffeeHouse, sep = "_"))
df = cbind(df, dummy(df$direction_same, sep = "_"))
df = df[ , -which(names(df) %in% c("destination", "passanger", "weather", "time", "expiration", "gender",
                                   "age", "education", "income", "CoffeeHouse", "direction_same"))]

colnames(df)[3] = "NoUrgentPlace"
colnames(df)[6] = "Friends"
colnames(df)[7] = "Kids"
colnames(df)[29:34] = c("Associates", "Bachelors", "Graduate", "HighSchool", "SomeCollege", "SomeHighSchool")
colnames(df)[35:43] = c("100000plus", "12500", "25000", "37500", "50000", "62500", "75000", "87500", "lessthan12500")
colnames(df)[44:48] = c("Coffee1", "Coffee4", "Coffee8plus", "Coffeelessthan1", "CoffeeNever")
colnames(df)[49:50] = c("sameDirection", "oppoDirection")


df <- dplyr::select(df, NoUrgentPlace, Friends, df_Sunny, df_Male, df_1d, CoffeeNever, sameDirection, Y)
colnames(df) = c("noUrgentPlace", "friends", "sunny", "male", "expOneDay", "noCoffee", "sameDirection", "Y")
variable.names(df)

###########################
# load functions

# compute total loss of a logistic predictor
# input: a dataset (X,y) and a logistic model beta
# output: total loss

totalloss = function(beta, X, y){
  sum(log(1 + exp(-y * (beta[1] + X %*% beta[-1]))))
}

# find rashomon set
# input: 
#   center of the rashomon set: beta
#   the size of the parameter space in which to search for good models: range
#   dataset: X, y
#   the cutoff of loss - a model is in the rashomon set if its loss < loss of the reference model 
#       (computed by the center) * ratio: ratio
# output: 
#   returns the rashomon set: rashomon

rash = function(beta, range, size = 11, X, y, ratio = 1.01){
  # create grid
  betagrid = expand.grid(rep(list(seq(from = -1, to = 1, length.out = size)), length(beta)))
  betagrid = sweep(betagrid, 2, range, "*")
  betagrid = sweep(betagrid, 2, beta, "+")
  # compute loss for each model in the grid
  l = apply(betagrid, 1, totalloss, X, y)
  # construct the rashomon set
  loss0 = totalloss(beta, X, y)
  # sum(l < loss0 * ratio)/nrow(betagrid)
  rashomon = betagrid[l < loss0 * ratio,]
  return(rashomon)
}

# model reliance of a single model on a specific variable
# input: 
#   vname - the variable to analyze
#   beta - the model to analyze
#   X, y - the dataset
# output: model reliance of a single model on a specific variable (ratio)

modelreliance = function(vname, beta, X, y){
  id = which(colnames(X) == vname)
  p = sum(X[,id] == 1)/nrow(X)
  X0 = X
  # replaced by 1 with prob p
  X0[,id] = 1
  loss = totalloss(beta, X0, y) * p
  X0[,id] = 0
  loss = loss + totalloss(beta, X0, y) * (1-p) # this is the loss after shuffle the obs
  mr = loss/totalloss(beta, X, y)
  return(mr)
}

# model reliance of a single model on all variables
# input: 
#   beta - the model to analyze
#   X, y - the dataset
# output: model reliance of a single model on all variables (ratio)

modelreliance.full = function(vlist, beta, X, y){
  mr = rep(0, length(vlist))
  for(i in 1:length(vlist)){
    vname = vlist[i]
    mr[i] = modelreliance(vname, beta, X, y)
  }
  return(mr)
}


#####################################
# data
X = as.matrix(df[,colnames(df) != "Y"])
yy = df$Y
y = yy
y[y == 0] = -1
# variables of interest
vlist = variable.names(X)
# parameter
cutoff = 1.05
NUM = 500
NUM_PCA = 3
size = 10 # for initial sampling
set.seed(2018)
epsScale = 1.25

# logistic regression
model = glm(yy ~ X, family = "binomial")
summary(model)
beta = model$coefficients
se = summary(model)$coefficients[, 2]
# prediction error
yhat = rep(-1, nrow(df))
yhat[beta[1] + X %*% beta[-1] > 0] = 1
prederror = sum(y != yhat)
prederror/nrow(X)

# total loss
loss0 = totalloss(beta, X, y)
bound = loss0*cutoff

# initial sampling
# sampling in a box
points = matrix(runif(length(beta) * NUM, -1, 1), nrow = NUM, ncol = length(beta))
points = points %*% diag(size*se)
points = sweep(points, 2, beta, "+")
# keep good models only
myLoss = rep(0, nrow(points))
for(i in 1:nrow(points)){
  myLoss[i] = totalloss(points[i,],X, y)
}
points = points[myLoss < bound, ]

source("ellipsoid.R")

for(t in 1:NUM_PCA){
  pca = prcomp(points, center = TRUE, scale = FALSE)
  cen = pca$center
  len = apply(abs(pca$x), 2, max) * epsScale
  rot = pca$rotation
  
  points = ellipsoid(cen, len, rot, NUM)
  myLoss = rep(0, nrow(points))
  for(i in 1:nrow(points)){
    myLoss[i] = totalloss(points[i,],X, y)
  }
  points = points[myLoss < bound, ]
}

# myLoss = rep(0, nrow(points))
# for(i in 1:nrow(points)){
#   myLoss[i] = totalloss(points[i,],X, y)
# }
# hist(myLoss)

# apply(points, 2, max)
# apply(points, 2, min)
# beta

rashomon = as.data.frame(points)

# model reliance for a single model - the center of the rashomon set
mr = modelreliance.full(vlist, beta, X, y)
mr = matrix(mr, nrow = 1)
colnames(mr) = vlist

# model class reliance for all models in the rashomon set
mcr = matrix(0, nrow = nrow(rashomon), ncol = length(vlist))
for(i in 1:nrow(rashomon)){
  mcr[i, ] = modelreliance.full(vlist, as.double(rashomon[i, ]), X, y)
}
mcr = data.frame(mcr)
colnames(mcr) = vlist

mcrMax = apply(mcr,2,max)
mcrMin = apply(mcr,2,min)
mcrMax
mcrMin


# model selection
score = 0.5 * mcr[,1] + 0.5 * mcr[,3]
mcr[which(score == min(score)),]
beta0 = rashomon[which(score == min(score)),]
beta0 = as.matrix(beta0, nrow = 1)
colnames(beta0) = c("intercept", colnames(X))
