
library(dplyr)
library(plotly)
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering visualization
library(dendextend) # for comparing two dendrograms
library(glmnet)

setwd("/Users/jiayun/Documents/coding/vgg")

# load data
rm(list = ls())

####################################

# compute total loss of a logistic predictor
# input: a dataset (X,y) and a logistic model beta
# output: total loss

totalloss = function(beta, X, y, lambda){
  sum(log(1 + exp(-y * (beta[1] + X %*% beta[-1])))) + lambda * sum(abs(beta))
}

modelreliance = function(index, beta, X, y){
  perm = sample(1:nrow(X), nrow(X), replace = FALSE) #permutation
  X0 = X
  loss_vector = rep(0, M)
  for(iter in 1:M){
    X0[,index] = X[perm, index]
    loss_vector[iter] = totalloss(beta, X0, y, lambda)
  }
  return(mean(loss_vector)/totalloss(beta, X, y, lambda))
}

# 
# # model reliance of a single model on all variables
# # input: 
# #   beta - the model to analyze
# #   X, y - the dataset
# # output: model reliance of a single model on all variables (ratio)
# 
modelreliance.full = function(beta, X, y){
  mr = rep(0, length(beta)-1)
  for(i in 1:length(mr)){
    mr[i] = modelreliance(i, beta, X, y)
  }
  return(mr)
}


####################################

feature = read.csv(file="features_small.csv", head = FALSE, sep=",")
outcome = read.csv(file="outcomes_small.csv", head = FALSE, sep=",")
feature = as.matrix(feature)
outcome = as.matrix(outcome)


cv.out <- cv.glmnet(feature, outcome, alpha = 1, family = "binomial", type.measure = "mse")
# best value of lambda
lambda <- cv.out$lambda.1se
# regression coefficients
beta = coef(cv.out, s = lambda)
# prediction accuracy
lasso_prob <- predict(cv.out, newx = feature, s = lambda, type = "response")
lasso_predict = ifelse(lasso_prob > 0.5, 1, 0)
sum(lasso_predict == outcome)/length(outcome)

# feature selection
id = beta@i[-1]
X = feature[,id]
beta_star = beta@x


##############################
# vid analysis

# X = feature
y = outcome
y[y==0] = -1

# parameter
epsilon = 0.05
NUM = 500 # number of points in each round of sampling
NUM_PCA = 10 # 10 for lower dim
size = 0.001 # for initial sampling
set.seed(2018)
epsScale = 3.75 # over-sample a bit # 3.75 for lower dim
M = 20

loss0 = totalloss(beta_star, X, y, lambda)
bound = loss0 * (1 + epsilon)

# initial sampling
# get NUM models in a box around beta*
points = matrix(runif(length(beta_star) * NUM, -1, 1), nrow = NUM, ncol = length(beta_star))
points = points * max(beta_star) * 0.1
points = sweep(points, 2, beta_star, "+")
# keep good models only
myLoss = rep(0, nrow(points))
for(i in 1:nrow(points)){
  myLoss[i] = totalloss(points[i,], X, y, lambda)
}
points = points[myLoss < bound, ]

source("ellipsoid.R")

for(t in 1:NUM_PCA){
  # pca for points from last sampling process
  pca = prcomp(points, center = TRUE, scale = FALSE, rank. = 513)
  cen = pca$center
  len = apply(abs(pca$x), 2, max) * epsScale
  rot = pca$rotation
  # next round of samping
  points = ellipsoid(cen, len, rot, NUM)
  # eliminates those points that exceed the bound
  myLoss = rep(0, nrow(points))
  for(i in 1:nrow(points)){
    myLoss[i] = totalloss(points[i,],X, y, lambda)
  }
  points = points[myLoss < bound, ]
}

rashomon = as.data.frame(points)

# model reliance for a single model - the center of the rashomon set
mr = modelreliance.full(beta_star, X, y)

# model class reliance for all models in the rashomon set
mcr = matrix(0, nrow = nrow(rashomon), ncol = length(beta_star)-1)
for(i in 1:nrow(rashomon)){
  mcr[i, ] = modelreliance.full(as.double(rashomon[i, ]), X, y)
  
}
mcr = data.frame(mcr)


# characteristics of the rachomon set (from the mcr perspective)
cor(mcr)

apply(mcr, 2, max)
order(apply(mcr, 2, max), decreasing = TRUE)

apply(mcr, 2, min)
order(apply(mcr, 2, min), decreasing = TRUE)

apply(mcr, 2, mean)
myOrder = order(apply(mcr, 2, mean), decreasing = TRUE)


# clustering

k2 = kmeans(mcr, centers = 4, nstart = 25)
sub_grp = k2$cluster

  # clustering results
  cat = sub_grp 
  sum(cat==1)
  sum(cat==2)
  sum(cat==3)
  sum(cat==4)

  # find centers
  red_center = order(rowSums(sweep(mcr, 2, k2$centers[1,])^2), decreasing = FALSE)[1]
  green_center = order(rowSums(sweep(mcr, 2, k2$centers[2,])^2), decreasing = FALSE)[1]
  blue_center = order(rowSums(sweep(mcr, 2, k2$centers[3,])^2), decreasing = FALSE)[1]
  purple_center = order(rowSums(sweep(mcr, 2, k2$centers[4,])^2), decreasing = FALSE)[1]
  model_id = c(red_center, green_center, blue_center, purple_center)


  # most important variable to a cluster
  mcr_red = apply(mcr[cat == 1,], 2, mean)
  mcr_green = apply(mcr[cat == 2,], 2, mean)
  mcr_blue = apply(mcr[cat == 3,], 2, mean)
  mcr_purple = apply(mcr[cat == 4,], 2, mean)
  mcr_avg = apply(mcr, 2, mean)
  
  red = order(mcr_red - mcr_avg, decreasing = TRUE)
  (mcr_red - mcr_avg)[red]
  green = order(mcr_green - mcr_avg, decreasing = TRUE)
  (mcr_green - mcr_avg)[green]
  blue = order(mcr_blue - mcr_avg, decreasing = TRUE)
  (mcr_blue - mcr_avg)[blue]
  purple = order(mcr_purple - mcr_avg, decreasing = TRUE)
  (mcr_purple - mcr_avg)[purple]

  vlist = c(red[1], green[1], blue[1], purple[1])
  mcr_small = mcr[,vlist]
  vlist = id[vlist]

  vname = NULL
  for(i in 1:length(vlist)){
    vname = c(vname, sprintf("feature_%i", vlist[i]))
  }


# VID
color = c("brown2", "chartreuse3", "aquamarine2", "darkorchid1")
m = max(mcr_small) * 1.02
plot_by_cat = function(x, y , name1, name2, m, cat){
  index = 1:length(x)
  i = 1
  plot(x[cat == i], y[cat == i], xlab = name1, ylab = name2, pch = ".", xlim = c(1,m), ylim = c(1,m), col = color[i],
       cex.lab=1.25, cex.axis=1.25, cex.main=1.5, cex.sub=1)
  for(i in 2:4){
    points(x[cat == i], y[cat == i], pch = "*", col = color[i])
  }
}

par(mfrow=c(4,4))
plot(0,type='n',axes=FALSE,ann=FALSE)
plot_by_cat(mcr_small[,2], mcr_small[,1], vname[2], vname[1], m, cat)
plot_by_cat(mcr_small[,3], mcr_small[,1], vname[3], vname[1], m, cat)
plot_by_cat(mcr_small[,4], mcr_small[,1], vname[4], vname[1], m, cat)

plot_by_cat(mcr_small[,1], mcr_small[,2], vname[1], vname[2], m, cat)
plot(0,type='n',axes=FALSE,ann=FALSE)
plot_by_cat(mcr_small[,3], mcr_small[,2], vname[3], vname[2], m, cat)
plot_by_cat(mcr_small[,4], mcr_small[,2], vname[4], vname[2], m, cat)

plot_by_cat(mcr_small[,1], mcr_small[,3], vname[1], vname[3], m, cat)
plot_by_cat(mcr_small[,2], mcr_small[,3], vname[2], vname[3], m, cat)
plot(0,type='n',axes=FALSE,ann=FALSE)
plot_by_cat(mcr_small[,4], mcr_small[,3], vname[4], vname[3], m, cat)

plot_by_cat(mcr_small[,1], mcr_small[,4], vname[1], vname[4], m, cat)
plot_by_cat(mcr_small[,2], mcr_small[,4], vname[2], vname[4], m, cat)
plot_by_cat(mcr_small[,3], mcr_small[,4], vname[3], vname[4], m, cat)
plot(0,type='n',axes=FALSE,ann=FALSE)

# visualization models
# model_id = c(36, 285, 50, 9)
model_coefficients = rashomon[model_id,]
filter_id = id
write.csv(model_coefficients, "model.csv")
write.csv(filter_id, "filter.csv")
