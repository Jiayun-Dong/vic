

library(plotly)
library(plyr)
library(dplyr)

setwd("/Users/jiayun/Dropbox/VID_local/compas")

# load data
rm(list = ls())

raw_data = read.csv("./compas.csv")
df <- dplyr::select(raw_data, Age.18.20, Age..45, Gender.Male, Race.African.American, Race.Caucasian,
                    Race.Asian, Race.Hispanic, Race.Native.American, Race.Other, Juvenile.Felonies.0,
                    Juvenile.Felonies.1.3, Juvenile.Felonies.3, Juvenile.Crimes.0, Juvenile.Crimes.1.3,
                    Juvenile.Crimes.3, Juvenile.Crimes.5, Prior.Crimes.0, Prior.Crimes.1.3, 
                    Prior.Crimes.3, Prior.Crimes.5, Current.Charge.Degree.Misdemeanor, 
                    Recidivate.Within.Two.Years)

###########################
# load functions

# compute total loss of a logistic predictor
# input: a dataset (X,y) and a logistic model beta
# output: total loss

totalloss = function(beta, X, y){
  sum(log(1 + exp(-y * (beta[1] + X %*% beta[-1]))))
}

# # find rashomon set
# # input: 
# #   center of the rashomon set: beta
# #   the size of the parameter space in which to search for good models: range
# #   dataset: X, y
# #   epsilon - a model is in the rashomon set if its loss < loss of the reference model 
# #       (computed by the center) * epsilon: epsilon
# # output: 
# #   returns the rashomon set: rashomon
# 
# rash = function(beta, range, size = 11, X, y, epsilon = 1.01){
#   # create grid
#   betagrid = expand.grid(rep(list(seq(from = -1, to = 1, length.out = size)), length(beta)))
#   betagrid = sweep(betagrid, 2, range, "*")
#   betagrid = sweep(betagrid, 2, beta, "+")
#   # compute loss for each model in the grid
#   l = apply(betagrid, 1, totalloss, X, y)
#   # construct the rashomon set
#   loss0 = totalloss(beta, X, y)
#   # sum(l < loss0 * epsilon)/nrow(betagrid)
#   rashomon = betagrid[l < loss0 * epsilon,]
#   return(rashomon)
# }

# model reliance of a single model on a specific variable
# input: 
#   vname - the variable to analyze
#   beta - the model to analyze
#   X, y - the dataset
# output: model reliance of a single model on a specific variable (ratio definition)

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


####################################
# data

X = cbind(df$Age.18.20, df$Race.African.American, df$Prior.Crimes.0, df$Gender.Male, df$Juvenile.Crimes.0, df$Current.Charge.Degree.Misdemeanor)
colnames(X) = c("age", "race", "prior", "gender", "juvenilecrime", "currentcharge")
y = df$Recidivate.Within.Two.Years
y[y==0] = -1

# variables of interest
vlist = c("age", "race", "prior", "gender")

# parameter
epsilon = 0.05
NUM = 500 # number of points in each round of sampling
NUM_PCA = 5
size = 10 # for initial sampling
set.seed(2018)
epsScale = 1.5 # over-sample a bit

# run logistic regression, get beta* and its s.e., compute its loss and prediction errors
  # logistic regression / center of the Rashomon set
  yy = y
  yy[y == -1] = 0
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
  bound = loss0*(1+epsilon)

# initial sampling
  # get NUM models in a box around beta*
  points = matrix(runif(length(beta) * NUM, -1, 1), nrow = NUM, ncol = length(beta))
  points = points %*% diag(size*se)
  points = sweep(points, 2, beta, "+")
  # keep good models only
  myLoss = rep(0, nrow(points))
  for(i in 1:nrow(points)){
    myLoss[i] = totalloss(points[i,],X, y)
  }
  points = points[myLoss < bound, ]
  
  # tuning the parameters here
  # choose epsScale = 1.5
  # NUM_PCA = 5
  
  source("ellipsoid.R")

  for(t in 1:NUM_PCA){
    # pca for points from last sampling process
    pca = prcomp(points, center = TRUE, scale = FALSE)
    cen = pca$center
    len = apply(abs(pca$x), 2, max) * epsScale
    rot = pca$rotation
    # next round of samping
    points = ellipsoid(cen, len, rot, NUM)
    # eliminates those points that exceed the bound
    myLoss = rep(0, nrow(points))
    for(i in 1:nrow(points)){
      myLoss[i] = totalloss(points[i,],X, y)
    }
    points = points[myLoss < bound, ]
  }
  
  
  rashomon = as.data.frame(points)
  
  # visualize rashomon set
  colnames(rashomon) = c("int", "v1", "v2", "v3")
  rashomonplot <- plot_ly(rashomon[,2:4], x = ~v1, y = ~v2, z = ~v3,
                          marker = list(size = 2)) %>%
    layout(title = "Rashomon Set",
           scene = list(xaxis = list(title = colnames(X)[1]),
                        yaxis = list(title = colnames(X)[2]),
                        zaxis = list(title = colnames(X)[3])))
  rashomonplot
  
  # model reliance for a single model - the center of the rashomon set
  mr = modelreliance.full(vlist, beta, X, y)
  
  # model class reliance for all models in the rashomon set
  mcr = matrix(0, nrow = nrow(rashomon), ncol = length(vlist))
  for(i in 1:nrow(rashomon)){
    mcr[i, ] = modelreliance.full(vlist, as.double(rashomon[i, ]), X, y)
  }
  mcr = data.frame(mcr)
  colnames(mcr) = c("v1", "v2", "v3", "v4")
  
  # visualize VID
  VIC <- plot_ly(mcr, x = ~v1, y = ~v2, z = ~v3,
                 marker = list(size = 2)) %>%
    layout(title = "Variable Importance Cloud",
           scene = list(xaxis = list(title = colnames(X)[1]),
                        yaxis = list(title = colnames(X)[2]),
                        zaxis = list(title = colnames(X)[3])))
  VIC
  
  # plot VID
    par(mfrow=c(4,4))
    plot(0,type='n',axes=FALSE,ann=FALSE)
    plot(mcr$v2, mcr$v1, xlab = "race", ylab = "age", pch = ".", xlim = c(1,1.15), ylim = c(1,1.15),
         cex.lab=1.25, cex.axis=1.25, cex.main=1.5, cex.sub=1)
    plot(mcr$v3, mcr$v1, xlab = "prior", ylab = "age", pch = ".", xlim = c(1,1.15), ylim = c(1,1.15),
         cex.lab=1.25, cex.axis=1.25, cex.main=1.5, cex.sub=1)
    plot(mcr$v4, mcr$v1, xlab = "gender", ylab = "age", pch = ".", xlim = c(1,1.15), ylim = c(1,1.15),
         cex.lab=1.25, cex.axis=1.25, cex.main=1.5, cex.sub=1)
    
    plot(mcr$v1, mcr$v2, xlab = "age", ylab = "race", pch = ".", xlim = c(1,1.15), ylim = c(1,1.15),
         cex.lab=1.25, cex.axis=1.25, cex.main=1.5, cex.sub=1)
    plot(0,type='n',axes=FALSE,ann=FALSE)
    plot(mcr$v3, mcr$v2, xlab = "prior", ylab = "race", pch = ".", xlim = c(1,1.15), ylim = c(1,1.15),
         cex.lab=1.25, cex.axis=1.25, cex.main=1.5, cex.sub=1)
    plot(mcr$v4, mcr$v2, xlab = "gender", ylab = "race", pch = ".", xlim = c(1,1.15), ylim = c(1,1.15),
         cex.lab=1.25, cex.axis=1.25, cex.main=1.5, cex.sub=1)

    plot(mcr$v1, mcr$v3, xlab = "age", ylab = "prior", pch = ".", xlim = c(1,1.15), ylim = c(1,1.15),
         cex.lab=1.25, cex.axis=1.25, cex.main=1.5, cex.sub=1)
    plot(mcr$v2, mcr$v3, xlab = "race", ylab = "prior", pch = ".", xlim = c(1,1.15), ylim = c(1,1.15),
         cex.lab=1.25, cex.axis=1.25, cex.main=1.5, cex.sub=1)
    plot(0,type='n',axes=FALSE,ann=FALSE)
    plot(mcr$v4, mcr$v3, xlab = "gender", ylab = "prior", pch = ".", xlim = c(1,1.15), ylim = c(1,1.15),
         cex.lab=1.25, cex.axis=1.25, cex.main=1.5, cex.sub=1)

    plot(mcr$v1, mcr$v4, xlab = "age", ylab = "gender", pch = ".", xlim = c(1,1.15), ylim = c(1,1.15),
         cex.lab=1.25, cex.axis=1.25, cex.main=1.5, cex.sub=1)
    plot(mcr$v2, mcr$v4, xlab = "race", ylab = "gender", pch = ".", xlim = c(1,1.15), ylim = c(1,1.15),
         cex.lab=1.25, cex.axis=1.25, cex.main=1.5, cex.sub=1)
    plot(mcr$v3, mcr$v4, xlab = "prior", ylab = "gender", pch = ".", xlim = c(1,1.15), ylim = c(1,1.15),
         cex.lab=1.25, cex.axis=1.25, cex.main=1.5, cex.sub=1)
    plot(0,type='n',axes=FALSE,ann=FALSE)
    
    