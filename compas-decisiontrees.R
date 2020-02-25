
library(dplyr)
library(plyr)
library(plotly)

setwd("/Users/jiayun/Dropbox/VID_local/compas")

# load data
rm(list = ls())

raw_data = read.csv("./compas.csv")
df <- dplyr::select(raw_data, Age.18.20, Gender.Male, Race.African.American, Juvenile.Felonies.0,
                    Juvenile.Crimes.0, Prior.Crimes.0, Current.Charge.Degree.Misdemeanor, 
                    Recidivate.Within.Two.Years)

###################################

# find candidates for a selection of variables
source("findCandidates.R")

# Compute prediction error of a prediction model
# inputs: 
#   classification rule: (uniqX, rule)
#   data: (X, y)
prederror = function(RashX, rule, X, y){
  pmf = count(dplyr::select(as.data.frame(cbind(X,y)), vars = c(variable.names(RashX)),"y"))
  freq = matrix(0, nrow = nrow(RashX), ncol = 2)
  for(i in 1:nrow(RashX)){
    for(j in 1:nrow(pmf)){
      if(prod(RashX[i,] == pmf[j,1:(ncol(pmf)-2)]) == 1){
        if(pmf[j,ncol(pmf)-1] == 0) freq[i,1] = pmf[j,ncol(pmf)]
        if(pmf[j,ncol(pmf)-1] == 1) freq[i,2] = pmf[j,ncol(pmf)]
      }
    }
  }
  errors = sum(freq[as.logical(rule), 1]) + sum(freq[as.logical(1-rule), 2])
  return(errors)
}

# correctness check
# myError1 = allTrees[[1]]$errors
# tic()
# myError2 = rep(0,length(myError1))
# for(i in 1:length(myError2)){
#   myError2[i] = prederror(RashX, allTrees[[1]]$candidates[,i],X,y)
# }
# toc()
# myError1 == myError2

# Compute loss after shuffling a variable for a given model

lossShuffled = function(vname, RashX, rule, X, y){
  id = which(colnames(X) == vname)
  p = sum(X[,id] == 1)/nrow(X)
  X0 = X
  # replaced by 1 with prob p
  X0[,id] = 1
  loss0 = prederror(RashX, rule, X0, y) * p
  X0[,id] = 0
  loss0 = loss0 + prederror(RashX, rule, X0, y) * (1-p)
  return(loss0)
}

# Compute reliance on vlist for a single model

modelreliance = function(vlist, RashX, rule, X, y){
  loss = prederror(RashX, rule, X, y)
  mr = rep(0, length(vlist))
  for(i in 1:length(vlist)){
    vname = vlist[i]
    mr[i] = lossShuffled(vname, RashX, rule, X, y)
  }
  mr = mr/loss
  return(mr)
}


#######################################

# this experiment select 6 variables
  y = df$Recidivate.Within.Two.Years
  X = cbind(df$Age.18.20, df$Race.African.American, df$Gender.Male, df$Prior.Crimes.0, df$Juvenile.Crimes.0 ,df$Current.Charge.Degree.Misdemeanor) 
  colnames(X) = c("age", "race", "gender", "prior", "juvenile", "charge") # MUST define variable names
  mydf = data.frame(X, y) 
  
# study VID of the following 
  vlist = colnames(X) # MUST match the predefined variable names

# cutoff: rashomon includes models with errors <= logistic loss * cutoff
  CUTOFF = 1.05

# look at decision trees with binary splits and can only split according to a fixed NUM of variables.
NUM = 4

# find all combinations of variables (excluding age and race, which must be included)
comb = combn(variable.names(X), NUM)

# loop over all combinations of variables, find all candidates, and get minError
allTrees = list()
minError = 2^30
numOfCandidates = 0

for(i in 1:ncol(comb)){
  # thisComb = c(vlist, comb[,i])
  thisComb = c(comb[,i])
  # get the candidates for rashomon set
  allTrees[[i]] = findCandidates(thisComb, CUTOFF)
  minError = min(c(minError, allTrees[[i]]$errors))
  numOfCandidates = numOfCandidates + length(allTrees[[i]]$errors)
}
numOfCandidates

# find the best tree, include any tree with loss <= cutoff*best
numOfModels = 0
del = NULL
for(i in 1:ncol(comb)){
  keep = allTrees[[i]]$errors < CUTOFF * minError
  if(sum(keep) == 0) {
    del = c(del,i)
  } else {
    theseTrees = allTrees[[i]]
    theseTrees$candidates = theseTrees$candidates[,keep]
    theseTrees$errors = theseTrees$errors[keep]
    numOfModels = numOfModels + sum(keep)
    allTrees[[i]] = theseTrees
  }
}
if(!is.null(del)){
  allTrees = allTrees[-del]
}

# this gives the Rashomon set
numOfModels

# compute model class reliance for the whole rashomon set, grouped by the selection of other variables
# this step takes about 40mins

mcrList = list()
for(i in 1:length(allTrees)){
  RashX = allTrees[[i]]$uniqX
  RashY = allTrees[[i]]$candidates
  RashY = as.matrix(RashY)
  mcr = matrix(0, nrow = ncol(RashY), ncol = length(vlist))
  for(j in 1:nrow(mcr)){
    mcr[j, ] = modelreliance(vlist, RashX, RashY[,j], X, y)
  }
  mcr = data.frame(mcr)
  mcrList[[i]] = mcr
}

# combine the trees
mcr = NULL
err = NULL
trees = matrix(0, nrow = numOfModels, ncol = 2)
colnames(trees) = c("X_index", "y_index")
id = 0
for(i in 1:length(allTrees)){
 mcr = rbind(mcr, mcrList[[i]])
 trees[(id+1):(id+length(allTrees[[i]]$errors)), 1] = i
 trees[(id+1):(id+length(allTrees[[i]]$errors)), 2] = 1:length(allTrees[[i]]$errors)
 err = c(err, allTrees[[i]]$errors)
 id = id + length(allTrees[[i]]$errors)
}
mcr = as.data.frame(mcr)
colnames(mcr) = colnames(X)


par(mfrow=c(4,4))
par(mar = c(4.5,4.5,2.5,2.5))


plot(0,type='n',axes=FALSE,ann=FALSE)
plot(mcr$race, mcr$age, xlab = "race", ylab = "age", pch = ".", xlim = c(1,1.25), ylim = c(1,1.25))
plot(mcr$prior, mcr$age, xlab = "prior", ylab = "age", pch = ".", xlim = c(1,1.25), ylim = c(1,1.25))
plot(mcr$gender, mcr$age, xlab = "gender", ylab = "age", pch = ".", xlim = c(1,1.25), ylim = c(1,1.25))

plot(mcr$age, mcr$race, xlab = "age", ylab = "race", pch = ".", xlim = c(1,1.25), ylim = c(1,1.25))
plot(0,type='n',axes=FALSE,ann=FALSE)
plot(mcr$prior, mcr$race, xlab = "prior", ylab = "race", pch = ".", xlim = c(1,1.25), ylim = c(1,1.25))
plot(mcr$gender, mcr$race, xlab = "gender", ylab = "race", pch = ".", xlim = c(1,1.25), ylim = c(1,1.25))

plot(mcr$age, mcr$prior, xlab = "age", ylab = "prior", pch = ".", xlim = c(1,1.25), ylim = c(1,1.25))
plot(mcr$race, mcr$prior, xlab = "race", ylab = "prior", pch = ".", xlim = c(1,1.25), ylim = c(1,1.25))
plot(0,type='n',axes=FALSE,ann=FALSE)
plot(mcr$gender, mcr$prior, xlab = "gender", ylab = "prior", pch = ".", xlim = c(1,1.25), ylim = c(1,1.25))

plot(mcr$age, mcr$gender, xlab = "age", ylab = "gender", pch = ".", xlim = c(1,1.25), ylim = c(1,1.25))
plot(mcr$race, mcr$gender, xlab = "race", ylab = "gender", pch = ".", xlim = c(1,1.25), ylim = c(1,1.25))
plot(mcr$prior, mcr$gender, xlab = "prior", ylab = "gender", pch = ".", xlim = c(1,1.25), ylim = c(1,1.25))
plot(0,type='n',axes=FALSE,ann=FALSE)


modelreliance(vlist, thisTree$uniqX, rash[,876], X, y)

prederror(thisTree$uniqX, rash[,350], X, y)
