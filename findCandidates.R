findCandidates = function(thisComb, CUTOFF){
  
  # get the frequency table
  uniqX = unique(dplyr::select(mydf, thisComb))
  pmf = count(mydf, vars = c(thisComb, "y"))
  freq = matrix(0, nrow = nrow(uniqX), ncol = 2)
  for(i in 1:nrow(uniqX)){
    for(j in 1:nrow(pmf)){
      if(prod(uniqX[i,] == pmf[j,1:(ncol(pmf)-2)]) == 1){
        if(pmf[j,ncol(pmf)-1] == 0) freq[i,1] = pmf[j,ncol(pmf)]
        if(pmf[j,ncol(pmf)-1] == 1) freq[i,2] = pmf[j,ncol(pmf)]
      }
    }
  }
  
  marg = freq[,1] - freq[,2]
  best = rep(0, nrow(uniqX))
  best[marg < 0] = 1
  
  # compute errors from freq and rule
  findErrors = function(freq, rule){
    errors = sum(freq[as.logical(rule), 1]) + sum(freq[as.logical(1-rule), 2])
    return(errors)
  }
  minErrors = findErrors(freq, best)
  
  marg = abs(marg)
  myOrder = order(marg)
  marg = marg[myOrder]
  bound = minErrors * (CUTOFF-1)
  
  source("findRash.R")
  NUM_MODELS = 5000
  selection = findRashomon(A = marg, m = NUM_MODELS, bound, msg = FALSE)
  
  candidates = matrix(0, nrow = nrow(uniqX), ncol = length(selection))
  for(i in 1:length(selection)){
    choice = selection[[i]]
    myRow = myOrder[choice]
    myPrediction = best
    myPrediction[myRow] = 1 - best[myRow]
    candidates[,i] = myPrediction
  }
  candidates = cbind(best, candidates)
  
  errors = rep(0,ncol(candidates))
  for(i in 1:ncol(candidates)){
    errors[i] = sum(freq[candidates[,i]==1,1]) + sum(freq[candidates[,i]==0,2])
  }

  outcome = list(uniqX = uniqX, candidates = candidates, errors = errors)
  
  return(outcome)
}