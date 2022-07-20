
##################################

printSelection = function(S){
  output = NULL
  for(i in 1:length(S)){
    sum = S[[i]]$sum
    choice = S[[i]]$choice
    output[i] = paste(c("(", sum, ", [", paste(choice, collapse = ","), "])"), collapse = "")
  }
  cat(paste0(output, sep = "\n"))
}

printTuple = function(tuple){
  sum = tuple$sum
  choice = tuple$choice
  output = paste(c("(", sum, ", [", paste(choice, collapse = ","), "])"), collapse = "")
}

##################################

findRashomon = function(A = c(7, 46, 74, 115, 167, 204, 294, 509), m = 10, bound = 150, msg = FALSE){
  # numRemoved = sum(A==0)
  # A = A[!A==0]
  n = length(A)
  
  S = list()
  for(i in 1:min(m, n)){
    tuple = list(A[i], i)
    names(tuple) = c("sum", "choice")
    S[[i]] = tuple
  }
  if(m > n){
    tuple = list(2^31, -1)
    names(tuple) = c("sum", "choice")
    aux = list()
    for(i in 1:(m-n)){
      aux[[i]] = tuple
    }
    S = append(S,aux)
  }
  
  noReplacement = FALSE
  exceedBound = FALSE
  for(i in 1:m){
    myTuple = S[[i]]
    if(msg == TRUE){
      display = sprintf("Add elements in A to the %ith tuple %s", i, printTuple(myTuple))
      print(display)
    }
    if(myTuple$sum > bound){
      exceedBound = TRUE
    }
    if(max(myTuple$choice) < n){
      for(j in (max(myTuple$choice)+1):n){
        tuple = list(myTuple$sum + A[j], c(myTuple$choice, j))
        names(tuple) = c("sum", "choice")
        maxSum = S[[length(S)]]$sum
        if(tuple$sum < maxSum){
          if(msg == TRUE){
            display = sprintf("  max tuple = %s is replaced by new tuple = %s", 
                              printTuple(S[[length(S)]]), printTuple(tuple))
            print(display)
          } 
          for(k in (i+1):m){
            if(S[[k]]$sum >= tuple$sum){
              myHead = S[1:(k-1)]
              myHead[[k]] = tuple
              myTail = S[k:(length(S)-1)]
              S = append(myHead, myTail)
              break
            }
          }
        } else {
          noReplacement = TRUE
          break
        }
      }
    }
    if(exceedBound == TRUE){
      if(msg == TRUE) print("  exceed bound! End of search!")
      while(S[[length(S)]]$sum > bound){
        S = S[-length(S)]
      }
      break
    }
    if(noReplacement == TRUE){
      if(msg == TRUE) print("  no replacement! End of search!")
      break
    }
    if(msg == TRUE) print('Current selection S is \n')
    if(msg == TRUE) printSelection(S)
    if(msg == TRUE) print('======= \n')
  }
  if(msg == TRUE) printSelection(S)
  
  # get a list of m choices
  selection = list()
  for(i in 1:length(S)){
    selection[[i]] = S[[i]]$choice
    # selection[[i]] = S[[i]]$choice + numRemoved
  }
  return(selection)
  
}


