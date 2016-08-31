##################
### aprioriGen ###
##################

aprioriGen <- function(thisSets, step, cl) {

  # this function generates candidate matrices. first, compute all possible
  # unions of two frequent sets of the last step

  cand <- combn(1:ncol(thisSets), 2)
  cand <- apply(cand, 2, function(x) { thisSets[, x[1]] | thisSets[, x[2]] })

  # drop candidates, that do not have enough items, and duplicate candidates

  cand <- cand[, apply(cand, 2, function(x) { sum(x) == step }), drop = F]
  cand <- unique(cand, MARGIN = 2)

  # multiply transposed candidate matrix with frequent sets of the last step.
  # if a candidate-set combination is equal to the last step count, then the
  # set is a subset of the candidate. this implies that a subset of the
  # candidate with cardinality step - 1 is equal to the set. if the number of
  # such candidate-set combinations for the candidate is equal to the current
  # step count, then every subset of the candidate with cardinality step - 1
  # is among the frequent sets of the last step

  if (hasArg(cl)) {
    if (ncol(cand)) {

      # parMM sometimes returns different results than the standard matrix
      # multiplication, e.g. it drops dimensions. working around this here

      tmp <- parMM(cl, t(cand), thisSets) == step - 1

      if (is.matrix(tmp)) {
        subsCheck <- rowSums(tmp) == step
      } else {
        subsCheck <- sum(tmp) == step
      }

      rm(tmp)

    } else {
      subsCheck <- logical(0)
    }
  } else {
    subsCheck <- rowSums(t(cand) %*% thisSets == step - 1) == step
  }

  # multiplication with one gives a numeric matrix

  cand <- cand[, subsCheck, drop = F] * 1

  return(cand)

}


###################
### aprioriSets ###
###################

aprioriSets <- function(data, minSupp, cl) {

  # store metadata for later use

  itemInfo <- data@itemInfo
  dataName <- deparse(substitute(data))
  data <- data@data

  # the identity matrix is the candidate matrix for single item sets. the row
  # sums of the data matrix are the absolute support of the single item sets

  cand <- diag(nrow(data))
  candSupp <- rowSums(data) / ncol(data)
  # candSupp <- rowSums(t(cand) %*% data == 1) / ncol(data)
  allSets <- thisSets <- cand[, candSupp >= minSupp, drop = F]
  allSupp <- candSupp[candSupp >= minSupp]
  step <- 2

  while (ncol(thisSets) >= 2) {

    cand <- aprioriGen(thisSets, step, cl)

    # multiply transposed candidate matrix with data matrix. in the resulting
    # matrix, each element represents a candidate-transaction combination. if
    # an element is equal to the step count, then the candidate is a subset of
    # the transaction. summing up the rows of this boolean matrix gives the
    # absolute support of the candidates

    if (hasArg(cl)) {
      if (ncol(cand)) {

        tmp <- parMM(cl, t(cand), as.matrix(data) * 1) == step

        if (is.matrix(tmp)) {
          candSupp <- rowSums(tmp) / ncol(data)
        } else {
          candSupp <- sum(tmp) / ncol(data)
        }

        rm(tmp)

      } else {
        candSupp <- numeric(0)
      }
    } else {
      candSupp <- rowSums(t(cand) %*% data == step) / ncol(data)
    }

    allSupp <- c(allSupp, candSupp[candSupp >= minSupp])
    thisSets <- cand[, candSupp >= minSupp, drop = F]
    allSets <- cbind(allSets, thisSets)
    step <- step + 1

  }

  res <- new("itemsets",
             items   = new("itemMatrix",
                           data     = as(allSets, "ngCMatrix"),
                           itemInfo = itemInfo),
             quality = data.frame(support = allSupp),
             info    = list(data          = dataName,
                            ntransactions = ncol(data),
                            support       = minSupp,
                            confidence    = 1))

  return(res)

}


####################
### aprioriRules ###
####################

aprioriRules <- function(data, allSets, minConf, cl) {

  data <- data@data
  nTransact <- ncol(data)
  nItems <- nrow(data)

  sets <- allSets@items@data
  quality <- allSets@quality
  itemInfo <- allSets@items@itemInfo
  info <- allSets@info
  rm(allSets)

  # uncomment the following lines if you do not want to allow for empty lhs.
  # following arules in allowing for these

  # setsIndex <- colSums(sets) >= 2
  # sets <- sets[, setsIndex]
  nSets <- ncol(sets)

  allLHS <- allRHS <- matrix(nrow = nItems, ncol = 0)
  allConv <- allLift <- allConf <- allSupp <- c()

  for (i in 1:nSets) {

    setItems <- which(sets[, i])
    setSupp <- quality$support[i]
    # setSupp <- quality$support[setsIndex][i]

    # initialize lhs and rhs matrices for single item rhs as zero matrices and
    # set items to one element by element

    thisLHS <- thisRHS <- matrix(0, nrow = nItems, ncol = length(setItems))

    for (j in 1:ncol(thisRHS)) {

      rhsItems <- setItems[j]
      lhsItems <- setdiff(setItems, rhsItems)
      thisRHS[rhsItems, j] <- 1
      thisLHS[lhsItems, j] <- 1

    }

    rhsSupp <- rowSums(t(thisRHS) %*% data == 1) / nTransact
    lhsSupp <- rowSums(t(thisLHS) %*% data == length(lhsItems)) / nTransact

    thisConf <- setSupp / lhsSupp
    thisConv <- (1 - rhsSupp) / (1 - thisConf)
    thisLift <- setSupp / (lhsSupp * rhsSupp)

    thisRHS <- thisRHS[, thisConf >= minConf, drop = F]
    thisLHS <- thisLHS[, thisConf >= minConf, drop = F]

    allRHS <- cbind(allRHS, thisRHS)
    allLHS <- cbind(allLHS, thisLHS)

    allSupp <- c(allSupp, rep(setSupp, ncol(thisRHS)))
    allConf <- c(allConf, thisConf[thisConf >= minConf])
    allLift <- c(allLift, thisLift[thisConf >= minConf])
    allConv <- c(allConv, thisConv[thisConf >= minConf])

    step <- 2

    # the same for multiple item rhs. if you want to allow for these, comment
    # the next line. following arules in not allowing for multiple item rhs

    next

    while (ncol(thisRHS) >= 2) {

      thisRHS <- aprioriGen(thisRHS, step, cl)
      thisLHS <- matrix(0, nrow = nItems, ncol = ncol(thisRHS))

      for (j in 1:ncol(thisRHS)) {

        rhsItems <- which(thisRHS[, j] == 1)
        lhsItems <- setdiff(setItems, rhsItems)
        thisLHS[lhsItems, j] <- 1

      }

      rhsSupp <- rowSums(t(thisRHS) %*% data == step) / nTransact
      lhsSupp <- rowSums(t(thisLHS) %*% data == length(lhsItems)) / nTransact
      thisLift <- setSupp / (lhsSupp * rhsSupp)
      thisConf <- setSupp / lhsSupp

      thisRHS <- thisRHS[, thisConf >= minConf, drop = F]
      thisLHS <- thisLHS[, thisConf >= minConf, drop = F]

      allRHS <- cbind(allRHS, thisRHS)
      allLHS <- cbind(allLHS, thisLHS)

      allSupp <- c(allSupp, rep(setSupp, ncol(thisRHS)))
      allConf <- c(allConf, thisConf[thisConf >= minConf])
      allLift <- c(allLift, thisLift[thisConf >= minConf])
      step <- step + 1

    }
  }

  res <- new("rules",
             lhs     = new("itemMatrix",
                           data     = as(allLHS, "ngCMatrix"),
                           itemInfo = itemInfo),
             rhs     = new("itemMatrix",
                           data     = as(allRHS, "ngCMatrix"),
                           itemInfo = itemInfo),
             quality = data.frame(support    = allSupp,
                                  confidence = allConf,
                                  lift       = allLift,
                                  conviction = allConv),
             info    = info)

  res@info$confidence <- minConf

  return(res)

}


###############
### apriori ###
###############

apriori <- function(data, minSupp, minConf, cl) {

  dataName <- deparse(substitute(data))
  res <- aprioriRules(data, aprioriSets(data, minSupp, cl), minConf, cl)
  res@info$data <- dataName
  return(res)

}


###############
### inspect ###
###############

setGeneric("inspect",
           function(data, search, minSupp, minConf, minLift, minConv, sortBy) {

  standardGeneric("inspect")

})

setMethod("inspect", signature(data = "itemsets"),
          function(data, search, minSupp, minConf, minLift, minConv, sortBy) {

  indSets <- 1:ncol(data@items@data)
  labels <- data@items@itemInfo$labels
  supp <- data@quality$support

  if (hasArg(search)) {
    indItems <- grep(search, labels)
    indSets <- apply(data@items@data[indItems, , drop = F], 1, which)
    indSets <- unlist(indSets)
  }

  if (hasArg(minSupp)) {
    indSets <- indSets[supp[indSets] >= minSupp]
  }

  if (! length(indSets)) {
    return(invisible(NULL))
  }

  if (hasArg(sortBy) && sortBy == "support") {
    indSets <- indSets[order(supp[indSets], decreasing = T)]
  }

  res <- data.frame(matrix(nrow = length(indSets), ncol = 2))
  names(res) <- c("Set", "Support")

  for (i in 1:length(indSets)) {
    tmpSets <- labels[data@items@data[, indSets[i]]]
    res[i, 1] <- if (length(tmpSets)) paste(tmpSets, collapse = ", ") else ""
    res[i, 2] <- supp[indSets[i]]
  }

  return(res)

})

setMethod("inspect", signature(data = "rules"),
          function(data, search, minSupp, minConf, minLift, minConv, sortBy) {

  labels <- data@lhs@itemInfo$labels

  supp <- data@quality$support
  conf <- data@quality$confidence
  conv <- data@quality$conviction
  lift <- data@quality$lift

  indRules <- 1:ncol(data@lhs@data)

  if (hasArg(search)) {
    indItems <- grep(search, labels)
    indLHS <- apply(data@lhs@data[indItems, , drop = F], 1, which)
    indRHS <- apply(data@rhs@data[indItems, , drop = F], 1, which)
    indRules <- sort(unique(c(unlist(indLHS), unlist(indRHS))))
  }

  if (hasArg(minSupp)) indRules <- indRules[supp[indRules] >= minSupp]
  if (hasArg(minConf)) indRules <- indRules[conf[indRules] >= minConf]
  if (hasArg(minConv)) indRules <- indRules[conv[indRules] >= minConv]
  if (hasArg(minLift)) indRules <- indRules[lift[indRules] >= minLift]

  if (! length(indRules)) return(invisible(NULL))

  if (hasArg(sortBy)) {
    if (sortBy == "support") {
      indRules <- indRules[order(supp[indRules], decreasing = T)]
    } else if (sortBy == "confidence") {
      indRules <- indRules[order(conf[indRules], decreasing = T)]
    } else if (sortBy == "conviction") {
      indRules <- indRules[order(conv[indRules], decreasing = T)]
    } else if (sortBy == "lift") {
      indRules <- indRules[order(lift[indRules], decreasing = T)]
    }
  }

  res <- data.frame(matrix(nrow = length(indRules), ncol = 6))
  names(res) <- c("LHS", "RHS", "Support", "Confidence", "Lift",
                  "Conviction")

  for (i in 1:length(indRules)) {
    tmpLHS <- labels[data@lhs@data[, indRules[i]]]
    tmpRHS <- labels[data@rhs@data[, indRules[i]]]
    res[i, 1] <- if (length(tmpLHS)) paste(tmpLHS, collapse = ", ") else ""
    res[i, 2] <- if (length(tmpRHS)) paste(tmpRHS, collapse = ", ") else ""
    res[i, 3] <- supp[indRules[i]]
    res[i, 4] <- conf[indRules[i]]
    res[i, 5] <- lift[indRules[i]]
    res[i, 6] <- conv[indRules[i]]
  }

  return(res)

})
