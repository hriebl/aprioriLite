# this file contains (almost) all generics and methods. most are required for
# the arules classes, especially for validity checking

##################
### itemMatrix ###
##################

setMethod("dim", signature(x = "itemMatrix"), function(x) {
  rev(dim(x@data))
})

setGeneric("itemInfo", function(x) {
  standardGeneric("itemInfo")
})

setMethod("itemInfo", signature(x = "itemMatrix"), function(x) {
  x@itemInfo
})

setGeneric("itemsetInfo", function(x) {
  standardGeneric("itemsetInfo")
})

setMethod("itemsetInfo", signature(x = "itemMatrix"), function(x) {
  x@itemsetInfo
})

setGeneric("itemUnion", function(x, y) {
  standardGeneric("itemUnion")
})

setMethod("itemUnion", signature(x = "itemMatrix", y = "itemMatrix"),
          function(x, y) {

  if (length(x) != length(y)) stop("x and y do not have the same length!")
  x@data <- as(x@data + y@data, "ngCMatrix")
  x

})

setMethod("length", signature(x = "itemMatrix"), function(x) {
  nrow(x)
})

setGeneric("nitems", function(x) {
  standardGeneric("nitems")
})

setMethod("nitems", signature(x = "itemMatrix"), function(x) {
  ncol(x)
})

setGeneric("size", function(x) {
  standardGeneric("size")
})

setMethod("size", signature(x = "itemMatrix"), function(x) {
  colSums(x@data)
})


#################
### itemsets ####
#################

setMethod("length", signature(x = "itemsets"), function(x) {
  length(x@items)
})

setMethod("plot", signature(x = "itemsets"), function(x) {

  par(mfrow = c(1, 2))

  sizes <- size(x@items)
  hist(sizes, breaks = seq(min(sizes) - 0.5, max(sizes) + 0.5, by = 1),
       main = "Histogram of Set Sizes", xlab = "Set Sizes")

  hist(quality(x)$support, main = "Histogram of Support Values",
       xlab = "Support Values")

  abline(v = x@info$support, col = "red", lwd = 2)

  invisible(NULL)

})

setGeneric("quality", function(x) {
  standardGeneric("quality")
})

setMethod("quality", signature(x = "itemsets"), function(x) {
  x@quality
})

setMethod("show", signature(object = "itemsets"), function(object) {

  nSets <- length(object)
  cat(nSets, "itemsets. Use the inspect function to find out more\n")
  invisible(NULL)

})

setMethod("summary", signature(object = "itemsets"), function(object) {

  cat("\n")
  show(object)
  cat("\n")

  cat("Set Sizes:\n")
  print(c(table(size(object@items))))
  cat("\n")

  cat("Summary of the Support Measure:\n")
  print(summary(quality(object)$support))
  cat("\n")

  cat("Other Info:\n")
  print(data.frame(object@info), row.names = F)
  cat("\n")

  invisible(NULL)

})

setGeneric("tidLists", function(x) {
  standardGeneric("tidLists")
})

setMethod("tidLists", signature(x = "itemsets"), function(x) {
  x@tidLists
})


#############
### rules ###
#############

setGeneric("items", function(x) {
  standardGeneric("items")
})

setMethod("items", signature(x = "rules"), function(x) {
  itemUnion(x@lhs, x@rhs)
})

setMethod("plot", signature(x = "rules"), function(x) {

  par(mfrow = c(3, 2))

  sizes <- size(x@lhs)
  hist(sizes, breaks = seq(min(sizes) - 0.5, max(sizes) + 0.5, by = 1),
       main = "Histogram of LHS Set Sizes", xlab = "LHS Set Sizes")

  main <- "Histogram of Support Values"
  hist(x@quality$support, main = main, xlab = "Support Values")
  abline(v = x@info$support, col = "red", lwd = 2)

  main <- "Histogram of Confidence Values"
  hist(x@quality$confidence, main = main, xlab = "Confidence Values")
  abline(v = x@info$confidence, col = "red", lwd = 2)

  main <- "Histogram of Lift Values"
  hist(x@quality$lift, main = main, xlab = "Lift Values")
  # no minimum option here

  main <- "Histogram of Conviction Values"
  hist(x@quality$conviction, main = main, xlab = "Conviction Values")
  # no minimum option here

  invisible(NULL)

})

setMethod("show", signature(object = "rules"), function(object) {

  nRules <- length(object@lhs)
  cat(nRules, "rules. Use the inspect function to find out more\n")
  invisible(NULL)

})

setMethod("summary", signature(object = "rules"), function(object) {

  cat("\n")
  show(object)
  cat("\n")

  cat("LHS Set Sizes:\n")
  print(c(table(size(object@lhs))))
  cat("\n")

  cat("Summary of the Quality Measures:\n")
  print(summary(object@quality))
  cat("\n")

  cat("Other Info:\n")
  print(data.frame(object@info), row.names = F)
  cat("\n")

  invisible(NULL)

})
