library(partykit)

runAll <- F

find_CART_split <-  function(response, data, weights, mtry = ncol(data), speed_control = min(125, nrow(data))){
  # speed control - randomly subset number of points to evaluate the CART criterion, for computational efficency. 
  y <- data[[response]]
  # Calculating within node variance
  var_Parent <- mean(y - sum(y*weights/sum(weights))^2)
  #var_Parent <- sum((y - sum(y*weights/sum(weights)))^2*weights/sum(weights))
  
  # Defining a function to evaluate the CART criterion at a point
  CART_evaluate <- function(j, z){
    left_ind <- which(data[,j] < z)
    right_ind <- which(data[,j] >= z)
    weights_l <- weights[left_ind]/sum(weights[left_ind])
    weights_r <- weights[right_ind]/sum(weights[right_ind])
    ytilde_l <- sum(weights_l*y[left_ind])
    ytilde_r <- sum(weights_r*y[right_ind])
    #CART_crit <- var_Parent - mean(y - ytilde_l*(data[,j] < z) - ytilde_r*(data[,j] >= z)^2)
    CART_crit <- var_Parent - sum((y - ytilde_l*(data[,j] < z) - ytilde_r*(data[,j] >= z))^2*weights/sum(weights))
    #CART_crit <- var_Parent - sum((y - ytilde_l)^2*weights_l/sum(weights_l)) - sum((y - ytilde_r)^2*weights_r/sum(weights_r))
    return(CART_crit)
  }
  
  # Randomly selecting indices upon which to make splits
  ind_splits <- sample(which(names(data) != response), mtry, replace = F)
  
  CART_matrix <- matrix(ncol = 3, nrow = 0)
  # Evaluating the criterion at each point
  for(j in ind_splits){
    split_pts <- sample(unique(data[,j]), min(length(unique(data[,j])), speed_control), replace = F)
    for(pt in split_pts){
      CART_val <- CART_evaluate(j, pt)
      CART_matrix <- rbind(CART_matrix, c(j, pt, CART_val))
    }
  }
  split_pt <- CART_matrix[which.max(CART_matrix[,3]),]
  ## split into two groups 
  #splitindex <- !(levels(data[,j]) %in% splitpoint)
  splitindex <- !(data[,j] < split_pt[2])
  #splitindex[!(levels(data[[xselect]]) %in% lev)] <- NA_integer_
  splitindex <- splitindex - min(splitindex, na.rm = TRUE) + 1L
  #print(splitindex)
  
  return(partysplit(varid = as.integer(split_pt[1]),
                    breaks = c(split_pt[2]),
                    index = 2:1,
                    info = list(CART = split_pt[3])))
}

## Testing the CART evaluation function
if(runAll){
weights <- rgamma(10000, 3)
weights <- weights/sum(weights)
data <- data.frame(matrix(rnorm(100000), ncol = 10), "y" = runif(10000))
CART_evaluate(3, .5)

(split1 <- find_CART_split("y", data, weights, mtry = 1, speed_control = 1000))
kidids_split(split1, data = data)
}


## Growing a weighted tree

grow_wtd_tree <- function(id = 1L, response, data, weights, 
                          minbucket = round(minsplit/3), minsplit = 40, mtry = ncol(data)-1,
                          maxnodes = 45){
  
  #if(nrow(data) < minsplit){ cat("Splitting stopped (minsplit)\n"); return(partynode(id = id))}
  if(nrow(data) < minsplit){  
    return(partynode(id = id))
  }
  ## find best split
  sp <- find_CART_split(response, data, weights, mtry = mtry)
  ## no split found, stop here
  if (is.null(sp)){
    return(partynode(id = id))
  }
  #if (sp$info$CART < cp) { cat("Splitting stopped\n"); return(partynode(id = id))} <- this was meant to add in a control parameter
  if (id > maxnodes){
    return(partynode(id = id))
  }
  
  ##  split the data
  kidids <- kidids_split(sp, data = data)
  min_node <- min(table(kidids))
  #if(min_node < minbucket){ cat("Splitting stopped (minbucket)\n"); return(partynode(id = id))}
  if(min_node < minbucket){ return(partynode(id = id))}
  ## set up all children nodes
  kids <- vector(mode = "list", length = max(kidids, na.rm = TRUE))
  for (kidid in 1:length(kids)){
    ## select observations for current node
    w <- weights[kidids == kidid]
    #print(sum(kidids == kidid))
    data_temp <- data[kidids == kidid,]
    #resp_temp <- response[kidids==kidid]
    ## get next node id
    if (kidid > 1) {
      myid <- max(nodeids(kids[[kidid - 1]]))
    } else {
      myid <- id
    }
    ## start recursion on this child node
    kids[[kidid]] <- grow_wtd_tree(id = as.integer(myid + 1), 
                                   response = response, data = data_temp, 
                                   weights = w, mtry = mtry, maxnodes = maxnodes)
  }
  ## return nodes
  return(partynode(id = as.integer(id), split = sp, kids = kids,
                   info = list(CART = min(info_split(sp)$CART, na.rm = TRUE))))
}

if(runAll){
tree1 <- grow_wtd_tree(response = "y", data = data, weights = rep(1, nrow(data)), minbucket = 30, minsplit = 30,
                       maxnodes = 15)

fitted_node(tree1, data = data)
}
## Formula interface for fitting a tree
fit_wtd_tree <- function(formula, data, weights = NULL, minsplit = 10, mtry = ncol(data) - 1, maxnodes = 550,
                         minbucket = round(minsplit/3)) {
  if(minbucket < 2) stop("Minbucket > 1 is necessary to ensure proper weighting calculations")
  old_d <- data
  ## name of the response variable
  respname <- all.vars(formula)[1]
  ## data without missing values, response comes last
  data <- data[complete.cases(data), c(attr(terms(formula, data = data), "term.labels"), respname)]
  #data <- data[complete.cases(data), -which(names(data) == respname)]
  if (is.null(weights)) weights <- rep(1L, nrow(data))
  weights <- weights/sum(weights)
  ## grow tree
  nodes <- grow_wtd_tree(id = 1L, response = respname, data = data, 
                         weights = weights, minbucket = minbucket, mtry = mtry,
                         maxnodes = maxnodes)
  ## compute terminal node number for each observation
  fitted <- fitted_node(nodes, data = old_d)
  ## return  constparty object
  ret <- party(nodes, data = old_d,
               fitted = data.frame("(fitted)" = fitted,
                                   "(response)" = data[[respname]],
                                   "(weights)" = weights,
                                   check.names = FALSE),
               terms = terms(formula, data = data))
  as.constparty(ret)
}

if(runAll){
library(dplyr)
N <- 1550
Ncol <- 10
d1 <- data.frame(matrix(rnorm(Ncol*N), ncol = Ncol)) %>% mutate(y = sin(pi*X1*X2) + rnorm(N, sd = .05))
weights <- rgamma(N, shape = 5)
(m1 <- fit_wtd_tree(y~., data = d1, mtry = Ncol/3, weights = weights, minbucket = 2, maxnodes = 150, minsplit = 1))
plot(m1)
mrpart <- as.party(rpart(y~., data = d1, control = rpart.control(minbucket = 150, cp = 0, minsplit =0, maxdepth = 25)))

plot(predict(m1, data = d1), d1$y)
plot(predict(mrpart, data = d1), d1$y)
}
