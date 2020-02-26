library(rpart)
library(dplyr)
library(party)
library(ranger)
library(tidyr)
library(ggplot2)
library(densratio)
library(cowplot)

######################################################################
# NECESSARY FUNCTION SECTION
######################################################################

rdirichlet <- function(alpha){
  x <- sapply(alpha, FUN = function(x) rgamma(n = 1, shape = x ))
  x <- x/sum(x)
  x
}

library(randomForest)

run_wtd_rpart <- function(n, ntest, l, nsim = 100){
  run_f <- function(null){
    #X_tr <- data.frame(cbind(t(replicate(n, rdirichlet(l^(1:5)))), replicate(1,runif(n)))) %>% 
      #mutate(Y = 5*exp(2*X1*X2 + X3)+ 5*X4*ifelse(X5 > X6, 1, -.5) + rnorm(n, sd = sqrt(5)))
   # X_te <- data.frame(cbind(t(replicate(ntest, rdirichlet(l^(5:1)))), replicate(1,runif(ntest)))) %>% 
      #mutate(Y = 5*exp(2*X1*X2 + X3) + 5*X4*ifelse(X5 > X6, 1, -.5) + rnorm(ntest, sd = sqrt(5)))
    tr_dist <- list(mean = -4, sd = 3.5)
    te_dist <- list(mean = 2.5, sd = 1.5)
    X_tr <- data.frame(X = rnorm(n, mean = -4, sd = 3.5)) %>% mutate(Y = sin(X/pi) + rnorm(n, sd = .5))
    X_te <- data.frame(X = rnorm(ntest, mean = 2.5, sd = 1.5)) %>% mutate(Y = sin(X/pi) + rnorm(ntest, sd = .5))
    
    resp <- as.factor(c(rep(0, n), rep(1, ntest)))
    design_mat <- rbind(data.frame(X = X_tr[,-which(names(X_tr) == "Y")]), 
                        data.frame(X = X_te[,-which(names(X_te) == "Y")]))
    likelihood_rf <- randomForest(x = design_mat, y = resp, ntree = 250, mtry = 5, nodesize = 25, do.trace =F,
                                  proximity = F)
    p_1 <- predict(likelihood_rf, data.frame(X = X_tr[,-which(names(X_tr) == "Y")]), type = "prob") 
    print(p_1)
    pi <- p_1[,2]
    # Convert to probabilities
    prox_max <- (pi + 1/n)/(1-pi + 1/n)
    plot(data.frame(x = X_tr$X, y = prox_max) %>% arrange(x), type = 'l')
    
    t_w <- as.party(rpart(Y~., data = X_tr, control = rpart.control(minbucket = 1, cp = 0, 
                                                                    maxcompete = 0, maxsurrogate = 0, xval = 0), weights = prox_max))
    t_u <- as.party(rpart(Y~., data = X_tr, control = rpart.control(minbucket = 1, cp = 0, 
                                                                    maxcompete = 0, maxsurrogate = 0, xval = 0)))
    RMSE_w <- sqrt(mean((predict(t_w, newdata = X_te) - X_te$Y)^2))
    RMSE_u <- sqrt(mean((predict(t_u, newdata = X_te) - X_te$Y)^2))
    print(c("RMSE_w" = RMSE_w, "RMSE_u" = RMSE_u))
    return(c("RMSE_w" = RMSE_w, "RMSE_u" = RMSE_u))}
  
  out <- do.call(rbind, lapply(1:nsim, run_f)) %>% data.frame()
  print(out %>% summarise(mRMSE_w = mean(RMSE_w), mRMSE_u = mean(RMSE_u), w_better_prop = mean(RMSE_u > RMSE_w)))
  return(out)
}
r1 <- run_wtd_rpart(250, 250, 1, 1)
#r1 <- run_wtd_rpart(250, 250, 12.5, 500)

## Quick function to bag a rpart model
bag_rpart <- function(formula, data, weights, control, ntree, replace, k, verbose, weighted_samp, wtd_struct, vtry){
  resp <- get_all_vars(formula, data = data)[1] %>% unlist()
  out_name <- all.vars(formula)[1]
  # Storing weights for out of bag measurements
  wts_oob <-weights
  if(weighted_samp){
    samp_weights <- weights
  }
  else{
    samp_weights <- rep(1, nrow(data))
  }
  if(!wtd_struct){
    weights <- rep(1, nrow(data))
  }
  # Updating the formula for technical purposes
  form_base <- Reduce(paste, deparse(formula))
  form <- paste(form_base, "-w") %>% as.formula()
  #form <- as.formula(form_base)
  # Fitting a resampled tree
  get_bagged_tree <- function(iter){
    temp_ind <- sample(nrow(data), size = k, replace = replace, prob = samp_weights)
    out_of_bag <- !(1:nrow(data) %in% temp_ind)
    w <- weights[temp_ind]
    if(vtry < ncol(data) -1){
      features_exclude <- sample(names(data[,-which(names(data) == out_name)]), ncol(data)-vtry-1)
      f_temp <- Reduce(paste, deparse(form)) %>% paste(paste(features_exclude, collapse = "-"), sep = "-") %>% as.formula()
    }
    else{
      f_temp <- Reduce(paste, deparse(form)) %>% as.formula()
    }
    #bag <- data.frame(data[temp_ind,], w = w)
    bag <- data[temp_ind,]
    bag[["w"]] <- w
    tree <- rpart(formula = f_temp, data = bag, control = control, weights = w, model = F)
    #out_of_bag_err <- sum(weights[-temp_ind]*(out_of_bag[[resp_name]] - predict(t_w, out_of_bag))^2)/sum(weights[-temp_ind])
    if(verbose & iter %% (ntree/10) == 0 | iter == 1 & verbose ){
      cat("Building tree number", iter, "...\n")
    }
    return(list("model" = tree, "oob_inds" = out_of_bag))
  }
  # Running the model
  outputs <- lapply(1:ntree, FUN = get_bagged_tree)
  #outputs <- vector(mode = "list", length = ntree)
  #or(i in 1:ntree){
    #outputs[[i]] <- get_bagged_tree(iter = i)
  #}
  models <- lapply(outputs, FUN = function(x) x[["model"]])
  oob_matrix <- do.call(rbind, lapply(outputs, FUN = function(x) x[["oob_inds"]]))
  #print(oob_matrix)
  #print(colSums(oob_matrix))
  data <- data.frame(data, "w" = weights)
  predictions <- do.call(rbind, lapply(models, FUN = function(x) predict(x, newdata = data)))
  print(dim(predictions))
  
  oob_preds <- ((predictions*oob_matrix) %>% colSums())/colSums(oob_matrix)
  #print(predictions*oob_matrix)
   
   par(mfrow = c(1, 2))
   plot(predictions %>% colMeans(), resp, ylim = range(resp), xlim = range(resp), main = "Training Predictions")
   abline(a = 0, b = 1, col = 'deepskyblue', lty = 'dashed', lwd = 2.5)
   plot(oob_preds, resp, ylim = range(resp), xlim = range(resp), main = "Out of Bag Predictions")
   abline(a = 0, b = 1, col = 'deepskyblue', lty = 'dashed', lwd = 2.5)
   
  
  oob_err <- sqrt(sum(weights*(oob_preds - resp)^2)/sum(weights))
  #print(sum(wts_oob))
  #print(sqrt(mean((oob_preds - resp)^2)))
  #print(dim(oob_matrix))
  #print(colSums(oob_matrix))
  
  if(verbose){
    cat("Full model out of bag error =", oob_err, "\n")
  }
  return(list("models" = models, "oob_err" = oob_err))
}


# Main Fitting Function ---------------------------------------------------
wtd_rpart_forest <- function(formula, df.train, df.test = NULL, fit_likelihood = "RF",  n_min_ess = nrow(df.train)/4, tune_lambda = T,
                             lambda = if(!tune_lambda) 0.5 else 1, weights = rep(1, nrow(df.train)),
                             control, ntree, replace, k, verbose, weighted_samp, wtd_struct, vtry = ncol(df.train) -1){
  n <- nrow(df.train)
  n.test <- nrow(df.test)
  out_name <- all.vars(formula)[1]
  # Calculating the weights
  if(fit_likelihood == "RF"){
    if(verbose) cat("Fitting likelihood forest ... \n")
    resp <- as.factor(c(rep(0, n), rep(1, n.test)))
    #print(length(resp))
    design_mat <- rbind(df.train[,-which(names(df.train) == out_name)], 
                        df.test[,-which(names(df.test) == out_name)])
    likelihood_rf <- randomForest(x = design_mat, y = resp, ntree = ntree, proximity = F, nodesize = 1)
    p_1 <- predict(likelihood_rf, df.train[,-which(names(df.train) == out_name)], type = "prob") 
    pi <- p_1[,2]
    # Convert to probabilities
    weights <- (pi + 1/n)/(1-pi + 1/n)
    par(mfrow = c(1,1))
    #plot(density.default(weights^lambda/sum(weights^lambda)))


    
    ## Calculating effective sample sizes, and tuning lambda
    if(tune_lambda){
      n_ess_0 <- sum(weights)^2/sum(weights^2)
      #if(n_ess_0 >= n_min_ess) lambda <- 1
      if(FALSE) lambda <- 3
      else{
        ess_diff <- function(l) sum(weights^l)^2/sum(weights^(2*l)) - n_min_ess
        # Finding the 0's of this equation
        lambda <- uniroot(ess_diff, interval = c(0,10), extendInt = "downX")[["root"]]
      }
      n_ess <- sum(weights^lambda)^2/sum(weights^(2*lambda))
    }
  }
  else if(fit_likelihood == "DR"){
    densratio_obj <- densratio(x1 = df.test[,-which(names(df.test) == out_name)], 
                               x2 = df.train[,-which(names(df.train) == out_name)], verbose = verbose, kernel_num = 100,
                               sigma = exp(seq(-4, ncol(df.train)/3, length.out = 15)))
    weights <- densratio_obj$compute_density_ratio(df.train[,-which(names(df.train) == out_name)])
    
    if(tune_lambda){
      n_ess_0 <- sum(weights)^2/sum(weights^2)
      #print(n_ess_0)
      if(FALSE) lambda <- 3
      else{
        ess_diff <- function(l) sum(weights^l)^2/sum(weights^(2*l)) - n_min_ess
        # Finding the 0's of this equation
        lambda <- uniroot(ess_diff, interval = c(0, 10), extendInt = "downX")[["root"]]
      }
      n_ess <- sum(weights^lambda)^2/sum(weights^(2*lambda))
    }
  }
  else{
    weights <- weights
    n_ess <- sum(weights^lambda)^2/sum(weights^(2*lambda))
  }
  # Fitting the model
  if(verbose & fit_likelihood == "RF")  cat("Fitting likelihood bagged rpart ... \n")
  else if(verbose) cat("Fitting bagged rpart ... \n")
  rpart_mods <- bag_rpart(formula = formula, data = df.train, weights = weights^lambda, ntree = ntree, replace = replace, k = k, verbose = verbose,
                          weighted_samp = weighted_samp, control = control, wtd_struct = wtd_struct, vtry = vtry)
  # Recording the test error
  if(!is.null(df.test)){
    predictions <- do.call(rbind, lapply(rpart_mods[["models"]], FUN = function(x) predict(x, newdata = data.frame(df.test, w = rep(0, nrow(df.test))))))
    MSE_test <- mean((colMeans(predictions) - df.test[,which(names(df.test) == out_name)])^2)
  }Vi
  else{
    MSE_test <- NULL
    predictions <- NULL
  }
  return(list("Models" = rpart_mods, "Test MSE" = MSE_test, "Test Preds" = predictions, 
              "response" = out_name, "n_ess" = n_ess, "lambda" = lambda, "weights" = weights))
}


#run_wtd_rpart(n = n, ntest = n.test, l = l, nsim = 1)

# Weighted Quantile Regression --------------------------------------------

QRF_Predict <- function(obj, newdata, tr_data, quantiles = c(0.1, .5, .9), returnMean = F){
  
  # Making sure the model used can actually be used for quantile regression
  if(is.null(obj$Models)) stop("Use an object that has a Models attribute")
  obj$Models <- lapply(obj$Models$models, FUN = partykit::as.party)
  pred_quants <- data.frame(matrix(nrow = nrow(newdata), ncol = length(quantiles)))
  names(pred_quants) <- paste("Quantile_", as.character(quantiles), sep = "")
  B <- length(obj$Models)
  
  # Getting the terminal nodes of all training data
  #tr_data <- obj$OriginalData
  nodes_tr <- do.call(rbind,lapply(obj$Models, FUN = function(mod) predict(mod, newdata = tr_data, type = "node")))
  # isolating the repsonse
  y <- tr_data[,which(names(tr_data) == obj$response)]
  
  # Getting the terminal nodes of the testing data
  nodes_te <- do.call(rbind,lapply(obj$Models, FUN = function(mod) predict(mod, newdata = newdata, type = "node")))
  
  # For test observation X_idx, calculates tree weights for tree j
  get_tree_weights <- function(idx, tree){
    idx_node <- nodes_te[tree,idx]
    tree_nodes <- nodes_tr[tree,]
    wts <- (tree_nodes == idx_node)/sum(tree_nodes == idx_node)
    return(wts)
  }
  
  # For test observation X_idx, calculates random forest weights (average of tree weights)
  get_forest_weights <- function(idx){
    all_tree_weights <- do.call(rbind, lapply(as.list(1:B), FUN = function(j) get_tree_weights(idx, j)))
    rf_weights <- colMeans(all_tree_weights, na.rm = T)
    return(rf_weights)
  }
  
  get_quantiles <- function(idx){
    wts <- get_forest_weights(idx)
    temp <- data.frame(resp = y, wts = wts) %>% arrange(y) %>% mutate(sumwt = cumsum(wts))
    out <- unlist(sapply(quantiles, FUN = function(q) temp %>% filter(sumwt >= q) %>% head(1) %>% select(resp)))
    return(out)
  }
  for(i in 1:nrow(newdata)) pred_quants[i,] <- get_quantiles(i)
  
  if(returnMean){
    preds_te <- do.call(rbind,lapply(obj$Models, FUN = function(mod) predict(mod, newdata = newdata)))
    pred_quants$mean <- preds_te
  }
  return(pred_quants)
}


Q_ex <- QRF_Predict(obj = debug_wtd_rpart, newdata = X_te[1,], tr_data = X_tr)

tune_wtd_RF <- function(formula, df.train, df.test,  ntree = 500, replace = T, 
                        k = if(replace) nrow(df.train) else .63*nrow(df.train), verbose = T,
                        mtry_grid = expand.grid(mtry = seq(ncol(df.train)/3 + 1, ncol(df.train) -1, length.out = 5)), ...){
  models_wtd <- vector(mode = "list", length = nrow(mtry_grid))
  models_unwtd <- vector(mode = "list", length = nrow(mtry_grid))
  
  # Putting in only whole number mtry parameters
  mtry_grid <- ceiling(mtry_grid)
  
  # Isolating the response
  resp <- get_all_vars(formula, data = df.test)[1] %>% unlist()
  
  # fits models 
  f_wtd <- function(mtry){
    cat("Fitting weighted model, mtry = ", mtry, "\n")
    resp_name <- all.vars(formula)[1]
    rf_wtd <- wtd_rpart_forest(formula = formula, df.train = df.train, df.test = df.test, fit_likelihood = "DR", wtd_struct = T, tune_lambda = T, 
                               n_min_ess = nrow(df.train)/2, vtry = mtry,
                               control = rpart.control( maxsurrogate = 0, cp = 0.00, maxdepth = 30, xval = 0, maxcompete = 0),
                               k = k, replace = replace, ntree = ntree, verbose = verbose, weighted_samp = F)
    return(list("model" = rf_wtd, "OOB_err" = rf_wtd$Models$oob_err, "TestRMSE" = sqrt(rf_wtd$`Test MSE`)))
  }
  f_unwtd <- function(mtry){
    cat("Fitting unweighted model, mtry = ", mtry, "\n")
    rf_unwtd <- ranger(formula, data = df.train, mtry = mtry, num.trees = ntree, 
                       replace = replace, sample.fraction = k/nrow(df.train), quantreg = T)
    RMSE <- sqrt(mean((predict(rf_unwtd, data = df.test)[["predictions"]] - resp)^2))
    return(list("model" = rf_unwtd, "OOB_err" = sqrt(rf_unwtd$prediction.error), "TestRMSE" = RMSE))
  }
  
  get_mod <- function(obj) obj$mod
  get_err <- function(obj) obj$OOB_err
  get_RMSE <- function(obj) obj$TestRMSE
  
  models_wtd <- lapply(mtry_grid$mtry, f_wtd)
  models_unwtd <- lapply(mtry_grid$mtry, f_unwtd)
  
  names(models_wtd) <- paste("mtry_", mtry_grid$mtry, sep = "")
  names(models_unwtd) <- paste("mtry_", mtry_grid$mtry, sep = "")
  
  errs_wtd <- sapply(models_wtd, get_err)
  errs_unwtd <- sapply(models_unwtd, get_err)
  RMSE_wtd <- sapply(models_wtd, get_RMSE)
  RMSE_unwtd <- sapply(models_unwtd, get_RMSE)
  
  OOB_table <- data.frame(mtry_grid, "WtdOOB" = errs_wtd, "UnwtdOOB" = errs_unwtd, 
                          "WtdRMSE" = RMSE_wtd, "UnwtdRMSE" = RMSE_unwtd)
  mtry_wtd <- OOB_table %>% arrange(errs_wtd) %>% head(1) %>% select(mtry) %>% unlist()
  mtry_unwtd <- OOB_table %>% arrange(errs_unwtd) %>% head(1) %>% select(mtry) %>% unlist()
  
  wtd_ind <- which(mtry_grid$mtry == mtry_wtd)
  unwtd_ind <- which(mtry_grid$mtry == mtry_unwtd)
  
  out <- list("Table" = OOB_table, "WtdModel" = get_mod(models_wtd[[wtd_ind]]), "UnwtdModel" = get_mod(models_unwtd[[unwtd_ind]]))
  return(out)  
}
n <- 500
n.test <- 250
l <- 1.25
sigma <- 0.5
X_tr <- data.frame(cbind(t(replicate(n, rdirichlet(l^(1:5)))), replicate(25,runif(n)))) %>%
  mutate(Y = 5*exp(2*X1*X2 - X3) + rnorm(n, sd = sigma))
X_te <- data.frame(cbind(t(replicate(n.test, rdirichlet(l^(5:1)))), replicate(25,runif(n.test)))) %>%
  mutate(Y = 5*exp(2*X1*X2 - X3) + rnorm(n.test, sd = sigma))

mod0 <- wtd_rpart_forest(Y~., df.train = X_tr, df.test = X_te, ntree = 500, verbose = F, weighted_samp = F, wtd_struct = T,
                         replace = F, k = 0.6*n, control = rpart.control(cp = 0, maxsurrogate = 0))
tune1 <- tune_wtd_RF(Y~., df.train = X_tr, df.test = X_te, ntree = 1500)

quant_tunes <- QRF_Predict(tune1$WtdModel, newdata = X_te, tr_data = X_tr)




######################################################################
# SIMULATION SECTION
######################################################################

# A High Dimensional Simulation -------------------------------------------

HighDimSim <- function(n = 2000, n.test = 500, l = 1, sigma = .25, p_additional = 10, Yx = function(X) 5*sqrt(X[,1]) + 5*sqrt(X[,3]) ){
  
  X_tr <- data.frame(cbind(t(replicate(n, rdirichlet(l^(1:6)))), replicate(p_additional,runif(n)))) 
  Y_tr <- Yx(X_tr)
  X_tr <- X_tr %>% mutate(Y = Y_tr + rnorm(n, sd = sigma))
  X_te <- data.frame(cbind(t(replicate(n.test, rdirichlet(l^(6:1)))), replicate(p_additional,runif(n.test)))) 
  Y_te <- Yx(X_te)
  X_te <- X_te %>% mutate(Y =  Y_te + rnorm(n.test, sd = sigma))
  debug_wtd_rpart <- wtd_rpart_forest(Y~., df.train = X_tr, df.test = X_te, fit_likelihood = "DR", tune_lambda = T,
                                      control = rpart.control(maxsurrogate = 0, cp = 0.0, maxdepth = 30, xval = 0, maxcompete = 0), 
                                      ntree = 500, k = .6*n, replace = F, verbose = F, weighted_samp = F, wtd_struct = T, vtry = ncol(X_tr)-2)
  ranger_out <- ranger(Y~., data = X_tr, num.trees = 500, sample.fraction = .6, replace = F, splitrule = "variance", quantreg = T)
  preds_ranger <- predict(ranger_out, X_te)[["predictions"]]
  quantiles_wtd <- QRF_Predict(debug_wtd_rpart, newdata = X_te, tr_data = X_tr, quantiles = c(.1,.5,.9))
  quantiles_unwtd <- predict(ranger_out, X_te, type = "quantiles", quantiles = c(.1, .5, .9))[["predictions"]]
  
  
  ### RMSE Calculations
  RMSE_wtd <-  sqrt(debug_wtd_rpart$`Test MSE`)
  RMSE_unwtd <-  sqrt(mean((preds_ranger - X_te$Y)^2))
  ### MAE Calculations - make sure 0.5 quantile is in column 2
  MAE_wtd <- mean(abs(quantiles_wtd[,2] - X_te$Y))
  MAE_unwtd <- mean(abs(quantiles_unwtd[,2] - X_te$Y))
  ### Calculating coverage percentages
  covg_wtd <- mean(quantiles_wtd[,1] < X_te$Y & quantiles_wtd[,ncol(quantiles_wtd)] > X_te$Y)
  covg_unwtd <- mean(quantiles_unwtd[,1] < X_te$Y & quantiles_unwtd[,ncol(quantiles_unwtd)] > X_te$Y)
  ### Interval Widths 
  width_wtd <- mean(quantiles_wtd[,ncol(quantiles_wtd)] - quantiles_wtd[,1])
  width_unwtd <- mean(quantiles_unwtd[,ncol(quantiles_unwtd)] - quantiles_unwtd[,1])
  
  
  score_table <- data.frame("Model" = c("Weighted", "Unweighted"),
                            "RMSE" = c(RMSE_wtd, RMSE_unwtd),
                            "MAE" = c(MAE_wtd, MAE_unwtd),
                            "CovgPCT" = c(covg_wtd, covg_unwtd),
                            "IntervalWidth" = c(width_wtd, width_unwtd)) %>% 
    mutate(score = (1/MAE + 1/RMSE + 4/IntervalWidth)*(CovgPCT)/.9)
  #print(score_table)
  
  preds_wtd <- debug_wtd_rpart$`Test Preds` %>% colMeans()
  preds_unwtd <- preds_ranger
  par(mfrow = c(1,2))
  # Weighted Predictions
  plot(1, type = "n", xlim = c(min(X_te$Y), max(X_te$Y)), 
       ylim = c(min(X_te$Y), max(X_te$Y)), main = paste("Proximity Weighted Forest"),
       xlab = "Mean Estimate", ylab = "Truth")
  grid()
  arrows(quantiles_wtd[,1], X_te$Y, quantiles_wtd[,ncol(quantiles_wtd)], X_te$Y, length = .05,
         angle = 90, code = 3, col = 'gray')
  try(points(preds_wtd, X_te$Y,
             pch = 20, col = rgb(1, 0, 0, alpha = .5)))
  abline(a = 0, b = 1, col = 'deepskyblue', lty = 'dashed', lwd = 2.5)
  
  # Unweighted predictions
  plot(1, type = "n", xlim = c(min(X_te$Y), max(X_te$Y)), 
       ylim = c(min(X_te$Y), max(X_te$Y)), main = paste("Unweighted Forest"),
       xlab = "Mean Estimate", ylab = "Truth")
  grid()
  arrows(quantiles_unwtd[,1], X_te$Y, quantiles_unwtd[,ncol(quantiles_unwtd)], X_te$Y, length = .05,
         angle = 90, code = 3, col = 'gray')
  points(preds_unwtd, X_te$Y,
         pch = 20, col = rgb(1, 0, 0, alpha = .5))
  abline(a = 0, b = 1, col = 'deepskyblue', lty = 'dashed', lwd = 2.5)
  
  plot <- recordPlot()
  return(list("wtdpreds" = preds_wtd, "wtd_model" = debug_wtd_rpart, "unwtdpreds" = preds_unwtd, 
              "scores" = score_table, "plot" = plot))
}

## Example
H1 <- HighDimSim(l = 2, p_additional = 5, sigma = .5, n = 1000, n.test = 200, Yx = function(X) 5*sin(pi*X[,3]) + 5*sqrt(X[,8]*X[,9]))


### Running this for several l values and getting out plots 
set.seed(1994)
l_values <- as.list(c(1, 1.5, 2, 3.5, 5, 10))
sim_apply_f <- function(l, plot = F, nsim = 1){
  d0 <- data.frame(Model = rep(0, 0), RMSE = rep(0, 0), MAE = rep(0, 0), CovgPCT = rep(0, 0),
                   IntervalWidth = rep(0, 0), score = rep(0, 0))
  for(i in 1:nsim){
  obj <- HighDimSim(l = l, p_additional = 25, sigma = .5, n = 1000, n.test = 200)
  d0 <- rbind(d0, obj$scores)
  }
  if(plot & nsim == 1){
  f_name <- paste("~/Documents/LANL/Figures/HighDimSim_l=", l, ".pdf", sep = "")
  print(f_name)
  pdf(file = f_name, width = 12, height = 6)
  replayPlot(obj$plot)
  dev.off()
  }
  d0_wtd <- d0 %>% filter(Model == "Weighted") %>% dplyr::select(-Model) %>% colMeans()
  d0_unwtd <- d0 %>% filter(Model == "Unweighted") %>% dplyr::select(-Model) %>% colMeans()
  return(data.frame(Model = c("Weighted", "Unweighted"), rbind(d0_wtd, d0_unwtd)))
}
HighDimSim_multiple_l <- lapply(l_values, sim_apply_f, nsim = 50)
## Getting out some xtables 
xtable_gen <- function(obj){
  tab <- obj %>% as.data.frame() %>% select(-score)
  print(xtable(tab, digits = 3), include.rownames = F)
}
lapply(HighDimSim_multiple_l, xtable_gen)

library(xtable)
print(xtable(data.frame("Model" = names(out), "RMSE" = out), digits = 3), include.rownames = F)
f_name <- paste("C:/Users/drain/Box Sync/LANL/Figures/PairsPlot_l=", l, ".pdf", sep = "")
print(f_name)
pdf(file = f_name, width = 6, height = 4)
pairs(d0, xlim = c(min(d0$Truth), max(d0$Truth)), ylim = c(min(d0$Truth), max(d0$Truth)))
dev.off()


# Multiple Regression Functions -------------------------------------------

# Our outcome models
Yx_list <- list(
  "Model_1" = function(X) 5*X[,1], 
  "Model_2" = function(X) 5*sin(pi*X[,1]),
  "Model_3" = function(X) 10*sin(pi*X[,1]*X[,2]) + 20*(X[,3] - 0.5)^2 + 10*X[,4] + 5*X[,5], 
  "Model_4" = function(X) 5*exp(2*sqrt(X[,1]*X[,2]) + X[,6]),
  "Model_5" = function(X) 5*(X[,1]^2 + X[,2]^2 + X[,3]^2 + X[,4]^2 + X[,5]^2)
)
TestX <- replicate(5, replicate(10, runif(10)), simplify = F)
mapply_test <- mapply(Yx_list, TestX, FUN = function(f, x) f(x), SIMPLIFY = F)

# Our possible lambda values
lambdas <- list(1, 1.5, 2, 3, 5, 10)

# Function for running the simulation
run_sims <- function(mods = Yx_list, lambdas = list(1, 1.5, 2, 3, 5, 10), nsim = 10, ...){
  # d0 <- data.frame(Model = rep("0", 0), RMSE = rep(0, 0), MAE = rep(0, 0), CovgPCT = rep(0, 0),
  #                  IntervalWidth = rep(0, 0), score = rep(0, 0))
  # Declaring where the results will live
  out_df <- expand.grid(Yx = rep("0", 0), lambdas =rep(0, 0), Model = rep("0", 0), 
                        RMSE = rep(0, 0), MAE = rep(0, 0), CovgPCT = rep(0, 0), 
                        IntervalWidth = rep(0, 0), score = rep(0, 0))
  # Defining a simulation function for a given model/lambda value
  gen_sim <- function(Yx, lambda){
    cat("Running Simulation For Lambda =", lambda, "\n")
    d0 <- data.frame(Model = rep("0", 0), RMSE = rep(0, 0), MAE = rep(0, 0), CovgPCT = rep(0, 0),
                                      IntervalWidth = rep(0, 0), score = rep(0, 0))
    for(i in seq_len(nsim)){
      obj <- HighDimSim(l = lambda, Yx = Yx, ...)
      d0 <- rbind(d0, obj$scores)
    }
    d0_wtd <- d0 %>% filter(Model == "Weighted") %>% dplyr::select(-Model) %>% colMeans()
    d0_unwtd <- d0 %>% filter(Model == "Unweighted") %>% dplyr::select(-Model) %>% colMeans()
    return(data.frame(Model = c("Weighted", "Unweighted"), rbind(d0_wtd, d0_unwtd)))
  }
  # Outer Loop: Yx_list
  # Inner Loop: Lambdas
  # Running this with the above structure
  for(m in names(Yx_list)){
    cat("Model:", m, "\n")
    f <- Yx_list[[m]]
    outs <- lapply(lambdas, FUN = function(l){
      obj <- gen_sim(lambda = l, Yx = f)
      obj$lambdas <- l
      obj$Yx <- m
      return(obj)
      })
    outs <- do.call(rbind, outs)
    out_df <- rbind(out_df, outs)
    print(outs)
  }
  return(out_df)
}

sim_small <- run_sims(nsim = 2, n = 150, n.test = 50, lambdas = list(1.05, 1.35), sigma = .5, p_additional = 25)




# Running the Big Simulation ----------------------------------------------
set.seed(1995)
sim_main <- run_sims(nsim = 150, n = 1000, n.test = 200, lambdas = as.list(seq(1, 1.5, length.out = 8)),
                     sigma = .5, p_additional = 25)
write.csv(sim_main, file = "~/Documents/LANL/Simulations/SimRedux.csv")




### 1-D vis example
library(ggplot2)
library(cowplot)
n <- 500
nt <- 250
tr_dist <- list(mean = -4, sd = 3)
te_dist <- list(mean = 3.5, sd = 1.5)
Yx <- function(X) ifelse(exp(X)/(1 + exp(X))*sin(X) > exp(- X)/(1 + exp(1-X))*sin(-X), exp(X)/(1 + exp(X))*sin(X),  exp(- X)/(1 + exp(1-X))*sin(-X))
#X <- data.frame(X = rnorm(n, mean = 0, sd = 1.5), Xt = rnorm(nt, mean = 1, sd = .75)) %>% mutate(Y = sin(pi*X) + rnorm(n, sd = 0.35),
                                                                                             #    Yt = sin(pi*Xt) + rnorm(nt, sd = 0.35))
gen_MSE_comparison <- function(iter, returnPlot = F){
  X <- data.frame(X = c(rnorm(n, mean = tr_dist$mean, sd = tr_dist$sd), rnorm(nt, mean = te_dist$mean, sd = te_dist$sd)), source = c(rep("Training", n), rep("Test", nt))) %>% 
    mutate(Y = Yx(X) + rnorm(n + nt, sd = 0.25))
  #ggplot(data = X, mapping = aes(col = source)) + stat_function(fun = function(x) exp(x)/(1 + exp(x))*sin(x), col = 'black', alpha = .5) + geom_point(aes(x = X, y = Y), alpha = .7) 
  
  X_tr <- X %>% filter(source == "Training") %>% dplyr::select(-source)
  X_te <- X %>% filter(source == "Test") %>% dplyr::select(-source)
  rf0 <- ranger(Y~X, data = X_tr, num.trees = 500, min.node.size = 50)
  rf_wtd <- wtd_rpart_forest(Y~X, df.train = X_tr, df.test = X_te, fit_likelihood = "DR", tune_lambda = T, n_min_ess = n*.1,
                              control = rpart.control(minsplit = 1, maxsurrogate = 0, cp = 0.0, maxdepth = 30, xval = 0, maxcompete = 0), 
                              ntree = 500, k = .8*n, replace = F, verbose = F, weighted_samp = F, wtd_struct = T)
  rf_oracle <- wtd_rpart_forest(Y~X, df.train = X_tr, df.test = X_te, fit_likelihood = F, tune_lambda = T, n_min_ess = n*.1, 
                                weights = dnorm(X_tr$X, mean = te_dist$mean, sd = te_dist$sd)/ dnorm(X_tr$X, mean = tr_dist$mean, sd = tr_dist$sd),
                                control = rpart.control(minsplit = 1, maxsurrogate = 0, cp = 0.0, maxdepth = 30, xval = 0, maxcompete = 0), 
                                ntree = 500, k = .8*n, replace = F, verbose = F, weighted_samp = F, wtd_struct = T)
  n_grid <- max(n*.1, 100)
  x_grid <- data.frame(X = seq(min(X$X), max(X$X), length.out = n_grid))
  rf0_preds <- data.frame(X = c(rep(x_grid$X,2), rep(X_te$X, 2)),  Y = c(Yx(x_grid$X), unlist(predict(rf0, x_grid)[["predictions"]]), 
                                                      rf_wtd$`Test Preds` %>% colMeans(), rf_oracle$`Test Preds` %>% colMeans()), 
                          fun = c(rep("Truth", n_grid), rep("RF", n_grid), rep("Wtd RF", nt), rep("Oracle RF", nt)))
  print(c("Unwtd MSE" = mean((predict(rf0, X_te)[["predictions"]] - X_te$Y )^2) , 
          "Wtd MSE" = rf_wtd$`Test MSE`, "Oracle MSE" = rf_oracle$`Test MSE`))
  #print(c("Unwtd OOB" = rf0$prediction.error , "Wtd OOB" = rf_wtd$Models$oob_err^2, "Oracle OOB" = rf_oracle$Models$oob_err^2))
  
  theme_textsize <- theme(legend.title = element_text(size = 14), legend.text = element_text(size = 12.5),
                          axis.text = element_text(size = 12.5))
  g_reg <- ggplot() +  geom_point(data = X, aes(x = X, y = Y, fill = source), alpha = .25, size = 1.55, shape = 21) + 
    scale_color_manual(name = "E(Y|X = x)", values = c("slateblue3", "orchid3", "grey65", "firebrick1")) + 
    #labels = c("RF\nEstimate", "Truth", "Wtd\nRF", "Oracle\nWtd RF")) + 
    scale_shape_manual(name = "Data Source", values = c(4, 20)) + xlim(c(min(X$X), max(X$X))) +
    scale_fill_discrete(name = "Data\nSource")+
    #scale_linetype_manual(name = "E(Y|X = x)", values = c("dashed", "solid")) +
    #scale_size_manual(name = "E(Y|X = x)", values = c(1.15, .5))+
    geom_line(data = rf0_preds, aes(x = X, y = Y, col = fun), size = 1, alpha = .8) + theme_bw() +
    guides(fill = guide_legend(override.aes = list(size=3.5))) + theme_textsize
  #plot(g_reg)
  density_grid <- x_grid %>% mutate(Training = dnorm(X, mean = tr_dist$mean, sd = tr_dist$sd), 
                                    Test = dnorm(X, mean = te_dist$mean, sd = te_dist$sd)) %>% 
    gather(Training, Test, key = "Source", value = "Density") 
  g_density <- ggplot(density_grid) +#geom_line(aes(x = X, y = Density, col = Source)) +
    geom_ribbon(aes(x = X, ymin = 0, ymax = Density, fill = Source, col = Source), alpha = .5)+ xlim(c(min(X$X), max(X$X))) +
    theme_bw() + scale_color_discrete(name = "P(X)") + scale_fill_discrete(name = "P(X)") + ylab("") +
    theme_textsize
  shift_grid <- data.frame(X = X_tr$X, Estimated = rf_wtd$weights/sum(rf_wtd$weights)) %>% 
    mutate(Oracle = dnorm(X, mean = te_dist$mean, sd = te_dist$sd)/ dnorm(X, mean = tr_dist$mean, sd = tr_dist$sd)) %>% 
    mutate(Oracle = Oracle/sum(Oracle)) %>% gather(Estimated, Oracle, value = "Ratio", key = "WeightType")
  g_shift <- ggplot(shift_grid) + xlim(c(min(X$X), max(X$X)))+
    geom_ribbon(aes(x = X, ymin = 0, ymax = Ratio, fill= WeightType, col = WeightType), alpha = .5) + theme_bw() +
    geom_hline(aes(yintercept = 0), size = 1.15, col = 'grey20') + ylab("") + scale_color_discrete(name = "W(X)") +
    scale_fill_discrete(name = "W(X)") + theme_textsize

  #cat("lambda" = rf_wtd$lambda)
#plot_grid(g_reg, g_density, g_shift, ncol = 1, rel_heights = c(2, 1, 1)) %>% plot()

if(returnPlot){
return(list("metrics" = c("Unwtd MSE" = mean((predict(rf0, X_te)[["predictions"]] - X_te$Y )^2) , 
               "Wtd MSE" = rf_wtd$`Test MSE`, "Oracle MSE" = rf_oracle$`Test MSE`, 
         c("Unwtd OOB" = rf0$prediction.error , "Wtd OOB" = rf_wtd$Models$oob_err^2, "Oracle OOB" = rf_oracle$Models$oob_err^2)),
       "plot" = plot_grid(g_reg, g_density, g_shift, ncol = 1, rel_heights  = c(2, 1, 1)),
       "reg_plot" = g_reg,
       "weights_plot" = plot_grid( g_density, g_shift, ncol = 1, rel_widths = c( 1, 1))))
}
else{
  return(c("Unwtd MSE" = mean((predict(rf0, X_te)[["predictions"]] - X_te$Y )^2) , 
                            "Wtd MSE" = rf_wtd$`Test MSE`, "Oracle MSE" = rf_oracle$`Test MSE`, 
                            c("Unwtd OOB" = rf0$prediction.error , "Wtd OOB" = rf_wtd$Models$oob_err^2, "Oracle OOB" = rf_oracle$Models$oob_err^2)))
}
}

### Running the simulation for real
set.seed(1994)
outs <- sapply(1:150, gen_MSE_comparison, returnPlot = F)
out_metrics <- outs %>% rowMeans()
print(out_metrics)

#X <- rbind(X, rf0_preds)
g_reg <- ggplot() +  geom_point(data = X, aes(x = X, y = Y, fill = source), alpha = .25, size = 1.55, shape = 21) + 
  scale_color_manual(name = "E(Y|X = x)", values = c("slateblue3", "orchid3", "black", "firebrick1")) + 
                     #labels = c("RF\nEstimate", "Truth", "Wtd\nRF", "Oracle\nWtd RF")) + 
  scale_shape_manual(name = "Data Source", values = c(4, 20)) + xlim(c(min(X$X), max(X$X))) +
  scale_fill_discrete(name = "Data\nSource")+
  #scale_linetype_manual(name = "E(Y|X = x)", values = c("dashed", "solid")) +
  #scale_size_manual(name = "E(Y|X = x)", values = c(1.15, .5))+
  geom_line(data = rf0_preds, aes(x = X, y = Y, col = fun), size = 1, alpha = .8) + theme_bw() +
  guides(fill = guide_legend(override.aes = list(size=3.5)))
#plot(g_reg)
density_grid <- x_grid %>% mutate(Training = dnorm(X, mean = tr_dist$mean, sd = tr_dist$sd), 
                                  Test = dnorm(X, mean = te_dist$mean, sd = te_dist$sd)) %>% 
  gather(Training, Test, key = "Source", value = "Density") 
g_density <- ggplot(density_grid) +#geom_line(aes(x = X, y = Density, col = Source)) +
  geom_ribbon(aes(x = X, ymin = 0, ymax = Density, fill = Source, col = Source), alpha = .5)+ xlim(c(min(X$X), max(X$X))) +
  theme_bw() + scale_color_discrete(name = "P(X)") + scale_fill_discrete(name = "P(X)") + ylab("")
shift_grid <- data.frame(X = X_tr$X, Estimated = rf_wtd$weights/sum(rf_wtd$weights)) %>% 
  mutate(Oracle = dnorm(X, mean = te_dist$mean, sd = te_dist$sd)/ dnorm(X, mean = tr_dist$mean, sd = tr_dist$sd)) %>% 
  mutate(Oracle = Oracle/sum(Oracle)) %>% gather(Estimated, Oracle, value = "Ratio", key = "WeightType")
g_shift <- ggplot(shift_grid) + 
  geom_ribbon(aes(x = X, ymin = 0, ymax = Ratio, fill = WeightType, col = WeightType), alpha = .5) + theme_bw() +
  geom_hline(aes(yintercept = 0), size = 1.15, col = 'grey20') + ylab("W(X)")



set.seed(1995)
out_plot <- lapply(1, gen_MSE_comparison, returnPlot = T)

## All 3 plots
pdf("~/Documents/LANL/Hurricane_Codes/Figures/CovShift1.pdf", width = 10, height = 7)
#plot_grid(g_reg, g_density, g_shift, ncol = 1, rel_heights = c(2, 1, 1))
out_plot[[1]]$plot
dev.off()

## Just the regression plot
pdf("~/Documents/LANL/Hurricane_Codes/Figures/CovShift_reg_plot.pdf", width = 8, height = 6)
#plot_grid(g_reg, g_density, g_shift, ncol = 1, rel_heights = c(2, 1, 1))
out_plot[[1]]$reg_plot + xlim(c(-3, 8))
dev.off()

## Just the shift plot
pdf("~/Documents/LANL/Hurricane_Codes/Figures/CovShift_shift_plot.pdf", width = 12, height = 6)
#plot_grid(g_reg, g_density, g_shift, ncol = 1, rel_heights = c(2, 1, 1))
out_plot[[1]]$weights_plot
dev.off()

## Making a table
out <- out_plot[[1]]$metrics[1:3]
library(xtable)
print(xtable(data.frame("Model" = names(out), "RMSE" = out), digits = 4), include.rownames = F)

#print(c("Unwtd MSE" = mean((predict(rf0, X_te)[["predictions"]] - X_te$Y )^2) , "Wtd MSE" = rf_wtd$`Test MSE`, "Oracle MSE" = rf_oracle$`Test MSE`))


# Quick plot example
pdf("~/Documents/LANL/Hurricane_Codes/Figures/VPhiX.pdf", width = 10, height = 4)
plot(Yx, xlim = c(-15, 12), ylab = expression(varphi(X)), xlab = "X", n = 1000, col = 'darkslateblue', lwd = 2.5)
dev.off()

## Making an example ecdf plot
qrf_plot_ex <- function(){
  qrf <- quantregForest(x = replicate(3, rnorm(100)), y = rnorm(100), keep.inbag = T)
  plot_inds <- sample(100, 1, replace = F)
  cols <- c("slateblue", "forestgreen", "tomato")
  plot_0 <- predict(qrf, what = ecdf)[[1]] %>% plot(main = "Example QRF ECDF's", xlab = "y", ylab = as.expression("P(Y <= y | X = x)"))
  k <- 1
  for(p in plot_inds){
    predict(qrf, what = ecdf)[[p]] %>% plot(add = T, col = cols[k])
    k <- k +1
  }
}
set.seed(192)
qrf_plot_ex()

pdf("~/Documents/LANL/Hurricane_Codes/Figures/ECDFs.pdf", width = 6, height = 6)
qrf_plot_ex()
dev.off()

######################################################################
# APPLICATION TO HURRICANES 
######################################################################
load("C:/Users/drain/Box Sync/LANL/Hurricane_Codes/WeatherData/2dfFullHD_weather_noNA.rda")
#load("~/Documents/LANL/Hurricane_Codes/WeatherData/2dfFullHD_weather_noNA.rda")
fit_wtd_forest <- function(storm, prop_storm_to_keep = 0, tune_only = F){
  ### ARGUMENTS:
  # storm : A string corresponding to the hurricane to be held out
  # prop_storm_to_keep : A numeric specifying how much of the storm can be made into training data - 0 is none, 1 is no split.
  
  
  # Creating a training data set for later on
  FHD_train <- FHD2_no_NA_w %>% filter(hurr != storm) %>% 
    mutate(logcust = log(customer, 10)) %>% select(-c(customer, FIPS, hurr))
  # Creating a test set
  FHD_test <- FHD2_no_NA_w %>% filter(hurr == storm) %>% mutate(logcust = log(customer, 10))  %>% 
    select(-c(customer, FIPS, hurr))
  FHD_log_out <- FHD2_no_NA_w %>% filter(hurr == storm) %>% 
    mutate(logcust = log(customer, 10)) %>% select(logcust) %>% unlist()
  
  if(prop_storm_to_keep > 0){
    train_from_storm_inds <- sample(nrow(FHD_test), prop_storm_to_keep*nrow(FHD_test))
    cat("Keeping", ceiling(prop_storm_to_keep*nrow(FHD_test)), "from ", storm, "in training data\n")
    FHD_train_from_storm <- FHD_test[train_from_storm_inds,]
    FHD_test <- FHD_test[-train_from_storm_inds,]
    FHD_train <- rbind(FHD_train, FHD_train_from_storm)
    FHD_log_out <- FHD_test %>% select(logcust) %>% unlist()
  }
  print(dim(FHD_train))
  print(dim(FHD_test))
  
  if(!tune_only){
  
  rf_baseline <- ranger(logcust~., data = FHD_train, splitrule = "variance", num.trees = 500,
                                replace = F, sample.fraction = .6, mtry = 50, min.node.size = 5, quantreg = T)
  cat("Unweighted model out of bag error = ", sqrt(rf_baseline$prediction.error), "\n")
  # rf_baseline <- wtd_rpart_forest(logcust~., df.train = FHD_train, df.test = FHD_test, fit_likelihood = F, lambda = .45,
  #                            control = rpart.control(minsplit = 1, maxsurrogate = 0, cp = 0.000, maxdepth = 30, xval = 0, maxcompete = 0), 
  #                            ntree = 250, k = .6*nrow(FHD_train), replace = T, verbose = T, weighted_samp = F, wtd_struct = T, vtry = 25)
  
   rf_wtd <- wtd_rpart_forest(logcust~., df.train = FHD_train, df.test = FHD_test, fit_likelihood = "DR", n_min_ess = .75*nrow(FHD_train),
            control = rpart.control(maxsurrogate = 0, cp = 0.000, maxdepth = 30, xval = 0, maxcompete = 0),
           ntree = 500, k = .6*nrow(FHD_train), replace = F, verbose = T, weighted_samp = F, wtd_struct = T, vtry = 73)
  }
  else{
  tune_wtd <- tune_wtd_RF(logcust~., df.train = FHD_train, df.test = FHD_test, fit_likelihood = "DR", 
                          control = rpart.control(maxsurrogate = 0, cp = 0.000, xval = 0, maxcompete = 0), 
                         ntree = 500, k = .6*nrow(FHD_train), replace = F, verbose = T, weighted_samp = F, wtd_struct = T)
  rf_wtd <- tune_wtd$WtdModel
  rf_baseline <- tune_wtd$UnwtdModel
    print(tune_wtd$Table)
  }
   #print(c("Baseline:" = RMSE_baseline, "Wtd" = rf_wtd$RMSE))
  
  preds_wtd <- rf_wtd[["Test Preds"]] %>% colMeans() 
  preds_unwtd <- predict(rf_baseline, FHD_test)[["predictions"]]
  #preds_unwtd <- rf_baseline[["Test Preds"]] %>% colMeans()
  #preds_unwtdstruct_wtdpred <- unwtdtree_weighted_preds(tree_unwtd, weights = prox_max, newdata = FHD_test)
  
  ### RMSE Calculations
  RMSE_wtd <- sqrt(rf_wtd[["Test MSE"]])
  RMSE_unwtd <- sqrt(mean((preds_unwtd - FHD_test$logcust)^2))
  if(!tune_only){
  ### Getting quantile estimates
  quantiles_wtd <- QRF_Predict(rf_wtd, newdata = FHD_test, tr_data = FHD_train)
  quantiles_unwtd <- predict(rf_baseline, FHD_test, type = "quantiles")[["predictions"]]
  ### MAE Calculations - make sure 0.5 quantile is in column 2
  MAE_wtd <- mean(abs(quantiles_wtd[,2] - FHD_log_out))
  MAE_unwtd <- mean(abs(quantiles_unwtd[,2] - FHD_log_out))
  ### Calculating coverage percentages
  covg_wtd <- mean(quantiles_wtd[,1] < FHD_log_out & quantiles_wtd[,ncol(quantiles_wtd)] > FHD_log_out)
  covg_unwtd <- mean(quantiles_unwtd[,1] < FHD_log_out & quantiles_unwtd[,ncol(quantiles_unwtd)] > FHD_log_out)
  ### Interval Widths 
  width_wtd <- mean(quantiles_wtd[,ncol(quantiles_wtd)] - quantiles_wtd[,1])
  width_unwtd <- mean(quantiles_unwtd[,ncol(quantiles_unwtd)] - quantiles_unwtd[,1])
  
  
  score_table <- data.frame("Model" = c("Weighted", "Unweighted"),
                            "RMSE" = c(RMSE_wtd, RMSE_unwtd),
                            "MAE" = c(MAE_wtd, MAE_unwtd),
                            "CovgPCT" = c(covg_wtd, covg_unwtd),
                            "IntervalWidth" = c(width_wtd, width_unwtd)) %>% 
    mutate(score = (1/MAE + 1/RMSE + 4/IntervalWidth)*(CovgPCT)/.9)
  print(score_table)
  
  par(mfrow = c(1,2))
  # Weighted Predictions
  plot(1, type = "n", xlim = c(min(FHD_log_out), max(FHD_log_out)), 
       ylim = c(min(FHD_log_out), max(FHD_log_out)), main = paste("Weighted Forest,\n", storm),
       xlab = "Mean Estimate", ylab = "Truth")
  grid()
  arrows(quantiles_wtd[,1], FHD_log_out, quantiles_wtd[,ncol(quantiles_wtd)], FHD_log_out, length = .05,
         angle = 90, code = 3, col = 'gray')
  try(points(preds_wtd, FHD_log_out,
             pch = 20, col = rgb(1, 0, 0, alpha = .5)))
  abline(a = 0, b = 1, col = 'deepskyblue', lty = 'dashed', lwd = 2.5)
  
  # Unweighted predictions
  plot(1, type = "n", xlim = c(min(FHD_log_out), max(FHD_log_out)), 
       ylim = c(min(FHD_log_out), max(FHD_log_out)), main = paste("Unweighted Forest,\n", storm),
       xlab = "Mean Estimate", ylab = "Truth")
  grid()
  arrows(quantiles_unwtd[,1], FHD_log_out, quantiles_unwtd[,ncol(quantiles_unwtd)], FHD_log_out, length = .05,
         angle = 90, code = 3, col = 'gray')
  points(preds_unwtd, FHD_log_out,
         pch = 20, col = rgb(1, 0, 0, alpha = .5))
  abline(a = 0, b = 1, col = 'deepskyblue', lty = 'dashed', lwd = 2.5)
  return(list("wtdpreds" = preds_wtd, "wtd_model" = rf_wtd, "unwtdpreds" = preds_unwtd, 
              "scores" = score_table, "storm" = storm))
  }
  else{
    return(list("wtdpreds" = preds_wtd, "wtd_model" = rf_wtd, "unwtdpreds" = preds_unwtd, 
                 "storm" = storm, "TuneTable"  = tune_wtd$Table))
  }
}

# Running for real
set.seed(1995)
fitted_Harvey <- fit_wtd_forest("Harvey-2017", 0)
fitted_Irma <- fit_wtd_forest("Irma-2017", 0)
fitted_Nate <- fit_wtd_forest("Nate-2017", 0)
fitted_Sandy <- fit_wtd_forest("Sandy-2012")
fitted_Matthew <- fit_wtd_forest("Matthew-2016", 0)
fitted_Arthur <- fit_wtd_forest("Arthur-2014", 0)

HaIrSa <- list(fitted_Harvey, fitted_Irma, fitted_Sandy, fitted_Nate, fitted_Matthew, fitted_Arthur)
names(HaIrSa) <- c("Harvey-2017", "Irma-2017", "Sandy-2012", "Nate-2017", "Matthew-2016", "Arthur-2014")

xtable_gen <- function(obj){
  tab <- obj$scores %>% as.data.frame()
  tab$lambda <- obj$wtd_model$lambda
  tab$storm <- obj$storm
  #print(xtable(tab, digits = 3), include.rownames = F)
  return(tab)
}
tabs <- do.call(rbind, lapply(HaIrSa, xtable_gen))
#tabs$storm <- rep(c("Harvey-2017", "Irma-2017", "Sandy-2012", "Nate-2017", "Matthew-2016", "Arthur-2014"), times = rep(2, 6))

tabs <- tabs %>% dplyr::select(storm, Model, lambda, RMSE, MAE, CovgPCT, IntervalWidth, score)
names(tabs) <- c("Storm", "Model", "$\\lambda$", "RMSE", "MAE", "Covg", "Interval Width", "Score")

# Outputting this as a table
xtabs <- xtable(tabs, digits = 4, align = "cc|c|c|ccccc")
caption(xtabs) <- "Model performance by storm, with weighted and unweighted storms fitted. Bolded values represent the better of the two by storm and loss function."
label(xtabs) <- "tab:WeightedStorms"
print(xtabs, include.rownames = F, 
      hline.after = c(-1,0, seq(2, nrow(tabs), by = 2)),
      size = "footnotesize", table.placement = "H",
      sanitize.text.function=function(x){x})


# Running the oob section
set.seed(1995)
fitted_Harvey_tune <-  fit_wtd_forest("Harvey-2017", 0, tune_only = T)
fitted_Irma_tune <- fit_wtd_forest("Irma-2017", 0, tune_only = T)
fitted_Sandy_tune <- fit_wtd_forest("Sandy-2012", 0, tune_only = T)
fitted_Nate_tune <- fit_wtd_forest("Nate-2017", 0, tune_only = T)
fitted_Matthew_tune <- fit_wtd_forest("Matthew-2016", 0, tune_only = T)
fitted_Arthur_tune <- fit_wtd_forest("Arthur-2014", 0, tune_only = T)

StormTunes <- list("Harvey-2017" = fitted_Harvey_tune$TuneTable, "Irma-2017" = fitted_Irma_tune$TuneTable,
                   "Sandy-2012" = fitted_Sandy_tune$TuneTable, "Nate-2017" = fitted_Nate_tune$TuneTable,
                   "Matthew-2016" = fitted_Matthew_tune$TuneTable, "Arthur-2014" = fitted_Arthur_tune$TuneTable)
for(i in 1:length(StormTunes)){
  StormTunes[[i]]$storm <- names(StormTunes)[i]
}
#save(StormTunes, file = "C:/Users/drain/Box Sync/LANL/Simulations/StormTuneDFs.rda")
load("C:/Users/drain/Box Sync/LANL/Simulations/StormTuneDFs.rda")
# Getting out the tables
storm_OOB_table <- do.call(rbind, StormTunes) %>% gather(WtdOOB, UnwtdOOB, key = "Model_OOB", value = "OOB") %>% 
  gather(WtdRMSE, UnwtdRMSE, key = "Model_RMSE", value = "RMSE") %>% 
  filter((Model_OOB == "WtdOOB" & Model_RMSE == "WtdRMSE") | (Model_OOB == "UnwtdOOB" & Model_RMSE == "UnwtdRMSE"))

Wtd_OOB_table <- storm_OOB_table %>% filter(Model_OOB == "WtdOOB")
Unwtd_OOB_table <- storm_OOB_table %>% filter(Model_OOB == "UnwtdOOB")

storm_OOB_table <- do.call(rbind, StormTunes) 
names(storm_OOB_table) <- c("mtry", "Weighted OOB", "Unweighted OOB", "Weighted RMSE", "Unweighted RMSE", "storm")
storm_OOB_mtry_tab <- storm_OOB_table %>% gather(`Weighted OOB`:`Unweighted OOB`, key = "Estimate", value = "Error") %>% 
  mutate(ErrorType = rep(c("OOB", "RMSE"), times = c(60, 60)))

## Calculating rank correlations
Wtd_cor <- Wtd_OOB_table %>% group_by(storm) %>% summarise(RankCor = cor(RMSE, OOB, method = "spearman"))
Unwtd_cor <-Unwtd_OOB_table %>% group_by(storm) %>% summarise(RankCor = cor(RMSE, OOB, method = "spearman")) 

# Making a plot
## mtry vs error plot
ggplot(data = storm_OOB_mtry_tab) + geom_line(aes(x = mtry, y = Error, col = Estimate, lty = ErrorType), size = 1.15) +
  geom_point(aes(x = mtry, y = Error, col = Estimate), size = 2.5) +   facet_wrap(storm~., scales = "free") + theme_bw() +
  theme(strip.text = element_text(size = 12), legend.text = element_text(size = 12),
                        legend.title = element_text(size = 13), legend.spacing.y = unit(.01, "npc"),
                        legend.background = element_rect(colour = 'black', size = .5))


## Scaterplot of the OOB/RMSE
pdf("C:/Users/drain/Box Sync/LANL/Hurricane_Codes/Figures/Weighted_OOB_scatter.pdf", width = 8, height = 4)
ggplot(data = Wtd_OOB_table) + geom_point(aes(x = OOB, y = RMSE)) + 
  facet_wrap(storm~., scales = "free") + ggtitle("Weighted Out of Bag Error versus Weighted RMSE")  + theme_bw() + 
  theme(strip.text = element_text(size = 12), legend.text = element_text(size = 12),
        legend.title = element_text(size = 13), legend.spacing.y = unit(.01, "npc"),
        legend.background = element_rect(colour = 'black', size = .5), plot.title = element_text(hjust = 0.5))
dev.off()
pdf("C:/Users/drain/Box Sync/LANL/Hurricane_Codes/Figures/Unweighted_OOB_scatter.pdf", width = 8, height = 4)
ggplot(data = Unwtd_OOB_table) + geom_point(aes(x = OOB, y = RMSE)) + facet_wrap(storm~., scales = "free") + 
  ggtitle("Unweighted Out of Bag Error versus Unweighted RMSE")  + theme_bw() + 
  theme(strip.text = element_text(size = 12), legend.text = element_text(size = 12),
        legend.title = element_text(size = 13), legend.spacing.y = unit(.01, "npc"),
        legend.background = element_rect(colour = 'black', size = .5), plot.title = element_text(hjust = 0.5))
dev.off()


### Analyzing some other aspects of this study - weights, lambda, and effective sample size
library(tidyr)
library(tibble)
wts_df <- do.call(rbind, lapply(HaIrSa, FUN = function(obj) 
  return(data.frame("Storm" = obj$storm, "Weights" = obj$wtd_model$weights^obj$wtd_model$lambda, "Lambda" = obj$wtd_model$lambda))))

  ggplot() + geom_density(data = wts_df, aes(x = Weights, y = ..density..), fill = 'tomato', alpha = .5) +
    facet_wrap(Storm~., scales = "free") + theme_bw()
  
  
### Running for various proportions

  
  prop_list <- c(.2, .5, .8)
  MC_storms <- function(nsim = 10){
    fit_f <- function(storm, p){
      tabs <- data.frame(matrix(nrow = 0, ncol = 8))
      names(tabs) <- c("Storm", "Model", "lambda", "RMSE", "MAE", "Covg", "Interval Width", "Score")
      for(iter in 1:nsim){
        cat("Fitting iteration ", iter, "for storm", storm, "; prop =", p, "\n")
        out <-  fit_wtd_forest(storm = storm, p = p)
        tab <- out$scores
        tab$lambda <- out$wtd_model$lambda
        tab$storm <- storm
        tab$p <- p
        tabs <- rbind(tabs, tab)
        print(tabs)
      }
      out <- tabs %>% group_by(Model) %>% summarise(p = mean(p), lambda = mean(lambda), RMSE = mean(RMSE), MAE = mean(MAE), Covg = mean(CovgPCT), 
                                                    `Interval Width` = mean(IntervalWidth), Score = mean(score)) %>% mutate(Storm = storm)
      return(out)
    }
    grids <- expand.grid(storm = c("Harvey-2017", "Irma-2017", "Sandy-2012", "Nate-2017", "Matthew-2016", "Arthur-2014"), p = prop_list)
    #grids <- expand.grid(storm = c("Harvey-2017", "Irma-2017"), p = prop_list)
    grids <- grids %>% split(seq(1, nrow(grids), by = 1))
    out <- lapply(grids, FUN = function(l) fit_f(l$storm, l$p))
    return(out)
  }

  set.seed(1995)
  prop_sim <- MC_storms()

# Costly to run this, so now we save and load it
save(prop_sim, file = "C:/Users/drain/Box Sync/LANL/Hurricane_Codes/StormPropResults.rda")
load("C:/Users/drain/Box Sync/LANL/Hurricane_Codes/StormPropResults.rda")

## Plotting these results
prop_table <- do.call(rbind, prop_sim)

prop_plot <- ggplot(prop_table) + geom_line(aes(x = p, y = Score, col = Model)) + geom_point(aes(x = p, y = Score, col = Model)) +
  xlab("Proportion of test set included") +
  facet_wrap(Storm~.) + theme_bw() +  theme(strip.text = element_text(size = 12), legend.text = element_text(size = 12),
                                            legend.title = element_text(size = 13), legend.spacing.y = unit(.01, "npc"),
                                            legend.background = element_rect(colour = 'black', size = .5))
prop_plot

pdf("C:/Users/drain/Box Sync/LANL/Simulations/Score_Plot_StormProps.pdf", width = 8, height = 4)
prop_plot
dev.off()
