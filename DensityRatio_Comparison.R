######
library(densratio)
library(randomForest)
library(dplyr)
library(tidyr)
library(ggplot2)


gen_shift_data <- function(p1 = list(mean = 0.65, sd = 1), p2 = list(mean = 0.5, sd = .95), n = 1000){
  D1 <- rnorm(n, mean = p1$mean, sd = p1$sd)
  D2 <- rnorm(n, mean = p2$mean, sd = p2$sd)
  #oracle_weights <- dnorm(seq(min(c(D1, D2)), max(c(D1, D2)), length.out = 750), mean = p2$mean, sd = p2$sd)/dnorm(seq(min(c(D1, D2)), max(c(D1, D2)), length.out = 750), mean = p1$mean, sd = p1$sd)
  oracle_weights <- dnorm(D1, mean = p2$mean, sd = p2$sd)/dnorm(D1, mean = p1$mean, sd = p1$sd)
  #return(list("D1" = D1, "D2" = D2, "OW" = data.frame("X" = seq(min(c(D1, D2)), max(c(D1, D2)), length.out = 750), "weights" = oracle_weights)))
  return(list("D1" = D1, "D2" = D2, "OW" = data.frame("X" = D1, "Oracle" = oracle_weights)))
  }
# Example
g_example <- gen_shift_data()
plot(g_example$OW %>% arrange(X), type = 'l')

# Quick DR/RF example
DR_RF_Compare <- function(...){
  # Generating the Data
    sim <- gen_shift_data(...)
  
  ## Random Forest Example
    resp <- c(rep(0, length(sim$D1)), rep(1, length(sim$D2)))
    design_mat <- data.frame(X = c(sim$D1, sim$D2))
    likelihood_rf <- randomForest(x = design_mat, y = resp, ntree = 500, nodesize = length(sim$D1)/3)
    #p_1 <- predict(likelihood_rf, data.frame(X = design_mat[1:length(sim$D1),]), type = "prob") 
    pi <- likelihood_rf$predicted[1:length(sim$D1)]
    weights_RF <- pi/(1-pi)
    
  ## Logistic Regression Example
    # likelihood_GLM <- glm(Y~., data = data.frame(design_mat, Y = resp), family = 'binomial')
    # weights_GLM <- fitted(likelihood_GLM)[1:length(sim$D1)]/(1 - fitted(likelihood_GLM)[1:length(sim$D1)])
    
  ## Density Ratio Example
    densratio_obj <- densratio(x = sim$D2, y = sim$D1, verbose = F, kernel_num = 100)
    weights_DR <- densratio_obj$compute_density_ratio(sim$D1)
    #pairs(data.frame(sim$OW, "RF_Wts" = weights_RF, "DR_Wts" = weights_DR, "GLM_Wts" = weights_GLM) %>% arrange(X)) 
    
  ## Plotting
    # df_plot <- data.frame(sim$OW,  "RF_Wts" = weights_RF, "DR_Wts" = weights_DR, "GLM_Wts" = weights_GLM) %>% 
    #   gather(Oracle:GLM_Wts, key = "Model", value = "weights")
    # plot0 <- ggplot(data = df_plot) + geom_line(aes(x = X, y = weights, col = Model, size = Model), alpha = .5) +
    #   scale_size_manual(values = c(1, 1, 1.5, .5)) + theme_bw()
    
    df_plot <- data.frame(sim$OW,  "RF_Weights" = weights_RF, "uLSIF_Weights" = weights_DR) %>% 
      gather(Oracle:uLSIF_Weights, key = "Model", value = "weights")
    plot0 <- ggplot(data = df_plot) + geom_line(aes(x = X, y = weights), alpha = 1, size = .85) + facet_wrap(Model~.) + theme_bw() +
      theme(strip.text = element_text(size = 13)) + ylab(expression(dP[2](X)/dP[1](X)))
  
    plot(plot0)
  ## Getting out the error metrics
    calc_RMSE <- function(wts){
      sqrt(mean((wts - sim$OW$Oracle)^2))
    }
    #RMSE_out <- lapply(list("RF" = weights_RF, "GLM"= weights_GLM, "DR" = weights_DR), calc_RMSE) %>% as.data.frame()
    RMSE_out <- lapply(list("RF" = weights_RF,  "DR" = weights_DR), calc_RMSE) %>% as.data.frame()
    print(RMSE_out)
    return(list("RMSEs" = RMSE_out, "plot" = plot0))
}

set.seed(1994)
DR_Compare_1 <- DR_RF_Compare(n = 1500, p1 = list(mean = 0, sd = 2.5))


pdf("C:/Users/drain/Box Sync/LANL/Simulations/dens_ratio_compare.pdf", width = 9, height = 3)
DR_Compare_1$plot
dev.off()
