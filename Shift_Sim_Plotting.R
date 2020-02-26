######################################################
# Purpose of this file is to generate tables and plots
######################################################
library(ggplot2)
library(xtable)
library(tidyr)
library(dplyr)

# Loading in the simulation results
SimRedux <- read.csv("~/Documents/LANL/Simulations/SimRedux.csv")



######################################################
# PLOTS
######################################################


### Defining a theme for all the plots
theme_facets <- theme(strip.text = element_text(size = 12), legend.text = element_text(size = 12),
                      legend.title = element_text(size = 13), legend.spacing.y = unit(.01, "npc"),
                      legend.background = element_rect(colour = 'black', size = .5),
                      #legend.box.margin = ggplot2::margin(.05, .02, .02, .02, "npc"),
                      legend.position = c(.93, .17), legend.justification = c(1, 0))

### Plotting the score variable
plot_score <- ggplot(SimRedux) + geom_line(aes(x = lambdas, y = score, col = Model)) +
  geom_point(aes(x = lambdas, y = score, col = Model)) + xlab(expression(lambda)) + ylab("Score") +
  scale_color_discrete(name = "Forest Type") + 
  facet_wrap(Yx~.) + theme_bw() + theme_facets 
pdf("~/Documents/LANL/Simulations/Score_Plot_highdim.pdf", width = 8, height = 4)
plot_score
dev.off()

### Plotting RMSE
plot_RMSE <- ggplot(SimRedux) + geom_line(aes(x = lambdas, y = RMSE, col = Model)) +
  geom_point(aes(x = lambdas, y = RMSE, col = Model)) + xlab(expression(lambda)) + ylab("RMSE") +
  scale_color_discrete(name = "Forest Type") + 
  facet_wrap(Yx~.) + theme_bw() + theme_facets
pdf("~/Documents/LANL/Simulations/RMSE_Plot_highdim.pdf", width = 8, height = 4)
plot_RMSE
dev.off()

### Plotting MAE
plot_MAE <- ggplot(SimRedux) + geom_line(aes(x = lambdas, y = MAE, col = Model)) +
  geom_point(aes(x = lambdas, y = MAE, col = Model)) + xlab(expression(lambda)) + ylab("MAE") +
  scale_color_discrete(name = "Forest Type") + 
  facet_wrap(Yx~.) + theme_bw() + theme_facets
pdf("~/Documents/LANL/Simulations/MAE_Plot_highdim.pdf", width = 8, height = 4)
plot_MAE
dev.off()

### Plotting CovgPCT
## Creating a fake data frame needed for the legend
SimPCT_nominal <- data.frame(X = NA, Model = c("Nominal level"), RMSE = NA, MAE = NA, CovgPCT = .8, 
                     IntervalWidth = NA, score = NA, lambdas = sort(rep(seq(min(SimRedux$lambdas), max(SimRedux$lambdas), length.out = 5), 5)), 
                     Yx = unique(SimRedux$Yx))
SimPCT <- rbind(SimRedux, SimPCT_nominal)
# Creating the color scheme
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols_base <- gg_color_hue(2)
plot_PCT <- ggplot(SimPCT) + geom_line(aes(x = lambdas, y = CovgPCT, col = Model, lty = Model)) +
  geom_point(aes(x = lambdas, y = CovgPCT, col = Model, size = Model)) + xlab(expression(lambda)) + ylab("Covg. Probability") +
  scale_color_manual(name = "Forest Type", values = c(cols_base, "slateblue3")) + scale_size_manual(name = "Forest Type", values = c(1, 1, -1)) +
  scale_linetype_manual(name = "Forest Type", values = c("solid", "solid", "dashed")) + 
  #geom_hline(aes(yintercept = 0.8), col = 'deepskyblue', lty = 'dashed') + 
  facet_wrap(Yx~.) + theme_bw() + theme_facets + theme(legend.position = c(.95, .15))
pdf("~/Documents/LANL/Simulations/PCT_Plot_highdim.pdf", width = 8, height = 4)
plot_PCT
dev.off()

### Plotting IntervalWidth
plot_IW <- ggplot(SimRedux) + geom_line(aes(x = lambdas, y = IntervalWidth, col = Model)) +
  geom_point(aes(x = lambdas, y = IntervalWidth, col = Model)) + xlab(expression(lambda)) + ylab("Interval Width") +
  scale_color_discrete(name = "Forest Type") + 
  facet_wrap(Yx~.) + theme_bw() + theme_facets
pdf("~/Documents/LANL/Simulations/IW_Plot_highdim.pdf", width = 8, height = 4)
plot_IW
dev.off()


## Saving these guys


######################################################
# TABLES
######################################################

# Making the table a bit easier to read
SimTable <- SimRedux %>% dplyr::select(-X) %>% 
  dplyr::select(Yx, lambdas, Model, RMSE, MAE, CovgPCT, IntervalWidth, score) %>% 
  rename("Y|bX" = Yx, "lambda" = lambdas, "Covg" = CovgPCT, "Interval Width" = IntervalWidth, "Score" = score)

# A big unwieldly table
print(xtable(SimTable), digits = 3, include.rownames = F)

# A table using only 3 lambda values
lambda_small <- unlist(unique(SimTable$lambda))[c(1, 4, 8)]
SimSmall <- SimTable %>% filter(lambda %in% lambda_small)
x_small <- xtable(SimSmall)
caption(x_small) <- "Simulation results for three lambda values. Bolded values represent the better value between the weighted and unweighted forest."
label(x_small) <- "tab:SimTab_HighDim"
print(x_small, digits = 3, include.rownames = F, 
      hline.after = c(-1,0, seq(2, nrow(SimSmall), by = 2), seq(0, nrow(SimSmall), by = 6)),
      size = "footnotesize", table.placement = "H")

# Printing xtable by model
print_model_table <- function(model){
  SimModel <- SimTable %>% filter(`Y|bX` == model) %>% select(-`Y|bX`)
  x_model <- xtable(SimModel, align = "ccc|ccccc", digits = 3)
  model_no_ <- paste(unlist(strsplit(as.character(model), split =  "_")), collapse = " ")
  caption(x_model) <- paste("Simulation results for ", model_no_, ". Bolded values represent the better for a given $lambda$ setting.", sep = "")
  label(x_model) <- paste("tab:", model, sep = "")
  print(x_model, include.rownames = F, 
        hline.after = c(-1,0, seq(2, nrow(SimModel), by = 2)),
        size = "footnotesize", table.placement = "H")
}

models <- unique(SimTable$`Y|bX`)
print_models <- sapply(models, print_model_table)
