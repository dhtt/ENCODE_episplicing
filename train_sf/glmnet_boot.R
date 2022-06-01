library(car)
library(glmnet)
library(caret)
library(boot)

all_df = readRDS("all_df.RDS")
all_df = all_df[names(all_df) %in% c("METRN", "NES", "OTX2", "RXRA", "SLAIN1") == FALSE]
n_boot = 1000L
#GLMNET MODEL
train_model_boot = function(df_lists, n_boots, remove_method, model){
  print(model)
  model_list = vector("list")
  for (i in 1:length(df_lists)){
    df = df_lists[[i]]
    print(names(df_lists)[i]) 
    
    if (remove_method == "linear") df = df[, -findLinearCombos(df)$remove]
    else if (remove_method == "both") {
      lin_comb = findLinearCombos(df)$remove
      if (length(lin_comb) != 0) df = df[, -lin_comb]
      cor_col = findCorrelation(cor(as.matrix(df)), cutoff = 0.8)
      if (length(cor_col) != 0) df = df[, -cor_col[cor_col!=1]]
    }
    
    
    if (nrow(df) >= 10){
      if (model == "lm"){
        fit_control = trainControl(method = "LOOCV", returnResamp = "all")
        bootSamples = train(form = AS ~ ., data = df,
                           method = "lm",
                           trControl = fit_control,
                           metric = c("RMSE"))
      }
      else {
        if (model == "glmnet") {
          fit_control = trainControl(method = "LOOCV")
          bootSamples <- boot(df, function(df_boot, idx) {
            bootstrapData <- df_boot[idx, ]
            
            bootstrapMod <- train(form = AS ~ ., data = bootstrapData, 
                                  method = "glmnet",
                                  trControl = fit_control,
                                  metric = c("RMSE"))
            return(as.vector(coef(bootstrapMod$finalModel, bootstrapMod$bestTune$lambda)))
          }, R=n_boots)
        }
        else {
          if (model == "ridge") glmnet_grid = expand.grid(alpha = 0, lambda = seq(0, 10, 0.1))
          else if (model == "lasso") glmnet_grid = expand.grid(alpha = 1, lambda = seq(0, 10, 0.1))
          else if (model == "enet") glmnet_grid = expand.grid(alpha = seq(0.1, 0.9, 0.1), lambda = seq(0, 10, 0.1))
          
          fit_control = trainControl(method = "LOOCV", search = "grid")
          
          bootSamples <- boot(df, function(df_boot, idx) {
            bootstrapData <- df_boot[idx, ]
            
            bootstrapMod <- train(form = AS ~ ., data = bootstrapData, 
                                  method = "glmnet",
                                  trControl = fit_control,
                                  tuneGrid = glmnet_grid,
                                  metric = c("RMSE"))
            return(as.vector(coef(bootstrapMod$finalModel, bootstrapMod$bestTune$lambda)))
          }, R=n_boots)
        }
      }
    }
    else {
      bootSamples = "FALSE"
    }
    model_list[[i]] = bootSamples
  }
  names(model_list) = names(df_lists)
  return(model_list)
}

all_model_boot_lm = train_model_boot(all_df, n_boot, "both", "lm")
saveRDS(all_model_boot_lm, "all_model_boot_lm.RDS")
all_model_boot_ridge = train_model_boot(all_df, n_boot, "both", "ridge")
saveRDS(all_model_boot_ridge, "all_model_boot_ridge.RDS")
all_model_boot_lasso = train_model_boot(all_df, n_boot, "both", "lasso" )
saveRDS(all_model_boot_lasso, "all_model_boot_lasso.RDS")
all_model_boot_enet = train_model_boot(all_df, n_boot, "both", "enet")
saveRDS(all_model_boot_enet, "all_model_boot_enet.RDS")
all_model_boot_glmnet = train_model_boot(all_df, n_boot, "both", "glmnet")
saveRDS(all_model_boot_glmnet, "all_model_boot_glmnet.RDS")
# all_model_boot = readRDS("all_model_boot.RDS")
# 
# all_model_boot_true = all_model_boot[all_model_boot != "FALSE"]
# names(all_model_boot_true[1])
# all_model_boot_sd = lapply(all_model_boot_true, function(x){
#   std = as.data.frame(apply(x$t, 2, sd))
#   colnames(std) = "SD"
#   std$RBP = c("Intercept", colnames(x$data)[2:ncol(x$data)])
#   std$Mean = x$t0
#   return(std)})
# all_model_boot_sd_df = do.call("rbind", all_model_boot_sd)
# all_model_boot_sd_df$gene = rep(names(all_model_boot_sd), sapply(all_model_boot_sd, nrow))
# all_model_boot_sd_df = all_model_boot_sd_df[all_model_boot_sd_df$Mean != 0, ]
# 
# head(all_model_boot_sd_df)
# dev.off()
# ggplot(all_model_boot_sd_df, aes(y = RBP, x = Mean, col = gene)) +
#   geom_point() +
#   geom_errorbar(aes(xmin = Mean - SD, xmax = Mean + SD, y = RBP), size=0.5, width = 0.25) +
#   geom_vline(xintercept = 0, linetype='dashed', size=0.25) +
#   facet_wrap(vars(gene), ncol = 5, scales = "free") +
#   theme(legend.position = "bottom", axis.title.y = element_blank())
