# modelling with Pterocles alchata observations
# cross-validation technique: Spatial Leave-One-Out
# each presence point is used alone as test and all points in a radius of 20km 
# around the test point are erased in each model
# this way to correct the effect of spatial autocorrelation basing on ecological criteria
# and training points per model are maximized

rm(list=ls())
source('scripts/utils.R') # all used libraries and my custom functions

library(doParallel) 
# n.cores <- strtoi(Sys.getenv("SLURM_CPUS_PER_TASK")) # slurm cores
n.cores <- parallel::detectCores() - 1 # silence if running in slurm
my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()

dataset <- na.omit(read.csv2("data/iberica_PA_present.csv"))
dataset$PA <- as.factor(dataset$PA)
dataset_points <- st_as_sf(dataset, coords = c("X", "Y"), crs = 4326)

# load env_df of every scenario

env_df <- read.csv2("data/envdf_present.csv")
env_df_ssp126 <- read.csv2("data/envdf_ssp126.csv")
env_df_ssp585 <- read.csv2("data/envdf_ssp585.csv")

# run models iterating through presence points only (375 models)
tic("RF-SLOO and projections")
foreach(i = 1:nrow(dataset[dataset$PA==1,]),
        .verbose = T,
        .packages=c("randomForest", "sf", "dplyr", "pROC",
                    "caret", "modEvA", "stats")) %dopar% {
  
  # set unique seed for each model
  set.seed(i)
  
  # generate 20km buffer around test point
  buffer <- st_buffer(dataset_points[i,], dist = 20000)
  
  # store all points inside the buffer
  to_erase <- st_intersection(dataset_points, buffer)
  coordinates <- data.frame(X = st_coordinates(to_erase)[,'X'],
                            Y = st_coordinates(to_erase)[,'Y'])
  to_erase <- to_erase %>% as.data.frame() %>% dplyr::select(-geometry) %>% cbind(coordinates)
  
  # select all points except to_erase to train
  train <- dataset %>%
    anti_join(to_erase, by = colnames(dataset))
  test <- dataset[i, ]
  
  # sample size
  prNum <- as.numeric(table(train$PA)["1"]) # number of presences
  
  # sampsize vector
  samsize <- c("0" = prNum, "1" = prNum)
  
  # train the model
  model <- randomForest(PA ~ . -X -Y,
                        data = train,
                        sampsize = samsize,
                        do.classif = TRUE,
                        importance = TRUE,
                        ntree = 1000,
                        mtry = 6,
                        nodesize = 90,
                        keep.inbag = T) 

  # save the model
  saveRDS(model, paste0("results/alchata/model",i,".rds"))
  
  # extract inbag samples to calculate training metrics
  has_inbag <- apply(model$inbag, 1, function(x) any(x > 0)) #boolean ib vs oob
  train_inbag <- train[has_inbag, ]
  
  # generate training predictions (only select probs of presence)
  preds<-predict(model, data = trian_inbag, type="prob")
  train_inbag$pred<-preds[,2] 
  
  # save training preds
  write.csv2(train_inbag['pred'], paste0('results/alchata/train_preds',i,'.csv'), row.names = F)
  
  # create ROC curve to calculate thresholds 
  roc_curve <- roc(train_inbag$PA, train_inbag$pred)
  
  # calculate threshold (maxTSS) and binarize preds
  maxTSS <- coords(roc_curve, "best", maximize="tss")$threshold
  train_inbag$pred_maxTSS <- ifelse(train_inbag$pred > maxTSS, 1, 0)
  
  # AUC training
  auc <- pROC::auc(as.numeric(train_inbag$PA), as.numeric(train_inbag$pred_maxTSS))

  # bin preds as factors because needed by confusionmatrix
  train_inbag$pred_maxTSS <- factor(train_inbag$pred_maxTSS, levels = c("0", "1"))

  # confusion matrix of the training per thresholds
  cm <- caret::confusionMatrix(train_inbag$pred_maxTSS, train_inbag$PA, '1')

  # TSS training per thresholds
  tss <- cm$byClass[['Sensitivity']] + cm$byClass[['Specificity']] - 1
  
  # boyce index
  boyce <- Boyce(model=model, plot=F)$Boyce
  
  # COR (point biserial correlation, dicotomous vs continuous)
  COR <- cor(as.numeric(train_inbag$PA), train_inbag$pred, method = "pearson")
  
  # regression metrics: MAE, RMSE, MSE, R-squared
  # errors between predictions and observed (1)
  train_inbag$PA <- as.numeric(train_inbag$PA)
  mae <- MAE(train_inbag$pred, train_inbag$PA)
  mse <- mean((train_inbag$pred - train_inbag$PA)^2)
  rmse <- RMSE(pred = train_inbag$pred, obs = train_inbag$PA)
  Rsq <- 1 - (sum((train_inbag$PA - train_inbag$pred)^2) / sum((train_inbag$PA - mean(train_inbag$PA))^2))
  
  # generate testing predictions and binaries from maxSSS
  preds <- data.frame(pred = predict(model, newdata = test, type ="prob"))
  test$pred <- preds[,2]
  test$pred_maxTSS <- ifelse(test$pred > maxTSS, 1, 0)

  # save all the relevant metrics 
  metrics_row <- c(test[1,'pred'], test[1,'pred_maxTSS'], maxTSS,
                   cm$byClass[['Sensitivity']], cm$byClass[['F1']],
                   cm$byClass[['Balanced Accuracy']], tss, auc,
                   boyce, COR, mae, mse, rmse, Rsq)
  write.csv2(metrics_row, paste0("results/alchata/metrics",i,".csv"), row.names=F)
  
  # generate predictions for the whole territory in present and save
  preds <- data.frame(pred = predict(model, newdata = env_df, type ="prob"))
  write.csv2(preds[,2], paste0("results/alchata/preds",i,".csv"), row.names=F)
  # same for ssp126
  preds <- data.frame(pred = predict(model, newdata = env_df_ssp126, type ="prob"))
  write.csv2(preds[,2], paste0("results/alchata/ssp126_preds",i,".csv"), row.names=F)
  # same for ssp585
  preds <- data.frame(pred = predict(model, newdata = env_df_ssp585, type ="prob"))
  write.csv2(preds[,2], paste0("results/alchata/ssp585_preds",i,".csv"), row.names=F)
}
toc() # 11543.58 sec, 11 i5 cores

# join all metrics rows into one dataframe
column_names <- c('cont_pred','bin_pred','threshold','Sens_train',
                  'F1_train','B.Accuracy_train','TSS_train','AUC_train',
                  'Boyce_train', 'COR_train', 'MAE_train', 'MSE_train',
                  'RMSE_train', 'Rsq_train')
metrics <- data.frame(matrix(ncol = length(column_names), nrow = 0))
file_list <- list.files("results/alchata", pattern = 'metrics', full.names = T)
for (file in file_list) {
  metrics_row <- t(read.csv2(file))
  metrics <- rbind(metrics, metrics_row)
}
colnames(metrics) <- column_names
write.csv2(metrics, "results/metrics_iberica.csv", row.names=F)

# join all preds columns into one dataframe
# present
preds <- data.frame(matrix(ncol = 0, nrow = nrow(env_df)))
file_list <- list.files("results/alchata", pattern = 'preds', full.names = T)
file_list <- file_list[!grepl('126', file_list)]
file_list <- file_list[!grepl('585', file_list)]
for (file in file_list) {
  preds_col <- read.csv2(file)
  preds <- cbind(preds, preds_col)
}
write.csv2(preds,"results/pres_preds_iberica.csv", row.names=F)
# ssp126
preds <- data.frame(matrix(ncol = 0, nrow = nrow(env_df)))
file_list <- list.files("results/alchata", pattern = 'ssp126_preds', full.names = T)
for (file in file_list) {
  preds_col <- read.csv2(file)
  preds <- cbind(preds, preds_col)
}
write.csv2(preds,"results/ssp126_preds_iberica.csv", row.names=F)
# ssp585
preds <- data.frame(matrix(ncol = 0, nrow = nrow(env_df)))
file_list <- list.files("results/alchata", pattern = 'ssp585_preds', full.names = T)
for (file in file_list) {
  preds_col <- read.csv2(file)
  preds <- cbind(preds, preds_col)
}
write.csv2(preds,"results/ssp585_preds_iberica.csv", row.names=F)

# metrics for test predictions (only positive class)
# overall sens
mean(metrics$bin_pred == 1) # 0.93
# overall cont pred
mean(metrics$cont_pred) # 0.81
# MAE
MAE(metrics$cont_pred, rep(1,nrow(metrics))) # 0.19
# MSE
mean((metrics$cont_pred - rep(1,nrow(metrics)))^2) # 0.09
# RMSE
RMSE(pred = metrics$cont_pred, obs = rep(1,nrow(metrics))) # 0.30

# metrics boxplots
par(mfrow = c(3, 4))
boxplot(metrics$cont_pred, main='Cont. Pred.')
boxplot(metrics$Sens_train, main='Sens.')
boxplot(metrics$F1_train, main='F1')
boxplot(metrics$B.Accuracy_train, main='B.Accuracy')
boxplot(metrics$TSS_train, main='TSS')
boxplot(metrics$AUC_train, main='AUC')
boxplot(metrics$Boyce_train, main='Boyce')
boxplot(metrics$COR_train, main='COR')
boxplot(metrics$MAE_train, main='MAE')
boxplot(metrics$MSE_train, main='MSE')
boxplot(metrics$RMSE_train, main='RMSE')
boxplot(metrics$Rsq_train, main='R-squared')

# variable importance calculation

# load models filenames
file_list <- list.files("results/alchata", pattern = 'model', full.names = T)

# empty dfs to store accuracy and gini values for each model
accu <- data.frame(matrix(nrow = 27, ncol = 0))
gini <- data.frame(matrix(nrow = 27, ncol = 0))
# iterate through the model list to extract importances
for (i in 1:length(file_list)) {
  
  model <- readRDS(file_list[[i]])
  col_name <- paste0('fold', i)
  
  accu <- cbind(accu, model$importance[, 3])
  names(accu)[ncol(accu)] <- col_name
  
  gini <- cbind(gini, model$importance[, 4])
  names(gini)[ncol(gini)] <- col_name
}

### coloured boxplots of variable importance

# first with accuracy

# transpose df
accut <- as.data.frame(t(accu))
# define variable categories (for color legend)
categorias <- list(climaticas = c("mean_temp", "mean_temp_seas", "prec", "prec_seas"),
                   topograficas = c("dem", "slope"),
                   perturbacion = c("linear_infrastr"))
# vector with all variable names
todas_variables <- colnames(accut)
# category for land uses
usos <- setdiff(todas_variables, unlist(categorias))
# calculate absolute means of each variable
accut_mean_abs <- accut %>%
  summarise_all(~ mean(abs(.))) %>%
  gather(key = "Variable", value = "Mean_Abs") %>%
  arrange(Mean_Abs) # order ascending
# define category of each variable
accut_mean_abs$categoria <- case_when(
  accut_mean_abs$Variable %in% categorias$climaticas ~ "Climate",
  accut_mean_abs$Variable %in% categorias$topograficas ~ "Topography",
  accut_mean_abs$Variable %in% categorias$perturbacion ~ "Anthropic",
  accut_mean_abs$Variable %in% usos ~ "Land Uses"
)
# vector for variables
ordered_vars <- accut_mean_abs$Variable
# transform to long format
accu_long_ordered <- accut %>%
  gather(key = "Variable", value = "Value") %>%
  mutate(Variable = factor(Variable, levels = ordered_vars))
# create category variable in long format df
accu_long_ordered <- accu_long_ordered %>%
  mutate(categoria = case_when(
    Variable %in% categorias$climaticas ~ "Climate",
    Variable %in% categorias$topograficas ~ "Topography",
    Variable %in% categorias$perturbacion ~ "Anthropic",
    TRUE ~ "Land Uses"  
  ))
# define colors
colores <- c("Climate" = "#95b8f6", 
             "Topography" = "#f9d99a", 
             "Anthropic" = "#dcd9f8",
             "Land Uses" = "#fa5f49")

# create boxplots
plot1 <- ggplot(accu_long_ordered, aes(x = Variable, y = abs(Value), fill = categoria)) +
  geom_boxplot() +
  coord_flip() +  # turn labels in x axis
  theme_minimal() +  
  labs(title = "Abs values of Mean Decrease in Accuracy",
       x = "Variable", y = "MeanDecreaseAccuracy", fill = "Category") +
  scale_fill_manual(values = colores)  # assign defined colors

# same with gini

ginit <- as.data.frame(t(gini))
categorias <- list(climaticas = c("mean_temp", "mean_temp_seas", "prec", "prec_seas"),
                   topograficas = c("dem", "slope"),
                   perturbacion = c("linear_infrastr"))
todas_variables <- colnames(ginit)
usos <- setdiff(todas_variables, unlist(categorias))
ginit_mean_abs <- ginit %>%
  summarise_all(~ mean(abs(.))) %>%
  gather(key = "Variable", value = "Mean_Abs") %>%
  arrange(Mean_Abs)
ginit_mean_abs$categoria <- case_when(
  ginit_mean_abs$Variable %in% categorias$climaticas ~ "Climate",
  ginit_mean_abs$Variable %in% categorias$topograficas ~ "Topography",
  ginit_mean_abs$Variable %in% categorias$perturbacion ~ "Anthropic",
  ginit_mean_abs$Variable %in% usos ~ "Land Uses"
)
ordered_vars <- ginit_mean_abs$Variable
gini_long_ordered <- ginit %>%
  gather(key = "Variable", value = "Value") %>%
  mutate(Variable = factor(Variable, levels = ordered_vars))
gini_long_ordered <- gini_long_ordered %>%
  mutate(categoria = case_when(
    Variable %in% categorias$climaticas ~ "Climate",
    Variable %in% categorias$topograficas ~ "Topography",
    Variable %in% categorias$perturbacion ~ "Anthropic",
    TRUE ~ "Land Uses"  
  ))
colores <- c("Climate" = "#95b8f6", 
             "Topography" = "#f9d99a", 
             "Anthropic" = "#dcd9f8",
             "Land Uses" = "#fa5f49")
plot2 <- ggplot(gini_long_ordered, aes(x = Variable, y = abs(Value), fill = categoria)) +
  geom_boxplot() +
  coord_flip() +  
  theme_minimal() +  
  labs(title = "Abs values of Mean Decrease in GINI",
       x = "Variable", y = "MeanDecreaseGini", fill = "Category") +
  scale_fill_manual(values = colores)  

plot3 <- grid.arrange(plot1, plot2, nrow = 1, ncol = 2)
ggsave("results/iberica_varimp.jpg", plot3, width = 10, height = 5, units = "in")

# partial plots

file_list <- list.files("results/alchata", pattern = 'model', full.names = T)

tic("PDPs")
foreach(i = 1:length(file_list),
        .verbose = T,
        .packages=c("randomForest", "pdp")) %dopar% {
  # load corresponding model from disk
  model <- readRDS(file_list[[i]])
  # prec
  pred <- partial(model, pred.var = "prec", which.class = "1",
                  prob = T, train = dataset, plot = FALSE)
  write.csv2(pred, paste0("results/alchata/pdp_prec",i,".csv"), row.names=F)
  # mean_temp
  pred <- partial(model, pred.var = "mean_temp", which.class = "1",
                  prob = T, train = dataset, plot = FALSE)
  write.csv2(pred, paste0("results/alchata/pdp_mean_temp",i,".csv"), row.names=F)
  # mean_temp_seas
  pred <- partial(model, pred.var = "mean_temp_seas", which.class = "1",
                  prob = T, train = dataset, plot = FALSE)
  write.csv2(pred, paste0("results/alchata/pdp_mean_temp_seas",i,".csv"), row.names=F)
  # prec_seas
  pred <- partial(model, pred.var = "prec_seas", which.class = "1",
                  prob = T, train = dataset, plot = FALSE)
  write.csv2(pred, paste0("results/alchata/pdp_prec_seas",i,".csv"), row.names=F)
  # linear_infrastr
  pred <- partial(model, pred.var = "linear_infrastr", which.class = "1",
                  prob = T, train = dataset, plot = FALSE)
  write.csv2(pred, paste0("results/alchata/pdp_linear_infrastr",i,".csv"), row.names=F)
  # elev
  pred <- partial(model, pred.var = "elev", which.class = "1",
                  prob = T, train = dataset, plot = FALSE)
  write.csv2(pred, paste0("results/alchata/pdp_elev",i,".csv"), row.names=F)
  # slope
  pred <- partial(model, pred.var = "slope", which.class = "1",
                  prob = T, train = dataset, plot = FALSE)
  write.csv2(pred, paste0("results/alchata/pdp_slope",i,".csv"), row.names=F)
  # forest (binary)
  pred <- partial(model, pred.var = "forest", which.class = c("0", "1"),
                  train = dataset, plot = FALSE)
  write.csv2(pred, paste0("results/alchata/pdp_forest",i,".csv"), row.names=F)
  # crop_low (binary)
  pred <- partial(model, pred.var = "crop_low", which.class = c("0", "1"),
                  train = dataset, plot = FALSE)
  write.csv2(pred, paste0("results/alchata/pdp_crop_low",i,".csv"), row.names=F)
  # ext_perm_crop (binary)
  pred <- partial(model, pred.var = "ext_perm_crop", which.class = c("0", "1"),
                  train = dataset, plot = FALSE)
  write.csv2(pred, paste0("results/alchata/pdp_ext_perm_crop",i,".csv"), row.names=F)
}
toc() # 15094 sec

# join all preds per variable and plot 

# prec
predictions_df <- data.frame(matrix(ncol = 2, nrow = 0))
file_list <- list.files("results/alchata", pattern = 'prec', full.names = T)
file_list <- file_list[!grepl('seas', file_list)]
for (file in file_list) {
  prec_file <- read.csv2(file)
  predictions_df <- rbind(predictions_df, prec_file)
}
summary_stats <- predictions_df %>%
  group_by(prec) %>%
  summarise(mean_yhat = mean(yhat), sd_yhat = sd(yhat))
prec <- ggplot(summary_stats, aes(x = prec, y = mean_yhat)) +
  geom_line() + # mean
  geom_ribbon(aes(ymin = mean_yhat - sd_yhat, ymax = mean_yhat + sd_yhat), alpha = 0.3) + # variability
  theme_bw() +
  labs(x = "prec", y = "Positive Prob")
print(prec)

# mean_temp
predictions_df <- data.frame(matrix(ncol = 2, nrow = 0))
file_list <- list.files("results/alchata", pattern = 'mean_temp', full.names = T)
file_list <- file_list[!grepl('seas', file_list)]
for (file in file_list) {
  mean_temp_file <- read.csv2(file)
  predictions_df <- rbind(predictions_df, mean_temp_file)
}
summary_stats <- predictions_df %>%
  group_by(mean_temp) %>%
  summarise(mean_yhat = mean(yhat), sd_yhat = sd(yhat))
mean_temp <- ggplot(summary_stats, aes(x = mean_temp, y = mean_yhat)) +
  geom_line() + # mean
  geom_ribbon(aes(ymin = mean_yhat - sd_yhat, ymax = mean_yhat + sd_yhat), alpha = 0.3) + # variability
  theme_bw() +
  labs(x = "mean_temp", y = "Positive Prob")
print(mean_temp)

# prec_seas
predictions_df <- data.frame(matrix(ncol = 2, nrow = 0))
file_list <- list.files("results/alchata", pattern = 'prec_seas', full.names = T)
for (file in file_list) {
  prec_seas_file <- read.csv2(file)
  predictions_df <- rbind(predictions_df, prec_seas_file)
}
summary_stats <- predictions_df %>%
  group_by(prec_seas) %>%
  summarise(mean_yhat = mean(yhat), sd_yhat = sd(yhat))
prec_seas <- ggplot(summary_stats, aes(x = prec_seas, y = mean_yhat)) +
  geom_line() + # mean
  geom_ribbon(aes(ymin = mean_yhat - sd_yhat, ymax = mean_yhat + sd_yhat), alpha = 0.3) + # variability
  theme_bw() +
  labs(x = "prec_seas", y = "Positive Prob")
print(prec_seas)

# mean_temp_seas
predictions_df <- data.frame(matrix(ncol = 2, nrow = 0))
file_list <- list.files("results/alchata", pattern = 'mean_temp_seas', full.names = T)
for (file in file_list) {
  mean_temp_seas_file <- read.csv2(file)
  predictions_df <- rbind(predictions_df, mean_temp_seas_file)
}
summary_stats <- predictions_df %>%
  group_by(mean_temp_seas) %>%
  summarise(mean_yhat = mean(yhat), sd_yhat = sd(yhat))
mean_temp_seas <- ggplot(summary_stats, aes(x = mean_temp_seas, y = mean_yhat)) +
  geom_line() + # mean
  geom_ribbon(aes(ymin = mean_yhat - sd_yhat, ymax = mean_yhat + sd_yhat), alpha = 0.3) + # variability
  theme_bw() +
  labs(x = "mean_temp_seas", y = "Positive Prob")
print(mean_temp_seas)

# elev
predictions_df <- data.frame(matrix(ncol = 2, nrow = 0))
file_list <- list.files("results/alchata", pattern = 'elev', full.names = T)
for (file in file_list) {
  elev_file <- read.csv2(file)
  predictions_df <- rbind(predictions_df, elev_file)
}
summary_stats <- predictions_df %>%
  group_by(elev) %>%
  summarise(mean_yhat = mean(yhat), sd_yhat = sd(yhat))
elev <- ggplot(summary_stats, aes(x = elev, y = mean_yhat)) +
  geom_line() + # mean
  geom_ribbon(aes(ymin = mean_yhat - sd_yhat, ymax = mean_yhat + sd_yhat), alpha = 0.3) + # variability
  theme_bw() +
  labs(x = "elev", y = "Positive Prob")
print(elev)

# slope
predictions_df <- data.frame(matrix(ncol = 2, nrow = 0))
file_list <- list.files("results/alchata", pattern = 'slope', full.names = T)
for (file in file_list) {
  slope_file <- read.csv2(file)
  predictions_df <- rbind(predictions_df, slope_file)
}
summary_stats <- predictions_df %>%
  group_by(slope) %>%
  summarise(mean_yhat = mean(yhat), sd_yhat = sd(yhat))
slope <- ggplot(summary_stats, aes(x = slope, y = mean_yhat)) +
  geom_line() + # mean
  geom_ribbon(aes(ymin = mean_yhat - sd_yhat, ymax = mean_yhat + sd_yhat), alpha = 0.3) + # variability
  theme_bw() +
  labs(x = "slope", y = "Positive Prob")
print(slope)

# linear_infrastr
predictions_df <- data.frame(matrix(ncol = 2, nrow = 0))
file_list <- list.files("results/alchata", pattern = 'linear_infrastr', full.names = T)
for (file in file_list) {
  linear_infrastr_file <- read.csv2(file)
  predictions_df <- rbind(predictions_df, linear_infrastr_file)
}
summary_stats <- predictions_df %>%
  group_by(linear_infrastr) %>%
  summarise(mean_yhat = mean(yhat), sd_yhat = sd(yhat))
linear_infrastr <- ggplot(summary_stats, aes(x = linear_infrastr, y = mean_yhat)) +
  geom_line() + # mean
  geom_ribbon(aes(ymin = mean_yhat - sd_yhat, ymax = mean_yhat + sd_yhat), alpha = 0.3) + # variability
  theme_bw() +
  labs(x = "linear_infrastr", y = "Positive Prob")
print(linear_infrastr)

# forest
predictions_df <- data.frame(matrix(ncol = 2, nrow = 0))
file_list <- list.files("results/alchata", pattern = 'forest', full.names = T)
for (file in file_list) {
  forest_file <- read.csv2(file)
  predictions_df <- rbind(predictions_df, forest_file)
}
# statistics for each level
summary_stats <- predictions_df %>%
  group_by(forest) %>%
  summarise(median_yhat = median(yhat), q1_yhat = quantile(yhat, 0.25), q3_yhat = quantile(yhat, 0.75))
# boxplots for both levels
forest <- ggplot(summary_stats, aes(x = factor(forest), y = median_yhat)) +
  geom_point() +
  geom_errorbar(aes(ymin = q1_yhat, ymax = q3_yhat), width = 0.1) + # error bars
  theme_bw() +
  labs(x = "forest", y = "Positive Prob")
print(forest)

# crop_low
predictions_df <- data.frame(matrix(ncol = 2, nrow = 0))
file_list <- list.files("results/alchata", pattern = 'crop_low', full.names = T)
for (file in file_list) {
  crop_low_file <- read.csv2(file)
  predictions_df <- rbind(predictions_df, crop_low_file)
}
# statistics for each level
summary_stats <- predictions_df %>%
  group_by(crop_low) %>%
  summarise(median_yhat = median(yhat), q1_yhat = quantile(yhat, 0.25), q3_yhat = quantile(yhat, 0.75))
# boxplots for both levels
crop_low <- ggplot(summary_stats, aes(x = factor(crop_low), y = median_yhat)) +
  geom_point() +
  geom_errorbar(aes(ymin = q1_yhat, ymax = q3_yhat), width = 0.1) + # error bars
  theme_bw() +
  labs(x = "crop_low", y = "Positive Prob")
print(crop_low)

# ext_perm_crop
predictions_df <- data.frame(matrix(ncol = 2, nrow = 0))
file_list <- list.files("results/alchata", pattern = 'ext_perm_crop', full.names = T)
for (file in file_list) {
  ext_perm_crop_file <- read.csv2(file)
  predictions_df <- rbind(predictions_df, ext_perm_crop_file)
}
# statistics for each level
summary_stats <- predictions_df %>%
  group_by(ext_perm_crop) %>%
  summarise(median_yhat = median(yhat), q1_yhat = quantile(yhat, 0.25), q3_yhat = quantile(yhat, 0.75))
# boxplots for both levels
ext_perm_crop <- ggplot(summary_stats, aes(x = factor(ext_perm_crop), y = median_yhat)) +
  geom_point() +
  geom_errorbar(aes(ymin = q1_yhat, ymax = q3_yhat), width = 0.1) + # error bars
  theme_bw() +
  labs(x = "ext_perm_crop", y = "Positive Prob")
print(ext_perm_crop)

pdps <- grid.arrange(mean_temp, mean_temp_seas, prec, prec_seas, elev, slope,
                     linear_infrastr, forest, crop_low, ext_perm_crop, nrow = 2, ncol = 5)
ggsave("results/pdps_iberica.jpg", pdps, width = 15, height = 5, units = "in")

registerDoSEQ()