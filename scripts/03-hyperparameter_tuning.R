# hyperparameter tuning on parallel

# we need to tune mtry, nodesize and ntree
# we need to cross validate to ensure all data is being used in test
# we're using random cross-validation in 5 folds

rm(list=ls())
source('scripts/utils.R') # all used libraries and my custom functions

# n_cores <- strtoi(Sys.getenv('SLURM_CPUS_PER_TASK')) # slurm cores
n.cores <- parallel::detectCores() - 1
my.cluster <- parallel::makeCluster(n.cores, type = 'PSOCK')
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()

x <- foreach(
  i = 1:10, 
  .combine = 'c'
) %dopar% {
  sqrt(i)
}
x # check parallel execution is on

dataset <- na.omit(read.csv('data/iberica_PA_present.csv'))
dataset$PA <- as.factor(dataset$PA)

# create metrics df to append results in the loop
column_names <- c('combination', 'fold', 'mtry', 'ntree', 'nodesize','thr_value',
                  'Sens_train', 'F1_train','B.Accuracy_train','TSS_train','AUC_train',
                  'Sens_test', 'F1_test', 'B.Accuracy_test','TSS_test','AUC_test')
metrics <- data.frame(matrix(ncol = length(column_names), nrow = 0))
colnames(metrics) <- column_names

# set random grid to search hyperparameters
# define ranges
set.seed(123)
grid <- expand.grid(mtry = 2:15, # variables in each tree
                    ntree = seq(1000, 10000, by = 1000), # trees in each forest
                    nodesize = seq(10, nrow(dataset[dataset$PA==1,]), by = 10)) # minimum size of terminal nodes
write.csv(grid, 'data/grid_iberica.csv', row.names=F)
# grid <- read.csv('data/grid_iberica.csv')

# run models iterating through each hyperpar combination
tic('random hyperparameters')
metrics <- foreach(i = 1:nrow(grid),
                   .packages=c('randomForest', 'sf', 'dplyr', 'pROC', 'caret'),
                   .combine='rbind') %:% # iterations through hyperpar combinations
  
  foreach(j = 1:5,
          .combine='rbind') %dopar% { # iterations through random cv folds (5)
            
            # set unique seed for each fold
            set.seed(j**2)
            
            # use caret to randomly split train/test data (80-20)
            index <- createDataPartition(dataset$PA, p = 0.8, list = FALSE)
            train <- dataset[index, ]
            test <- dataset[-index, ]  
            
            # extracting sample size
            prNum <- as.numeric(table(train$PA)['1']) # number of presences
            
            # set sample size as vector, same size for 0 and for 1
            samsize <- c('0' = prNum, '1' = prNum)
            
            # train the model
            model <- randomForest(PA ~ . -X -Y -Species,
                                  data = train,
                                  sampsize = samsize,
                                  do.classif = TRUE,
                                  mtry = grid[i,1],
                                  ntree = grid[i,2],
                                  nodesize = grid[i,3])
            
            # generate training predictions (only select probs of presence)
            preds<-predict(model, type='prob')
            train$pred<-preds[,2]
            
            # create ROC curve to calculate thresholds
            roc_curve <- roc(train$PA, train$pred)
            
            # calculate thresholds and binarize preds
            # maxTSS
            maxTSS <- coords(roc_curve, 'best', maximize='tss')$threshold
            train$pred <- ifelse(train$pred > maxTSS, 1, 0)
            
            # AUC training per thresholds
            auc_train <- pROC::auc(as.numeric(train$PA), as.numeric(train$pred))
            
            # bin preds as factors because needed by confusionmatrix
            train$pred <- factor(train$pred, levels = c('0', '1'))
            
            # confusion matrix of the training per thresholds
            cm_train <- confusionMatrix(train$pred, train$PA, '1')
            
            # TSS training per thresholds
            tss_train <- cm_train$byClass[['Sensitivity']] + cm_train$byClass[['Specificity']] - 1
            
            # generate testing predictions and binaries from maxSSS
            preds <- data.frame(pred = predict(model, newdata = test, type ='prob'))
            test$pred <- preds[,2]
            test$pred <- ifelse(test$pred > maxTSS, 1, 0)
            
            # AUC test
            auc_test <- pROC::auc(as.numeric(test$PA), as.numeric(test$pred))
            
            # bin preds as factors because needed by confusionmatrix
            test$pred <- factor(test$pred, levels = c('0', '1'))
            
            # confusion matrix of the testing per thresholds
            cm_test <- confusionMatrix(test$pred, test$PA, '1')
            
            # TSS testing per thresholds
            tss_test <- cm_test$byClass[['Sensitivity']] + cm_test$byClass[['Specificity']] - 1
            
            # store all the relevant values in the metrics dataframe
            # column_names <- c('mtry', 'ntree', 'nodesize','thr_value',
            # 'Sens_train', 'F1_train','B.Accuracy_train','TSS_train','AUC_train',
            # 'Sens_test', 'F1_test', 'B.Accuracy_test','TSS_test','AUC_test')
            
            metrics <- c(i, j, grid[i,1], grid[i,2], grid[i,3], maxTSS,
                         cm_train$byClass[['Sensitivity']],
                         cm_train$byClass[['F1']],
                         cm_train$byClass[['Balanced Accuracy']],
                         tss_train, auc_train, # training
                         cm_test$byClass[['Sensitivity']],
                         cm_test$byClass[['F1']],
                         cm_test$byClass[['Balanced Accuracy']],
                         tss_test, auc_test # testing
            )
            return(metrics)
          }
toc() 

# assign col names
metrics <- as.data.frame(metrics)
colnames(metrics) <- column_names
write.csv(metrics, 'results/hyperparameters_iberica.csv', row.names = F)
# metrics <- read.csv('results/hyperparameters_iberica.csv')

metrics_mean <- metrics %>%
  group_by(combination) %>%
  summarise(across(-fold, mean))
write.csv(metrics_mean, 'results/hyperparameters_iberica_mean.csv', row.names = F)
# metrics_mean <- read.csv('results/hyperparameters_iberica_mean.csv')

# best combination should be the one that gets higher sensitivity
# and lower difference between training and testing sensitivity
# because if test sens < train sens = overfitting
# if test sens > train sens = underfitting
# equal metrics mean generalization, which is our goal

metrics_mean$sum <- rowSums(metrics_mean[,5:15])
metrics_mean[which.max(metrics_mean$sum),]

metrics_mean$sens_diff <- abs(metrics_mean$Sens_test-metrics_mean$Sens_train)
metrics_mean[which.min(metrics_mean$sens_diff),]

# we choose nodesize=90 and mtry=6, since it's the combination with lowest sens diff (0)
# and still high metrics sum (7.62, max=7.96)
# this combination forms the most complex trees in terms of depth and still maximizes generalization
# better combinations obtain the same metrics from 1000 to 10000 trees, so we choose 1000

### now with ortega
rm(list=ls())
source('scripts/utils.R') # all used libraries and my custom functions

# n_cores <- strtoi(Sys.getenv('SLURM_CPUS_PER_TASK')) # slurm cores
n.cores <- parallel::detectCores() - 1
my.cluster <- parallel::makeCluster(n.cores, type = 'PSOCK')
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()
dataset <- na.omit(read.csv('data/ortega_PA_present.csv'))
dataset$PA <- as.factor(dataset$PA)

# create metrics df to append results in the loop
column_names <- c('combination', 'fold', 'mtry', 'ntree', 'nodesize','thr_value',
                  'Sens_train', 'F1_train','B.Accuracy_train','TSS_train','AUC_train',
                  'Sens_test', 'F1_test', 'B.Accuracy_test','TSS_test','AUC_test')
metrics <- data.frame(matrix(ncol = length(column_names), nrow = 0))
colnames(metrics) <- column_names

# set random grid to search hyperparameters
# define ranges
set.seed(123)
grid <- expand.grid(mtry = 2:15, # variables in each tree
                    ntree = seq(1000, 10000, by = 1000), # trees in each forest
                    nodesize = seq(10, nrow(dataset[dataset$PA==1,]), by = 10)) # minimum size of terminal nodes
write.csv(grid, 'data/grid_ortega.csv', row.names=F)
# grid <- read.csv('data/grid_ortega.csv')

# run models iterating through each hyperpar. combination
tic('random hyperparameters')
metrics <- foreach(i = 1:nrow(grid),
                   .packages=c('randomForest', 'sf', 'dplyr', 'pROC', 'caret'),
                   .combine='rbind') %:% # iterations through hyperpar combinations
  
  foreach(j = 1:5,
          .combine='rbind') %dopar% { # iterations through random cv folds (5)
            
            # set unique seed for each fold
            set.seed(j**2)
            
            # use caret to randomly split train/test data (80-20)
            index <- createDataPartition(dataset$PA, p = 0.8, list = FALSE)
            train <- dataset[index, ]
            test <- dataset[-index, ]  
            
            # extracting sample size
            prNum <- as.numeric(table(train$PA)['1']) # number of presences
            
            # set sample size as vector, same size for 0 and for 1
            samsize <- c('0' = prNum, '1' = prNum)
            
            # train the model
            model <- randomForest(PA ~ . -X -Y -Species,
                                  data = train,
                                  sampsize = samsize,
                                  do.classif = TRUE,
                                  mtry = grid[i,1],
                                  ntree = grid[i,2],
                                  nodesize = grid[i,3])
            
            # generate training predictions (only select probs of presence)
            preds<-predict(model, type='prob')
            train$pred<-preds[,2]
            
            # create ROC curve to calculate thresholds
            roc_curve <- roc(train$PA, train$pred)
            
            # calculate thresholds and binarize preds
            # maxTSS
            maxTSS <- coords(roc_curve, 'best', maximize='tss')$threshold
            train$pred <- ifelse(train$pred > maxTSS, 1, 0)
            
            # AUC training per thresholds
            auc_train <- pROC::auc(as.numeric(train$PA), as.numeric(train$pred))
            
            # bin preds as factors because needed by confusionmatrix
            train$pred <- factor(train$pred, levels = c('0', '1'))
            
            # confusion matrix of the training per thresholds
            cm_train <- confusionMatrix(train$pred, train$PA, '1')
            
            # TSS training per thresholds
            tss_train <- cm_train$byClass[['Sensitivity']] + cm_train$byClass[['Specificity']] - 1
            
            # generate testing predictions and binaries from maxSSS
            preds <- data.frame(pred = predict(model, newdata = test, type ='prob'))
            test$pred <- preds[,2]
            test$pred <- ifelse(test$pred > maxTSS, 1, 0)
            
            # AUC test
            auc_test <- pROC::auc(as.numeric(test$PA), as.numeric(test$pred))
            
            # bin preds as factors because needed by confusionmatrix
            test$pred <- factor(test$pred, levels = c('0', '1'))
            
            # confusion matrix of the testing per thresholds
            cm_test <- confusionMatrix(test$pred, test$PA, '1')
            
            # TSS testing per thresholds
            tss_test <- cm_test$byClass[['Sensitivity']] + cm_test$byClass[['Specificity']] - 1
            
            # store all the relevant values in the metrics dataframe
            # column_names <- c('mtry', 'ntree', 'nodesize','thr_value',
            # 'Sens_train', 'F1_train','B.Accuracy_train','TSS_train','AUC_train',
            # 'Sens_test', 'F1_test', 'B.Accuracy_test','TSS_test','AUC_test')
            
            metrics <- c(i, j, grid[i,1], grid[i,2], grid[i,3], maxTSS,
                         cm_train$byClass[['Sensitivity']],
                         cm_train$byClass[['F1']],
                         cm_train$byClass[['Balanced Accuracy']],
                         tss_train, auc_train, # training
                         cm_test$byClass[['Sensitivity']],
                         cm_test$byClass[['F1']],
                         cm_test$byClass[['Balanced Accuracy']],
                         tss_test, auc_test # testing
            )
            return(metrics)
          }
toc() 

# assign col names
metrics <- as.data.frame(metrics)
colnames(metrics) <- column_names
write.csv(metrics, 'results/hyperparameters_ortega.csv', row.names = F)
# metrics <- read.csv('results/hyperparameters_ortega.csv')

metrics_mean <- metrics %>%
  group_by(combination) %>%
  summarise(across(-fold, mean))
write.csv(metrics_mean, 'results/hyperparameters_ortega_mean.csv', row.names = F)
# metrics_mean <- read.csv('results/hyperparameters_ortega_mean.csv')

metrics_mean$sum <- rowSums(metrics_mean[,5:15])
metrics_mean[which.max(metrics_mean$sum),]

metrics_mean$sens_diff <- abs(metrics_mean$Sens_test-metrics_mean$Sens_train)
metrics_mean[which.min(metrics_mean$sens_diff),]

# min sens diff and maximum sum reached at: nodesize=70, mtry=4; ntree=2000

parallel::stopCluster(cl = my.cluster)
