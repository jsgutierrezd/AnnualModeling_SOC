#============================================================================
# Proc03b_DataModelingP2D1 ------------------------------------------------
#============================================================================
rm(list = ls())
Sys.setenv(language="EN")


# 1) Working directory ----------------------------------------------------


setwd("~/AARHUS_PhD/DSMactivities/MonitoringGridData/SOCModeling/AnnualModeling_SOC")


# 2) Libraries ------------------------------------------------------------

pckg <- c('caret',     
          'magrittr',
          'tidyr',
          'readr',
          'dplyr',
          'MASS',
          'parallel',
          'doParallel',
          'e1071',
          'ranger',
          'hydroGOF',
          'Boruta',
          'randomForest',
          'xgboost',
          'soilassessment',
          'reshape',
          'caretEnsemble',
          'terra',
          'Hmisc',
          'snow',
          'quantregForest',
          'raster',
          'prospectr',
          'solitude'
)

usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE)
  require(p, character.only = TRUE)
}
lapply(pckg,usePackage)

goof <- function(observed,predicted, plot.it = FALSE, type = "DSM"){
  # Coefficient of determination
  rLM <- lm(predicted ~ observed)
  R2 <- as.matrix(summary(rLM)$adj.r.squared)
  # Standard error of prediction ˆ2 or MSE
  SEP2 <- mean((observed - predicted)^2)
  # Standard error of prediction or RMSE
  SEP <- sqrt(SEP2)
  #Bias
  bias <- mean(predicted) - mean(observed)
  # residual variance
  SEP2c <- sum(((predicted - bias - observed)^2) / length(observed))
  SEPc <- sqrt(SEP2c)
  # ratio of performance to deviation
  RPD <- sd(observed) / SEP
  # Ratio of performance to interquartile distance
  IQ <- c(quantile(observed))[3] - c(quantile(observed))[2]
  RPIQ <- IQ / SEP
  # Concordance
  mx <- mean(observed)
  my <- mean(predicted)
  s2x <- var(observed)
  s2y <- var(predicted)
  sxy <- mean((observed-mx) * (predicted-my))
  ccc <- 2 * sxy / (s2x + s2y + (mx - my)^2)
  if (plot.it==TRUE){
    eqscplot(observed, predicted)
    abline(a = 0, b = 1, col = "brown4")
  }
  if (type == "DSM"){
    gf <- data.frame(R2 = R2,
                     concordance = ccc,
                     MSE = SEP2,
                     RMSE = SEP,
                     bias = bias,
                     row.names = NULL
    )
  }
  else if (type == "spec"){
    gf <- data.frame(R2 = R2,
                     concordance = ccc,
                     MSE = SEP2,
                     RMSE = SEP,
                     bias = bias,
                     MSEc = SEP2c,
                     RMSEc = SEPc,
                     RPD = RPD,
                     RPIQ = RPIQ,
                     row.names = NULL
    )
  }
  else {
    stop("ERROR: Revise the type of output you require. Select from either DSM or spec")
  }
  return(gf)
}

# 3) Data loading ---------------------------------------------------------

dataP2D1 <- read_delim("RegMat_P2D1.csv",
                       delim = ",") %>% 
  mutate(SOC_1997t = log(SOC_1997)) %>% na.omit %>% data.frame
summary(dataP2D1)
# dataP2D1$bsiP10_1985 <- ifelse(dataP2D1$bsiP10_1985=="Inf"|dataP2D1$bsiP10_1985=="-Inf",NA,dataP2D1$bsiP10_1985)
# dataP2D1$strP10_1985 <- ifelse(dataP2D1$strP10_1985=="Inf"|dataP2D1$strP10_1985=="-Inf",NA,dataP2D1$strP10_1985)
# dataP2D1$strP10_1986 <- ifelse(dataP2D1$strP10_1986=="Inf"|dataP2D1$strP10_1986=="-Inf",NA,dataP2D1$strP10_1986)
dataP2D1$geology_9 <- NULL
dataP2D1$geology_3 <- NULL
dataP2D1$geology_11 <- NULL
dataP2D1 <- dataP2D1 %>% na.omit
summary(dataP2D1)
names(dataP2D1)



# 4) Correlation matrix ---------------------------------------------------
# 
# corMat <- cor(as.matrix(dataP2D1[,-c(6:23)])) %>% na.omit %>% as.data.frame
# corMat <- dataP2D1[,-c(6:23)] %>% na.omit %>% cor %>% data.frame
# 
# write_csv(corMat,"CorMatRegMat_P2D1.csv")
# 
# names(dataP2D1)
# hist(log(dataP2D1$SOC_1997))

# 5) Data splitting -------------------------------------------------------

set.seed(1000)
inTrain <- createDataPartition(y = dataP2D1$SOC_1997t, p = .70, list = FALSE) # Random
# inTrain <- kenStone(dataP2D1, k = nrow(dataP2D1)*0.70, metric = "mahal") # Kennard Stone
# data1 <- dataP2D1[,c("SOC_1997t")]
# set.seed(58)
# indx <- clhs(data1, size = round(nrow(data1)*0.7),
#              progress = T, iter = 1000,use.cpp = F,  simple = FALSE) # CLHS

train_data <- dataP2D1[ inTrain,] %>% data.frame #Random
# train_data <- dataP2D1[ inTrain$model,] %>% data.frame #Kennard Stone
# train_data <- dataP2D1[indx$index_samples,] %>% data.frame  # CLHS

names(train_data)
y_train <- train_data[,143]
x_train <- train_data[,c(2:44,105:142)]

max_train <- apply(x_train, 2, max)
min_train <- apply(x_train, 2, min)
x_train <- scale(x_train, center = min_train, scale = max_train-min_train)
x_train <- data.frame(SOC_1997t=y_train,x_train)

summary(x_train)
# PCA on spectral indices ----------------------------------------------
names(train_data)
pca1996Med<-prcomp(train_data[,c(45:54)], scale=TRUE) 
summary(pca1996Med)
(corvar <- pca1996Med$rotation %*% diag(pca1996Med$sdev))
Pred.pcs<-predict(pca1996Med,train_data[,c(45:54)])
x_train$PCA1_1996Med=Pred.pcs[,1] 
x_train$PCA2_1996Med=Pred.pcs[,2]
x_train$PCA3_1996Med=Pred.pcs[,3] 

pca1996P10<-prcomp(train_data[,c(55:64)], scale=TRUE) 
summary(pca1996P10)
(corvar <- pca1996P10$rotation %*% diag(pca1996P10$sdev))
Pred.pcs<-predict(pca1996P10,train_data[,c(55:64)])
x_train$PCA1_1996P10=Pred.pcs[,1] 
x_train$PCA2_1996P10=Pred.pcs[,2]
x_train$PCA3_1996P10=Pred.pcs[,3]

pca1996P90<-prcomp(train_data[,c(65:74)], scale=TRUE) 
summary(pca1996P90)
(corvar <- pca1996P90$rotation %*% diag(pca1996P90$sdev))
Pred.pcs<-predict(pca1996P90,train_data[,c(65:74)])
x_train$PCA1_1996P90=Pred.pcs[,1] 
x_train$PCA2_1996P90=Pred.pcs[,2]
x_train$PCA3_1996P90=Pred.pcs[,3]

pca1997Med<-prcomp(train_data[,c(75:84)], scale=TRUE) 
summary(pca1997Med)
(corvar <- pca1997Med$rotation %*% diag(pca1997Med$sdev))
Pred.pcs<-predict(pca1997Med,train_data[,c(75:84)])
x_train$PCA1_1997Med=Pred.pcs[,1] 
x_train$PCA2_1997Med=Pred.pcs[,2]
x_train$PCA3_1997Med=Pred.pcs[,3]

pca1997P10<-prcomp(train_data[,c(85:94)], scale=TRUE) 
summary(pca1997P10)
(corvar <- pca1997P10$rotation %*% diag(pca1997P10$sdev))
Pred.pcs<-predict(pca1997P10,train_data[,c(85:94)])
x_train$PCA1_1997P10=Pred.pcs[,1] 
x_train$PCA2_1997P10=Pred.pcs[,2]
x_train$PCA3_1997P10=Pred.pcs[,3]

pca1997P90<-prcomp(train_data[,c(95:104)], scale=TRUE) 
summary(pca1997P90)
(corvar <- pca1997P90$rotation %*% diag(pca1997P90$sdev))
Pred.pcs<-predict(pca1997P90,train_data[,c(95:104)])
x_train$PCA1_1997P90=Pred.pcs[,1] 
x_train$PCA2_1997P90=Pred.pcs[,2]
x_train$PCA3_1997P90=Pred.pcs[,3]

x_train

test_data <- dataP2D1[-inTrain,] #Random
# test_data <- dataP2D1[inTrain$test,] # Kennard Stone
# test_data <- dataP2D1[-indx$index_samples,] # CLHS

y_test <- test_data[,143]
x_test <- test_data[c(2:44,105:142)]
x_test <- scale(x_test, center = min_train, scale = max_train-min_train)
x_test <- data.frame(SOC_1997t=y_test,x_test)


Pred.pcs<-predict(pca1996Med,test_data[,c(45:54)])
x_test$PCA1_1996Med=Pred.pcs[,1] 
x_test$PCA2_1996Med=Pred.pcs[,2]
x_test$PCA3_1996Med=Pred.pcs[,3] 

Pred.pcs<-predict(pca1996P10,test_data[,c(55:64)])
x_test$PCA1_1996P10=Pred.pcs[,1] 
x_test$PCA2_1996P10=Pred.pcs[,2]
x_test$PCA3_1996P10=Pred.pcs[,3]

Pred.pcs<-predict(pca1996P90,test_data[,c(65:74)])
x_test$PCA1_1996P90=Pred.pcs[,1] 
x_test$PCA2_1996P90=Pred.pcs[,2]
x_test$PCA3_1996P90=Pred.pcs[,3]

Pred.pcs<-predict(pca1997Med,test_data[,c(75:84)])
x_test$PCA1_1997Med=Pred.pcs[,1] 
x_test$PCA2_1997Med=Pred.pcs[,2]
x_test$PCA3_1997Med=Pred.pcs[,3]

Pred.pcs<-predict(pca1997P10,test_data[,c(85:94)])
x_test$PCA1_1997P10=Pred.pcs[,1] 
x_test$PCA2_1997P10=Pred.pcs[,2]
x_test$PCA3_1997P10=Pred.pcs[,3]

Pred.pcs<-predict(pca1997P90,test_data[,c(95:104)])
x_test$PCA1_1997P90=Pred.pcs[,1] 
x_test$PCA2_1997P90=Pred.pcs[,2]
x_test$PCA3_1997P90=Pred.pcs[,3]

summary(x_train)



summary(x_test)



train_data <- x_train
test_data <- x_test

summary(train_data)
train_data <- train_data[complete.cases(train_data),]
summary(test_data)
test_data <- test_data[complete.cases(test_data),]


# Novelty detection -------------------------------------------------------

# svm.model<-svm(train_data,y=NULL,
#                type='one-classification',
#                nu=0.2,
#                scale=F,
#                kernel="polynomial")
# 
# 
# svm.predtrain<-predict(svm.model,train_data)
# table(svm.predtrain)
# 
# # svm.predtest<-predict(svm.model,test_data)
# # table(svm.predtest)
# 
# 
# test_data <- test_data[svm.predtest=="TRUE",]
# library(solitude)
# 
# iso = isolationForest$new()
# iso$fit(train_data)
# 
# library(h2o)
# h2o.init()
# train_data <- as.h2o(train_data)
# test_data <- as.h2o(test_data)
# model <- h2o.isolationForest(training_frame = train_data,
#                              sample_rate = 0.1,
#                              max_depth = 20,
#                              ntrees = 50)
# score <- h2o.predict(model, test_data)
# result_pred <- score$predict
# ln_pred <- h2o.predict_leaf_node_assignment(model, test_data)
# score <- as.data.frame(result_pred)
# res <- data.frame(Instance=as.character(paste0("I",1:68)),score=score)
# res$label <- ifelse(res$predict>=0.5,"Outlier","Normal")
# res
# table(res$label)
# 
# 
# test_data <- as.data.frame(test_data)
# train_data <- as.data.frame(train_data)
# 
# test_data <- test_data[res$label=="Normal",]




# 6) Features selection ---------------------------------------------------

# 6.1) Boruta algorithm ---------------------------------------------------
names(train_data)
train_data <- train_data %>% na.omit %>% data.frame
{
  start <- Sys.time()
  set.seed(1841)
  (bor <- Boruta(x = train_data[,c(2:100)],
                 y = train_data[,1], 
                 #data = train_data, 
                 doTrace = 0, 
                 ntree = 500,
                 maxRuns=500))
  plot(bor, xlab = "", xaxt = "n")
  lz<-lapply(1:ncol(bor$ImpHistory),function(i)
    bor$ImpHistory[is.finite(bor$ImpHistory[,i]),i])
  names(lz) <- colnames(bor$ImpHistory)
  Labels <- sort(sapply(lz,median))
  axis(side = 1,las=2,labels = names(Labels),
       at = 1:ncol(bor$ImpHistory), cex.axis = 0.7)
  print(Sys.time() - start)
}

print(bor)
names(bor$finalDecision[bor$finalDecision %in% c("Confirmed")])

bor <- TentativeRoughFix(bor)
print(bor)
names(bor$finalDecision[bor$finalDecision %in% c("Confirmed")])

preds <- names(bor$finalDecision[bor$finalDecision %in% c("Confirmed")])

# saveRDS(preds,"Outputs/NamesPreds/P2D1/PredictorsP2D1_01082022.rds")

# 7) Model fitting --------------------------------------------------------

fm <- as.formula(paste("SOC_1997t ~", paste0(preds,
                                             collapse = "+")))
fm

# 7.1) Randon forest - Ranger ---------------------------------------------

rctrlG <- trainControl(method = "repeatedcv",
                       number = 20,
                       repeats = 20,
                       returnResamp = "all",
                       search = "grid"
)

grid <- expand.grid(mtry = c(3,5,7),
                    splitrule = c("variance", "extratrees"),
                    min.node.size = c(2,3,4,5)
)

set.seed(1629)
names(train_data)

model_rf <- caret::train(fm,
                  data=train_data,
                  method = "ranger",
                  trControl = rctrlG,
                  tuneGrid = grid,
                  num.trees = 500,
                  importance = "impurity"
)
model_rf
model_rf$finalModel$r.squared
model_rf$bestTune

mod_imp <- varImp(model_rf)
plot(mod_imp, top = 20)


pred_rf <- predict(model_rf, newdata = test_data[,preds])

(rf.goof <- goof(observed = exp(test_data$SOC_1997t), predicted = exp(pred_rf)))


# 7.2) SVM ----------------------------------------------------------------

tuneResult <- tune(svm, fm, data = train_data,
                   ranges = list(epsilon = seq(0.1,0.3,0.02),
                                 cost = c(5,7,15)))

model_svm <- tuneResult$best.model
print(model_svm)
pred_svm <- predict(model_svm, newdata = test_data[,names(bor$finalDecision[bor$finalDecision %in% c("Confirmed")])])
(svm.goof <- goof(observed = test_data$SOC_1997t, predicted = pred_svm))


# 7.3) Multiple linear regression -----------------------------------------

model_lm <- lm(fm,
               data=train_data)
summary(model_lm)       
model_lm_up <- stepAIC(model_lm,direction="both")

summary(model_lm_up)

pred_lm <- predict(model_lm_up, newdata = test_data[,names(bor$finalDecision[bor$finalDecision %in% c("Confirmed")])])
(lm.goof <- goof(observed = test_data$SOC_1997t, predicted = pred_lm))




# 7.4) Cubist -------------------------------------------------------------

ctrl <- trainControl(method = "boot",
                     summaryFunction = defaultSummary,
                     selectionFunction = "best"
)
grid <- expand.grid(committees = c(1, 20, 50), 
                    neighbors = c(1, 5, 9))
set.seed(49)
model_cubist <- train(fm,
                      data=train_data,
                      method = "cubist",
                      tuneGrid = grid,
                      trControl=ctrl
)

model_cubist
summary(model_cubist)

pred_cub <- predict(model_cubist, newdata = test_data[,names(bor$finalDecision[bor$finalDecision %in% c("Confirmed")])])
(cub.goof <- goof(observed = exp(test_data$SOC_1997t), predicted = exp(pred_cub)))



# 7.5) Extreme gradient boosting ------------------------------------------

xgb_trcontrol = trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 10,
  allowParallel = TRUE,
  verboseIter = FALSE,
  returnData = FALSE
)

xgbGrid <- expand.grid(nrounds = c(100,200),  
                       max_depth = c(2,4,6,8),
                       colsample_bytree = seq(0.5, 0.9, length.out = 5),
                       eta = 0.3,
                       gamma=1,
                       min_child_weight = 1,
                       subsample = 0.5
)


set.seed(1105) 
xgb_model = train(fm,
                  data=train_data,
                  method = "xgbTree",
                  tuneGrid = xgbGrid,
                  trControl=xgb_trcontrol
)


xgb_model
summary(xgb_model)

pred_xgb <- predict(xgb_model, newdata = test_data[,names(bor$finalDecision[bor$finalDecision %in% c("Confirmed")])])
(xgb.goof <- goof(observed = exp(test_data$SOC_1997t), predicted = exp(pred_xgb)))



# 7.6) Quantile random forest ---------------------------------------------

set.seed(1841)
model_qrf <- quantregForest(y = train_data[,"SOC_1997t"],
                            x = train_data[,preds],
                            data=train_data,
                            keep.inbag=TRUE,
                            mtry = as.numeric(model_rf$bestTune))
model_qrf
importance(model_qrf,type = 2)

pred_qrf <- predict(model_qrf, newdata = test_data[,preds])
(qrf.goof <- goof(observed = exp(test_data$SOC_1997t), predicted = exp(pred_qrf[,2])))

# saveRDS(model_qrf,"Outputs/Models/P2D1/ModelP2D1qrf_010822.rds")


# 7.7) Caret model ensemble -----------------------------------------------

data.frame(pred_cubist,pred_lm,pred_svm,pred_qrf[,2],pred_rf,pred_xgb) %>% cor
ls(getModelInfo())
my_control <- trainControl(method = "repeatedcv",
                           number = 20,
                           repeats = 20,
                           returnResamp = "all",
                           search = "grid"
)

model_list <- caretList(
  fm, data=train_data,
  trControl=my_control,
  methodList=c("xgbTree", "svmLinear")
)

glm_ensemble <- caretStack(
  model_list,
  method="glm",
  metric="RMSE")

pred_ens <- predict(glm_ensemble, newdata=test_data)
ens.goof <- gof(sim = pred_ens,obs = test_data$SOC_1986t)

ens.goof



# 8) Spatial prediction ---------------------------------------------------

fm

covP2 <- c(rast("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/StatPreds.tif"),
           rast("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/DynPredsP1996_1997.tif"))
names(covP2) <- c(readRDS("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/NamesStatPreds.rds"),
                  readRDS("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/NamesDynPredsP1996_1997.rds"))

preds %in% names(covP2)
fm

# Not mandatory) Missing covariates - i.e. PCA layers ---------------------

Pred.pcs.layers1 <- predict(covP2[[c(41:50)]],pca1996Med,cores=15)
# Pred.pcs.layers2 <- predict(covP2[[c(71:80)]],pca1997Med,cores=15)
# Pred.pcs.layers3 <- predict(covP2[[c(91:100)]],pca1997P90,cores=15)
#  
 
writeRaster(raster(Pred.pcs.layers1[[1]]),"ExtraCovariates/P2D1/PCA1_1996Med.tif",overwrite=T)


# 8.1) Maps generation ----------------------------------------------------

model_qrf <- readRDS("Outputs/Models/P2D1/ModelP2D1qrf_010822.rds")
print(model_qrf)
model_qrf$importance
covP2 <- stack(stack("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/StatPreds.tif"),
           stack("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/DynPredsP1996_1997.tif"))
names(covP2) <- c(readRDS("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/NamesStatPreds.rds"),
                  readRDS("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/NamesDynPredsP1996_1997.rds"))

fm
PCA1_1996Med <- raster("ExtraCovariates/P2D1/PCA1_1996Med.tif")

covP2 <- stack(covP2,PCA1_1996Med)
names(covP2[[-25]])

covP2 <- covP2[[preds]]

beginCluster(n=detectCores()-2,type='SOCK')

covP2sc <- clusterR(covP2[[-25]], scale, 
                    args=list(center = min_train[preds[-25]],
                              scale = max_train[preds[-25]]-
                                min_train[preds[-25]]))

covP2sc <- stack(covP2sc,PCA1_1996Med)
names(covP2sc) <- preds

median <- clusterR(covP2sc, predict,
                   args=list(model=model_qrf, what=0.5))
median <- exp(median)


UppL <- clusterR(covP2sc, predict,
                 args=list(model=model_qrf, what=0.95))
UppL <- exp(UppL)


LowL <- clusterR(covP2sc, predict,
                 args=list(model=model_qrf, what=0.05))
LowL <- exp(LowL)


mean <- clusterR(covP2sc, predict,
                 args=list(model=model_qrf, what=mean))
mean <- exp(mean)

sd <- clusterR(covP2sc, predict,
               args=list(model=model_qrf, what=sd))
sd <- exp(sd)

writeRaster(median,"Outputs/Layers/P2D1/ModelP2D1QrfMedian_010822.tif",overwrite=T)
writeRaster(UppL,"Outputs/Layers/P2D1/ModelP2D1QrfUppL_010822.tif",overwrite=T)
writeRaster(LowL,"Outputs/Layers/P2D1/ModelP2D1QrfLowL_010822.tif",overwrite=T)
writeRaster(mean,"Outputs/Layers/P2D1/ModelP2D1QrfMean_010822.tif",overwrite=T)
writeRaster(sd,"Outputs/Layers/P2D1/ModelP2D1QrfSd_010822.tif",overwrite=T)

endCluster()


