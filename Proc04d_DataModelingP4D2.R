#============================================================================
# Proc04c_DataModelingP3D2 ------------------------------------------------
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
          'MVN',
          'rsample',
          'h2o'
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

h2o.no_progress() # turn off h2o progress bars
h2o.init() # launch h2o

# 3) Data loading ---------------------------------------------------------

dataP4D2 <- read_delim("RegMat_P4D2.csv",
                       delim = ",") %>% 
  mutate(SOC_2019t = log(SOC_2019)) %>% na.omit
summary(dataP4D2)
dataP4D2$strP10_2019 <- ifelse(dataP4D2$strP10_2019=="Inf"|dataP4D2$strP10_2019=="-Inf",NA,dataP4D2$strP10_2019)
dataP4D2$strP10_2018 <- ifelse(dataP4D2$strP10_2018=="Inf"|dataP4D2$strP10_2018=="-Inf",NA,dataP4D2$strP10_2018)
# dataP4D2$ndviP10_2018 <- ifelse(dataP4D2$ndviP10_2018=="Inf"|dataP4D2$ndviP10_2018=="-Inf",NA,dataP4D2$ndviP10_2018)

dataP4D2 <- dataP4D2 %>% na.omit
summary(dataP4D2)
names(dataP4D2)



# 4) Correlation matrix ---------------------------------------------------

# corMat <- cor(as.matrix(dataP4D2[,-c(6:23,121)])) %>% na.omit %>% as.data.frame
# corMat <- dataP4D2[,-c(6:23,121:127)] %>% na.omit %>% cor %>% data.frame
# 
# write_csv(corMat,"CorMatRegMat_P2D2.csv")
# 
# names(dataP4D2)
# hist(log(dataP4D2$SOC_2019))


# 5) Data splitting -------------------------------------------------------

set.seed(1604)
# inTrain <- createDataPartition(y = dataP4D2$SOC_2019t, p = .70, list = FALSE) # Random
# inTrain <- kenStone(dataP4D2, k = nrow(dataP4D2)*0.70, metric = "mahal") # Kennard Stone
# data1 <- dataP4D2[,c("SOC_2019t")]
# set.seed(58)
# indx <- clhs(data1, size = round(nrow(data1)*0.7),
#              progress = T, iter = 1000,use.cpp = F,  simple = FALSE) # CLHS

dataP4D2.hex <- as.h2o(dataP4D2,destination_frame ="dataP4D2.hex")

dataP4D2.2 <- h2o.splitFrame(dataP4D2, ratios = 0.7,
                          seed = 1604)

iris.gbm <- h2o.gbm(y = 1, x = 2:5, training_frame =
                      iris.hex, ntrees = 10,
                    6 max_depth = 3,min_rows = 2, learn_rate = 0.2,
                    distribution= "gaussian")

train <- dataP4D2.2[[1]]

train[,1]
test <- dataP4D2.2[[2]]

str(train)
train <- "Train"
test_data$label <- "Test"
df <- rbind(train_data,test_data)

ggplot(df, aes(x=SOC_2019t, color=label)) +
  geom_density()


# train_data <- dataP4D2[ inTrain,] %>% data.frame #Random
train_data <- dataP4D2[ inTrain$model,] %>% data.frame #Kennard Stone
# train_data <- dataP4D2[indx$index_samples,] %>% data.frame  # CLHS

names(train_data)
y_train <- train_data[,128]
x_train <- train_data[,c(2:41,102:127)]
max_train <- apply(x_train, 2, max)
min_train <- apply(x_train, 2, min)
x_train <- scale(x_train, center = min_train, scale = max_train-min_train)
x_train <- data.frame(SOC_2019t=y_train,x_train)

# 6) PCA on spectral indices ----------------------------------------------
names(train_data)
pca2018Med<-prcomp(train_data[,c(42:51)], scale=TRUE) 
summary(pca2018Med)
(corvar <- pca2018Med$rotation %*% diag(pca2018Med$sdev))
Pred.pcs<-predict(pca2018Med,train_data[,c(42:51)])
x_train$PCA1_2018Med=Pred.pcs[,1] 
x_train$PCA2_2018Med=Pred.pcs[,2]
x_train$PCA3_2018Med=Pred.pcs[,3] 

pca2018P10<-prcomp(train_data[,c(52:61)], scale=TRUE) 
summary(pca2018P10)
(corvar <- pca2018P10$rotation %*% diag(pca2018P10$sdev))
Pred.pcs<-predict(pca2018P10,train_data[,c(52:61)])
x_train$PCA1_2018P10=Pred.pcs[,1] 
x_train$PCA2_2018P10=Pred.pcs[,2]
x_train$PCA3_2018P10=Pred.pcs[,3]

pca2018P90<-prcomp(train_data[,c(62:71)], scale=TRUE) 
summary(pca2018P90)
(corvar <- pca2018P90$rotation %*% diag(pca2018P90$sdev))
Pred.pcs<-predict(pca2018P90,train_data[,c(62:71)])
x_train$PCA1_2018P90=Pred.pcs[,1] 
x_train$PCA2_2018P90=Pred.pcs[,2]
x_train$PCA3_2018P90=Pred.pcs[,3]

pca2019Med<-prcomp(train_data[,c(72:81)], scale=TRUE) 
summary(pca2019Med)
(corvar <- pca2019Med$rotation %*% diag(pca2019Med$sdev))
Pred.pcs<-predict(pca2019Med,train_data[,c(72:81)])
x_train$PCA1_2019Med=Pred.pcs[,1] 
x_train$PCA2_2019Med=Pred.pcs[,2]
x_train$PCA3_2019Med=Pred.pcs[,3]

pca2019P10<-prcomp(train_data[,c(82:91)], scale=TRUE) 
summary(pca2019P10)
(corvar <- pca2019P10$rotation %*% diag(pca2019P10$sdev))
Pred.pcs<-predict(pca2019P10,train_data[,c(82:91)])
x_train$PCA1_2019P10=Pred.pcs[,1] 
x_train$PCA2_2019P10=Pred.pcs[,2]
x_train$PCA3_2019P10=Pred.pcs[,3]

pca2019P90<-prcomp(train_data[,c(92:101)], scale=TRUE) 
summary(pca2019P90)
(corvar <- pca2019P90$rotation %*% diag(pca2019P90$sdev))
Pred.pcs<-predict(pca2019P90,train_data[,c(92:101)])
x_train$PCA1_2019P90=Pred.pcs[,1] 
x_train$PCA2_2019P90=Pred.pcs[,2]
x_train$PCA3_2019P90=Pred.pcs[,3]

x_train

# test_data <- dataP4D2[-inTrain,] #Random
test_data <- dataP4D2[inTrain$test,] # Kennard Stone
# test_data <- dataP4D2[-indx$index_samples,] # CLHS

y_test <- test_data[,128]
x_test <- test_data[c(2:41,102:127)]
x_test <- scale(x_test, center = min_train, scale = max_train-min_train)
x_test <- data.frame(SOC_2019t=y_test,x_test)


Pred.pcs<-predict(pca2018Med,test_data[,c(42:51)])
x_test$PCA1_2018Med=Pred.pcs[,1] 
x_test$PCA2_2018Med=Pred.pcs[,2]
x_test$PCA3_2018Med=Pred.pcs[,3] 

Pred.pcs<-predict(pca2018P10,test_data[,c(52:61)])
x_test$PCA1_2018P10=Pred.pcs[,1] 
x_test$PCA2_2018P10=Pred.pcs[,2]
x_test$PCA3_2018P10=Pred.pcs[,3]

Pred.pcs<-predict(pca2018P90,test_data[,c(62:71)])
x_test$PCA1_2018P90=Pred.pcs[,1] 
x_test$PCA2_2018P90=Pred.pcs[,2]
x_test$PCA3_2018P90=Pred.pcs[,3]

Pred.pcs<-predict(pca2019Med,test_data[,c(72:81)])
x_test$PCA1_2019Med=Pred.pcs[,1] 
x_test$PCA2_2019Med=Pred.pcs[,2]
x_test$PCA3_2019Med=Pred.pcs[,3]

Pred.pcs<-predict(pca2019P10,test_data[,c(82:91)])
x_test$PCA1_2019P10=Pred.pcs[,1] 
x_test$PCA2_2019P10=Pred.pcs[,2]
x_test$PCA3_2019P10=Pred.pcs[,3]

Pred.pcs<-predict(pca2019P90,test_data[,c(92:101)])
x_test$PCA1_2019P90=Pred.pcs[,1] 
x_test$PCA2_2019P90=Pred.pcs[,2]
x_test$PCA3_2019P90=Pred.pcs[,3]

x_train$geology_9 <- NULL
x_train$geology_11 <- NULL
x_train$BareSoil18 <- NULL
x_train$bluespot <- NULL

x_test$geology_9 <- NULL
x_test$geology_11 <- NULL
x_test$BareSoil18 <- NULL
x_test$bluespot <- NULL



train_data <- x_train
test_data <- x_test

summary(train_data)
train_data <- train_data[complete.cases(train_data),]
summary(test_data)
test_data <- test_data[complete.cases(test_data),]


train_data$label <- "Train"
test_data$label <- "Test"
df <- rbind(train_data,test_data)

ggplot(df, aes(x=SOC_2019t, color=label)) +
  geom_density()



# 6) Features selection ---------------------------------------------------


# 6.1) Recursive feature elimination --------------------------------------

# names(train_data)
# train_data <- train_data %>% na.omit %>% data.frame
# {
#   set.seed(2023)
#   start <- Sys.time()
#   cl <- parallel::makeCluster(detectCores(), type='PSOCK')
#   registerDoParallel(cl)
#   control2 <- rfeControl(functions=rfFuncs, method="repeatedcv", number=10, repeats=10,allowParallel = TRUE)
#   (rfe <- rfe(x=train_data[,c(2:41,102:139,141:158)], y=train_data[,140], sizes=c(1:150), rfeControl=control2))
#   print(Sys.time() - start)
# }
# 
# #names(train_data[,c(2:25,150:187,189:206)])
# plot(rfe, type=c("g", "o"))
# predictors(rfe)

# 6.2) Boruta algorithm ---------------------------------------------------
names(train_data)
train_data <- train_data %>% na.omit %>% data.frame
{
  start <- Sys.time()
  set.seed(1610)
  (bor <- Boruta(x = train_data[,c(2:81)],
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

saveRDS(preds,"Outputs/NamesPreds/P4D2/PredictorsP4D2_01082022.rds")

# 7) Model fitting --------------------------------------------------------

fm <- as.formula(paste("SOC_2019t ~", paste0(preds,
                                             collapse = "+")))
fm

# 7.1) Randon forest - Ranger ---------------------------------------------

rctrlG <- trainControl(method = "repeatedcv",
                       number = 10,
                       repeats = 10,
                       returnResamp = "all",
                       search = "grid"
)

grid <- expand.grid(mtry = c(3,4,5),
                    splitrule = c("variance", "extratrees"),
                    min.node.size = c(1,2,4)
)

set.seed(1611)

model_rf <- train(fm,
                  data=train_data,
                  method = "ranger",
                  trControl = rctrlG,
                  tuneGrid = grid,
                  num.trees = 500,
                  importance = "impurity"
)
model_rf
model_rf$finalModel
model_rf$bestTune

mod_imp <- varImp(model_rf)
plot(mod_imp, top = 20)


pred_rf <- predict(model_rf, newdata = test_data[,preds])

(rf.goof <- goof(observed = exp(test_data$SOC_2019t), predicted = exp(pred_rf)))


# 7.2) SVM ----------------------------------------------------------------

tuneResult <- tune(svm, fm, data = train_data,
                   ranges = list(epsilon = seq(0.1,0.3,0.02),
                                 cost = c(5,7,15)))

model_svm <- tuneResult$best.model
print(model_svm)
pred_svm <- predict(model_svm, newdata = test_data[,names(bor$finalDecision[bor$finalDecision %in% c("Confirmed")])])

(svm.goof <- goof(observed = exp(test_data$SOC_2019t), predicted = exp(pred_svm)))



# 7.3) Multiple linear regression -----------------------------------------

model_lm <- lm(fm,
               data=train_data)
summary(model_lm)       
model_lm_up <- stepAIC(model_lm,direction="both")

summary(model_lm_up)

pred_lm <- predict(model_lm, newdata = test_data[,names(bor$finalDecision[bor$finalDecision %in% c("Confirmed")])])

(lm.goof <- goof(observed = exp(test_data$SOC_2019t), predicted = exp(pred_lm)))



# 7.4) Cubist -------------------------------------------------------------

ctrl <- trainControl(method = "boot",
                     summaryFunction = defaultSummary,
                     selectionFunction = "best"
)

set.seed(49)
model_cubist <- train(fm,
                      data=train_data,
                      method = "cubist",
                      # tuneGrid = grid,
                      trControl=ctrl
)

model_cubist
summary(model_cubist)

pred_cubist <- predict(model_cubist, newdata = test_data[,names(bor$finalDecision[bor$finalDecision %in% c("Confirmed")])])

(cubist.goof <- goof(observed = exp(test_data$SOC_2019t), predicted = exp(pred_cubist)))



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

(xgb.goof <- goof(observed = exp(test_data$SOC_2019t), predicted = exp(pred_xgb)))



# 7.6) Quantile random forest ---------------------------------------------

set.seed(415)
model_qrf <- quantregForest(y = train_data[,"SOC_2019t"],
                            x = train_data[,preds],
                            data=train_data,
                            keep.inbag=TRUE,
                            mtry = model_rf$bestTune$mtry)
model_qrf
importance(model_qrf,type = 2)
pred_qrf <- predict(model_qrf, newdata = test_data[,preds])

(qrf.goof <- goof(observed = exp(test_data$SOC_2019t), predicted = exp(pred_qrf[,2])))


saveRDS(model_qrf,"Outputs/Models/P4D2/ModelP4D2qrf_010822.rds")


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
  methodList=c("ranger", "lm")
)

glm_ensemble <- caretStack(
  model_list,
  method="glm",
  metric="RMSE")

pred_ens <- predict(glm_ensemble, newdata = test_data[,preds])

(ens.goof <- goof(observed = exp(test_data$SOC_2019t), predicted = exp(pred_ens)))



# 8) Spatial prediction ---------------------------------------------------

fm

covP2 <- c(rast("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/StatPreds.tif"),
           rast("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/DynPredsP2018_2019.tif"))
names(covP2) <- c(readRDS("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/NamesStatPreds.rds"),
                  readRDS("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/NamesDynPredsP2018_2019.rds"))

preds %in% names(covP2)
fm

# Not mandatory) Missing covariates - i.e. PCA layers ---------------------
names(covP2)
Pred.pcs.layers1 <- predict(covP4[[c(41:50)]],pca2018Med,cores=15)
Pred.pcs.layers2 <- predict(covP4[[c(61:70)]],pca2018P90,cores=15)



writeRaster(raster(Pred.pcs.layers1[[3]]),"ExtraCovariates/P4D1/PCA3_2019Med.tif",overwrite=T)
writeRaster(raster(Pred.pcs.layers2[[1]]),"ExtraCovariates/P4D1/PCA1_2019P90.tif",overwrite=T)
writeRaster(raster(Pred.pcs.layers2[[2]]),"ExtraCovariates/P4D1/PCA2_2019P90.tif",overwrite=T)

# 8.1) Maps generation ----------------------------------------------------

model_qrf <- readRDS("Outputs/Models/ModelP3D1qrf_280722.rds")
print(model_qrf)

covP4 <- stack(stack("O:/Tech_AGRO/Jord/Sebastian/Multiannual2019_2019/YearbyYear/StatPreds.tif"),
               stack("O:/Tech_AGRO/Jord/Sebastian/Multiannual2019_2019/YearbyYear/DynPredsP2018_2019.tif"))
names(covP4) <- c(readRDS("O:/Tech_AGRO/Jord/Sebastian/Multiannual2019_2019/YearbyYear/NamesStatPreds.rds"),
                  readRDS("O:/Tech_AGRO/Jord/Sebastian/Multiannual2019_2019/YearbyYear/NamesDynPredsP2018_2019.rds"))


PCA3_2019Med <- raster("ExtraCovariates/P4D1/PCA3_2019Med.tif")
PCA1_2019P90 <- raster("ExtraCovariates/P4D1/PCA1_2019P90.tif")
PCA2_2019P90 <- raster("ExtraCovariates/P4D1/PCA2_2019P90.tif")

covP4 <- stack(covP4,PCA3_2019Med,PCA1_2019P90,PCA2_2019P90)

covP4 <- covP4[[names(bor$finalDecision[bor$finalDecision %in% c("Confirmed")])]]

names(covP4)
beginCluster(n=detectCores()-2,type='SOCK')

covP4sc <- clusterR(covP4[[-c(22:24)]], scale, 
                    args=list(center = min_train[names(bor$finalDecision[bor$finalDecision %in% c("Confirmed")])][-c(22:24)],
                              scale = max_train[names(bor$finalDecision[bor$finalDecision %in% c("Confirmed")])][-c(22:24)]-
                                min_train[names(bor$finalDecision[bor$finalDecision %in% c("Confirmed")])][-c(22:24)]))

covP4sc <- stack(covP4sc,PCA3_2019Med,PCA1_2019P90,PCA2_2019P90)
names(covP4sc) <- names(bor$finalDecision[bor$finalDecision %in% c("Confirmed")])

median <- clusterR(covP4sc, predict,
                   args=list(model=model_qrf, what=0.5))
median <- exp(median)


UppL <- clusterR(covP4sc, predict,
                 args=list(model=model_qrf, what=0.95))
UppL <- exp(UppL)


LowL <- clusterR(covP4sc, predict,
                 args=list(model=model_qrf, what=0.05))
LowL <- exp(LowL)


mean <- clusterR(covP4sc, predict,
                 args=list(model=model_qrf, what=mean))
mean <- exp(mean)



# #quantile probs=0.5
# plot(median)
# plot(UppL)
# plot(LowL)

writeRaster(median,"Outputs/Layers/P4D1/ModelP4D1QrfMedian_290722.tif",overwrite=T)
writeRaster(UppL,"Outputs/Layers/P4D1/ModelP4D1QrfUppL_290722.tif",overwrite=T)
writeRaster(LowL,"Outputs/Layers/P4D1/ModelP4D1QrfLowL_290722.tif",overwrite=T)
writeRaster(mean,"Outputs/Layers/P4D1/ModelP4D1QrfMean_290722.tif",overwrite=T)

endCluster()
