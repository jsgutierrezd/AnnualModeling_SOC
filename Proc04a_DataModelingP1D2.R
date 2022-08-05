#============================================================================
# Proc04a_DataModelingP1D2 ------------------------------------------------
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
          'clhs',
          'prospectr',
          'factoextra'
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

dataP1D2 <- read_delim("RegMat_P1D2.csv",
                       delim = ",") %>% 
  mutate(SOC_1986t = log(SOC_1986)) %>% na.omit
summary(dataP1D2)
dataP1D2$bsiP10_1985 <- ifelse(dataP1D2$bsiP10_1985=="Inf"|dataP1D2$bsiP10_1985=="-Inf",NA,dataP1D2$bsiP10_1985)
dataP1D2$strP10_1985 <- ifelse(dataP1D2$strP10_1985=="Inf"|dataP1D2$strP10_1985=="-Inf",NA,dataP1D2$strP10_1985)
dataP1D2$strP10_1986 <- ifelse(dataP1D2$strP10_1986=="Inf"|dataP1D2$strP10_1986=="-Inf",NA,dataP1D2$strP10_1986)
dataP1D2 <- dataP1D2 %>% na.omit
summary(dataP1D2)
names(dataP1D2)


# 4) Correlation matrix ---------------------------------------------------

# corMat <- cor(as.matrix(dataP1D2[,-c(6:23,121:127)])) %>% na.omit %>% as.data.frame
# corMat <- dataP1D2[,-c(6:23,121:127)] %>% na.omit %>% cor %>% data.frame
# 
# write_csv(corMat,"CorMatRegMat_P12.csv")
# 
# names(dataP1D2)
# hist(log(dataP1D2$SOC_1986))

# 5) Data splitting -------------------------------------------------------

set.seed(840)
#inTrain <- createDataPartition(y = dataP1D2$SOC_1986t, p = .70, list = FALSE) # Random
inTrain <- kenStone(dataP1D2, k = nrow(dataP1D2)*0.70, metric = "mahal") # Kennard Stone
# data1 <- dataP1D2[,c("SOC_1986t")]
# set.seed(58)
# indx <- clhs(data1, size = round(nrow(data1)*0.7),
#              progress = T, iter = 1000,use.cpp = F,  simple = FALSE) # CLHS




# train_data <- dataP1D2[ inTrain,] %>% data.frame #Random
train_data <- dataP1D2[ inTrain$model,] %>% data.frame #Kennard Stone
# train_data <- dataP1D2[indx$index_samples,] # CLHS
y_train <- train_data[,140]
x_train <- train_data[,c(2:41,102:139)]
max_train <- apply(x_train, 2, max)
min_train <- apply(x_train, 2, min)
x_train <- scale(x_train, center = min_train, scale = max_train-min_train)
x_train <- data.frame(SOC_1986t=y_train,x_train)

# 6) PCA on spectral indices ----------------------------------------------
names(train_data)
pca1985Med<-prcomp(train_data[,c(42:51)], scale=TRUE) 
summary(pca1985Med)
(corvar <- pca1985Med$rotation %*% diag(pca1985Med$sdev))
Pred.pcs<-predict(pca1985Med,train_data[,c(42:51)])
x_train$PCA1_1985Med=Pred.pcs[,1] 
x_train$PCA2_1985Med=Pred.pcs[,2]
x_train$PCA3_1985Med=Pred.pcs[,3] 

pca1985P10<-prcomp(train_data[,c(52:61)], scale=TRUE) 
summary(pca1985P10)
(corvar <- pca1985P10$rotation %*% diag(pca1985P10$sdev))
Pred.pcs<-predict(pca1985P10,train_data[,c(52:61)])
x_train$PCA1_1985P10=Pred.pcs[,1] 
x_train$PCA2_1985P10=Pred.pcs[,2]
x_train$PCA3_1985P10=Pred.pcs[,3]

pca1985P90<-prcomp(train_data[,c(62:71)], scale=TRUE) 
summary(pca1985P90)
(corvar <- pca1985P90$rotation %*% diag(pca1985P90$sdev))
Pred.pcs<-predict(pca1985P90,train_data[,c(62:71)])
x_train$PCA1_1985P90=Pred.pcs[,1] 
x_train$PCA2_1985P90=Pred.pcs[,2]
x_train$PCA3_1985P90=Pred.pcs[,3]

pca1986Med<-prcomp(train_data[,c(72:81)], scale=TRUE) 
summary(pca1986Med)
(corvar <- pca1986Med$rotation %*% diag(pca1986Med$sdev))
Pred.pcs<-predict(pca1986Med,train_data[,c(72:81)])
x_train$PCA1_1986Med=Pred.pcs[,1] 
x_train$PCA2_1986Med=Pred.pcs[,2]
x_train$PCA3_1986Med=Pred.pcs[,3]

pca1986P10<-prcomp(train_data[,c(82:91)], scale=TRUE) 
summary(pca1986P10)
(corvar <- pca1986P10$rotation %*% diag(pca1986P10$sdev))
Pred.pcs<-predict(pca1986P10,train_data[,c(82:91)])
x_train$PCA1_1986P10=Pred.pcs[,1] 
x_train$PCA2_1986P10=Pred.pcs[,2]
x_train$PCA3_1986P10=Pred.pcs[,3]

pca1986P90<-prcomp(train_data[,c(92:101)], scale=TRUE) 
summary(pca1986P90)
(corvar <- pca1986P90$rotation %*% diag(pca1986P90$sdev))
Pred.pcs<-predict(pca1986P90,train_data[,c(92:101)])
x_train$PCA1_1986P90=Pred.pcs[,1] 
x_train$PCA2_1986P90=Pred.pcs[,2]
x_train$PCA3_1986P90=Pred.pcs[,3]

x_train

# test_data <- dataP1D2[-inTrain,] #Random
test_data <- dataP1D2[inTrain$test,] # Kennard Stone
# test_data <- dataP1D2[-indx$index_samples,] # CLHS
y_test <- test_data[,140]
x_test <- test_data[c(2:41,102:139)]
x_test <- scale(x_test, center = min_train, scale = max_train-min_train)
x_test <- data.frame(SOC_1986t=y_test,x_test)


Pred.pcs<-predict(pca1985Med,test_data[,c(42:51)])
x_test$PCA1_1985Med=Pred.pcs[,1] 
x_test$PCA2_1985Med=Pred.pcs[,2]
x_test$PCA3_1985Med=Pred.pcs[,3] 

Pred.pcs<-predict(pca1985P10,test_data[,c(52:61)])
x_test$PCA1_1985P10=Pred.pcs[,1] 
x_test$PCA2_1985P10=Pred.pcs[,2]
x_test$PCA3_1985P10=Pred.pcs[,3]

Pred.pcs<-predict(pca1985P90,test_data[,c(62:71)])
x_test$PCA1_1985P90=Pred.pcs[,1] 
x_test$PCA2_1985P90=Pred.pcs[,2]
x_test$PCA3_1985P90=Pred.pcs[,3]

Pred.pcs<-predict(pca1986Med,test_data[,c(72:81)])
x_test$PCA1_1986Med=Pred.pcs[,1] 
x_test$PCA2_1986Med=Pred.pcs[,2]
x_test$PCA3_1986Med=Pred.pcs[,3]

Pred.pcs<-predict(pca1986P10,test_data[,c(82:91)])
x_test$PCA1_1986P10=Pred.pcs[,1] 
x_test$PCA2_1986P10=Pred.pcs[,2]
x_test$PCA3_1986P10=Pred.pcs[,3]

Pred.pcs<-predict(pca1986P90,test_data[,c(92:101)])
x_test$PCA1_1986P90=Pred.pcs[,1] 
x_test$PCA2_1986P90=Pred.pcs[,2]
x_test$PCA3_1986P90=Pred.pcs[,3]

x_train$geology_9 <- NULL
x_train$geology_11 <- NULL
x_train$geology_2 <- NULL
x_train$geology_3 <- NULL

x_test$geology_9 <- NULL
x_test$geology_11 <- NULL
x_test$geology_2 <- NULL
x_test$geology_3 <- NULL

train_data <- x_train
test_data <- x_test

summary(train_data)
train_data <- train_data[complete.cases(train_data),]
summary(test_data)
test_data <- test_data[complete.cases(test_data),]

train_data$label <- "Train"
test_data$label <- "Test"
df <- rbind(train_data,test_data)
ggplot(df, aes(x=SOC_1986t, color=label)) +
  geom_density()
train_data$label <- NULL
test_data$label <- NULL

# Novelty detection -------------------------------------------------------

svm.model<-svm(train_data,y=NULL,
               type='one-classification',
               nu=0.1,
               scale=F,
               kernel="radial")


svm.predtrain<-predict(svm.model,train_data)
table(svm.predtrain)

svm.predtest<-predict(svm.model,test_data)
table(svm.predtest)

train_data <- train_data[svm.predtrain=="TRUE",]
test_data <- test_data[svm.predtest=="TRUE",]



# 6) Features selection ---------------------------------------------------


# 6.1) Boruta algorithm ---------------------------------------------------
names(train_data)
train_data <- train_data %>% na.omit %>% data.frame
{
  start <- Sys.time()
  set.seed(845)
  (bor <- Boruta(x = train_data[,c(2:93)],
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

#saveRDS(preds,"Outputs/NamesPreds/P1D2/PredictorsP1D2_01082022.rds")

# 7) Model fitting --------------------------------------------------------
names(train_data)
fm <- as.formula(paste("SOC_1986t ~", paste0(preds,
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

set.seed(850)

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

(rf.goof <- goof(observed = exp(test_data$SOC_1986t), predicted = exp(pred_rf)))


# 7.2) SVM ----------------------------------------------------------------

tuneResult <- tune(svm, fm, data = train_data,
                   ranges = list(epsilon = seq(0.1,0.3,0.02),
                                 cost = c(5,7,15)))

model_svm <- tuneResult$best.model
print(model_svm)
pred_svm <- predict(model_svm, newdata = test_data[,names(bor$finalDecision[bor$finalDecision %in% c("Confirmed")])])

(svm.goof <- goof(observed = exp(test_data$SOC_1986t), predicted = exp(pred_svm)))



# 7.3) Multiple linear regression -----------------------------------------

model_lm <- lm(fm,
               data=train_data)
summary(model_lm)       
model_lm_up <- stepAIC(model_lm,direction="both")

summary(model_lm_up)

pred_lm <- predict(model_lm, newdata = test_data[,names(bor$finalDecision[bor$finalDecision %in% c("Confirmed")])])

(lm.goof <- goof(observed = exp(test_data$SOC_1986t), predicted = exp(pred_lm)))




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

(cubist.goof <- goof(observed = exp(test_data$SOC_1986t), predicted = exp(pred_cubist)))



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

(xgb.goof <- goof(observed = exp(test_data$SOC_1986t), predicted = exp(pred_xgb)))



# 7.6) Quantile random forest ---------------------------------------------

set.seed(915)
model_qrf <- quantregForest(y = train_data[,"SOC_1986t"],
                            x = train_data[,preds],
                            data=train_data,
                            keep.inbag=TRUE,
                            mtry = model_rf$bestTune$mtry)
model_qrf
importance(model_qrf,type = 2)
pred_qrf <- predict(model_qrf, newdata = test_data[,preds])

(qrf.goof <- goof(observed = exp(test_data$SOC_1986t), predicted = exp(pred_qrf[,2])))


#saveRDS(model_qrf,"Outputs/Models/P4D1/ModelP4D1qrf_290722.rds")


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

(ens.goof <- goof(observed = exp(test_data$SOC_1986t), predicted = exp(pred_ens)))



# 8) Spatial prediction ---------------------------------------------------

fm

covP4 <- c(rast("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_1986/YearbyYear/StatPreds.tif"),
           rast("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_1986/YearbyYear/DynPredsP1985_1986.tif"),
           rast("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_1986/YearbyYear/FieldBlock18.tif"),
           rast("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_1986/YearbyYear/FieldBlock19.tif"))
names(covP4) <- c(readRDS("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_1986/YearbyYear/NamesStatPreds.rds"),
                  readRDS("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_1986/YearbyYear/NamesDynPredsP1985_1986.rds"),
                  readRDS("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_1986/YearbyYear/NamesFieldBlock18.rds"),
                  readRDS("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_1986/YearbyYear/NamesFieldBlock19.rds"))

names(bor$finalDecision[bor$finalDecision %in% c("Confirmed")]) %in% names(covP4)
fm

# Not mandatory) Missing covariates - i.e. PCA layers ---------------------
names(covP4)
Pred.pcs.layers1 <- predict(covP4[[c(71:80)]],pca1986Med,cores=15)
Pred.pcs.layers2 <- predict(covP4[[c(91:100)]],pca1986P90,cores=15)



writeRaster(raster(Pred.pcs.layers1[[3]]),"ExtraCovariates/P4D1/PCA3_1986Med.tif",overwrite=T)
writeRaster(raster(Pred.pcs.layers2[[1]]),"ExtraCovariates/P4D1/PCA1_1986P90.tif",overwrite=T)
writeRaster(raster(Pred.pcs.layers2[[2]]),"ExtraCovariates/P4D1/PCA2_1986P90.tif",overwrite=T)

# 8.1) Maps generation ----------------------------------------------------

model_qrf <- readRDS("Outputs/Models/ModelP3D1qrf_280722.rds")
print(model_qrf)

covP4 <- stack(stack("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_1986/YearbyYear/StatPreds.tif"),
               stack("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_1986/YearbyYear/DynPredsP1985_1986.tif"))
names(covP4) <- c(readRDS("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_1986/YearbyYear/NamesStatPreds.rds"),
                  readRDS("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_1986/YearbyYear/NamesDynPredsP1985_1986.rds"))


PCA3_1986Med <- raster("ExtraCovariates/P4D1/PCA3_1986Med.tif")
PCA1_1986P90 <- raster("ExtraCovariates/P4D1/PCA1_1986P90.tif")
PCA2_1986P90 <- raster("ExtraCovariates/P4D1/PCA2_1986P90.tif")

covP4 <- stack(covP4,PCA3_1986Med,PCA1_1986P90,PCA2_1986P90)

covP4 <- covP4[[names(bor$finalDecision[bor$finalDecision %in% c("Confirmed")])]]

names(covP4)
beginCluster(n=detectCores()-2,type='SOCK')

covP4sc <- clusterR(covP4[[-c(22:24)]], scale, 
                    args=list(center = min_train[names(bor$finalDecision[bor$finalDecision %in% c("Confirmed")])][-c(22:24)],
                              scale = max_train[names(bor$finalDecision[bor$finalDecision %in% c("Confirmed")])][-c(22:24)]-
                                min_train[names(bor$finalDecision[bor$finalDecision %in% c("Confirmed")])][-c(22:24)]))

covP4sc <- stack(covP4sc,PCA3_1986Med,PCA1_1986P90,PCA2_1986P90)
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
