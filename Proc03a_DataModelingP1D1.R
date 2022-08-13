#============================================================================
# Proc03a_DataModelingP1D1 ------------------------------------------------
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
          'prospectr'
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

dataP1D1 <- read_delim("RegMat_P1D1.csv",
                       delim = ",") %>% 
  mutate(SOC_1986t = log(SOC_1986)) %>% na.omit
summary(dataP1D1)
dataP1D1$bsiP10_1985 <- ifelse(dataP1D1$bsiP10_1985=="Inf"|dataP1D1$bsiP10_1985=="-Inf",NA,dataP1D1$bsiP10_1985)
dataP1D1$strP10_1985 <- ifelse(dataP1D1$strP10_1985=="Inf"|dataP1D1$strP10_1985=="-Inf",NA,dataP1D1$strP10_1985)
dataP1D1$strP10_1986 <- ifelse(dataP1D1$strP10_1986=="Inf"|dataP1D1$strP10_1986=="-Inf",NA,dataP1D1$strP10_1986)
dataP1D1$geology_9 <- NULL
dataP1D1$geology_11 <- NULL
dataP1D1 <- dataP1D1 %>% na.omit
names(dataP1D1)




# 4) Correlation matrix ---------------------------------------------------

# corMat <- cor(as.matrix(dataP1D1[,-c(6:23)])) %>% na.omit %>% as.data.frame
# corMat <- dataP1D1[,-c(6:23)] %>% na.omit %>% cor %>% data.frame
# 
# write_csv(corMat,"CorMatRegMat_P1D1.csv")
# 
# names(dataP1D1)
# hist(log(dataP1D1$SOC_1986))

# 5) Data splitting -------------------------------------------------------

set.seed(840)
inTrain <- createDataPartition(y =dataP1D1$SOC_1986t, p = .70, list = FALSE) #Random
# inTrain <- kenStone(dataP1D1, k = nrow(dataP1D1)*0.70, metric = "mahal") # Kennard Stone

names(dataP1D1)
train_data <- dataP1D1[ inTrain,] %>% data.frame #Random
# train_data <- dataP1D1[ inTrain$model,] %>% data.frame #Kennard Stone
y_train <- train_data[,148]
x_train <- train_data[,c(2:45,106:147)]
max_train <- apply(x_train, 2, max)
min_train <- apply(x_train, 2, min)
x_train <- scale(x_train, center = min_train, scale = max_train-min_train)
x_train <- data.frame(SOC_1986t=y_train,x_train)

# 6) PCA on spectral indices ----------------------------------------------
names(train_data)
pca1985Med<-prcomp(train_data[,c(46:55)], scale=TRUE) 
summary(pca1985Med)
(corvar <- pca1985Med$rotation %*% diag(pca1985Med$sdev))
Pred.pcs<-predict(pca1985Med,train_data[,c(46:55)])
x_train$PCA1_1985Med=Pred.pcs[,1] 
x_train$PCA2_1985Med=Pred.pcs[,2]
x_train$PCA3_1985Med=Pred.pcs[,3] 

pca1985P10<-prcomp(train_data[,c(56:65)], scale=TRUE) 
summary(pca1985P10)
(corvar <- pca1985P10$rotation %*% diag(pca1985P10$sdev))
Pred.pcs<-predict(pca1985P10,train_data[,c(56:65)])
x_train$PCA1_1985P10=Pred.pcs[,1] 
x_train$PCA2_1985P10=Pred.pcs[,2]
x_train$PCA3_1985P10=Pred.pcs[,3]

pca1985P90<-prcomp(train_data[,c(66:75)], scale=TRUE) 
summary(pca1985P90)
(corvar <- pca1985P90$rotation %*% diag(pca1985P90$sdev))
Pred.pcs<-predict(pca1985P90,train_data[,c(66:75)])
x_train$PCA1_1985P90=Pred.pcs[,1] 
x_train$PCA2_1985P90=Pred.pcs[,2]
x_train$PCA3_1985P90=Pred.pcs[,3]

pca1986Med<-prcomp(train_data[,c(76:85)], scale=TRUE) 
summary(pca1986Med)
(corvar <- pca1986Med$rotation %*% diag(pca1986Med$sdev))
Pred.pcs<-predict(pca1986Med,train_data[,c(76:85)])
x_train$PCA1_1986Med=Pred.pcs[,1] 
x_train$PCA2_1986Med=Pred.pcs[,2]
x_train$PCA3_1986Med=Pred.pcs[,3]

pca1986P10<-prcomp(train_data[,c(86:95)], scale=TRUE) 
summary(pca1986P10)
(corvar <- pca1986P10$rotation %*% diag(pca1986P10$sdev))
Pred.pcs<-predict(pca1986P10,train_data[,c(86:95)])
x_train$PCA1_1986P10=Pred.pcs[,1] 
x_train$PCA2_1986P10=Pred.pcs[,2]
x_train$PCA3_1986P10=Pred.pcs[,3]

pca1986P90<-prcomp(train_data[,c(96:105)], scale=TRUE) 
summary(pca1986P90)
(corvar <- pca1986P90$rotation %*% diag(pca1986P90$sdev))
Pred.pcs<-predict(pca1986P90,train_data[,c(96:105)])
x_train$PCA1_1986P90=Pred.pcs[,1] 
x_train$PCA2_1986P90=Pred.pcs[,2]
x_train$PCA3_1986P90=Pred.pcs[,3]

x_train

test_data <- dataP1D1[-inTrain,] #Random
# test_data <- dataP1D1[inTrain$test,] # Kennard Stone
# test_data <- dataP1D1[-indx$index_samples,] # CLHS


y_test <- test_data[,148]
x_test <- test_data[c(2:45,106:147)]
x_test <- scale(x_test, center = min_train, scale = max_train-min_train)
x_test <- data.frame(SOC_1986t=y_test,x_test)


Pred.pcs<-predict(pca1985Med,test_data[,c(46:55)])
x_test$PCA1_1985Med=Pred.pcs[,1] 
x_test$PCA2_1985Med=Pred.pcs[,2]
x_test$PCA3_1985Med=Pred.pcs[,3] 

Pred.pcs<-predict(pca1985P10,test_data[,c(56:65)])
x_test$PCA1_1985P10=Pred.pcs[,1] 
x_test$PCA2_1985P10=Pred.pcs[,2]
x_test$PCA3_1985P10=Pred.pcs[,3]

Pred.pcs<-predict(pca1985P90,test_data[,c(66:75)])
x_test$PCA1_1985P90=Pred.pcs[,1] 
x_test$PCA2_1985P90=Pred.pcs[,2]
x_test$PCA3_1985P90=Pred.pcs[,3]

Pred.pcs<-predict(pca1986Med,test_data[,c(76:85)])
x_test$PCA1_1986Med=Pred.pcs[,1] 
x_test$PCA2_1986Med=Pred.pcs[,2]
x_test$PCA3_1986Med=Pred.pcs[,3]

Pred.pcs<-predict(pca1986P10,test_data[,c(86:95)])
x_test$PCA1_1986P10=Pred.pcs[,1] 
x_test$PCA2_1986P10=Pred.pcs[,2]
x_test$PCA3_1986P10=Pred.pcs[,3]

Pred.pcs<-predict(pca1986P90,test_data[,c(96:105)])
x_test$PCA1_1986P90=Pred.pcs[,1] 
x_test$PCA2_1986P90=Pred.pcs[,2]
x_test$PCA3_1986P90=Pred.pcs[,3]

train_data <- x_train
test_data <- x_test

summary(train_data)
train_data <- train_data[complete.cases(train_data),]
summary(test_data)
test_data <- test_data[complete.cases(test_data),]


# 6) Features selection ---------------------------------------------------

# 6.2) Boruta algorithm ---------------------------------------------------
names(train_data)
train_data <- train_data %>% na.omit %>% data.frame
{
  start <- Sys.time()
  set.seed(1523)
  (bor <- Boruta(x = train_data[,c(2:101)],
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

saveRDS(preds,"Outputs/NamesPreds/P1D1/PredictorsP1D1_08082022.rds")

# 7) Model fitting --------------------------------------------------------

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

grid <- expand.grid(mtry = c(3,4,5,6),
                    splitrule = c("variance", "extratrees"),
                    min.node.size = c(3,4,5,6)
)

set.seed(1759)
model_rf <- train(fm,
                  data=train_data,
                  method = "ranger",
                  trControl = rctrlG,
                  tuneGrid = grid,
                  num.trees = 500,
                  importance = "impurity"
)
model_rf$finalModel$variable.importance
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
pred_svm <- predict(model_svm, newdata = test_data[,preds])
(svm.goof <- goof(observed = exp(test_data$SOC_1986t), predicted = exp(pred_svm)))




# 7.3) Multiple linear regression -----------------------------------------

model_lm <- lm(fm,
               data=train_data)
summary(model_lm)       
model_lm_up <- stepAIC(model_lm,direction="both")

summary(model_lm_up)

pred_lm <- predict(model_lm_up, newdata = test_data[,preds])
(lm.goof <- goof(observed = exp(test_data$SOC_1986t), predicted = exp(pred_lm)))



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

pred_cub <- predict(model_cubist, newdata = test_data[,preds])
(cub.goof <- goof(observed = exp(test_data$SOC_1986t), predicted = exp(pred_cub)))



# 7.5) Gradient boosting regression ---------------------------------------

xgb_trcontrol = trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 10,
  allowParallel = TRUE
)

gbmGrid <-  expand.grid(interaction.depth = 3, 
                        n.trees = 30000, 
                        shrinkage = 0.01,
                        n.minobsinnode = 5)

set.seed(1205) 
gb_model = train(fm,
                  data=train_data,
                  method = "gbm",
                  tuneGrid = gbmGrid,
                  trControl=xgb_trcontrol,
                  verbose=FALSE
)

gb_model
summary(gb_model)

pred_gb <- predict(gb_model, newdata = test_data[,preds])

(gb.goof <- goof(observed = exp(test_data$SOC_1986t), predicted = exp(pred_gb)))


# 7.6) Quantile random forest ---------------------------------------------

set.seed(1554)
model_qrf <- quantregForest(y = train_data[,"SOC_1986t"],
                            x = train_data[,preds],
                            data=train_data,
                            keep.inbag=TRUE,
                            mtry = as.numeric(model_rf$bestTune))
model_qrf$importance
importance(model_qrf,type = 2)


pred_qrf <- predict(model_qrf, newdata = test_data[,preds])
(qrf.goof <- goof(observed = exp(test_data$SOC_1986t), predicted = exp(pred_qrf[,2])))


saveRDS(model_qrf,"Outputs/Models/P1D1/ModelP1D1qrf_08082022.rds")


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
  methodList=c("lm", "qrf")
)

glm_ensemble <- caretStack(
  model_list,
  method="glm",
  metric="RMSE")

pred_ens <- predict(glm_ensemble, newdata = test_data[,preds])
(ens.goof <- goof(observed = test_data$SOC_1986t, predicted = pred_ens))



# 8) Spatial prediction ---------------------------------------------------

fm

covP1 <- c(rast("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/StatPreds.tif"),
           rast("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/DynPredsP1985_1986.tif"))
names(covP1) <- c(readRDS("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/NamesStatPreds.rds"),
                  readRDS("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/NamesDynPredsP1985_1986.rds"))
preds <- readRDS("Outputs/NamesPreds/P1D1/PredictorsP1D1_29072022.rds")

preds %in% names(covP1)
fm

# Not mandatory) Missing covariates - i.e. PCA layers ---------------------

Pred.pcs.layers1 <- predict(covP1[[c(67:76)]],pca1985P90,cores=15)
# 
# 
writeRaster(raster(Pred.pcs.layers1[[2]]),"ExtraCovariates/P1D1/PCA2_1985P90.tif",overwrite=T)

# 8.1) Maps generation ----------------------------------------------------

model_qrf <- readRDS("Outputs/Models/P1D1/ModelP1D1qrf_08082022.rds")
print(model_qrf)

covP1 <- stack(stack("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/StatPreds.tif"),
               stack("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/DynPredsP1985_1986.tif"))
names(covP1) <- c(readRDS("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/NamesStatPreds.rds"),
                  readRDS("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/NamesDynPredsP1985_1986.rds"))

names(covP1)
covP1 <- covP1[[-c(145:148)]]
PCA2_1985P90 <- raster("ExtraCovariates/P1D1/PCA2_1985P90.tif")
covP1 <- stack(covP1,PCA2_1985P90)
names(covP1)

covP1 <- covP1[[preds]]

names(covP1)

beginCluster(n=detectCores()-2,type='SOCK')

covP1sc <- clusterR(covP1[[-16]], scale, 
                    args=list(center = min_train[preds[-16]],
                              scale = max_train[preds[-16]]-
                                min_train[preds[-16]]))

covP1sc <- stack(covP1sc,PCA2_1985P90)
names(covP1sc) <- preds

median <- clusterR(covP1sc, predict,
                   args=list(model=model_qrf, what=0.5))
median <- exp(median)


UppL <- clusterR(covP1sc, predict,
                 args=list(model=model_qrf, what=0.95))
UppL <- exp(UppL)


LowL <- clusterR(covP1sc, predict,
                 args=list(model=model_qrf, what=0.05))
LowL <- exp(LowL)


mean <- clusterR(covP1sc, predict,
                 args=list(model=model_qrf, what=mean))
mean <- exp(mean)

sd <- clusterR(covP1sc, predict,
                 args=list(model=model_qrf, what=sd))
sd <- exp(sd)

writeRaster(median,"Outputs/Layers/P1D1/ModelP1D1QrfMedian_080822.tif",overwrite=T)
writeRaster(UppL,"Outputs/Layers/P1D1/ModelP1D1QrfUppL_080822.tif",overwrite=T)
writeRaster(LowL,"Outputs/Layers/P1D1/ModelP1D1QrfLowL_080822.tif",overwrite=T)
writeRaster(mean,"Outputs/Layers/P1D1/ModelP1D1QrfMean_080822.tif",overwrite=T)
writeRaster(sd,"Outputs/Layers/P1D1/ModelP1D1QrfSd_080822.tif",overwrite=T)

endCluster()


