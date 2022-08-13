#============================================================================
# Proc03d_DataModelingP4D1 ------------------------------------------------
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

dataP4D1 <- read_delim("RegMat_P4D1.csv",
                       delim = ",") %>% 
  mutate(SOC_2019t = log(SOC_2019)) %>% na.omit
summary(dataP4D1)
dataP4D1$strP10_2018 <- ifelse(dataP4D1$strP10_2018=="Inf"|dataP4D1$strP10_2018=="-Inf",NA,dataP4D1$strP10_2018)
dataP4D1$strP10_2019 <- ifelse(dataP4D1$strP10_2019=="Inf"|dataP4D1$strP10_2019=="-Inf",NA,dataP4D1$strP10_2019)
dataP4D1$geology_9 <- NULL
dataP4D1$geology_11 <- NULL
dataP4D1 <- dataP4D1 %>% na.omit
summary(dataP4D1)
names(dataP4D1)



# 4) Correlation matrix ---------------------------------------------------
# 
# corMat <- cor(as.matrix(dataP4D1[,-c(6:23,121:127)])) %>% na.omit %>% as.data.frame
# corMat <- dataP4D1[,-c(6:23,121:127)] %>% na.omit %>% cor %>% data.frame
# 
# write_csv(corMat,"CorMatRegMat_P4D1.csv")
# 
# names(dataP4D1)
# hist(log(dataP4D1$SOC_2019))

# 5) Data splitting -------------------------------------------------------

set.seed(934)
# inTrain <- createDataPartition(y = dataP2D1$SOC_1997t, p = .70, list = FALSE) # Random
inTrain <- kenStone(dataP4D1, k = nrow(dataP4D1)*0.70, metric = "mahal") # Kennard Stone
# data1 <- dataP4D1[,c("SOC_1997t")]
# set.seed(58)
# indx <- clhs(data1, size = round(nrow(data1)*0.7),
#              progress = T, iter = 1000,use.cpp = F,  simple = FALSE) # CLHS


# train_data <- dataP4D1[ inTrain,] %>% data.frame #Random
train_data <- dataP4D1[ inTrain$model,] %>% data.frame #Kennard Stone
# train_data <- dataP4D1[indx$index_samples,] %>% data.frame  # CLHS
names(train_data)
y_train <- train_data[,131]
x_train <- train_data[,c(2:45,106:130)]
max_train <- apply(x_train, 2, max)
min_train <- apply(x_train, 2, min)
x_train <- scale(x_train, center = min_train, scale = max_train-min_train)
x_train <- data.frame(SOC_2019t=y_train,x_train)

# 6) PCA on spectral indices ----------------------------------------------
names(train_data)
pca2018Med<-prcomp(train_data[,c(46:55)], scale=TRUE) 
summary(pca2018Med)
(corvar <- pca2018Med$rotation %*% diag(pca2018Med$sdev))
Pred.pcs<-predict(pca2018Med,train_data[,c(46:55)])
x_train$PCA1_2018Med=Pred.pcs[,1] 
x_train$PCA2_2018Med=Pred.pcs[,2]
x_train$PCA3_2018Med=Pred.pcs[,3] 

pca2018P10<-prcomp(train_data[,c(56:65)], scale=TRUE) 
summary(pca2018P10)
(corvar <- pca2018P10$rotation %*% diag(pca2018P10$sdev))
Pred.pcs<-predict(pca2018P10,train_data[,c(56:65)])
x_train$PCA1_2018P10=Pred.pcs[,1] 
x_train$PCA2_2018P10=Pred.pcs[,2]
x_train$PCA3_2018P10=Pred.pcs[,3]

pca2018P90<-prcomp(train_data[,c(66:75)], scale=TRUE) 
summary(pca2018P90)
(corvar <- pca2018P90$rotation %*% diag(pca2018P90$sdev))
Pred.pcs<-predict(pca2018P90,train_data[,c(66:75)])
x_train$PCA1_2018P90=Pred.pcs[,1] 
x_train$PCA2_2018P90=Pred.pcs[,2]
x_train$PCA3_2018P90=Pred.pcs[,3]

pca2019Med<-prcomp(train_data[,c(76:85)], scale=TRUE) 
summary(pca2019Med)
(corvar <- pca2019Med$rotation %*% diag(pca2019Med$sdev))
Pred.pcs<-predict(pca2019Med,train_data[,c(76:85)])
x_train$PCA1_2019Med=Pred.pcs[,1] 
x_train$PCA2_2019Med=Pred.pcs[,2]
x_train$PCA3_2019Med=Pred.pcs[,3]

pca2019P10<-prcomp(train_data[,c(86:95)], scale=TRUE) 
summary(pca2019P10)
(corvar <- pca2019P10$rotation %*% diag(pca2019P10$sdev))
Pred.pcs<-predict(pca2019P10,train_data[,c(86:95)])
x_train$PCA1_2019P10=Pred.pcs[,1] 
x_train$PCA2_2019P10=Pred.pcs[,2]
x_train$PCA3_2019P10=Pred.pcs[,3]

pca2019P90<-prcomp(train_data[,c(96:105)], scale=TRUE) 
summary(pca2019P90)
(corvar <- pca2019P90$rotation %*% diag(pca2019P90$sdev))
Pred.pcs<-predict(pca2019P90,train_data[,c(96:105)])
x_train$PCA1_2019P90=Pred.pcs[,1] 
x_train$PCA2_2019P90=Pred.pcs[,2]
x_train$PCA3_2019P90=Pred.pcs[,3]

x_train


# test_data <- dataP4D1[-inTrain,] #Random
test_data <- dataP4D1[inTrain$test,] # Kennard Stone
# test_data <- dataP4D1[-indx$index_samples,] # CLHS

y_test <- test_data[,131]
x_test <- test_data[c(2:45,106:130)]
x_test <- scale(x_test, center = min_train, scale = max_train-min_train)
x_test <- data.frame(SOC_2019t=y_test,x_test)


Pred.pcs<-predict(pca2018Med,test_data[,c(46:55)])
x_test$PCA1_2018Med=Pred.pcs[,1] 
x_test$PCA2_2018Med=Pred.pcs[,2]
x_test$PCA3_2018Med=Pred.pcs[,3] 

Pred.pcs<-predict(pca2018P10,test_data[,c(56:65)])
x_test$PCA1_2018P10=Pred.pcs[,1] 
x_test$PCA2_2018P10=Pred.pcs[,2]
x_test$PCA3_2018P10=Pred.pcs[,3]

Pred.pcs<-predict(pca2018P90,test_data[,c(66:75)])
x_test$PCA1_2018P90=Pred.pcs[,1] 
x_test$PCA2_2018P90=Pred.pcs[,2]
x_test$PCA3_2018P90=Pred.pcs[,3]

Pred.pcs<-predict(pca2019Med,test_data[,c(76:85)])
x_test$PCA1_2019Med=Pred.pcs[,1] 
x_test$PCA2_2019Med=Pred.pcs[,2]
x_test$PCA3_2019Med=Pred.pcs[,3]

Pred.pcs<-predict(pca2019P10,test_data[,c(86:95)])
x_test$PCA1_2019P10=Pred.pcs[,1] 
x_test$PCA2_2019P10=Pred.pcs[,2]
x_test$PCA3_2019P10=Pred.pcs[,3]

Pred.pcs<-predict(pca2019P90,test_data[,c(96:105)])
x_test$PCA1_2019P90=Pred.pcs[,1] 
x_test$PCA2_2019P90=Pred.pcs[,2]
x_test$PCA3_2019P90=Pred.pcs[,3]

train_data <- x_train
test_data <- x_test

summary(train_data)
train_data <- train_data[complete.cases(train_data),]
summary(test_data)
test_data <- test_data[complete.cases(test_data),]
# 6) Features selection ---------------------------------------------------


# 6.1) Boruta algorithm ---------------------------------------------------
names(train_data)
train_data <- train_data %>% na.omit %>% data.frame
{
  start <- Sys.time()
  set.seed(1910)
  (bor <- Boruta(x = train_data[,c(2:88)],
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

# saveRDS(preds,"Outputs/NamesPreds/P4D1/PredictorsP4D1_01082022.rds")

# 7) Model fitting --------------------------------------------------------

fm <- as.formula(paste("SOC_2019t ~", paste0(preds,
                                             collapse = "+")))
fm

# 7.1) Randon forest - Ranger ---------------------------------------------

rctrlG <- trainControl(method = "repeatedcv",
                       number = 20,
                       repeats = 20,
                       returnResamp = "all",
                       search = "grid"
)

grid <- expand.grid(mtry = c(5,7,9),
                    splitrule = c("variance", "extratrees"),
                    min.node.size = c(2,4,6)
)

set.seed(949)

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

(svm.goof <- goof(observed = test_data$SOC_2019t, predicted = pred_svm))



# 7.3) Multiple linear regression -----------------------------------------

model_lm <- lm(fm,
               data=train_data)
summary(model_lm)       
model_lm_up <- stepAIC(model_lm,direction="both")

summary(model_lm_up)

pred_lm <- predict(model_lm, newdata = test_data[,names(bor$finalDecision[bor$finalDecision %in% c("Confirmed")])])

(lm.goof <- goof(observed = test_data$SOC_2019t, predicted = pred_lm))



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




# 7.5) Gradient boosting regression ---------------------------------------

xgb_trcontrol = trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 10,
  allowParallel = TRUE
  )

gbmGrid <-  expand.grid(interaction.depth = c(1,2,4), 
                        n.trees = (1:30)*50, 
                        shrinkage = 0.1,
                        n.minobsinnode = 20)



set.seed(1105) 
xgb_model = train(fm,
                  data=train_data,
                  method = "gbm",
                  tuneGrid = gbmGrid,
                  trControl=xgb_trcontrol,
                  verbose=FALSE
)

gb_model
summary(gb_model)

pred_gb <- predict(gb_model, newdata = test_data[,preds])

(gb.goof <- goof(observed = exp(test_data$SOC_2019t), predicted = exp(pred_gb)))



# 7.6) Quantile random forest ---------------------------------------------

set.seed(955)
model_qrf <- quantregForest(y = train_data[,"SOC_2019t"],
                            x = train_data[,preds],
                            data=train_data,
                            keep.inbag=TRUE,
                            mtry = model_rf$bestTune$mtry)
model_qrf
importance(model_qrf,type = 2)
pred_qrf <- predict(model_qrf, newdata = test_data[,preds])

(qrf.goof <- goof(observed = exp(test_data$SOC_2019t), predicted = exp(pred_qrf[,2])))


# saveRDS(model_qrf,"Outputs/Models/P4D1/ModelP4D1qrf_010822.rds")


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

pred_ens <- predict(glm_ensemble, newdata=test_data)
ens.goof <- gof(sim = pred_ens,obs = test_data$SOC_2019t)

ens.goof



# 8) Spatial prediction ---------------------------------------------------
preds <- readRDS("Outputs/NamesPreds/P4D1/PredictorsP4D1_01082022.rds")
model_qrf <- readRDS("Outputs/Models/P4D1/ModelP4D1qrf_010822.rds")
model_qrf$importance
rownames(model_qrf$importance) %in% preds

covP4 <- c(rast("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/StatPreds.tif"),
           rast("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/DynPredsP2018_2019.tif"),
           rast("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/FieldBlock18.tif"),
           rast("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/FieldBlock19.tif"))
names(covP4) <- c(readRDS("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/NamesStatPreds.rds"),
                  readRDS("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/NamesDynPredsP2018_2019.rds"),
                  readRDS("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/NamesFieldBlock18.rds"),
                  readRDS("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/NamesFieldBlock19.rds"))

preds %in% names(covP4)
preds

# Not mandatory) Missing covariates - i.e. PCA layers ---------------------
names(covP4)
Pred.pcs.layers1 <- predict(covP4[[c(77:86)]],pca2019Med,cores=15)
Pred.pcs.layers2 <- predict(covP4[[c(97:106)]],pca2019P90,cores=15)
Pred.pcs.layers3 <- predict(covP4[[c(97:106)]],pca2019P90,cores=15)


writeRaster(raster(Pred.pcs.layers1[[3]]),"ExtraCovariates/P4D1/PCA3_2019Med.tif",overwrite=T)
writeRaster(raster(Pred.pcs.layers2[[1]]),"ExtraCovariates/P4D1/PCA1_2019P90.tif",overwrite=T)
writeRaster(raster(Pred.pcs.layers2[[2]]),"ExtraCovariates/P4D1/PCA2_2019P90.tif",overwrite=T)

# 8.1) Maps generation ----------------------------------------------------

model_qrf <- readRDS("Outputs/Models/P4D1/ModelP4D1qrf_010822.rds")
print(model_qrf)

covP4 <- stack(stack("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/StatPreds.tif"),
               stack("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/DynPredsP2018_2019.tif"))
names(covP4) <- c(readRDS("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/NamesStatPreds.rds"),
                  readRDS("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/NamesDynPredsP2018_2019.rds"))


PCA3_2019Med <- raster("ExtraCovariates/P4D1/PCA3_2019Med.tif")
PCA1_2019P90 <- raster("ExtraCovariates/P4D1/PCA1_2019P90.tif")
PCA2_2019P90 <- raster("ExtraCovariates/P4D1/PCA2_2019P90.tif")

covP4 <- stack(covP4,PCA3_2019Med,PCA1_2019P90,PCA2_2019P90)

covP4 <- covP4[[preds]]

names(covP4)
beginCluster(n=detectCores()-2,type='SOCK')

covP4sc <- clusterR(covP4[[-c(26:28)]], scale, 
                    args=list(center = min_train[preds[-c(26:28)]],
                              scale = max_train[preds[-c(26:28)]]-
                                min_train[preds[-c(26:28)]]))

covP4sc <- stack(covP4sc,PCA3_2019Med,PCA1_2019P90,PCA2_2019P90)
names(covP4sc) <- preds

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

sd <- clusterR(covP4sc, predict,
               args=list(model=model_qrf, what=sd))
sd <- exp(sd)


writeRaster(median,"Outputs/Layers/P4D1/ModelP4D1QrfMedian_010822.tif",overwrite=T)
writeRaster(UppL,"Outputs/Layers/P4D1/ModelP4D1QrfUppL_010822.tif",overwrite=T)
writeRaster(LowL,"Outputs/Layers/P4D1/ModelP4D1QrfLowL_010822.tif",overwrite=T)
writeRaster(mean,"Outputs/Layers/P4D1/ModelP4D1QrfMean_010822.tif",overwrite=T)
writeRaster(sd,"Outputs/Layers/P4D1/ModelP4D1QrfSd_010822.tif",overwrite=T)

endCluster()
