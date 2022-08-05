#============================================================================
# Proc03c_DataModelingP3D1 ------------------------------------------------
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
          'raster'
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

dataP3D1 <- read_delim("RegMat_P3D1.csv",
                       delim = ",") %>% 
  mutate(SOC_2009t = log(SOC_2009)) %>% na.omit
summary(dataP3D1)
dataP3D1$ndviP10_2008 <- ifelse(dataP3D1$ndviP10_2008=="Inf"|dataP3D1$ndviP10_2008=="-Inf",NA,dataP3D1$ndviP10_2008)
dataP3D1$strP10_2008 <- ifelse(dataP3D1$strP10_2008=="Inf"|dataP3D1$strP10_2008=="-Inf",NA,dataP3D1$strP10_2008)
dataP3D1$strP10_2009 <- ifelse(dataP3D1$strP10_2009=="Inf"|dataP3D1$strP10_2009=="-Inf",NA,dataP3D1$strP10_2009)
dataP3D1 <- dataP3D1 %>% na.omit
summary(dataP3D1)
names(dataP3D1)



# 4) Correlation matrix ---------------------------------------------------

# corMat <- cor(as.matrix(dataP3D1[,-c(6:23)])) %>% na.omit %>% as.data.frame
# corMat <- dataP3D1[,-c(6:23)] %>% na.omit %>% cor %>% data.frame
# 
# write_csv(corMat,"CorMatRegMat_P3D1.csv")
# 
# names(dataP3D1)
# hist(log(dataP3D1$SOC_2009))

# 5) Data splitting -------------------------------------------------------

set.seed(1638)
inTrain <- createDataPartition(y = dataP3D1$SOC_2009t, p = .70, list = FALSE) # Random
# inTrain <- kenStone(dataP3D1, k = nrow(dataP3D1)*0.70, metric = "mahal") # Kennard Stone
# data1 <- dataP3D1[,c("SOC_1997t")]
# set.seed(58)
# indx <- clhs(data1, size = round(nrow(data1)*0.7),
#              progress = T, iter = 1000,use.cpp = F,  simple = FALSE) # CLHS


train_data <- dataP3D1[ inTrain,] %>% data.frame #Random
# train_data <- dataP3D1[ inTrain$model,] %>% data.frame #Kennard Stone
# train_data <- dataP3D1[indx$index_samples,] %>% data.frame  # CLHS

y_train <- train_data[,140]
x_train <- train_data[,c(2:41,102:139)]
max_train <- apply(x_train, 2, max)
min_train <- apply(x_train, 2, min)
x_train <- scale(x_train, center = min_train, scale = max_train-min_train)
x_train <- data.frame(SOC_2009t=y_train,x_train)

# 6) PCA on spectral indices ----------------------------------------------
pca2008Med<-prcomp(train_data[,c(42:51)], scale=TRUE) 
summary(pca2008Med)
(corvar <- pca2008Med$rotation %*% diag(pca2008Med$sdev))
Pred.pcs<-predict(pca2008Med,train_data[,c(42:51)])
x_train$PCA1_2008Med=Pred.pcs[,1] 
x_train$PCA2_2008Med=Pred.pcs[,2]
x_train$PCA3_2008Med=Pred.pcs[,3] 

pca2008P10<-prcomp(train_data[,c(52:61)], scale=TRUE) 
summary(pca2008P10)
(corvar <- pca2008P10$rotation %*% diag(pca2008P10$sdev))
Pred.pcs<-predict(pca2008P10,train_data[,c(52:61)])
x_train$PCA1_2008P10=Pred.pcs[,1] 
x_train$PCA2_2008P10=Pred.pcs[,2]
x_train$PCA3_2008P10=Pred.pcs[,3]

pca2008P90<-prcomp(train_data[,c(62:71)], scale=TRUE) 
summary(pca2008P90)
(corvar <- pca2008P90$rotation %*% diag(pca2008P90$sdev))
Pred.pcs<-predict(pca2008P90,train_data[,c(62:71)])
x_train$PCA1_2008P90=Pred.pcs[,1] 
x_train$PCA2_2008P90=Pred.pcs[,2]
x_train$PCA3_2008P90=Pred.pcs[,3]

pca2009Med<-prcomp(train_data[,c(72:81)], scale=TRUE) 
summary(pca2009Med)
(corvar <- pca2009Med$rotation %*% diag(pca2009Med$sdev))
Pred.pcs<-predict(pca2009Med,train_data[,c(72:81)])
x_train$PCA1_2009Med=Pred.pcs[,1] 
x_train$PCA2_2009Med=Pred.pcs[,2]
x_train$PCA3_2009Med=Pred.pcs[,3]

pca2009P10<-prcomp(train_data[,c(82:91)], scale=TRUE) 
summary(pca2009P10)
(corvar <- pca2009P10$rotation %*% diag(pca2009P10$sdev))
Pred.pcs<-predict(pca2009P10,train_data[,c(82:91)])
x_train$PCA1_2009P10=Pred.pcs[,1] 
x_train$PCA2_2009P10=Pred.pcs[,2]
x_train$PCA3_2009P10=Pred.pcs[,3]

pca2009P90<-prcomp(train_data[,c(92:101)], scale=TRUE) 
summary(pca2009P90)
(corvar <- pca2009P90$rotation %*% diag(pca2009P90$sdev))
Pred.pcs<-predict(pca2009P90,train_data[,c(92:101)])
x_train$PCA1_2009P90=Pred.pcs[,1] 
x_train$PCA2_2009P90=Pred.pcs[,2]
x_train$PCA3_2009P90=Pred.pcs[,3]

x_train

test_data <- dataP3D1[-inTrain,] #Random
# test_data <- dataP3D1[inTrain$test,] # Kennard Stone
# test_data <- dataP3D1[-indx$index_samples,] # CLHS

y_test <- test_data[,140]
x_test <- test_data[c(2:41,102:139)]
x_test <- scale(x_test, center = min_train, scale = max_train-min_train)
x_test <- data.frame(SOC_2009t=y_test,x_test)


Pred.pcs<-predict(pca2008Med,test_data[,c(42:51)])
x_test$PCA1_2008Med=Pred.pcs[,1] 
x_test$PCA2_2008Med=Pred.pcs[,2]
x_test$PCA3_2008Med=Pred.pcs[,3] 

Pred.pcs<-predict(pca2008P10,test_data[,c(52:61)])
x_test$PCA1_2008P10=Pred.pcs[,1] 
x_test$PCA2_2008P10=Pred.pcs[,2]
x_test$PCA3_2008P10=Pred.pcs[,3]

Pred.pcs<-predict(pca2008P90,test_data[,c(62:71)])
x_test$PCA1_2008P90=Pred.pcs[,1] 
x_test$PCA2_2008P90=Pred.pcs[,2]
x_test$PCA3_2008P90=Pred.pcs[,3]

Pred.pcs<-predict(pca2009Med,test_data[,c(72:81)])
x_test$PCA1_2009Med=Pred.pcs[,1] 
x_test$PCA2_2009Med=Pred.pcs[,2]
x_test$PCA3_2009Med=Pred.pcs[,3]

Pred.pcs<-predict(pca2009P10,test_data[,c(82:91)])
x_test$PCA1_2009P10=Pred.pcs[,1] 
x_test$PCA2_2009P10=Pred.pcs[,2]
x_test$PCA3_2009P10=Pred.pcs[,3]

Pred.pcs<-predict(pca2009P90,test_data[,c(92:101)])
x_test$PCA1_2009P90=Pred.pcs[,1] 
x_test$PCA2_2009P90=Pred.pcs[,2]
x_test$PCA3_2009P90=Pred.pcs[,3]

x_train$geology_9 <- NULL
x_train$geology_11 <- NULL

x_test$geology_9 <- NULL
x_test$geology_11 <- NULL

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
  set.seed(1901)
  (bor <- Boruta(x = train_data[,c(2:95)],
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

saveRDS(preds,"Outputs/NamesPreds/P3D1/PredictorsP3D1_01082022.rds")

# 7) Model fitting --------------------------------------------------------

fm <- as.formula(paste("SOC_2009t ~", paste0(preds,
                                             collapse = "+")))
fm

# 7.1) Randon forest - Ranger ---------------------------------------------

rctrlG <- trainControl(method = "repeatedcv",
                       number = 20,
                       repeats = 20,
                       returnResamp = "all",
                       search = "grid"
                       )

grid <- expand.grid(mtry = c(5,7,8,10),
                    splitrule = c("variance", "extratrees"),
                    min.node.size = c(3,4,8)
)

set.seed(1650)

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

(rf.goof <- goof(observed = exp(test_data$SOC_2009t), predicted = exp(pred_rf)))


# 7.2) SVM ----------------------------------------------------------------

tuneResult <- tune(svm, fm, data = train_data,
                   ranges = list(epsilon = seq(0.1,0.3,0.02),
                                 cost = c(5,7,15)))

model_svm <- tuneResult$best.model
print(model_svm)
pred_svm <- predict(model_svm, newdata = test_data[,names(bor$finalDecision[bor$finalDecision %in% c("Confirmed")])])

(svm.goof <- goof(observed = test_data$SOC_2009t, predicted = pred_svm))



# 7.3) Multiple linear regression -----------------------------------------

model_lm <- lm(fm,
               data=train_data)
summary(model_lm)       
model_lm_up <- stepAIC(model_lm,direction="both")

summary(model_lm_up)

pred_lm <- predict(model_lm, newdata = test_data[,names(bor$finalDecision[bor$finalDecision %in% c("Confirmed")])])

(lm.goof <- goof(observed = test_data$SOC_2009t, predicted = pred_lm))



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

(cubist.goof <- goof(observed = test_data$SOC_2009t, predicted = pred_cubist))



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

(xgb.goof <- goof(observed = test_data$SOC_2009t, predicted = pred_xgb))



# 7.6) Quantile random forest ---------------------------------------------

set.seed(1702)
model_qrf <- quantregForest(y = train_data[,"SOC_2009t"],
                            x = train_data[,preds],
                            data=train_data,
                            keep.inbag=TRUE,
                            mtry = model_rf$bestTune$mtry)
model_qrf
importance(model_qrf,type = 2)
pred_qrf <- predict(model_qrf, newdata = test_data[,preds])

(qrf.goof <- goof(observed = exp(test_data$SOC_2009t), predicted = exp(pred_qrf[,2])))


saveRDS(model_qrf,"Outputs/Models/P3D1/ModelP3D1qrf_010822.rds")


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
ens.goof <- gof(sim = pred_ens,obs = test_data$SOC_2009t)

ens.goof



# 8) Spatial prediction ---------------------------------------------------

fm

covP3 <- c(rast("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/StatPreds.tif"),
           rast("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/DynPredsP2008_2009.tif"))
names(covP3) <- c(readRDS("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/NamesStatPreds.rds"),
                  readRDS("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/NamesDynPredsP2008_2009.rds"))

preds %in% names(covP3)
fm

# Not mandatory) Missing covariates - i.e. PCA layers ---------------------

Pred.pcs.layers1 <- predict(covP3[[c(91:100)]],pca2009P90,cores=15)
writeRaster(raster(Pred.pcs.layers1[[1]]),"ExtraCovariates/P3D1/PCA1_2009P90.tif",overwrite=T)



# 8.1) Maps generation ----------------------------------------------------

model_qrf <- readRDS("Outputs/Models/P3D1/ModelP3D1qrf_010822.rds")
print(model_qrf)

covP3 <- stack(stack("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/StatPreds.tif"),
               stack("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/DynPredsP2008_2009.tif"))
names(covP3) <- c(readRDS("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/NamesStatPreds.rds"),
                  readRDS("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/NamesDynPredsP2008_2009.rds"))


PCA1_2009P90 <- raster("ExtraCovariates/P3D1/PCA1_2009P90.tif")

covP3 <- stack(covP3,PCA1_2009P90)
fm
names(covP3[[-33]])

covP3 <- covP3[[preds]]

beginCluster(n=detectCores()-2,type='SOCK')

covP3sc <- clusterR(covP3[[-33]], scale, 
                    args=list(center = min_train[preds[-33]],
                              scale = max_train[preds[-33]]-
                                min_train[preds[-33]]))

covP3sc <- stack(covP3sc,PCA1_2009P90)
names(covP3sc) <- preds

median <- clusterR(covP3sc, predict,
                   args=list(model=model_qrf, what=0.5))
median <- exp(median)


UppL <- clusterR(covP3sc, predict,
                 args=list(model=model_qrf, what=0.95))
UppL <- exp(UppL)


LowL <- clusterR(covP3sc, predict,
                 args=list(model=model_qrf, what=0.05))
LowL <- exp(LowL)


mean <- clusterR(covP3sc, predict,
                 args=list(model=model_qrf, what=mean))
mean <- exp(mean)

sd <- clusterR(covP3sc, predict,
               args=list(model=model_qrf, what=sd))
sd <- exp(sd)

writeRaster(median,"Outputs/Layers/P3D1/ModelP3D1QrfMedian_010822.tif",overwrite=T)
writeRaster(UppL,"Outputs/Layers/P3D1/ModelP3D1QrfUppL_010822.tif",overwrite=T)
writeRaster(LowL,"Outputs/Layers/P3D1/ModelP3D1QrfLowL_010822.tif",overwrite=T)
writeRaster(mean,"Outputs/Layers/P3D1/ModelP3D1QrfMean_010822.tif",overwrite=T)
writeRaster(sd,"Outputs/Layers/P3D1/ModelP3D1QrfSd_010822.tif",overwrite=T)

endCluster()