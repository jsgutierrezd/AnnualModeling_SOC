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
          'raster'
)

usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE)
  require(p, character.only = TRUE)
}
lapply(pckg,usePackage)



# 3) Data loading ---------------------------------------------------------

dataP2D1 <- read_delim("RegMat_P2D1.csv",
                       delim = ",") %>% 
  mutate(SOC_1997t = log(SOC_1997)) %>% na.omit
summary(dataP2D1)
# dataP2D1$bsiP10_1985 <- ifelse(dataP2D1$bsiP10_1985=="Inf"|dataP2D1$bsiP10_1985=="-Inf",NA,dataP2D1$bsiP10_1985)
# dataP2D1$strP10_1985 <- ifelse(dataP2D1$strP10_1985=="Inf"|dataP2D1$strP10_1985=="-Inf",NA,dataP2D1$strP10_1985)
# dataP2D1$strP10_1986 <- ifelse(dataP2D1$strP10_1986=="Inf"|dataP2D1$strP10_1986=="-Inf",NA,dataP2D1$strP10_1986)
dataP2D1 <- dataP2D1 %>% na.omit
summary(dataP2D1)
names(dataP2D1)



# 4) Correlation matrix ---------------------------------------------------

corMat <- cor(as.matrix(dataP2D1[,-c(6:23)])) %>% na.omit %>% as.data.frame
corMat <- dataP2D1[,-c(6:23)] %>% na.omit %>% cor %>% data.frame

write_csv(corMat,"CorMatRegMat_P2D1.csv")

names(dataP2D1)
hist(log(dataP2D1$SOC_1997))

# 5) Data splitting -------------------------------------------------------

set.seed(1843)
inTrain <- createDataPartition(y = dataP2D1$SOC_1997t, p = .70, list = FALSE)

train_data <- dataP2D1[ inTrain,]
test_data <- dataP2D1[-inTrain,]

names(train_data)
train_data <- train_data %>% na.omit %>% data.frame
# train_data[,-c(1,6,7,188:206)] <- scale(train_data[,-c(1,6,7,188:206)])



# 6) PCA on spectral indices ----------------------------------------------

names(train_data)
names(test_data)

pca1996Med<-prcomp(train_data[,c(42:51)], scale=TRUE) 
summary(pca1996Med)
(corvar <- pca1996Med$rotation %*% diag(pca1996Med$sdev))
Pred.pcs<-predict(pca1996Med,train_data[,c(42:51)])
train_data$PCA1_1996Med=Pred.pcs[,1] 
train_data$PCA2_1996Med=Pred.pcs[,2]
train_data$PCA3_1996Med=Pred.pcs[,3] 
Pred.pcs<-predict(pca1996Med,test_data[,c(42:51)])
test_data$PCA1_1996Med=Pred.pcs[,1] 
test_data$PCA2_1996Med=Pred.pcs[,2]
test_data$PCA3_1996Med=Pred.pcs[,3] 

pca1996P10<-prcomp(train_data[,c(52:61)], scale=TRUE) 
summary(pca1996P10)
(corvar <- pca1996P10$rotation %*% diag(pca1996P10$sdev))
Pred.pcs<-predict(pca1996P10,train_data[,c(52:61)])
train_data$PCA1_1996P10=Pred.pcs[,1] 
train_data$PCA2_1996P10=Pred.pcs[,2]
train_data$PCA3_1996P10=Pred.pcs[,3]
Pred.pcs<-predict(pca1996P10,test_data[,c(52:61)])
test_data$PCA1_1996P10=Pred.pcs[,1] 
test_data$PCA2_1996P10=Pred.pcs[,2]
test_data$PCA3_1996P10=Pred.pcs[,3]

pca1996P90<-prcomp(train_data[,c(62:71)], scale=TRUE) 
summary(pca1996P90)
(corvar <- pca1996P90$rotation %*% diag(pca1996P90$sdev))
Pred.pcs<-predict(pca1996P90,train_data[,c(62:71)])
train_data$PCA1_1996P90=Pred.pcs[,1] 
train_data$PCA2_1996P90=Pred.pcs[,2]
train_data$PCA3_1996P90=Pred.pcs[,3]
Pred.pcs<-predict(pca1996P90,test_data[,c(62:71)])
test_data$PCA1_1996P90=Pred.pcs[,1] 
test_data$PCA2_1996P90=Pred.pcs[,2]
test_data$PCA3_1996P90=Pred.pcs[,3]

pca1997Med<-prcomp(train_data[,c(72:81)], scale=TRUE) 
summary(pca1997Med)
(corvar <- pca1997Med$rotation %*% diag(pca1997Med$sdev))
Pred.pcs<-predict(pca1997Med,train_data[,c(72:81)])
train_data$PCA1_1997Med=Pred.pcs[,1] 
train_data$PCA2_1997Med=Pred.pcs[,2]
train_data$PCA3_1997Med=Pred.pcs[,3]
Pred.pcs<-predict(pca1997Med,test_data[,c(72:81)])
test_data$PCA1_1997Med=Pred.pcs[,1] 
test_data$PCA2_1997Med=Pred.pcs[,2]
test_data$PCA3_1997Med=Pred.pcs[,3]

pca1997P10<-prcomp(train_data[,c(82:91)], scale=TRUE) 
summary(pca1997P10)
(corvar <- pca1997P10$rotation %*% diag(pca1997P10$sdev))
Pred.pcs<-predict(pca1997P10,train_data[,c(82:91)])
train_data$PCA1_1997P10=Pred.pcs[,1] 
train_data$PCA2_1997P10=Pred.pcs[,2]
train_data$PCA3_1997P10=Pred.pcs[,3]
Pred.pcs<-predict(pca1997P10,test_data[,c(82:91)])
test_data$PCA1_1997P10=Pred.pcs[,1] 
test_data$PCA2_1997P10=Pred.pcs[,2]
test_data$PCA3_1997P10=Pred.pcs[,3]

pca1997P90<-prcomp(train_data[,c(92:101)], scale=TRUE) 
summary(pca1997P90)
(corvar <- pca1997P90$rotation %*% diag(pca1997P90$sdev))
Pred.pcs<-predict(pca1997P90,train_data[,c(92:101)])
train_data$PCA1_1997P90=Pred.pcs[,1] 
train_data$PCA2_1997P90=Pred.pcs[,2]
train_data$PCA3_1997P90=Pred.pcs[,3]
Pred.pcs<-predict(pca1997P90,test_data[,c(92:101)])
test_data$PCA1_1997P90=Pred.pcs[,1] 
test_data$PCA2_1997P90=Pred.pcs[,2]
test_data$PCA3_1997P90=Pred.pcs[,3]


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
  (bor <- Boruta(x = train_data[,c(2:41,102:139,141:158)],
                 y = train_data[,140], 
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

# 7) Model fitting --------------------------------------------------------

fm <- as.formula(paste("SOC_1997t ~", paste0(names(bor$finalDecision[bor$finalDecision %in% c("Confirmed")]),
                                             collapse = "+")))
fm

# 7.1) Randon forest - Ranger ---------------------------------------------

rctrlG <- trainControl(method = "repeatedcv",
                       number = 20,
                       repeats = 20,
                       returnResamp = "all",
                       search = "grid"
)

grid <- expand.grid(mtry = c(2,3,4),
                    splitrule = c("variance", "extratrees"),
                    min.node.size = c(2,3,4)
)

set.seed(37)
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

pred_rf <- predict(model_rf, newdata = test_data)

rf.goof <- gof(sim = pred_rf,obs = test_data$SOC_1997t)

rf.goof


# 7.2) SVM ----------------------------------------------------------------

tuneResult <- tune(svm, fm, data = train_data,
                   ranges = list(epsilon = seq(0.1,0.3,0.02),
                                 cost = c(5,7,15)))

model_svm <- tuneResult$best.model
print(model_svm)
pred_svm <- predict(model_svm, newdata = test_data)
svm.goof <- gof(sim = pred_svm,obs = test_data$SOC_1997t)
svm.goof



# 7.3) Multiple linear regression -----------------------------------------

model_lm <- lm(fm,
               data=train_data)
summary(model_lm)       
model_lm_up <- stepAIC(model_lm,direction="both")

summary(model_lm_up)

pred_lm <- predict(model_lm_up, newdata = test_data)

lm.goof <- gof(sim = pred_lm,obs = test_data$SOC_1997t)
lm.goof



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

pred_cubist <- predict(model_cubist, newdata = test_data)

cubist.goof <- gof(sim = pred_cubist,obs = test_data$SOC_1997t)
cubist.goof



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

pred_xgb <- predict(xgb_model, newdata = test_data)

xgb.goof <- gof(sim = pred_xgb,obs = test_data$SOC_1997t)
xgb.goof



# 7.6) Quantile random forest ---------------------------------------------

set.seed(37)
model_qrf <- quantregForest(y = train_data[,"SOC_1997t"],
                            x = train_data[,names(bor$finalDecision[bor$finalDecision %in% c("Confirmed")])],
                            data=train_data,
                            keep.inbag=TRUE,
                            mtry = as.numeric(model_rf$bestTune))
model_qrf
importance(model_qrf,type = 2)
pred_qrf <- predict(model_qrf, newdata = test_data)

qrf.goof <- gof(sim = pred_qrf[,2],obs = test_data$SOC_1997t)

qrf.goof
  
saveRDS(model_qrf,"Outputs/Models/ModelP2D1qrf_260722.rds")


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

names(bor$finalDecision[bor$finalDecision %in% c("Confirmed")]) %in% names(covP2)
fm

# Not mandatory) Missing covariates - i.e. PCA layers ---------------------

Pred.pcs.layers1 <- predict(covP2[[c(41:50)]],pca1996Med,cores=15)
Pred.pcs.layers2 <- predict(covP2[[c(71:80)]],pca1997Med,cores=15)
Pred.pcs.layers3 <- predict(covP2[[c(91:100)]],pca1997P90,cores=15)
 
 
writeRaster(raster(Pred.pcs.layers1[[1]]),"ExtraCovariates/P2D1/PCA1_1996Med.tif",overwrite=T)
writeRaster(raster(Pred.pcs.layers2[[2]]),"ExtraCovariates/P2D1/PCA2_1997Med.tif",overwrite=T)
writeRaster(raster(Pred.pcs.layers3[[1]]),"ExtraCovariates/P2D1/PCA1_1997P90.tif",overwrite=T)


# 8.1) Maps generation ----------------------------------------------------

model_qrf <- readRDS("Outputs/Models/ModelP2D1qrf_260722.rds")
print(model_qrf)

covP2 <- stack(stack("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/StatPreds.tif"),
           stack("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/DynPredsP1996_1997.tif"))
names(covP2) <- c(readRDS("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/NamesStatPreds.rds"),
                  readRDS("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/NamesDynPredsP1996_1997.rds"))

fm
PCA1_1996Med <- raster("ExtraCovariates/P2D1/PCA1_1996Med.tif")
PCA2_1997Med <- raster("ExtraCovariates/P2D1/PCA2_1997Med.tif") 
PCA1_1997P90 <- raster("ExtraCovariates/P2D1/PCA1_1997P90.tif")

covP2 <- stack(covP2,PCA1_1996Med,PCA2_1997Med,PCA1_1997P90)
names(covP2)

beginCluster(n=detectCores()-2,type='SOCK')

unc <- clusterR(covP2[[names(bor$finalDecision[bor$finalDecision %in% c("Confirmed")])]], predict, 
                args=list(model=model_qrf,what=sd))
unc <- exp(unc)

mean <- clusterR(covP2[[names(bor$finalDecision[bor$finalDecision %in% c("Confirmed")])]], predict,
                 args=list(model=model_qrf, what=mean))
mean <- exp(mean)


plot(mean)
plot(unc)
writeRaster(mean,"ModelP2D1QrfMean_260722.tif")
writeRaster(unc,"ModelP2D1QrfUnc_260722.tif")

endCluster()
