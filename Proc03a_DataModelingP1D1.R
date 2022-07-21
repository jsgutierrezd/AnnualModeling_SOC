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
          'snow'
)

usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE)
  require(p, character.only = TRUE)
}
lapply(pckg,usePackage)



# 3) Data loading ---------------------------------------------------------

dataP1D1 <- read_delim("RegMat_P1D1.csv",
                       delim = ",") %>% 
  mutate(SOC_1986t = log(SOC_1986))
summary(dataP1D1)
dataP1D1$bsiP10_1985 <- NULL
dataP1D1$strP10_1985 <- NULL
dataP1D1$strP10_1986 <- NULL
dataP1D1$geology <- as.factor(dataP1D1$geology)
dataP1D1$georeg <- as.factor(dataP1D1$georeg)



# 4) Correlation matrix ---------------------------------------------------

corMat <- cor(as.matrix(dataP1D1[,-c(6,7)])) %>% as.data.frame
write_csv(corMat,"CorMatRegMat_P1D1.csv")

names(dataP1D1)
hist(log(dataP1D1$SOC_1986))

# 5) Data splitting -------------------------------------------------------


inTrain <- createDataPartition(y = dataP1D1$SOC_1986t, p = .70, list = FALSE)

train_data <- dataP1D1[ inTrain,]
test_data <- dataP1D1[-inTrain,]

names(train_data)
train_data <- train_data %>% na.omit %>% data.frame
train_data[,-c(1,6,7,188:206)] <- scale(train_data[,-c(1,6,7,188:206)])



# 6) PCA on spectral indices ----------------------------------------------

names(train_data)
names(test_data)

pca1985Med<-prcomp(train_data[,c(53:61)], scale=TRUE) 
summary(pca1985Med)
(corvar <- pca1985Med$rotation %*% diag(pca1985Med$sdev))
Pred.pcs<-predict(pca1985Med,train_data[,c(53:61)])
train_data$PCA1_1985Med=Pred.pcs[,1] 
train_data$PCA2_1985Med=Pred.pcs[,2]
train_data$PCA3_1985Med=Pred.pcs[,3] 
Pred.pcs<-predict(pca1985Med,test_data[,c(53:61)])
test_data$PCA1_1985Med=Pred.pcs[,1] 
test_data$PCA2_1985Med=Pred.pcs[,2]
test_data$PCA3_1985Med=Pred.pcs[,3] 

pca1985P10<-prcomp(train_data[,c(62:68)], scale=TRUE) 
summary(pca1985P10)
(corvar <- pca1985P10$rotation %*% diag(pca1985P10$sdev))
Pred.pcs<-predict(pca1985P10,train_data[,c(62:68)])
train_data$PCA1_1985P10=Pred.pcs[,1] 
train_data$PCA2_1985P10=Pred.pcs[,2]
train_data$PCA3_1985P10=Pred.pcs[,3]
Pred.pcs<-predict(pca1985P10,test_data[,c(62:68)])
test_data$PCA1_1985P10=Pred.pcs[,1] 
test_data$PCA2_1985P10=Pred.pcs[,2]
test_data$PCA3_1985P10=Pred.pcs[,3]

pca1985P90<-prcomp(train_data[,c(69:77)], scale=TRUE) 
summary(pca1985P90)
(corvar <- pca1985P90$rotation %*% diag(pca1985P90$sdev))
Pred.pcs<-predict(pca1985P90,train_data[,c(69:77)])
train_data$PCA1_1985P90=Pred.pcs[,1] 
train_data$PCA2_1985P90=Pred.pcs[,2]
train_data$PCA3_1985P90=Pred.pcs[,3]
Pred.pcs<-predict(pca1985P90,test_data[,c(69:77)])
test_data$PCA1_1985P90=Pred.pcs[,1] 
test_data$PCA2_1985P90=Pred.pcs[,2]
test_data$PCA3_1985P90=Pred.pcs[,3]

pca1986Med<-prcomp(train_data[,c(78:86)], scale=TRUE) 
summary(pca1986Med)
(corvar <- pca1986Med$rotation %*% diag(pca1986Med$sdev))
Pred.pcs<-predict(pca1986Med,train_data[,c(78:86)])
train_data$PCA1_1986Med=Pred.pcs[,1] 
train_data$PCA2_1986Med=Pred.pcs[,2]
train_data$PCA3_1986Med=Pred.pcs[,3]
Pred.pcs<-predict(pca1986Med,test_data[,c(78:86)])
test_data$PCA1_1986Med=Pred.pcs[,1] 
test_data$PCA2_1986Med=Pred.pcs[,2]
test_data$PCA3_1986Med=Pred.pcs[,3]

pca1986P10<-prcomp(train_data[,c(87:94)], scale=TRUE) 
summary(pca1986P10)
(corvar <- pca1986P10$rotation %*% diag(pca1986P10$sdev))
Pred.pcs<-predict(pca1986P10,train_data[,c(87:94)])
train_data$PCA1_1986P10=Pred.pcs[,1] 
train_data$PCA2_1986P10=Pred.pcs[,2]
train_data$PCA3_1986P10=Pred.pcs[,3]
Pred.pcs<-predict(pca1986P10,test_data[,c(87:94)])
test_data$PCA1_1986P10=Pred.pcs[,1] 
test_data$PCA2_1986P10=Pred.pcs[,2]
test_data$PCA3_1986P10=Pred.pcs[,3]

pca1986P90<-prcomp(train_data[,c(95:103)], scale=TRUE) 
summary(pca1986P90)
(corvar <- pca1986P90$rotation %*% diag(pca1986P90$sdev))
Pred.pcs<-predict(pca1986P90,train_data[,c(95:103)])
train_data$PCA1_1986P90=Pred.pcs[,1] 
train_data$PCA2_1986P90=Pred.pcs[,2]
train_data$PCA3_1986P90=Pred.pcs[,3]
Pred.pcs<-predict(pca1986P90,test_data[,c(95:103)])
test_data$PCA1_1986P90=Pred.pcs[,1] 
test_data$PCA2_1986P90=Pred.pcs[,2]
test_data$PCA3_1986P90=Pred.pcs[,3]


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
#   (rfe <- rfe(x=train_data[,c(2:25,150:187,189:206)], y=train_data[,188], sizes=c(1:180), rfeControl=control2))
#   print(Sys.time() - start)
# }
# 
# #names(train_data[,c(2:25,150:187,189:206)])
# plot(rfe, type=c("g", "o"))
# predictors(rfe)

# 6.2) Boruta algorithm ---------------------------------------------------
names(train_data)
train_data <- train_data %>% na.omit %>% data.frame
{start <- Sys.time()
(bor <- Boruta(x = train_data[,c(2:25,150:187,189:206)],
               y = train_data[,188], 
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
print(Sys.time() - start)}

print(bor)
names(bor$finalDecision[bor$finalDecision %in% c("Confirmed")])


# 7) Model fitting --------------------------------------------------------

fm <- as.formula(paste("SOC_1986t ~", paste0(names(bor$finalDecision[bor$finalDecision %in% c("Confirmed")]),
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

pred_rf <- predict(model_rf, newdata = test_data)

rf.goof <- gof(sim = pred_rf,obs = test_data$SOC_1986t)

rf.goof


# 7.2) SVM ----------------------------------------------------------------

tuneResult <- tune(svm, fm, data = train_data,
                   ranges = list(epsilon = seq(0.1,0.3,0.02),
                                 cost = c(5,7,15)))

model_svm <- tuneResult$best.model
print(model_svm)
pred_svm <- predict(model_svm, newdata = test_data)
svm.goof <- gof(sim = pred_svm,obs = test_data$SOC_1986t)
svm.goof



# 7.3) Multiple linear regression -----------------------------------------

model_lm <- lm(fm,
               data=train_data)
summary(model_lm)       
model_lm_up <- stepAIC(model_lm,direction="both")

summary(model_lm_up)

pred_lm <- predict(model_lm_up, newdata = test_data)

lm.goof <- gof(sim = pred_lm,obs = test_data$SOC_1986t)
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

cubist.goof <- gof(sim = pred_cubist,obs = test_data$SOC_1986t)
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

xgb.goof <- gof(sim = pred_xgb,obs = test_data$SOC_1986t)
xgb.goof



# 7.6) Quantile random forest ---------------------------------------------

set.seed(37)
model_qrf <- train(fm,
                  data=train_data,
                  method = "qrf")
model_qrf
varImpPlot(model_qrf[11][[1]])
plot(model_qrf[11][[1]])
model_qrf$finalModel

pred_qrf <- predict(model_qrf, newdata = test_data)

qrf.goof <- gof(sim = pred_qrf,obs = test_data$SOC_1986t)

qrf.goof




# 7.7) Caret model ensemble -----------------------------------------------

data.frame(pred_cubist,pred_lm,pred_qrf,pred_rf,pred_svm,pred_xgb) %>% cor
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

pred_ens <- predict(glm_ensemble, newdata=test_data)
ens.goof <- gof(sim = pred_ens,obs = test_data$SOC_1986t)

ens.goof



# 8) Spatial prediction ---------------------------------------------------

fm
library(raster)
covP1 <- rast("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/CovP1984_1986.tif")
names(covP1) <- readRDS("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/YearbyYear/NamesCovP1984_1986.rds")
names(covP1)

names(bor$finalDecision[bor$finalDecision %in% c("Confirmed")]) %in% names(covP1)

# 8.1) Missign covariates - i.e. principal component analysis layer -------

Pred.pcs.layers <- predict(covP1[[c(52:60)]],pca1985Med)

covP1$PCA1_1985Med <- Pred.pcs.layers[[1]]



georeg <- raster(covP1$georeg) %>% ratify
(catg <- levels(georeg)[[1]])
catg$index <- paste0("georeg", c(1,2,3,4,6,8,10))
(levels(georeg) <- catg)
georeg <- rast(georeg)


covP1$georeg <- georeg


beginCluster(n=detectCores()-2,type='SOCK')

PRED <- clusterR(covP1[[names(bor$finalDecision[bor$finalDecision %in% c("Confirmed")])]], predict,
                 args=list(model=glm_ensemble))

endCluster()
