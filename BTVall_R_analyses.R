------------------------------------------------------------------
  ################# Code to generate BTV niche models################
#
#------------------------------------------------------------------
 #Written by Moh Alkhamis and Nick FOuntain-Jones (nfj@umn.edu)

rm(list = ls())

#load packages
library("randomForest")
library("caret")
library("pROC")
library("ROCR")
library("dplyr")
library("missForest")
library("tidyverse")
library("gbm")
library("pdp")
library("ggplot2")
library("iml")
library(maps) 
library(raster)
library(rgdal)
library(rminer)
library(kernlab)
library(rJava)
library(spatstat)
library(maptools)
library(mapdata)
library(doParallel)
library(xgboost)
# load functions
source("functions.R")
source("functionsCV.R")


data(wrld_simpl)

#read BTV outbreaks locations##

btvall<-read.csv("BTVall.csv", head=TRUE,sep=",")
head(btvall)
dim(btvall)

geno1<-read.csv("Geno1.csv", head=TRUE,sep=",")
head(geno1)
dim(geno1)

geno4<-read.csv("Geno4.csv", head=TRUE,sep=",")
head(geno4)
dim(geno4)

geno8<-read.csv("Geno8.csv", head=TRUE,sep=",")
head(geno8)
dim(geno8)

#dont do genotype 16


#map 
#EU <- shapefile("EU.shp")
#### projection(EU)<-CRS('+proj=longlat +datum=WGS84')
#plot(EU, xlim=c(-11.71815,35.47866),ylim=c(34.7099,58.51793),border='black', col='light grey', lwd=0.2)
#points(geno1$long, geno1$lat, col='red', pch=4, cex=0.5)
#points(geno4$long, geno4$lat, col='blue', pch=5, cex=0.5)
#points(geno8$long, geno8$lat, col='dark green', pch=6, cex=0.5)
#### text(coordinates(EU),labels=EU$FAO,cex=0.6)
#SpatialPolygonsRescale(layout.north.arrow(),offset=locator(1),scale=3, plot.grid=FALSE)
#SpatialPolygonsRescale(layout.scale.bar(),offset=locator(1),scale=7, fill=c("transparent","black"),plot.grid=FALSE)
#text(locator(1),"0", cex=1)
#text(locator(1),"800 km", cex=1)
#legend(locator(1),c("Serotype 1","Serotype 4","Serotype 8"),pch=c(4,5,6), col=c('red','blue','dark green'), cex=1,bg='white')
#map.axes()

### Selected Variables
### Reading Bioclim rasters

bio1<-raster("Bioclime/wc2.0_bio_2.5m_01.tif") ##BIO1 = Annual Mean Temperature
bio2<-raster("Bioclime/wc2.0_bio_2.5m_02.tif") ##BIO2 = Mean Diurnal Range (Mean of monthly (max temp - min temp))
bio3<-raster("Bioclime/wc2.0_bio_2.5m_03.tif") ##BIO3 = Isothermality (BIO2/BIO7) (* 100)
bio4<-raster("Bioclime/wc2.0_bio_2.5m_04.tif") ##BIO4 = Temperature Seasonality (standard deviation *100)
bio5<-raster("Bioclime/wc2.0_bio_2.5m_05.tif") ##BIO5 = Max Temperature of Warmest Month
bio6<-raster("Bioclime/wc2.0_bio_2.5m_06.tif") ##BIO6 = Min Temperature of Coldest Month
bio7<-raster("Bioclime/wc2.0_bio_2.5m_07.tif") ##BIO7 = Temperature Annual Range (BIO5-BIO6)
bio10<-raster("Bioclime/wc2.0_bio_2.5m_10.tif") ##BIO10 = Mean Temperature of Warmest Quarter
bio11<-raster("Bioclime/wc2.0_bio_2.5m_11.tif") ##BIO11 = Mean Temperature of Coldest Quarter
bio12<-raster("Bioclime/wc2.0_bio_2.5m_12.tif") ##BIO12 = Annual Precipitation
bio13<-raster("Bioclime/wc2.0_bio_2.5m_13.tif") ##BIO13 = Precipitation of Wettest Month
bio14<-raster("Bioclime/wc2.0_bio_2.5m_14.tif") ##BIO14 = Precipitation of Driest Month
bio15<-raster("Bioclime/wc2.0_bio_2.5m_15.tif") ##BIO15 = Precipitation Seasonality (Coefficient of Variation)
bio16<-raster("Bioclime/wc2.0_bio_2.5m_16.tif") ##BIO16 = Precipitation of Wettest Quarter
bio17<-raster("Bioclime/wc2.0_bio_2.5m_17.tif") ##BIO17 = Precipitation of Driest Quarter

### Reading Other rasters
Buffalo<-raster("Buffalo.density/5_Bf_2010_Da.tif")
Cattle<-raster("Cattle.density/5_Ct_2010_Da.tif")
Goat<-raster("Goat.density/5_Gt_2010_Da.tif")
Sheep<-raster("Sheep.density/5_Sh_2010_Da.tif")
#missing this one
NDVI<-raster("NDVI_Europe/NDVI_EU2018.tif")
land<-raster("Landcover/newifpri30/w001001.adf")
wind<-raster("Wind.Speed/meanwind.tif")
CS<-raster("C_density.tif")
#glps<-raster<-("/Users/malkhamis/Desktop/OneDrive2/BTV Niche/glps.tif")
#z <- raster(glps)
#writeRaster(z, filename="glps.tif", overwrite=TRUE)
#glps1<-aggregate(glps, fact=4, fun=max) 
#glps<-projectRaster(glps, param, method='ngb')


bio1<-crop(bio1, extent(-13.60022,41.14927,34.14933,61.21744))
plot(bio1)
bio2<-crop(bio2, extent(-13.60022,41.14927,34.14933,61.21744))
bio3<-crop(bio3, extent(-13.60022,41.14927,34.14933,61.21744))
bio4<-crop(bio4, extent(-13.60022,41.14927,34.14933,61.21744))
bio5<-crop(bio5, extent(-13.60022,41.14927,34.14933,61.21744))
bio6<-crop(bio6, extent(-13.60022,41.14927,34.14933,61.21744))
bio7<-crop(bio7, extent(-13.60022,41.14927,34.14933,61.21744))
bio10<-crop(bio10, extent(-13.60022,41.14927,34.14933,61.21744))
bio11<-crop(bio11, extent(-13.60022,41.14927,34.14933,61.21744))
bio12<-crop(bio12, extent(-13.60022,41.14927,34.14933,61.21744))
bio13<-crop(bio13, extent(-13.60022,41.14927,34.14933,61.21744))
bio14<-crop(bio14, extent(-13.60022,41.14927,34.14933,61.21744))
bio15<-crop(bio15, extent(-13.60022,41.14927,34.14933,61.21744))
bio16<-crop(bio16, extent(-13.60022,41.14927,34.14933,61.21744))
bio17<-crop(bio17, extent(-13.60022,41.14927,34.14933,61.21744))
Buffalo<-crop(Buffalo,extent(-13.60022,41.14927,34.14933,61.21744))
Cattle<-crop(Cattle,extent(-13.60022,41.14927,34.14933,61.21744))
Goat<-crop(Goat,extent(-13.60022,41.14927,34.14933,61.21744))
Sheep<-crop(Sheep,extent(-13.60022,41.14927,34.14933,61.21744))
NDVI<-crop(NDVI,extent(-13.60022,41.14927,34.14933,61.21744))
land<-crop(land,extent(-13.60022,41.14927,34.14933,61.21744))
wind<-crop(wind,extent(-13.60022,41.14927,34.14933,61.21744))
CS<-crop(CS,extent(-13.60022,41.14927,34.14933,61.21744))
#glps<-crop(glps,extent(-13.60022,41.14927,34.14933,61.21744))


param<-raster(nrow=720, ncol=1200, xmn=-13.60022, xmx=41.14927, ymn=34.14933, ymx=61.21744)
bio1a<-aggregate(bio1, fact=4, fun=max) 
bio1b<-resample(bio1a, param, method='bilinear')
bio2a<-aggregate(bio2, fact=4, fun=max) 
bio2b<-resample(bio2a, param, method='bilinear')
bio3a<-aggregate(bio3, fact=4, fun=max) 
bio3b<-resample(bio3a, param, method='bilinear')
bio4a<-aggregate(bio4, fact=4, fun=max) 
bio4b<-resample(bio4a, param, method='bilinear')
bio5a<-aggregate(bio5, fact=4, fun=max) 
bio5b<-resample(bio5a, param, method='bilinear')
bio6a<-aggregate(bio6, fact=4, fun=max) 
bio6b<-resample(bio6a, param, method='bilinear')
bio7a<-aggregate(bio7, fact=4, fun=max) 
bio7b<-resample(bio7a, param, method='bilinear')
bio10a<-aggregate(bio10, fact=4, fun=max) 
bio10b<-resample(bio10a, param, method='bilinear')
bio11a<-aggregate(bio11, fact=4, fun=max) 
bio11b<-resample(bio11a, param, method='bilinear')
bio12a<-aggregate(bio12, fact=4, fun=max) 
bio12b<-resample(bio12a, param, method='bilinear')
bio13a<-aggregate(bio13, fact=4, fun=max) 
bio13b<-resample(bio13a, param, method='bilinear')
bio14a<-aggregate(bio14, fact=4, fun=max) 
bio14b<-resample(bio14a, param, method='bilinear')
bio15a<-aggregate(bio15, fact=4, fun=max) 
bio15b<-resample(bio15a, param, method='bilinear')
bio16a<-aggregate(bio16, fact=4, fun=max) 
bio16b<-resample(bio16a, param, method='bilinear')
bio17a<-aggregate(bio17, fact=4, fun=max) 
bio17b<-resample(bio17a, param, method='bilinear')
Buffaloa<-aggregate(Buffalo, fact=4, fun=max) 
Buffalob<-resample(Buffaloa, param, method='bilinear')
Cattlea<-aggregate(Cattle, fact=4, fun=max) 
Cattleb<-resample(Cattlea, param, method='bilinear')
Goata<-aggregate(Goat, fact=4, fun=max) 
Goatb<-resample(Goata, param, method='bilinear')
Sheepa<-aggregate(Sheep, fact=4, fun=max) 
Sheepb<-resample(Sheepa, param, method='bilinear')
NDVIa<-aggregate(NDVI, fact=4, fun=max) 
NDVIb<-resample(NDVIa, param, method='bilinear')
landa<-aggregate(land, fact=4, fun=max) 
landb<-resample(landa, param, method='ngb')
winda<-aggregate(wind, fact=4, fun=max) 
windb<-resample(winda, param, method='bilinear')
CSa<-aggregate(CS, fact=4, fun=max) 
CSb<-resample(CSa, param, method='bilinear')


# first layer of the RasterStack
predictors<-stack(bio1b,bio2b,bio3b,bio4b,bio5b,bio6b,bio7b,bio10b,bio11b,bio12b,bio13b,bio14b,bio15b,bio16b,bio17b,Buffalob,Cattleb,Goatb,Sheepb,NDVIb,landb,windb,CSb)
names(predictors)<-c("bio1","bio2","bio3","bio4","bio5","bio6","bio7","bio10","bio11","bio12","bio13","bio14","bio15","bio16","bio17","Buffalo","Cattle","Goat","Sheep","NDVI","land","wind","CS")

#clear memory
gc()
#################All BTV#################
btvp <- raster::extract(predictors, btvall)
head(btvp)
dim(btvp)
i <- which(is.na(btvp[,1])) 
i

plot(btvall, cex=0.5, col='blue')
plot(EU, add=TRUE)
points(btvall[i, ], pch=20, cex=3, col='red')
btvp1 <- btvp[-i, ]
dim(btvp1)



e <- extent(-11.71815,35.47866,34.7099,58.51793)
e
set.seed(0)

#null data
bg <- sampleRandom(predictors, 10000, ext=e)
dim(bg)

#join everything together
d <- rbind(cbind(pa=1, btvp1), cbind(pa=0, bg)) 
d <- data.frame(d)
d[,'land'] <- as.factor(d[,'land'])
dim(d)
data <- d; str(data)
colnames(data)[1] <- "Class"


#################Geno 1#################
btvGeno1 <- raster::extract(predictors, geno1)
str(geno1)
head(btvGeno1)
dim(btvGeno1)
i <- which(is.na(btvGeno1[,1])) 
i
#null data
bgGeno1 <- sampleRandom(predictors, 912, ext=e) #same extent to all data
dim(bgGeno1)

#join everything together
d1 <- rbind(cbind(pa=1, btvGeno1), cbind(pa=0, bgGeno1)) 
d1 <- data.frame(d1)
d1[,'land'] <- as.factor(d1[,'land'])
dim(d1)

data1 <- d1
colnames(data1)[1] <- "Class"

#################Geno 4#################
btvGeno4 <- raster::extract(predictors, geno4)

i <- which(is.na(btvGeno4[,1])) 
i
#null data
bgGeno4 <- sampleRandom(predictors, 6166, ext=e) #same extent to all data
dim(bgGeno4)

#join everything together
d2 <- rbind(cbind(pa=1, btvGeno4), cbind(pa=0, bgGeno4)) 
d2 <- data.frame(d2)
d2[,'land'] <- as.factor(d2[,'land'])
dim(d2)

data2 <- d2
colnames(data2)[1] <- "Class"

#################Geno 8#################
btvGeno8 <- raster::extract(predictors, geno8)

i <- which(is.na(btvGeno8[,1])) 
i
#null data
bgGeno8 <- sampleRandom(predictors, 1887, ext=e) #same extent to all data

#join everything together
d3 <- rbind(cbind(pa=1, btvGeno8), cbind(pa=0, bgGeno8)) 
d3 <- data.frame(d3)
d3[,'land'] <- as.factor(d3[,'land'])
dim(d3)

data3 <- d3
colnames(data3)[1] <- "Class"

------------------------------------------------------------------
  ################# A. Drop NAs#################
#------------------------------------------------------------------
# 
data <-data3[complete.cases(data3), ];str(data) #dropping NAs - only 18 outbreaks with NAs so didn't seem worth it to impute

#------------------------------------------------------------------
########Check correlation structure#################
#------------------------------------------------------------------

CorData <-  data[-which(names(data) == "Class")] #all predictors
CorDataNum <- select_if(CorData, is.numeric)
# # calculate correlation matrix
correlationMatrix <- cor(CorDataNum)

# # summarize the correlation matrix and plot
#library(ggcorrplot)
#ggcorrplot(correlationMatrix, hc.order = TRUE, type = "lower",
           #lab = TRUE)

#remove correlated variabeles
hc <-  findCorrelation(correlationMatrix , cutoff=0.9) # put any value as a "cutoff" 
hc <- sort(hc); hc
reduced_Data <- data[,-c(hc)]
dataR<- cbind(data[1], reduced_Data);str(data)
str(dataR)

#------------------------------------------------------------------
########Select only relevent features#################
#------------------------------------------------------------------
library("Boruta")

X <- dataR[-which(names(dataR) == "Class")] 
Y <- dataR$Class

set.seed(124)
ImpVar <- Boruta(X, Y, doTrace = 0, maxRuns = 500) #, getImp = getImpRfGini seems to be a bit more severe compared to rfZ
print(ImpVar) #looks like no predictor can be excluded.
#plot with x axis label vertical

#plot(ImpVar, xlab = "", xaxt = "n")
#lz<-lapply(1:ncol(ImpVar$ImpHistory),function(i)
  #ImpVar$ImpHistory[is.finite(ImpVar$ImpHistory[,i]),i])
#names(lz) <- colnames(ImpVar$ImpHistory)
#Labels <- sort(sapply(lz,median))
#axis(side = 1,las=2,labels = names(Labels),
     #at = 1:ncol(ImpVar$ImpHistory), cex.axis = 0.7)

#plotImpHistory(ImpVar) #optional - look at the runs in more detail

#set <- getSelectedAttributes(ImpVar, withTentative = TRUE)#keep tentative
#dataReduced <- data[which(names(data) %in% set)]
#data<- cbind(data[1], dataReduced);str(data)


-----------------------------------------------------
  ########create machine learning objects#################
#------------------------------------------------------------------

#have to have 1's as 'Positive and 0's as negatives
as.character(dataR$Class)
dataR$Class[dataR$Class== '0']='Negative'
dataR$Class[dataR$Class== '1']= 'Positive' 
str(dataR)
#for XGRboost we need all numeric
dataR$land <- as.numeric(dataR$land)
# create two partitions: training and testing data
inTrain <-  createDataPartition(y=dataR$Class, p = .8, list = FALSE)

data.train.descr <- dataR[inTrain,-1]
data.train.class <- dataR[inTrain,1]
data.test.descr <- dataR[-inTrain,-1]
data.test.class <- dataR[-inTrain,1]


# ## create folds for CV
set.seed(123)
myFolds <- createMultiFolds(y=data.train.class,k=10,times=10)


# for parallel model training
registerDoParallel(cores = 3)

#------------------------------------------------------------------
#################Set up cross validation#################
#------------------------------------------------------------------

#all data - no downsampling

myControl <- trainControl(## 10-fold CV)
  method = "repeatedcv",
  number = 10,
  repeats = 10,
  index = myFolds,
  savePredictions=TRUE,        
  classProbs = TRUE,
  summaryFunction = twoClassSummary,  
  allowParallel=TRUE,
  selectionFunction = "best")


#------------------------------------------------------------------
#################Random Forests #################
#------------------------------------------------------------------


# tune mtry parameters (number of predictors) and perform cross-validation on data
rf.grid <- expand.grid(.mtry=(3:6)) #number of predictors (mtry) to test

#Run model
set.seed(125)
rf.fit <- train(data.train.descr,data.train.class,
                method = "rf",
                metric = "ROC",
                verbose = FALSE,
                trControl = myControl,
                tuneGrid = rf.grid,
                verboseIter=TRUE
)

plot(rf.fit)

# estimate model performance in terms of a confusion matrix from repeated cross-validation
oof <- rf.fit$pred #01 fold cross validation -check if worked
head(oof)
tail(oof)#mtry = covariates used
# consider only the best tuned model - this checks performance
oof <- oof[oof$mtry==rf.fit$bestTune[,'mtry'],]
repeats<-rf.fit$control$repeats
rf.performance.cv <- EstimatePerformanceCV(oof=oof,repeats=repeats) #Matthew's correlation coefficient (MCC)
rf.performance.cv 
# estimate error rates from repeated cross-validation
rf.error.cv <- EstimateErrorRateCV(oof=oof,repeats=repeats)
rm(oof,repeats)


#save the model
save(rf.fit,file="Geno8RF.RData")
load("Geno4RF.RData")

## #------------------------------------------------------------------
## ################# Gradient Boosting #################
## #-----------------------------------------------------------------
#This didnt work weall with this data so it seems - did XGB instead

## # #set up GBM tuning paramters
gbm.grid <-  expand.grid(interaction.depth = c(1,3,5,7,9),
                          shrinkage = 0.1,
                        n.minobsinnode = c(10), n.trees = (0:10)*200) #n.tree tests 0,50...2500 # will stop when is 10 observation in terminal node
  nrow(gbm.grid)
 set.seed(123)
  gbm.fit <- train(data.train.descr, data.train.class,
                   method = "gbm",
                  metric = "ROC",
                   verbose = FALSE,
                   trControl = myControl,
                   ## Now specify the exact models
                   ## to evaludate:
                   tuneGrid = gbm.grid)
## #
 plot(gbm.fit)
 save(gbm.fit,file="CDV Age Sampled GBMAll.RData")
## # # estimate model performance in terms of a confusion matrix from repeated cross-validation
  oof <- gbm.fit$pred
## # # consider only the best tuned model
oof <- oof[intersect(which(oof$n.trees==gbm.fit$bestTune[,'n.trees']),which(oof$interaction.depth==gbm.fit$bestTune[,'interaction.depth'])),]
repeats <- gbm.fit$control$repeats
 gbm.performance.cv <- EstimatePerformanceCV(oof=oof,repeats=repeats)
 gbm.performance.cv
 # # estimate error rates from repeated cross-validation
 gbm.error.cv <- EstimateErrorRateCV(oof=oof,repeats=repeats)
 gbm.error.cv
 
 ## #------------------------------------------------------------------
 ## ################# XGBoost #################
 ## #------------------------------------------------------------------ 
 
 
 #XGBoost performs better than GBM with larger datasets 
 xgbGrid <-  expand.grid(nrounds = c(100,200),  
                         max_depth = c(10, 15, 20, 25),
                         colsample_bytree = seq(0.5, 0.9, length.out = 5),
                         ## The values below are default values in the sklearn-api. 
                         eta = 0.1,
                         gamma=0,
                         min_child_weight = 1,
                         subsample = 1)
nrow(xgbGrid)
set.seed(123) 
xgbfit <-  train( data.train.descr, data.train.class, trControl = myControl,
   tuneGrid = xgbGrid,
   method = "xgbTree")
 
 ## #
 plot(xgbGridfit)
 save(xgbfit,file="geno8_XGRboost.RData")
 ## # # estimate model performance in terms of a confusion matrix from repeated cross-validation
 oof <- xgbfit$pred
 ## # # consider only the best tuned model   -need to update oof with all variables in the grid
 oof <- oof[intersect(which(oof$nrounds==xgbfit$bestTune[,'nrounds']),which(oof$max_depth==xgbfit$bestTune[,'max_depth'])),]
 repeats <- xgbfit$control$repeats
 xgb.performance.cv <- EstimatePerformanceCV(oof=oof,repeats=repeats)
 xgb.performance.cv 
 # # estimate error rates from repeated cross-validation
 xgb.error.cv <- EstimateErrorRateCV(oof=oof,repeats=repeats)
 xgb.error.cv
 
## 
## #------------------------------------------------------------------
## ################# Support vector machine #################
## #------------------------------------------------------------------
## #
## #
set.seed(123)
svm.fit <- train(Class ~ .,data=cbind(Class=data.train.class,data.train.descr),
                 method = "svmRadial",
                 tuneLength = 9,
                 metric="ROC",
                 trControl = myControl)
plot(svm.fit)
## #
## # # estimate model performance in terms of a confusion matrix from repeated cross-validation
oof <- svm.fit$pred
## # # consider only the best tuned model
oof <- oof[intersect(which(oof$sigma==svm.fit$bestTune[,'sigma']),which(oof$C==svm.fit$bestTune[,'C'])),]
 repeats <- svm.fit$control$repeats
svm.performance.cv <- EstimatePerformanceCV(oof=oof,repeats=repeats)
 svm.performance.cv
# # estimate error rates from repeated cross-validation
 svm.error.cv <- EstimateErrorRateCV(oof=oof,repeats=repeats)
 rm(oof,repeats)
 # save the model
save(svm.fit, file="Geno8_SVM.RData")

## #------------------------------------------------------------------
## ################# Interpreting the best model #################
## #------------------------------------------------------------------

X <-dataR[-which(names(dataR) == "Class")] #load data again for the visualization
Y <- dataR$Class

# create the iml object
#create predictor object. Add your model object name here (GLM, RF, GBM or SVM)
mod <- Predictor$new(rf.fit, data = X, y = Y) #create predictor object
set.seed(123)
imp <-FeatureImp$new(mod, loss = "ce", compare='ratio', n.repetitions = 5) 
imp.dat<- imp$results
plot(imp)+ theme_bw()#plot results

#' 
#' ###4.2 Partial dependency (PD) plots and centered individual conditional expectation (cICE) plots.
#' 
#' PD plots and cICE plots provide a valuable way to visualize what effect each predictor in the model has in shaping exposure risk. Here we plot the two most important features in our RF model (but we provide code to plot all of the relevant features). This code plots the feature effect for both classes (negative and positive). The y-axis representing the probability of parvovirus exposure. The vertical lines on the x-axis (or 'rug plot') shows the distribution of observations across the feature. The highlighted line in the PD plot (i.e. model average) with the other lines (the ICE lines) reflect the predictions for one observation when we vary that feature.
#' 
#' The age exposed cICE plots demonstrate that, whilst exposure risk of parvovirus is higher in cubs ~2 years of age, some observations have a much more subtle peak. This could be evidence for an interaction effect (i.e., risk is lower for some cubs). CICE curves are diffiult to interpret for categorical predictors so we use standard ICE plots. The resultant box plot shows parvovirus exposure risk was higher in the 1992 epidemic compared to the others.
#' 
#' 
## ---- warning=F----------------------------------------------------------

#plot pd plots for the top predictors (any number appropriate for the data). Cateforical features don't plot properly

top<- imp.dat$feature[1:5]#n = number of predictors you want to display

ice_curves <- lapply(top, FUN = function(x) {
  cice <- partial(xgbfit, pred.var = x, center = TRUE, ice = TRUE, which.class="Positive",
                  prob = T) #note that we center values in the plotso these are centered ICE plots (cICE)
  autoplot(cice, rug = TRUE, train = dataR, alpha = 0.1) +
    theme_bw() +
    ylab("c-ICE")
})
#put them together
grid.arrange(grobs = c(ice_curves), ncol = 2)


#for plots of a particular feature

#rf.fit %>%
#partial(pred.var = "Age_exposed", ice=T, which.class="Positive", prob = T) %>%
#autoplot(rug = TRUE, center =TRUE, train = data, alpha = 0.1)+theme_bw() 

#' 
#' ###4.3 Interactions using Friedman's H index
#' 
#' Calculating Friedman's H index provides a robust way to assess the importance of interactions in shaping risk across models. The interactions identified can then be visualized using PD plots.
#' 
## ------------------------------------------------------------------------
set.seed(345)
mod <- Predictor$new(rf.fit, data = X, y=Y, type='prob', class='Positive') #we just want  positive class results now

interact <- Interaction$new(mod)
plot(interact)+ theme_bw()


##plot most important interactions. We could explore others though.
interact1 <- Interaction$new(mod, feature = "Cattle")
plot(interact1)+theme_bw()
#plot dominant interaction
pdp.obj <-  FeatureEffect$new(mod, feature = c("Cattle","Buffalo"), method='pdp')


plot(pdp.obj)+ scale_fill_gradient(low = "white", high = "red")+ theme_bw()


#' ###Explain single predictions using Shapely values
 
#' To better understand model predictions, the last step is to use a game theory and specifically Shapely values to understand how the model is applied to individual observations (see main text for more details).
#' 
## ------------------------------------------------------------------------

shapley <- Shapley$new(mod, x.interest = X[3000,])
shapley$plot()+ theme_bw()
results <- shapley$results #for each instance you can view these results as a table
head(results, n = 7L)
sum(results$phi)# <0 indicate model prediction was negative, > 0 model prediction was postive.

Y[3000] #see if observation was negative postive
## ------------------------------------------------------------------------
#shapley2 <- Shapley$new(mod, x.interest = X[131,]) #LUPINE (parvovirus positive)
#shapley2$plot()+ theme_bw()

