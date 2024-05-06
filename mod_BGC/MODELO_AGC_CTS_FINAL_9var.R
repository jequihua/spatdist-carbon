# Load packages.
library("here")
library("raster")
library("rgdal")
library("sp")
library("randomForest")
library("caret")
library("Metrics")
library("quantregForest")
library("snow")

DATAPATH = "mod_BGC"

# list feature rasters
raster_files = list.files(here(DATAPATH, "covariables_9var"),
                          pattern = "\\.tif$",
                          full.names = TRUE)

covs = brick()
for (j in 1:length(raster_files)){
  rast = raster(raster_files[j])
  covs = addLayer(covs,rast)
}

# Load train data points.
train_bgc1 <- readOGR(here(DATAPATH,"train_table_bgc.shp"))
extraction1 <- data.frame(extract(covs,train_bgc1))

# Add target variable.
extraction1$BGC <- train_bgc1@data$BGC

# Remove rows with missing values.
train_table = extraction1[complete.cases(extraction1),]
head(train_table)
# Categorical variables to factors.
train_table$Final_.veg <- factor(train_table$Final_.veg)
train_table$Final_.soils <- factor(train_table$Final_.soils)
train_table$Final_.geomor <- factor(train_table$Final_.geomor)


################################
# RF simple model with data frame (train_bgc1)
###############################
rf_bgc1 = randomForest(y = train_table$BGC,x = train_table[,-c(10)], ntree = 500, importance = TRUE, proximity=TRUE)

# Accuracy metrics.
plot(rf_bgc1)
cor(rf_bgc1$predicted,train_table$BGC)
sqrt(sum((rf_bgc1$predicted - train_bgc1$BGC)^2)/length(train_bgc1$BGC))
varImpPlot(rf_bgc1)
round(importance(rf_bgc1),2)

# Prediction.
pred_bgc1_simple <- predict(covs, rf_bgc1)
writeRaster(pred_bgc1_simple, here(DATAPATH, "rfmodel_bgc_9var_character.tif"),overwrite=TRUE)

#################################
# rf with Cross validation (10fold)
#################################
# As spatialpoints 
train_bgc2 <- readOGR(paste0(getwd(), "/datos_campo/","train_table_bgc.shp"))
str(train_bgc2)
class(train_bgc2)
train_bgc2@proj4string
train_bgc2 <- extract(covs, train_bgc2, sp=TRUE)
str(train_bgc2)
#train_bgc2$Final__veg <- factor(train_bgc2$Final__veg)    # Agregado
#train_bgc2$Final__soils <- factor(train_bgc2$Final__soils)
#train_bgc1$Final__geomor <- factor(train_bgc1$Final__geomor)
#train_bgc1$Final__veg <- as.character(train_bgc1$Final__veg)    # Agregado
#train_bgc1$Final__soils <- as.character(train_bgc1$Final__soils)
#train_bgc1$Final__geomor <- as.character(train_bgc1$Final__geomor)

fm = as.formula(paste("BGC ~", paste0(names(covs[[c(1:9)]]),
                                           collapse = "+")))

ctrl <- trainControl(method = "cv", savePred=T)
rfmodel <- train(fm, data=train_bgc2@data, method = "rf", trControl = ctrl,
                 importance=TRUE)
rfmodel
varImpPlot(rfmodel[11][[1]])
plot(rfmodel)

# Predicion
pred <- predict(covs, rfmodel)
plot(pred, main="BGC")

writeRaster(pred, file=paste0(getwd(), "/Modelo_raster/results_9var/","rfmodel_bgc_cv_9var_num.tif"),overwrite=TRUE)

#################################
#SENSITIVITY + UNCERTAINTY
#################################
dat <- train_bgc2
str(dat)
validation_bgc_cv_9var <- data.frame(rmse=numeric(), r2=numeric())
for (i in 1:10){
  # We will build 10 models using random samples of 25%
  smp_size <- floor(0.25 * nrow(dat))
  train_ind <- sample(seq_len(nrow(dat)), size = smp_size)
  train <- dat[train_ind, ]
  test <- dat[-train_ind, ]
  modn <- train(fm, data=train@data, method = "rf",
                trControl = ctrl)
  pred <- stack(pred, predict(covs, modn))
  test$pred <- extract(pred[[i+1]], test)
  # Store the results in a dataframe
  validation_bgc_cv_9var[i, 1] <- rmse(test$BGC, test$pred)     #cambiar ABGlog
  validation_bgc_cv_9var[i, 2] <- cor(test$BGC, test$pred)^2    #cambiar ABGlog
}
summary(validation_bgc_cv_9var)
#################

sensitivity <- calc(pred[[-1]], sd)
plot(sensitivity, col=rev(topo.colors(10)),
     main='Sensitivity based on 10 realizations using 25% samples')

writeRaster(sensitivity, file= paste0(getwd(), "/Modelo_raster/results_9var/","rfmodel_bgc_9var_sensitivity.tif"),
            overwrite=TRUE)

prediction75 <- (pred[[4]])  #quitar exp si no es LOG
plot(prediction75, main='AGC prediction based on 75% of data',
     col=rev(topo.colors(10)))
writeRaster(prediction75, file=paste0(getwd(), "/Modelo_raster/results_6var/","cv_rf_log_6var_predict_75.tif"),
            overwrite=TRUE)

########################
# Julian
########################
str(train_bgc1)
rf_quant <- quantregForest(x=train_bgc1[,-10],y=train_bgc1$BGC, 
                           ntree=1000) 

rf_quant1<- as.data.frame(predict(rf_quant,train_bgc1,what = c(0.05,0.95))) # aqu? estimas la predicci?n para el perenctil 0.05 y 0.95 y 
intervalo_de_prediccion_90 <- rf_quant1$`quantile= 0.95`-rf_quant1$`quantile= 0.05`
plot(intervalo_de_prediccion_90)
print(rf_quant_prediction)
plot(rf_quant1)
######
str(intervalo_de_prediccion_90)
#################################
#ESTIMATE THE FULL CONDITIONAL DISTRIBUTION OF BGC
#################################

model <- quantregForest(y=dat@data$BGC, x=dat@data[,3:11],   
                        ntree=1000, keep.inbag=TRUE,
                        mtry = as.numeric(rfmodel$bestTune))

beginCluster()
unc <- clusterR(covs, predict, args=list(model=model,what=sd))
mean <- clusterR(covs, predict, args=list(model=model, what=mean))
unc <- unc + sensitivity
Total_unc_Percent <- (unc)/(mean)   #QUITAR o poner el exp si hay o no el  log  exp(unc)/exp(mean)
endCluster()
plot(Total_unc_Percent, main='Total uncertainty')

plot((mean), main='AGC mean based in all data', col=rev(topo.colors(10))) #Quitar el exp si no hay log
writeRaster(exp(mean), file=paste0(getwd(), "/Modelo_raster/results_6var/","cv_log_6var_mean_quantrf.tif"),
            overwrite=TRUE)


plot(unc, zlim=c(0,10), main='Uncertainty')
summary(Total_unc_Percent)
summary(unc)

writeRaster(Total_unc_Percent, file=paste0(getwd(), "/Modelo_raster/results_9var/","rfmodel_bgc_9var_total_uncertainty.tif"),
            overwrite=TRUE)

writeRaster(unc, file='bgc_rf_unc.tif',
            overwrite=TRUE)


plot(unc)
plot(mean)
plot(Total_unc_Percent)
