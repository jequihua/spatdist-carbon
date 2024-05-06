# Load packages
library("here")
library("raster")
library("rgdal")
library("sp")
library("randomForest")
library("caret")
library("Metrics")
library("quantregForest")
#library("snow")

DATAPATH = "mod_BGC"

# Load vegetatopm raster
veg = raster("D://Documentos personales/Doctorado_2014/Cartografia/AGB/MODELO_FINAL/CovariablesModel_original/Final_Veget_conabio_90m.tif")

# list feature rasters
raster_files = list.files(here(DATAPATH, "covariables_9var"),
                          pattern = "\\.tif$",
                          full.names = TRUE)
here(DATAPATH, "/Cov_mod_6var/")
raster_files

covs = brick()
for (j in 1:length(raster_files))
{
  print(j)
  
  name = paste("Final_",name= basename(raster_files[j]))
  
  rast = raster(raster_files[j])

  covs = addLayer(covs,rast)
  
}

# Assign correct names to layers.
names(covs) <- c("CHM", "NDVI", "OPT_RED", "RATIO_vv_v", "DIFER_vh_v", "NDVI_vv_vh")


names(covs)
covs
plot(covs)

# load train data points
train1 = readOGR(paste0(getwd(), "/Cov_mod_6var/","train_table_6var.shp"))
train1$AGClog <- log(train1$AGC)
str(train1)

# As data frame 
train2 <- as.data.frame(train1)
train2 <- train2[,-(11:12)]
str(train2)
train_table = train2[complete.cases(train2),]
rm(train_table)
################################
# RF simple model with data frame (train2)
###############################
rf_agc1 = randomForest(y = train2$AGC, x = train2[,c(4:9)],ntree = 500,
                       importance = TRUE, proximity=TRUE)
rf_agc1

cor(rf_agc1$predicted,train2$AGC)
sqrt(sum((rf_agc1$predicted - train2$AGC)^2)/length(train2$AGC))
varImpPlot(rf_agc1)
plot(rf_agc1)
round(importance(rf_agc1), 2)

# Predicion
pred_simple <- predict(covs, rf_agc1)
#pred_simple <- exp(pred_simple)
rfmodel_6var <- pred_simple
plot(pred_simple)
writeRaster(pred_simple, file= paste0(getwd(), "/Modelo_raster/results_6var/","rfmodel_6var.tif"),overwrite=TRUE)


#################################
# rf with Cross validation (10fold)
#################################
fm = as.formula(paste("log(AGC) ~", paste0(names(covs),  #(covs[[-c(4,5,6,7,8,9,13)]])
                                        collapse = "+")))

ctrl <- trainControl(method = "cv", savePred=T)
rfmodel <- train(fm, data=train1@data, method = "rf", trControl = ctrl,
                 importance=TRUE)
rfmodel
varImpPlot(rfmodel[11][[1]])
plot(rfmodel)

# Predicion
pred <- predict(covs, rfmodel)
pred <- exp(pred)
plot(pred, main="AGC")

writeRaster(pred, file=paste0(getwd(), "/Modelo_raster/results_6var/","cv_rf_log_6var.tif"),overwrite=TRUE)

#################################
#SENSITIVITY + UNCERTAINTY
#################################
dat <- train1
str(dat)
validation_cv_log_6var <- data.frame(rmse=numeric(), r2=numeric())
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
  validation_cv_log_6var[i, 1] <- rmse(test$AGClog, test$pred)     #cambiar ABGlog
  validation_cv_log_6var[i, 2] <- cor(test$AGClog, test$pred)^2    #cambiar ABGlog
}

summary(validation_cv_log_6var)
#################

sensitivity <- calc(pred[[-1]], sd)
plot(sensitivity, col=rev(topo.colors(10)),
     main='Sensitivity based on 10 realizations using 25% samples')

writeRaster(sensitivity, file= paste0(getwd(), "/Modelo_raster/results_6var/","cv_rf_log_6var_sensit_25.tif"),
            overwrite=TRUE)

prediction75 <- exp(pred[[4]])  #quitar exp si no es LOG
plot(prediction75, main='AGC prediction based on 75% of data',
     col=rev(topo.colors(10)))
writeRaster(prediction75, file=paste0(getwd(), "/Modelo_raster/results_6var/","cv_rf_log_6var_predict_75.tif"),
            overwrite=TRUE)

#################################
#ESTIMATE THE FULL CONDITIONAL DISTRIBUTION OF AGClog
#################################
str(dat)
model <- quantregForest(y=dat@data$AGClog, x=dat@data[,4:9],   #CAMBIAR ABGlog
                        ntree=1000, keep.inbag=TRUE,
                        mtry = as.numeric(rfmodel$bestTune))

beginCluster()
unc <- clusterR(covs, predict, args=list(model=model,what=sd))
mean <- clusterR(covs, predict, args=list(model=model, what=mean))
unc <- unc + sensitivity
Total_unc_Percent <- exp(unc)/exp(mean)   #QUITAR el exp si no hay log
endCluster()
plot(Total_unc_Percent, zlim=c(0,0.2), main='Total uncertainty')

plot((mean), main='AGC mean based in all data', col=rev(topo.colors(10))) #Quitar el exp si no hay log
writeRaster(exp(mean), file=paste0(getwd(), "/Modelo_raster/results_6var/","cv_log_6var_mean_quantrf.tif"),
            overwrite=TRUE)


plot(unc, zlim=c(0,10), main='Uncertainty')
summary(Total_unc_Percent)
summary(unc)

writeRaster(Total_unc_Percent, file=paste0(getwd(), "/Modelo_raster/results_6var/","cv_rf_logo_6var_total_uncertainty.tif"),
            overwrite=TRUE)

writeRaster(unc, file='AGC_rf_unc.tif',
            overwrite=TRUE)

###########################
#Graphs with ggplot2
###########################
r_points = rasterToPoints(pred_log)
r_df = data.frame(r_points)
head(r_df) #breaks will be set to column "layer"
r_df$cuts=cut(r_df$layer,breaks=c(1,10,25,50,80,120,150,200,240)) #set breaks

ggplot(data=r_df) + 
  geom_tile(aes(x=x,y=y,fill=cuts)) + 
  scale_fill_brewer("AGClog (Mg/ha)",type = "seq", palette = "YlOrRd") +
  coord_equal() +
  theme_bw() +
  theme(panel.grid.major = element_blank()) +
  xlab("Longitude") + ylab("Latitude")

