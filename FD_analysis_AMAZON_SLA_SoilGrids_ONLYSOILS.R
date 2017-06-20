
setwd("D:/AMAZON/ATDN")

source("ATDN_functions_J.R")
source("Diversity_Functions.R")
source("ATDN_Read_Composition.R")
source("Diversity_Functions.R")
library(raster)
library(FD)
library(raster)
library(boot)
library(arm)
library(ggdendro)
library(lattice)
library(lme4)  
library(ggplot2)
library(sp)
library(gstat)
library(scales) 
library(INLA)
library(ggmap)
library(dismo)
library(usdm)
library(reshape2)
source("./HighstatLibV10.R")


###### 
#read input files for this script and make vegtab!
#      make vector with plots to be removed

load("D:/AMAZON/ATDN/Traits_VegTable_12June17/VegtabGen20170613.RData")
data_env = data_env
vegtab   = vegtab; vegtab[1:10,1:10]

######################################## Show map with all plots
###### 
show.map.plots(); add.geography()

tree.traits = read.table("D:/AMAZON/ATDN/Traits_VegTable_12June17/Traits4ATDN.csv", header = T, sep = ",", as.is = T)
tree.traits[tree.traits == -99] = NA
names(tree.traits)
head(tree.traits)

sp     = 0.25
mapres = 0.5
#######################################
###############Include UTM coordinates to Plots dataframe ##################

utmcoor <- SpatialPoints(cbind(data_env$Longitude, data_env$Latitude), proj4string = CRS("+proj=longlat"))
longlatcoor   <- spTransform(utmcoor, CRS("+init=epsg:29193"))
data_env$Longitude_m <- longlatcoor@coords[,1]
data_env$Latitude_m  <- longlatcoor@coords[,2]
head(data_env)

###ploting the locations on the map
range(data_env$Longitude)
range(data_env$Latitude)
glgmap   <- get_map(location = c( -100, -6, -20, -15),
                    zoom = 4,
                    maptype= "terrain")       
p <- ggmap(glgmap)
p <- p + geom_point(aes(Longitude, Latitude), 
                    data = data_env) 
p <- p + xlab("Longitude") 
p <- p + ylab("Latitude")  
p <- p + theme(text = element_text(size=15)) 
p
###############      END of maps Amazon and plot locations #################


###############################################################################################################
################################### FUNCTIONAL DIVERSITY ANALYSIS #############################################
###############################################################################################################


### Prepare data ###
# All the CW mean data was calculated in the script for Wood Density
# so that part of the script is not here. Also for a complete script about 
# the spatial data, resolution and VIF values for variable selection see the script
# "FD_analysis_AMAZON_WD_SoilGrids.R"

########################################################################
## Relation FD or CWmean and Environment #########################################
########################################################################
setwd('D:/AMAZON/ATDN')

#Environmental data
xx <- stack('Biovars_CHIRPS_and_MaxMinTempBioclim/brazil_CHIRPS_MaxMinTempBIOCLIM_biovars.img')
cc <- stack(list.files('./CloudCover_Jetz2017/', pattern="*.img$", full.names = T)) #cloud cover
sr <- raster('./PET_Aridity_rasters/SolarRadi_sa.img')                              #average solar radiation
ar <- raster('./PET_Aridity_rasters/AridityAnnual_sa.img')                          #aridity (+ values =  >humidity)
pet<- raster('./PET_Aridity_rasters/PET_sa.img')                                    #potential evapotranspiration    
soil<-stack(list.files('./Soils_SOILGRIDS/recent/USING/', pattern="*_cropped.img$",full.names = T))
soil1<-stack(list.files('./Soils_SOILGRIDS/recent/SoilTypesProbability/cropped/', pattern="*.img$",full.names = T))
#soil<-raster('./Soils_SOILGRIDS/recent/SoilTypesProbability/TAXNWRB_250m_ll.tif')


#----Mask to South America----#
#sa<-shapefile('SouthAmericaStudyArea_JA.shp')
#soil<-crop(soil,sa)
forestborder = readOGR("../ATDN/ESRI_WORLD","forestborder")
xs<-mask(xx,forestborder)

#set extent
extend_use <- extent(xs)                                                                     
extent(cc) <- extend_use; extent(sr) <- extend_use; extent(ar) <- extend_use; 
extent(pet)<- extend_use; extent(soil) <- extend_use; extent(soil1) <- extend_use

#Specify new resolution based on CHIRPS data
cc1<-projectRaster(cc,xs);sr1<-projectRaster(sr,xs);ar1<-projectRaster(ar,xs);pet1<-projectRaster(pet,xs);
soil_a<-projectRaster(soil,xs); soil_b<-projectRaster(soil1,xs)

#stack all predictors
env_var<- stack(cc1,ar1,sr1,pet1,soil_a,soil_b,xs); env_var<-mask(env_var,forestborder)

######################
# Collinearity #######
######################
library(usdm)
names(env_var)
out<-vifstep(env_var,th=3)
env_var1<-exclude(env_var,out)# or manual as: env_var1 <- env_var[[-c(4....)]]

names(env_var1)
names(env_var1)<- c("MODCF_meanannual","PET","BLDFIE_M_sl6",
                    "CECSOL_M_sl6","CLYPPT_M_sl4", "ORCDRC_M_sl4","PHIHOX_M_sl6", "AA","CaH","CG",                                        
                   "CrH","FA","GP","HA","HAC","HFA","HFC","HFD","HFE","HG","HGD","HGE","HP","MG","PA","SH","UG",
                    "CHIRPS_bio4","CHIRPS_bio8","CHIRPS_bio13","CHIRPS_bio18","CHIRPS_bio19")

#Extract env var info for plots
fd_env <- read.csv("D:/AMAZON/ATDN/ResultsFD/CWM_WDSMCSLA_Amazon_13June17.csv")#Remember only CW mean now
nrow((fd_env))

coordinates(fd_env) <-c("Longitude","Latitude") 
env_var_pl <- extract(env_var1, fd_env, all=T)
env_var_pl <- as.data.frame(env_var_pl)

fd_env<-as.data.frame(fd_env) #Merge plot data with environmental data
fd_env<-cbind(fd_env, env_var_pl)
head(fd_env)

fd_env1 <-na.omit(fd_env) #eliminate NA's & outliers from table
nrow(fd_env1)#1617 plots with full information

MyVar <- names(fd_env1[,c(11:42)])

MyMultipanel.ggp2(Z = fd_env1,  varx = MyVar, 
                  vary = "SLA", 
                  ylab = "SLA",addSmoother = F,addRegressionLine = T, addHorizontalLine = TRUE)

#Delete variables with no much variation
fd_env1<-fd_env1[,c(-20,-21,-25,-26,-35)]
names(fd_env1)

# We better standardize the covariates to avoid numerical problems.
fd_env1$MODCF_meanannual.std  <- MyStd(fd_env1$MODCF_meanannua)
fd_env1$PET.std               <- MyStd(fd_env1$PET)
fd_env1$BLDFIE_M_sl6.std      <- MyStd(fd_env1$BLDFIE_M_sl6)
fd_env1$CECSOL_M_sl6.std      <- MyStd(fd_env1$CECSOL_M_sl6)
fd_env1$CLYPPT_M_sl4.std      <- MyStd(fd_env1$CLYPPT_M_sl4)
fd_env1$ORCDRC_M_sl4.std      <- MyStd(fd_env1$ORCDRC_M_sl4)
fd_env1$PHIHOX_M_sl6.std      <- MyStd(fd_env1$PHIHOX_M_sl6)
fd_env1$AA.std      <- MyStd(fd_env1$AA)
fd_env1$CaH.std     <- MyStd(fd_env1$CaH)
fd_env1$FA.std      <- MyStd(fd_env1$FA)
fd_env1$GP.std      <- MyStd(fd_env1$GP)
fd_env1$HA.std     <- MyStd(fd_env1$HA)
fd_env1$HFC.std     <- MyStd(fd_env1$HFC)
fd_env1$HFD.std     <- MyStd(fd_env1$HFD)
fd_env1$HFE.std     <- MyStd(fd_env1$HFE)
fd_env1$HG.std      <- MyStd(fd_env1$HG)
fd_env1$HGD.std     <- MyStd(fd_env1$HGD)
fd_env1$HGE.std     <- MyStd(fd_env1$HGE)
fd_env1$HP.std      <- MyStd(fd_env1$HP)
fd_env1$MG.std      <- MyStd(fd_env1$MG)
fd_env1$SH.std      <- MyStd(fd_env1$SH)
fd_env1$UG.std      <- MyStd(fd_env1$UG)
fd_env1$CHIRPS_bio4.std       <- MyStd(fd_env1$CHIRPS_bio4)
fd_env1$CHIRPS_bio8.std      <- MyStd(fd_env1$CHIRPS_bio8)
fd_env1$CHIRPS_bio13.std      <- MyStd(fd_env1$CHIRPS_bio13)
fd_env1$CHIRPS_bio18.std      <- MyStd(fd_env1$CHIRPS_bio18)
fd_env1$CHIRPS_bio19.std      <- MyStd(fd_env1$CHIRPS_bio19)


###########
# outlier #  NOT DONE THIS TIME for pH
###########
MyVar <- names(fd_env1[,c(27:41)])

MyMultipanel.ggp2(Z = fd_env1,  varx = MyVar, 
                  vary = "SLA", 
                  ylab = "SLA",addSmoother = TRUE,addRegressionLine = T, addHorizontalLine = TRUE)

##############
# Models
##############

I1 <- inla(SLA ~ MODCF_meanannual.std + PET.std + BLDFIE_M_sl6.std + 
             CECSOL_M_sl6.std + CLYPPT_M_sl4.std + ORCDRC_M_sl4.std + PHIHOX_M_sl6.std + 
             AA.std + CaH.std + FA.std + GP.std + HFC.std + HA.std + HFD.std + HFE.std + HG.std + HGD.std + HGE.std + 
             HP.std + MG.std + SH.std + UG.std + CHIRPS_bio4.std +     
             CHIRPS_bio8.std + CHIRPS_bio13.std + CHIRPS_bio18.std + CHIRPS_bio19.std,
             control.predictor = list(compute = TRUE, quantiles = c(0.025, 0.975)),
             control.compute = list(dic = TRUE, waic = TRUE),data = fd_env1)

summary(I1)

# Assess overdispersion using frequentist way of thinking
Fit1 <- I1$summary.fitted.values[1:1617, "mean"]
E1   <- (fd_env1$SLA - Fit1) / sqrt(Fit1)
N    <- nrow(fd_env1)
p    <- nrow(I1$summary.fixed)
sum(E1^2) / (N - p)# much smaller than 1 so no problem

##################################
# Model validation
# A. Outliers?

par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
plot(x = Fit1, 
     y = E1,
     xlab = "Fitted values",
     ylab = "Pearson residuals")
abline(h = 0, lty = 2)     
# No clear outliers?

# Plot residuals vs each covariate in the model
# And vs each covariate not in the model.
MyVar <- names(fd_env1[,c(27:41)])
fd_env1$E1 <- E1
MyMultipanel.ggp2(Z = fd_env1, 
                  varx = MyVar, 
                  vary = "E1", 
                  ylab = "Pearson residuals",
                  addSmoother = FALSE,
                  addRegressionLine = FALSE,
                  addHorizontalLine = TRUE)
# Nothing major jumps out?


# Spatial dependency?
# The covariates do not differ within a year.
# Hence, any differences between observations made
# in the same year end up in the Pearson
# residuals.

# Adding a spatial correlation will capture missing
# spatial covariates and real spatial correlation.


# Check for spatial correlation
MyData1 <- data.frame(E1 = E1, 
                      Xkm = fd_env1$Longitude_m / 1000, 
                      Ykm = fd_env1$Latitude_m / 1000)
coordinates(MyData1) <- c("Xkm", "Ykm")
V1 <- variogram(E1 ~ 1, 
                MyData1, 
                cressie = TRUE)
plot(V1, 
     main = "", 
     xlab = list(label = "Distance", cex = 1.5), 
     ylab = list(label = "Semi-variogram", cex = 1.5),
     pch = 16,
     col = 1,
     cex = 1.5
)
# Inconclusive? or up to ~150k there is strong spatial autocorrelation ?


####################################################
# We will implement the following 8 steps.
# 1. Make a mesh.
# 2. Define the weighting factors a_ik (also called 
#    the projector matrix).
# 3. Define the SPDE.
# 4. Define the spatial field.
# 5. Make a stack. In this process we tell INLA 
#    at which points on the mesh we sampled the 
#    response variable and the covariates. We 
#    also need to inform INLA (via the stack) at 
#    which points of the mesh we have any other 
#    terms, e.g. the random effects as we know 
#    them from Chapter 8.
# 6. Specify the model formula in terms of the 
#    response variable, covariates and the 
#    spatial correlated term.
# 7. Run the spatial model in INLA.
# 8. Inspect the results.



########################
#1. Make a mesh.
#   Step 1 of making a mesh:  Get a 
#   sense for the distribution of 
#   distances between sampling locations. 

Loc   <- cbind(fd_env1$Longitude_m, fd_env1$Latitude_m)
Loc

mesh1 <- inla.mesh.2d(Loc, max.edge=c(100000, 200000), cutoff = 0)
par(mfrow = c(1,1), mar = c(0,0,2,0))
plot(mesh1)
points(Loc, col = 1, pch = 16, cex = 1)
mesh1$n #number of nodes in the mesh

#what are the distances between the points?
Locations <- fd_env1[,c("Latitude_m", "Longitude_m")]
D <- dist(Locations)
par(mfrow = c(1,2), mar = c(5,5,2,2), cex.lab = 1.5)
hist(D, 
     freq = TRUE,
     main = "", 
     xlab = "Distance between sites (m)",
     ylab = "Frequency")

plot(x = sort(D), 
     y = (1:length(D))/length(D), 
     type = "l",
     xlab = "Distance between sites (m)",
     ylab = "Cumulative proportion")


# 2. Define the weighting factors a_ik (also called 
#    the projector matrix).

A2 <- inla.spde.make.A(mesh1, loc = Loc)
dim(A2)  #1694 observations on a 14735 grid


# 3. Define the SPDE. (Stokastic partial differential equation)
spde   <- inla.spde2.matern(mesh1, alpha = 2)


# 4. Define the spatial field.
w.index <- inla.spde.make.index(
  name    = 'w', 
  n.spde  = spde$n.spde,
  n.group = 1,
  n.repl  = 1)


# 5. Make a stack. In this process we tell INLA 
#    at which points on the mesh we sampled the 
#    response variable and the covariates. We 
#    also need to inform INLA (via the stack) at 
#    which points of the mesh we have any other 
#    terms, e.g. the random effects as we know 
#    them from Chapter 8.

# Create a data frame with an intercept and covariates.
N <- nrow(fd_env1)                
X <- data.frame(Intercept   = rep(1, N), 
                MODCF_meanannual.std = fd_env1$MODCF_meanannual.std,
                PET.std              = fd_env1$PET.std,
                BLDFIE_M_sl6.std     = fd_env1$BLDFIE_M_sl6.std,
                CECSOL_M_sl6.std     = fd_env1$CECSOL_M_sl6.std,
                CLYPPT_M_sl4.std     = fd_env1$CLYPPT_M_sl4.std,
                ORCDRC_M_sl4.std     = fd_env1$ORCDRC_M_sl4.std,
                PHIHOX_M_sl6.std     = fd_env1$PHIHOX_M_sl6.std,
                AA.std               = fd_env1$AA.std,
                CaH.std              = fd_env1$CaH.std,
                FA.std               = fd_env1$FA.std,
                GP.std               = fd_env1$GP.std,
                HA.std               = fd_env1$HA.std,
                HFC.std              = fd_env1$HFC.std,
                HFD.std              = fd_env1$HFD.std,
                HFE.std              = fd_env1$HFE.std,
                HG.std               = fd_env1$HG.std,
                HGD.std              = fd_env1$HGD.std,
                HGE.std              = fd_env1$HGE.std,
                HP.std               = fd_env1$HP.std,
                MG.std               = fd_env1$MG.std,
                SH.std               = fd_env1$SH.std,
                UG.std               = fd_env1$UG.std,
                CHIRPS_bio4.std      = fd_env1$CHIRPS_bio4.std,
                CHIRPS_bio8.std      = fd_env1$CHIRPS_bio8.std,
                CHIRPS_bio13.std     = fd_env1$CHIRPS_bio13.std,
                CHIRPS_bio18.std     = fd_env1$CHIRPS_bio18.std,
                CHIRPS_bio19.std     = fd_env1$CHIRPS_bio19.std
)


Stack2 <- inla.stack(
  tag  = "Fit",
  data = list(y = fd_env1$SLA),  
  A    = list(1,A2),                      
  effects = list(
    X = X,                #Covariates
    w = w.index))         #Spatial field  


# 6. Specify the model formula in terms of the 
#    response variable, covariates and the 
#    spatial correlated term.


f1 <- y ~ -1 + Intercept + 
  MODCF_meanannual.std + PET.std + BLDFIE_M_sl6.std + 
  CECSOL_M_sl6.std + CLYPPT_M_sl4.std + ORCDRC_M_sl4.std + PHIHOX_M_sl6.std + 
  AA.std + CaH.std + FA.std + GP.std + HFC.std + HA.std + HFD.std + HFE.std + HG.std + HGD.std + HGE.std + 
  HP.std + MG.std + SH.std + UG.std + CHIRPS_bio4.std +     
  CHIRPS_bio8.std + CHIRPS_bio13.std + CHIRPS_bio18.std + CHIRPS_bio19.std

f2 <- y ~ -1 + Intercept + 
  MODCF_meanannual.std + PET.std + BLDFIE_M_sl6.std + 
  CECSOL_M_sl6.std + CLYPPT_M_sl4.std + ORCDRC_M_sl4.std + PHIHOX_M_sl6.std + 
  AA.std + CaH.std + FA.std + GP.std + HFC.std + HA.std + HFD.std + HFE.std + HG.std + HGD.std + HGE.std + 
  HP.std + MG.std + SH.std + UG.std + CHIRPS_bio4.std +     
  CHIRPS_bio8.std + CHIRPS_bio13.std + CHIRPS_bio18.std + CHIRPS_bio19.std + 
  f(w, model = spde)

f3 <- y ~ -1 + Intercept +
  f(w, model = spde)

# 7. Run the non spatial and spatial model in INLA.
I1 <- inla(f1,
           family = "gaussian", 
           data = inla.stack.data(Stack2),
           control.compute = list(dic = TRUE, waic = TRUE),
           control.predictor = list(A = inla.stack.A(Stack2)),
           control.inla = list(strategy = "gaussian"))  

I2 <- inla(f2,
           family = "gaussian", 
           data = inla.stack.data(Stack2),
           control.compute = list(dic = TRUE, waic = TRUE),
           control.predictor = list(A = inla.stack.A(Stack2)),
           control.inla = list(strategy = "gaussian"))  

I3 <- inla(f3,
           family = "gaussian", 
           data = inla.stack.data(Stack2),
           control.compute = list(dic = TRUE, waic = TRUE),
           control.predictor = list(A = inla.stack.A(Stack2)),
           control.inla = list(strategy = "gaussian"))  

#Inspect the results
summary(I2)

I1$summary.fixed[, c("mean", "0.025quant", "0.975quant")]
I2$summary.fixed[, c("mean", "0.025quant", "0.975quant")]
I3$summary.fixed[, c("mean", "0.025quant", "0.975quant")]

# And compare the two models with DICs and WAICs
dic2  <- c(I1$dic$dic, I2$dic$dic, I3$dic$dic)
waic2 <- c(I1$waic$waic, I2$waic$waic, I3$waic$waic)
Z     <- cbind(dic2, waic2)
rownames(Z) <- c("Gaussian GLM",  "Spatial-Soil GLM", "Spatial only")
Z
# Conclusion: 
# The model with soil and spatial correlation is better.


# Assess overdispersion using frequentist way of thinking
Fit1 <- I1$summary.fitted.values[1:1617, "mean"];E2<- (fd_env1$SLA - Fit1) / sqrt(Fit1);p<- nrow(I1$summary.fixed);N<- nrow(fd_env1)
sum(E2^2) / (N - p)#well below 1 which is good

Fit2 <- I2$summary.fitted.values[1:1617, "mean"];E2<- (fd_env1$SLA - Fit2) / sqrt(Fit2);p<- nrow(I2$summary.fixed);N<- nrow(fd_env1)
sum(E2^2) / (N - p)#well below 1 which is good

Fit3 <- I3$summary.fitted.values[1:1617, "mean"];E2<- (fd_env1$SLA - Fit3) / sqrt(Fit3);p<- nrow(I3$summary.fixed);N<- nrow(fd_env1)
sum(E2^2) / (N - p)#well below 1 which is good

# Coefficient of determination
SSTotal <- var( fd_env1$SLA ) * (nrow(fd_env1)-1)
SSE     <- sum((fd_env1$SLA - Fit2)^2)
r2      <- 1-(SSE/SSTotal)

#or the same below
R2_I1 <- 1 - (sum((fd_env1$SLA - Fit1 )^2)/sum((fd_env1$SLA-mean(fd_env1$SLA))^2))
R2_I2 <- 1 - (sum((fd_env1$SLA - Fit2 )^2)/sum((fd_env1$SLA-mean(fd_env1$SLA))^2))
R2_I3 <- 1 - (sum((fd_env1$SLA - Fit3 )^2)/sum((fd_env1$SLA-mean(fd_env1$SLA))^2))

##################################
# Model validation
# A. Outliers?
par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
plot(x = Fit2, 
     y = E2,
     xlab = "Fitted values",
     ylab = "Pearson residuals")
abline(h = 0, lty = 2)     
# No clear outliers?

# Plot residuals vs each covariate in the model
# And vs each covariate not in the model.
fd_env1$E2 <- E2
MyMultipanel.ggp2(Z = fd_env1, 
                  varx = MyVar, 
                  vary = "E2", 
                  ylab = "Pearson residuals",
                  addSmoother = FALSE,
                  addRegressionLine = FALSE,
                  addHorizontalLine = TRUE)

####If I want to extract the model Fit to the table and see
####how much my model fit the data use:
####e.g  mymodelfit<- cbind(fd_env1,Fit2); where Fit2 is defined above in line 532 

########################################
# Let's plot the results of the model, without and with
# the spatial correlation side by side.
# Better not try to understand all this R code.
# RUN FROM HERE......

NumberOfBetas <- nrow(I2$summary.fixed) 
Combined <- rbind(I1$summary.fixed[, c("mean", "0.025quant", "0.975quant")],
                  I2$summary.fixed[, c("mean", "0.025quant", "0.975quant")],
                  I3$summary.fixed[, c("mean", "0.025quant", "0.975quant")]
)
Combined$WhichModel <- rep(c("GLM","Soil-Spatial GLM", "Spatial GLM"), 
                           each = NumberOfBetas)
Combined$WhichVariable <- rep(rownames(I2$summary.fixed), 2)
colnames(Combined) <- c("Mean", "Lo", "Up", "WhichModel", "WhichVariable")
Combined


p <- ggplot()
p <- p + geom_point(data = Combined,
                    aes(x = WhichModel, 
                        y = Mean)
)

p <- p + geom_errorbar(data = Combined,
                       aes(x = WhichModel, 
                           ymax = Up, 
                           ymin = Lo), 
                       width=0.2)

p <- p + xlab("Parameters") + ylab("Values") + ggtitle("CWM SLA 14/06/17")
p <- p + theme(text = element_text(size = 15)) 
p <- p + facet_wrap( ~ WhichVariable, scales = "free_y")
p <- p + theme(legend.position="none") 
p

##############################
# Plot the spatial field
##############################
w.pm <- I2$summary.random$w$mean  
length(w.pm)

# There are various ways to plot the spatial 
# field. One option is to use the function 
# inla.mesh.projector; it creates a lattice 
# using the mesh and specified ranges. 
w.proj <- inla.mesh.projector(mesh1,
                              xlim = range(Loc[,1]),
                              ylim = range(Loc[,2])) 

w.pm100_100 <- inla.mesh.project(w.proj, w.pm)

# This w.pm100_100 is of dimension 100 by 100 
# and is a projection (interpolation and extrapolation) 
# of the random field w. We can use the levelplot 
# function from the lattice package to plot w.pm100_100

library(grid)
grid <- expand.grid(x = w.proj$x, 
                    y = w.proj$y)
grid$z <- as.vector(w.pm100_100)               
levelplot(z ~ x * y,
          data = grid, 
          scales = list(draw = TRUE),
          xlab = list("Longitude", cex = 1.5),
          ylab = list("Latitude", cex = 1.5),
          main = list("Posterior mean spatial random field", cex = 1.5),
          panel=function(...){
            panel.levelplot(...)
            grid.points(fd_env$Longitude_m,fd_env$Latitude_m, pch=1)  
          })


# Plot results in the map of Brazil
# The code fo this section works fine.
# it runs and it makes a nice picture.


# Add the values of u_i (the spatial random effects)
# - Big points for large u_i
# - Small dots for small u_i
# - Different colour for positive and negative u_i 
# The u.pm below is the spatial random effect at the 
# sampling locations.
xy     <- cbind(fd_env1$Longitude_m, fd_env1$Latitude_m)
u.proj <- inla.mesh.projector(mesh1, loc = Loc)
u.pm   <- inla.mesh.project(u.proj, I3$summary.random$w$mean)


# Use different font sizes and colours dependening on the values of g.mean
MyCex <- 2 * abs(u.pm) / max(u.pm)
SignResidual <- as.numeric(u.pm >=0) + 1
MyCol <- c("black", "yellow")[SignResidual] 
#Yellow point is a positive u_i
#black point is a negative u_i   


# We will gmap figure out for itself which part
# of the world to take the map from. To do this
# we will give ggmap the extrapolated spatial 
# field w of dimension 100 by 100.

# Take the coordinates of the 100 by 100 field
# and use expand.grid to make a grid of coordinates
# and match the corresponding w.pm100_100 values.
# gmap doesn;t want to have km, hence the * 1000.
# Change that for your own data!
xygrid <- expand.grid(w.proj$x, w.proj$y)
Data3D <- data.frame(x = xygrid[,1],   
                     y = xygrid[,2],
                     z = melt(w.pm100_100)[,3])
names(Data3D) <- c("x", "y", "z")
head(Data3D, 10)

coordinates(Data3D) <- c("x", "y") 
head(Data3D, 10)


# Need some spatial stuff now that we won't explain
r1 <- rasterFromXYZ(Data3D)
projection(r1) <- paste("+init=epsg:29193")

r2<-projectRaster(r1,crs="+proj=longlat +ellps=GRS80 +towgs84=0,0,0 +no_defs")


# Let gmap load the map of South America based on r1
migmap <- gmap(x = r2, type = "terrain", zoom = 4)
plot(migmap)


# Now we need some Mercator magic to convert coordinates
# from one system to another system. The reason for this
# is that we will project the spatial field on top
# of a gmap graph, which requires Mercator coordinates.
# The study area with the sampling locations
xy     <- cbind(fd_env1$Longitude, fd_env1$Latitude)

plot(migmap)
points(Mercator(xy), 
       col = MyCol, 
       pch = c(16,17)[SignResidual], 
       cex = MyCex)

# Note the correlation in the residuals! The spatial field
# picks up a NorthEast-SouthWest pattern for positive and West-central pattern? What may this be?

# And we also want to plot the entire surface on 
# top of the gmap.
# This means that we also have to convert the  
# spatial field to Mercator coordinates. 
rprobGM <- projectRaster(r2,crs= "+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +no_defs")

plot(migmap)
points(Mercator(xy), 
       col = MyCol, 
       pch = c(16,17)[SignResidual], 
       cex = MyCex)
plot(rprobGM, 
     add = T, 
     legend = F, 
     col = rev(cm.colors(5, alpha = 0.6)))

# This is the spatial field converted in to 
# Mercator coordinates
# We only need to do some cropping to remove
# the extrapollated field from the sea.
####################################

# Let's focus on the hyper-parameters.
SpatField.w <- inla.spde2.result(inla = I2,
                                 name = "w",
                                 spde = spde,
                                 do.transfer = TRUE)
names(SpatField.w)
Kappa <- inla.emarginal(function(x) x, 
                        SpatField.w$marginals.kappa[[1]] )
Sigma_u <- inla.emarginal(function(x) sqrt(x), 
                          SpatField.w$marginals.variance.nominal[[1]] )
r <- inla.emarginal(function(x) x, 
                    SpatField.w$marginals.range.nominal[[1]] )
Kappa
Sigma_u
r    #Distance at which the dependency diminishes

# This is perhaps a nicer graph to make and present.
# Show correlation structure
# First we obtain the locations of each point of the mesh.
LocMesh <- mesh1$loc[,1:2]
# And then we calculate the distance between each vertex.
D     <- as.matrix(dist(LocMesh))
# Using the estimated parameters from the model (see above)
# we can calculate the imposed Matern correlation values.
d.vec <- seq(0, max(D), length = 100)      
Cor.M <- (Kappa * d.vec) * besselK(Kappa * d.vec, 1) 
Cor.M[1] <- 1

# Which we plot here:
par(mfrow=c(1,1), mar = c(5,5,2,2))
plot(x = d.vec, 
     y = Cor.M, 
     pch = 16, 
     type = "l", 
     cex.lab = 1.5,
     xlab = "Distance (m)", 
     ylab = "Correlation",
     xlim = c(0, 100000))


##########################################
#PREDICT
##########################################

# For this we need a predict function, but there is
# no predict function here. Remember the two options
# that are available to us. We can either add a vector
# with NAs to MyData and glue that to the iph object,
# or use the lincomp option.
# But because of the stack we need to do this slightly different:
#  1. Use the stack for the original
#  2. Make a stack for the data for which predictions are needed.
#  3. Combine them.
#  4. Run INLA
#  5. Extract the relevant pieces and plot them.


# We already did step 1. Here is step 2.
# Recall that the stack does not like factors.

#########################
# Get climatic data from full study area for the data to predict to
# The original resolution (0.05 or ~5km is too fine for INLA and the full study area, 
# I need ~10k or 0.08 degrees)

#env_var_sa <- env_var[[c(-29,-21,-19,-28,-27,-30,-32,-31,-26,-33,-20,-17,-15,-13,-11,-9,-7,-23,-24,-22,-1,-4)]]#used variables in model

env_var_full   <-projectRaster(env_var1,res=0.5,                                            #higher resolution not computable here
                          crs="+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs")#Layers with 0.8 resolution for prediction

xx1    <-as.data.frame(env_var_full,xy=T); xx1<-na.omit(xx1)#create dataframe with the variables values for full SA

names(xx1)<- c("x","y","MODCF_meanannual","PET","BLDFIE_M_sl6",
               "CECSOL_M_sl6","CLYPPT_M_sl4", "ORCDRC_M_sl4","PHIHOX_M_sl6", "AA","CaH","CG",                                        
               "CrH","FA","GP","HA","HAC","HFA","HFC","HFD","HFE","HG","HGD","HGE","HP","MG","PA","SH","UG",
               "CHIRPS_bio4","CHIRPS_bio8","CHIRPS_bio13","CHIRPS_bio18","CHIRPS_bio19")


# We better standardize the covariates to avoid numerical problems.
xx1$MODCF_meanannual.std  <- MyStd(xx1$MODCF_meanannua)
xx1$PET.std               <- MyStd(xx1$PET)
xx1$BLDFIE_M_sl6.std      <- MyStd(xx1$BLDFIE_M_sl6)
xx1$CECSOL_M_sl6.std      <- MyStd(xx1$CECSOL_M_sl6)
xx1$CLYPPT_M_sl4.std      <- MyStd(xx1$CLYPPT_M_sl4)
xx1$ORCDRC_M_sl4.std      <- MyStd(xx1$ORCDRC_M_sl4)
xx1$PHIHOX_M_sl6.std      <- MyStd(xx1$PHIHOX_M_sl6)
xx1$AA.std      <- MyStd(xx1$AA)
xx1$CaH.std     <- MyStd(xx1$CaH)
xx1$FA.std      <- MyStd(xx1$FA)
xx1$GP.std      <- MyStd(xx1$GP)
xx1$HA.std      <- MyStd(xx1$HA)
xx1$HFC.std     <- MyStd(xx1$HFC)
xx1$HFD.std     <- MyStd(xx1$HFD)
xx1$HFE.std     <- MyStd(xx1$HFE)
xx1$HG.std      <- MyStd(xx1$HG)
xx1$HGD.std     <- MyStd(xx1$HGD)
xx1$HGE.std     <- MyStd(xx1$HGE)
xx1$HP.std      <- MyStd(xx1$HP)
xx1$MG.std      <- MyStd(xx1$MG)
xx1$SH.std      <- MyStd(xx1$SH)
xx1$UG.std      <- MyStd(xx1$UG)
xx1$CHIRPS_bio4.std       <- MyStd(xx1$CHIRPS_bio4)
xx1$CHIRPS_bio8.std       <- MyStd(xx1$CHIRPS_bio8)
xx1$CHIRPS_bio13.std      <- MyStd(xx1$CHIRPS_bio13)
xx1$CHIRPS_bio18.std      <- MyStd(xx1$CHIRPS_bio18)
xx1$CHIRPS_bio19.std      <- MyStd(xx1$CHIRPS_bio19)

MyData2<-xx1; nrow(MyData2)

# And this is the corresponding X matrix
Xmm <- model.matrix(~ MODCF_meanannual.std + PET.std + BLDFIE_M_sl6.std + 
                      CECSOL_M_sl6.std + CLYPPT_M_sl4.std + ORCDRC_M_sl4.std + PHIHOX_M_sl6.std + 
                      AA.std + CaH.std + FA.std + GP.std + HFC.std + HA.std + HFD.std + HFE.std + HG.std + HGD.std + HGE.std + 
                      HP.std + MG.std + SH.std + UG.std + CHIRPS_bio4.std +     
                      CHIRPS_bio8.std + CHIRPS_bio13.std + CHIRPS_bio18.std + CHIRPS_bio19.std, 
                      data = MyData2)
head(Xmm)

Xp <- data.frame(MODCF_meanannual.std    = Xmm[,2],
                 PET.std                 = Xmm[,3],
                 BLDFIE_M_sl4.std        = Xmm[,4],
                 CECSOL_M_sl4.std        = Xmm[,5],
                 CLYPPT_M_sl4.std        = Xmm[,6],
                 ORCDRC_M_sl4.std        = Xmm[,7],
                 PHIHOX_M_sl4.std        = Xmm[,8],
                 AA.std                  = Xmm[,9],
                 CaH.std                 = Xmm[,10],
                   FA.std                = Xmm[,11],
                   GP.std                = Xmm[,12],
                   HFC.std               = Xmm[,13],
                   HA.std                = Xmm[,14],
                   HFD.std               = Xmm[,15],
                   HFE.std               = Xmm[,16],
                   HG.std                = Xmm[,17],
                   HGD.std               = Xmm[,18],
                   HGE.std               = Xmm[,19],
                   HP.std                = Xmm[,20],
                   MG.std                = Xmm[,21],
                   SH.std                = Xmm[,22],
                   UG.std                = Xmm[,23],
                   CHIRPS_bio4.std       = Xmm[,24],  
                   CHIRPS_bio8.std       = Xmm[,25],
                   CHIRPS_bio13.std      = Xmm[,26],
                   CHIRPS_bio18.std      = Xmm[,27],
                   CHIRPS_bio19.std      = Xmm[,28]
)

#stack of full study area
StackCov <- inla.stack(
  tag = "Covariates",
  data = list(y = NA),  
  A = list(1, 1),                  
  effects = list(
    Intercept = rep(1, nrow(Xp)),
    Xp = Xp))

All.stacks <- inla.stack(StackCov,Stack2)	              

# And run the model with the combined stack
f2 <- y ~ -1 + Intercept +
  MODCF_meanannual.std + PET.std + BLDFIE_M_sl6.std + 
  CECSOL_M_sl6.std + CLYPPT_M_sl4.std + ORCDRC_M_sl4.std + PHIHOX_M_sl6.std + 
  AA.std + CaH.std + FA.std + GP.std + HFC.std + HA.std + HFD.std + HFE.std + HG.std + HGD.std + HGE.std + 
  HP.std + MG.std + SH.std + UG.std + CHIRPS_bio4.std +     
  CHIRPS_bio8.std + CHIRPS_bio13.std + CHIRPS_bio18.std + CHIRPS_bio19.std +
  f(w, model = spde)

# 7. Run the spatial model in INLA.
I2 <- inla(f2,
           family = "gaussian", 
           data = inla.stack.data(All.stacks),
           control.compute = list(dic = TRUE, waic = TRUE),
           control.predictor = list(A = inla.stack.A(All.stacks)),
           control.inla = list(strategy = "gaussian")) 

# This is the crucial code for extracting the correct rows.
# It provides an index for which rows in the combined stack
# belong to the observed data and which rows to the 
# artificial  covariate data. 
index.Fit <- inla.stack.index(All.stacks,
                              tag = "Fit")$data

index.Cov <- inla.stack.index(All.stacks,
                              tag = "Covariates")$data


#And we can extact the correct rows     
F7.fit  <- I2$summary.fitted.values[index.Fit, c(1,3,5)]  #1654  by 3
F7.pred <- I2$summary.fitted.values[index.Cov, c(1,3,5)]  #38907 by 3

########################################################
# Plot map predictions                                 #
########################################################
#Bind to full data and plot the spatially explicit model

library(RColorBrewer);library(rasterVis)
sa<-shapefile('SouthAmericaStudyArea_JA.shp')

MyData3 <- cbind(MyData2, F7.pred)
dim(MyData3)
colnames(MyData3)

cw_wd_pred_ras      <- rasterFromXYZ(MyData3)
cw_wd_pred_ras_mean <- cw_wd_pred_ras[[60]];   writeRaster(cw_wd_pred_ras_mean,"D:/AMAZON/ATDN/ResultsFD/Plots 20 June 2017/cw_SLA_pred_ras_mean1.img",overwrite=T)
cw_wd_pred_ras_025q <- cw_wd_pred_ras[[61]]
cw_wd_pred_ras_0975 <- cw_wd_pred_ras[[62]]

colrxx <- colorRampPalette((brewer.pal(n=5, 'BrBG')))#other option #colrxx <- colorRampPalette(rev(brewer.pal(n=5, 'YlGnBu')))

levelplot(cw_wd_pred_ras_mean, xlim=c(-85,-30),ylim=c(-25,15),col.regions=colrxx,
          margin=FALSE, main="CWM SLA 20/06/17")+layer(sp.polygons(sa, lwd=1, col="grey"))+ 
  #layer(sp.points(locat, col = MyCol2,pch = c(16,17)[SignResidual],cex = MyCex))+
  layer(sp.polygons(forestborder))

levelplot(cw_wd_pred_ras_025q, xlim=c(-85,-30),ylim=c(-25,15),col.regions=colrxx,
          margin=FALSE, main="CWM SLA 0.025q 20/06/17")+layer(sp.polygons(sa, lwd=1, col="grey"))+ 
  #layer(sp.points(locat, col = MyCol2,pch = c(16,17)[SignResidual],cex = MyCex))+
  layer(sp.polygons(forestborder))

levelplot(cw_wd_pred_ras_0975, xlim=c(-85,-30),ylim=c(-25,15),col.regions=colrxx,
          margin=FALSE, main="CWM SLA 0.0975q 13/06/17")+layer(sp.polygons(sa, lwd=1, col="grey"))+ 
  #layer(sp.points(locat, col = MyCol2,pch = c(16,17)[SignResidual],cex = MyCex))+
  layer(sp.polygons(forestborder))

#OR
colr <- colorRampPalette(brewer.pal(11, 'RdYlBu'))
MyCol2 <- c("black", "yellow")[SignResidual] #Plot prediction with residuals and see pattern catched by the randomfield
locat  <- SpatialPoints(cbind(fd_env$Longitude, fd_env$Latitude), proj4string = CRS("+proj=longlat"))


levelplot(cw_wd_pred_ras, xlim=c(-100,-30),col.regions=colr,
          margin=FALSE, main="CWM SLA 14/06/17")+layer(sp.polygons(sa, lwd=1, col="grey"))+ 
  layer(sp.points(locat, col = MyCol2,pch = c(16,17)[SignResidual],cex = MyCex))     


#OR
summary(mask_ras)
breakpoints <- c(12,12.57483,12.83683,13,14,15.18755)
colors <- c("dark blue","blue","green","yellow","brown")

levelplot(mask_ras, xlim=c(-85,-30),ylim=c(-25,15),at=breakpoints,col.regions=colors,#col.regions=colrxx,
          margin=FALSE, main="CWM SLA 14/06/17")+layer(sp.polygons(sa, lwd=1, col="grey"))+
  layer(sp.polygons(forestborder))

#OR
levelplot(mask_ras, xlim=c(-85,-30),ylim=c(-25,15),at=breakpoints,col.regions=colrxx,pretty=T,
          margin=FALSE, main="CWM SLA 14/06/17")+layer(sp.polygons(sa, lwd=1, col="grey"))+
  layer(sp.polygons(forestborder))
