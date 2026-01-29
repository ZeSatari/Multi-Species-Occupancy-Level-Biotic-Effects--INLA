
library(timechange) 
library(lubridate) 
library(INLA) 
library(tidyverse)
library(spOccupancy)
library(scico)
library(patchwork)
library(kableExtra)
library(sf)
library(terra)
library(tidyterra)

plot_inla_effects = function(effect)
{
  p1 = ggplot(data.frame(effect)) +

    geom_line(aes(ID, -mean)) + 

    geom_ribbon(aes(ID, ymin = -X0.025quant, ymax = -X0.975quant),

                alpha = 0.5)
  print(p1)
  
}


theme_maps = theme(axis.line=element_blank(),
                   axis.text.x=element_blank(),
                   axis.text.y=element_blank(),
                   axis.ticks=element_blank(),
                   axis.title.x=element_blank(),
                   axis.title.y=element_blank()
                   #legend.position="none",
                   #panel.background=element_blank(),
                   #panel.border=element_blank(),
                   #panel.grid.major=element_blank(),
                   #panel.grid.minor=element_blank(),
                   #plot.background=element_blank()
)

#### LOAD AND PREPARE THE DATA
data(hbefTrends)

revi.data <- hbefTrends
sp.names <- dimnames(hbefTrends$y)[[1]] 


revi.data$y <- revi.data$y[sp.names == 'REVI', , , ]
revi.data$coords = revi.data$coords/1000 #get this in km!


###===============Detection-based biotic information

data <- hbefTrends

OVEN<- data$y[sp.names == 'OVEN', , , ]
OVEN.Y = data.frame(OVEN[,1,])

for(i in 2:9)
{
OVEN.Y  = rbind(OVEN.Y , data.frame(OVEN[,i,]))
 }
head(OVEN.Y)
dim(OVEN.Y)
OVEN.y=OVEN.Y[,1]
head(OVEN.y)


Y = data.frame(revi.data$y[,1,])
Xdet.day = revi.data$det.covs$day[,1,]
Xdet.tod = revi.data$det.covs$tod[,1,]

for(i in 2:9)
{
  Y = rbind(Y, data.frame(revi.data$y[,i,]))
  Xdet.day = rbind(Xdet.day, revi.data$det.covs$day[,i,])
  Xdet.tod = rbind(Xdet.tod, revi.data$det.covs$tod[,i,])
}
Xdet.day = (Xdet.day - mean(Xdet.day, na.rm = T))/sd(Xdet.day,na.rm = T) 
Xdet.tod = (Xdet.tod - mean(Xdet.tod, na.rm = T))/sd(Xdet.tod,na.rm = T) 
 head(Xdet.day)
Xdet.O = cbind(
  1, Xdet.day[,1], Xdet.tod[,1], OVEN.Y[,1],
  1, Xdet.day[,2], Xdet.tod[,2], OVEN.Y[,2],
  1, Xdet.day[,3], Xdet.tod[,3], OVEN.Y[,3]
)

head(Xdet.O )


Xocc = data.frame(x = rep(revi.data$coords[,1],9),
                  y = rep(revi.data$coords[,2],9),
                  elev = rep(revi.data$occ.covs$elev,9),
                  scale_elev = scale(rep(revi.data$occ.covs$elev,9)),
                 scale_elev2 = scale(rep(revi.data$occ.covs$elev,9))^2,
                   site = rep(1:373,9),
                  time = rep(1:9,each = 373),
                  Int_occ = 1,
                 spatialfield = NA) %>%
  mutate(scale_time = scale(time))



elev_raster= rast(data.frame(x =  hbefElev$Easting/1000,
                                    y = hbefElev$Northing/1000,
                                    elev =  hbefElev$val))

elev_raster$scale_elev = (elev_raster$elev - mean(revi.data$occ.covs$elev)) / sd(revi.data$occ.covs$elev)
elev_raster$scale_elev2 = elev_raster$scale_elev^2 

all_elev = inla.group(c(Xocc$elev, values(elev_raster$elev)))
val_elev = sort(na.omit(unique(all_elev)))
Xocc$group_elev = all_elev[1:dim(Xocc)[1]]
elev_group_raster = elev_raster$elev
values(elev_group_raster) <- all_elev[-c(1:dim(Xocc)[1])]


pred_df = data.frame(x = rep(crds(elev_raster, na.rm = F)[,1],9),
                     y = rep(crds(elev_raster, na.rm = F)[,2],9),
                     scale_elev = rep(values(elev_raster$scale_elev),9),
                     scale_elev2 = rep(values(elev_raster$scale_elev2),9),
                     group_elev = rep(values(elev_group_raster),9),
                     time = rep(c(1:9), each= length(values(elev_raster$elev)))
                     ) %>%
  dplyr::filter(!is.na(scale_elev))




data_list = as.list(Xocc)
data_list$Y = Y
data_list$X =  Xdet.O


#### MODEL FIT

boundary = inla.nonconvex.hull(points = revi.data$coords, convex = .3)
mesh = inla.mesh.2d(boundary = boundary,
                    #   loc = cbind(data$X, data$Y),
                    max.edge = c(0.1,0.7),
                    min.angle = 20,
                    offset = c(.01, 1),
                    cutoff = 0.12,
)


spde <- inla.spde2.pcmatern(
  mesh = mesh, 
  prior.range = c(5, 0.01), # prior for range
  prior.sigma = c(1, 0.5))  # prior for sd parameter


 
h.spec <- list(rho = list(prior = 'pc.cor0', param = c(0.5, 0.1)))
spde <- inla.spde2.pcmatern(
  mesh = mesh, 
  prior.range = c(5, 0.7),
  prior.sigma = c(1, 0.5),
  constr = T) 

A_sp2 <- inla.spde.make.A(mesh = mesh, group = Xocc$time,
                         loc = cbind(Xocc$x, Xocc$y))

formula_2 <- inla.mdata(Y,X) ~ 
  -1   + Int_occ +              
  scale_elev +  scale_elev2 +    
  f(time, model = "ar1") +       
  f(spatialfield,                
    model=spde,                  
    A.local = A_sp2,            
    group = time,               
    control.group = list(model = 'ar1',hyper = h.spec))

model_2 <- inla(formula_2, #the formula
                                     data=data_list,  #the data stack
                                     family= 'occupancy',   #which family the data comes from
                                     control.fixed =  list(prec = 1, prec.intercept = 1),
                                     control.compute = list(waic = TRUE,config = TRUE),
                                   verbose = F,
                                     control.family = list(control.link = list(model = "logit"), 
                                                                                    link.simple = "logit",
                                                           hyper = list(beta1 = list(param = c(0,10), 
                                                                                   initial = 0),
                                                                      beta2 = list(param = c(0,10)),
                                                                      beta3 = list(param = c(0,10))
                                                           )))



summary(model_2)

loocv_m2_auto <- inla.group.cv(result = model_2, num.level.sets = 3)
UlogCV_sim <- mean(log(loocv_detect_auto$cv), na.rm=TRUE)

#####============Baseline detection model

revi.data <- hbefTrends
sp.names <- dimnames(hbefTrends$y)[[1]] 
revi.data$y <- revi.data$y[sp.names == 'REVI', , , ]
revi.data$coords = revi.data$coords/1000 #get this in km!


Y = data.frame(revi.data$y[,1,])
Xdet.day = revi.data$det.covs$day[,1,]
Xdet.tod = revi.data$det.covs$tod[,1,]

for(i in 2:9)
{
  Y = rbind(Y, data.frame(revi.data$y[,i,]))
  Xdet.day = rbind(Xdet.day, revi.data$det.covs$day[,i,])
  Xdet.tod = rbind(Xdet.tod, revi.data$det.covs$tod[,i,])
}
Xdet.day = (Xdet.day - mean(Xdet.day, na.rm = T))/sd(Xdet.day,na.rm = T) 
Xdet.tod = (Xdet.tod - mean(Xdet.tod, na.rm = T))/sd(Xdet.tod,na.rm = T) 
 head(Xdet.day)
Xdet= cbind(1, Xdet.day[,1], Xdet.tod[,1],
             1, Xdet.day[,2], Xdet.tod[,2],
             1, Xdet.day[,3], Xdet.tod[,3])


Xocc = data.frame(x = rep(revi.data$coords[,1],9),
                  y = rep(revi.data$coords[,2],9),
                  elev = rep(revi.data$occ.covs$elev,9),
                  scale_elev = scale(rep(revi.data$occ.covs$elev,9)),
                 scale_elev2 = scale(rep(revi.data$occ.covs$elev,9))^2,
                   site = rep(1:373,9),
                  time = rep(1:9,each = 373),
                  Int_occ = 1,
                 spatialfield = NA) %>%
  mutate(scale_time = scale(time))



data_list = as.list(Xocc)
data_list$Y = Y
data_list$X = Xdet



#####
h.spec <- list(rho = list(prior = 'pc.cor0', param = c(0.5, 0.1)))
spde <- inla.spde2.pcmatern(
  mesh = mesh, 
  prior.range = c(5, 0.7),
  prior.sigma = c(1, 0.5),
  constr = T) 

A_sp2 <- inla.spde.make.A(mesh = mesh, group = Xocc$time,
                         loc = cbind(Xocc$x, Xocc$y))

formula <- inla.mdata(Y,X) ~ 
  -1   + Int_occ +               
  scale_elev +  scale_elev2 +    
  f(time, model = "ar1") +       
  f(spatialfield,              
    model=spde,                 
    A.local = A_sp2,           
    group = time,               
    control.group = list(model = 'ar1',hyper = h.spec))

model <- inla(formula, #the formula
                                     data=data_list,  #the data stack
                                     family= 'occupancy',   #which family the data comes from
                                     control.fixed =  list(prec = 1, prec.intercept = 1),
                                     control.compute = list(waic = TRUE,config = TRUE),
                                     verbose = F,
                                     control.family = list(control.link = list(model = "logit"), 
                                                                                    link.simple = "logit",
                                                           hyper = list(beta1 = list(param = c(0,10), 
                                                                                   initial = 0),
                                                                      beta2 = list(param = c(0,10)),
                                                                      beta3 = list(param = c(0,10))
                                                           )))



summary(model)

loocv_m1_auto <- inla.group.cv(result =  model,num.level.sets = 3)
ULOOCV_auto  = mean(log(loocv_m1_auto$cv),na.rm=T)



