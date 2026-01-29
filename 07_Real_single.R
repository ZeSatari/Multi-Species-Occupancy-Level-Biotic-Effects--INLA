

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
    geom_ribbon(aes(ID, ymin = -X0.025quant, ymax = -X0.975quant),  alpha = 0.5)
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
#
data(hbefTrends)

data <- hbefTrends
sp.names <- dimnames(hbefTrends$y)[[1]]

####================AMRE

data$y <- data$y[sp.names == 'AMRE', , , ]
data$coords = data$coords/1000 #get this in km!


Y = data.frame(data$y[,1,])

###=====With identical visit-specific detection covariates(Xdet)

Xdet.day = data$det.covs$day[,1,]
Xdet.tod = data$det.covs$tod[,1,]

for(i in 2:9)
{
  Y = rbind(Y, data.frame(data$y[,i,]))
  Xdet.day = rbind(Xdet.day, data$det.covs$day[,i,])
  Xdet.tod = rbind(Xdet.tod, data$det.covs$tod[,i,])
}
Xdet.day = (Xdet.day - mean(Xdet.day, na.rm = T))/sd(Xdet.day,na.rm = T) 
Xdet.tod = (Xdet.tod - mean(Xdet.tod, na.rm = T))/sd(Xdet.tod,na.rm = T) 
 head(Xdet.day)
Xdet= cbind(1, Xdet.day[,1], Xdet.tod[,1],
             1, Xdet.day[,2], Xdet.tod[,2],
             1, Xdet.day[,3], Xdet.tod[,3])


Xocc = data.frame(x = rep(data$coords[,1],9),
                  y = rep(data$coords[,2],9),
                  elev = rep(data$occ.covs$elev,9),
                  scale_elev = scale(rep(data$occ.covs$elev,9)),
                 scale_elev2 = scale(rep(data$occ.covs$elev,9))^2,
                   site = rep(1:373,9),
                  time = rep(1:9,each = 373),
                  Int_occ = 1,
                 spatialfield = NA) %>%
  mutate(scale_time = scale(time))



elev_raster= rast(data.frame(x =  hbefElev$Easting/1000,
                                    y = hbefElev$Northing/1000,
                                    elev =  hbefElev$val))

elev_raster$scale_elev = (elev_raster$elev - mean(data$occ.covs$elev)) / sd(data$occ.covs$elev)
elev_raster$scale_elev2 = elev_raster$scale_elev^2 

all_elev = inla.group(c(Xocc$elev, values(elev_raster$elev)))

val_elev = sort(na.omit(unique(all_elev)))
Xocc$group_elev = all_elev[1:dim(Xocc)[1]]
elev_group_raster = elev_raster$elev
values(elev_group_raster) <- all_elev[-c(1:dim(Xocc)[1])]



data_list = as.list(Xocc)
data_list$Y = Y
data_list$X = Xdet




boundary = inla.nonconvex.hull(points = data$coords, convex = .3)
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

formula1 <- inla.mdata(Y,X) ~ 
  -1   + Int_occ +              
  scale_elev +  scale_elev2 +   
  f(time, model = "ar1") +        
  f(spatialfield,                
    model=spde,                  
    A.local = A_sp2,             
    group = time,                
    control.group = list(model = 'ar1',hyper = h.spec))

model1 <- inla(formula1, #the formula
                                     data=data_list,  #the data stack
                                     family= 'occupancy',   #which family the data comes from
                                     control.fixed =  list(prec = 1, prec.intercept = 1),
                                     control.compute = list(waic = TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
                                     verbose = F,
                                     control.family = list(control.link = list(model = "logit"), 
                                                                                    link.simple = "logit",
                                                           hyper = list(beta1 = list(param = c(0,10), 
                                                                                   initial = 0),
                                                                      beta2 = list(param = c(0,10)),
                                                                      beta3 = list(param = c(0,10))
                                                           )))



summary(model1)

loocv_m1_auto <- inla.group.cv(result =  model1,num.level.sets = 3)
ULOOCV_auto  = mean(log(loocv_m1_auto$cv),na.rm=T)

######============================BAWW

data <- hbefTrends

data$y <- data$y[sp.names == "BAWW", , , ]
data$coords = data$coords/1000 #get this in km!


Y = data.frame(data$y[,1,])



Xocc = data.frame(x = rep(data$coords[,1],9),
                  y = rep(data$coords[,2],9),
                  elev = rep(data$occ.covs$elev,9),
                  scale_elev = scale(rep(data$occ.covs$elev,9)),
                 scale_elev2 = scale(rep(data$occ.covs$elev,9))^2,
                   site = rep(1:373,9),
                  time = rep(1:9,each = 373),
                  Int_occ = 1,
                 spatialfield = NA) %>%
  mutate(scale_time = scale(time))


elev_raster= rast(data.frame(x =  hbefElev$Easting/1000,
                                    y = hbefElev$Northing/1000,
                                    elev =  hbefElev$val))

elev_raster$scale_elev = (elev_raster$elev - mean(data$occ.covs$elev)) / sd(data$occ.covs$elev)
elev_raster$scale_elev2 = elev_raster$scale_elev^2 

all_elev = inla.group(c(Xocc$elev, values(elev_raster$elev)))

val_elev = sort(na.omit(unique(all_elev)))
Xocc$group_elev = all_elev[1:dim(Xocc)[1]]
elev_group_raster = elev_raster$elev
values(elev_group_raster) <- all_elev[-c(1:dim(Xocc)[1])]

data_list = as.list(Xocc)
data_list$Y = Y
data_list$X = Xdet


boundary = inla.nonconvex.hull(points = data$coords, convex = .3)
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



#####
h.spec <- list(rho = list(prior = 'pc.cor0', param = c(0.5, 0.1)))
spde <- inla.spde2.pcmatern(
  mesh = mesh, 
  prior.range = c(5, 0.7),
  prior.sigma = c(1, 0.5),
  constr = T) 

A_sp2 <- inla.spde.make.A(mesh = mesh, group = Xocc$time,
                         loc = cbind(Xocc$x, Xocc$y))

formula2 <- inla.mdata(Y,X) ~ 
  -1   + Int_occ +              
  scale_elev +  scale_elev2 +    
  f(time, model = "ar1") +       
  f(spatialfield,               
    model=spde,                  
    A.local = A_sp2,            
    group = time,               
    control.group = list(model = 'ar1',hyper = h.spec))

model2 <- inla(formula2, #the formula
                                     data=data_list,  #the data stack
                                     family= 'occupancy',   #which family the data comes from
                                     control.fixed =  list(prec = 1, prec.intercept = 1),
                                     control.compute = list(waic = TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
                                     verbose = F,
                                     control.family = list(control.link = list(model = "logit"), 
                                                                                    link.simple = "logit",
                                                           hyper = list(beta1 = list(param = c(0,10), 
                                                                                   initial = 0),
                                                                      beta2 = list(param = c(0,10)),
                                                                      beta3 = list(param = c(0,10))
                                                           )))



summary(model2)


loocv_m1_auto <- inla.group.cv(result =  model2,num.level.sets = 3)
ULOOCV_auto  = mean(log(loocv_m1_auto$cv),na.rm=T)



########=================BHVI


data <- hbefTrends

data$y <- data$y[sp.names == "BHVI", , , ]
data$coords = data$coords/1000 #get this in km!


Y = data.frame(data$y[,1,])

Xocc = data.frame(x = rep(data$coords[,1],9),
                  y = rep(data$coords[,2],9),
                  elev = rep(data$occ.covs$elev,9),
                  scale_elev = scale(rep(data$occ.covs$elev,9)),
                 scale_elev2 = scale(rep(data$occ.covs$elev,9))^2,
                   site = rep(1:373,9),
                  time = rep(1:9,each = 373),
                  Int_occ = 1,
                 spatialfield = NA) %>%
  mutate(scale_time = scale(time))


elev_raster= rast(data.frame(x =  hbefElev$Easting/1000,
                                    y = hbefElev$Northing/1000,
                                    elev =  hbefElev$val))

elev_raster$scale_elev = (elev_raster$elev - mean(data$occ.covs$elev)) / sd(data$occ.covs$elev)
elev_raster$scale_elev2 = elev_raster$scale_elev^2 

all_elev = inla.group(c(Xocc$elev, values(elev_raster$elev)))

val_elev = sort(na.omit(unique(all_elev)))
Xocc$group_elev = all_elev[1:dim(Xocc)[1]]
elev_group_raster = elev_raster$elev
values(elev_group_raster) <- all_elev[-c(1:dim(Xocc)[1])]

data_list = as.list(Xocc)
data_list$Y = Y
data_list$X = Xdet


boundary = inla.nonconvex.hull(points = data$coords, convex = .3)
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

#####
h.spec <- list(rho = list(prior = 'pc.cor0', param = c(0.5, 0.1)))
spde <- inla.spde2.pcmatern(
  mesh = mesh, 
  prior.range = c(5, 0.7),
  prior.sigma = c(1, 0.5),
  constr = T) 

A_sp2 <- inla.spde.make.A(mesh = mesh, group = Xocc$time,#در مدل 3 گروهبندي زمان را اضافه مي کند
                         loc = cbind(Xocc$x, Xocc$y))

formula3<- inla.mdata(Y,X) ~ 
  -1   + Int_occ +               #-1: بدون عامل ثابت ----Int_occ: عامل کلی Int_occ
  scale_elev +  scale_elev2 +    #  متغیر ارتفاع
  f(time, model = "ar1") +       #اثر زمانی AR1 
  f(spatialfield,                #میدان فضایی
    model=spde,                  #استفاده از SPDE
    A.local = A_sp2,             # ماتریس جدید SPDE
    group = time,                #گروه بندي برا اساس زمان
    control.group = list(model = 'ar1',hyper = h.spec))

model3 <- inla(formula3, #the formula
                                     data=data_list,  #the data stack
                                     family= 'occupancy',   #which family the data comes from
                                     control.fixed =  list(prec = 1, prec.intercept = 1),
                                     control.compute = list(waic = TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
                                     verbose = F,
                                     control.family = list(control.link = list(model = "logit"), 
                                                                                    link.simple = "logit",
                                                           hyper = list(beta1 = list(param = c(0,10), 
                                                                                   initial = 0),
                                                                      beta2 = list(param = c(0,10)),
                                                                      beta3 = list(param = c(0,10))
                                                           )))



summary(model3)
loocv_m1_auto <- inla.group.cv(result =  model3,num.level.sets = 3)
ULOOCV_auto  = mean(log(loocv_m1_auto$cv),na.rm=T)

####================BLBW

data <- hbefTrends

data$y <- data$y[sp.names == "BLBW", , , ]
data$coords = data$coords/1000 #get this in km!


Y = data.frame(data$y[,1,])



Xocc = data.frame(x = rep(data$coords[,1],9),
                  y = rep(data$coords[,2],9),
                  elev = rep(data$occ.covs$elev,9),
                  scale_elev = scale(rep(data$occ.covs$elev,9)),
                 scale_elev2 = scale(rep(data$occ.covs$elev,9))^2,
                   site = rep(1:373,9),
                  time = rep(1:9,each = 373),
                  Int_occ = 1,
                 spatialfield = NA) %>%
  mutate(scale_time = scale(time))


elev_raster= rast(data.frame(x =  hbefElev$Easting/1000,
                                    y = hbefElev$Northing/1000,
                                    elev =  hbefElev$val))

elev_raster$scale_elev = (elev_raster$elev - mean(data$occ.covs$elev)) / sd(data$occ.covs$elev)
elev_raster$scale_elev2 = elev_raster$scale_elev^2 

all_elev = inla.group(c(Xocc$elev, values(elev_raster$elev)))

val_elev = sort(na.omit(unique(all_elev)))
Xocc$group_elev = all_elev[1:dim(Xocc)[1]]
elev_group_raster = elev_raster$elev
values(elev_group_raster) <- all_elev[-c(1:dim(Xocc)[1])]




data_list = as.list(Xocc)
data_list$Y = Y
data_list$X = Xdet



boundary = inla.nonconvex.hull(points = data$coords, convex = .3)
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

#####
h.spec <- list(rho = list(prior = 'pc.cor0', param = c(0.5, 0.1)))
spde <- inla.spde2.pcmatern(
  mesh = mesh, 
  prior.range = c(5, 0.7),
  prior.sigma = c(1, 0.5),
  constr = T) 

A_sp2 <- inla.spde.make.A(mesh = mesh, group = Xocc$time,
                         loc = cbind(Xocc$x, Xocc$y))

formula4<- inla.mdata(Y,X) ~ 
  -1   + Int_occ +              
  scale_elev +  scale_elev2 +    
  f(time, model = "ar1") +     
  f(spatialfield,                
    model=spde,                  
    A.local = A_sp2,           
    group = time,                
    control.group = list(model = 'ar1',hyper = h.spec))

model4 <- inla(formula4, #the formula
                                     data=data_list,  #the data stack
                                     family= 'occupancy',   #which family the data comes from
                                     control.fixed =  list(prec = 1, prec.intercept = 1),
                                     control.compute = list(waic = TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
                                     verbose = F,
                                     control.family = list(control.link = list(model = "logit"), 
                                                                                    link.simple = "logit",
                                                           hyper = list(beta1 = list(param = c(0,10), 
                                                                                   initial = 0),
                                                                      beta2 = list(param = c(0,10)),
                                                                      beta3 = list(param = c(0,10))
                                                           )))



summary(model4)

loocv_m1_auto <- inla.group.cv(result =  model4,num.level.sets = 3)
ULOOCV_auto  = mean(log(loocv_m1_auto$cv),na.rm=T)


####===============BLPW

data <- hbefTrends

data$y <- data$y[sp.names == "BLPW", , , ]
data$coords = data$coords/1000 #get this in km!


Y = data.frame(data$y[,1,])



Xocc = data.frame(x = rep(data$coords[,1],9),
                  y = rep(data$coords[,2],9),
                  elev = rep(data$occ.covs$elev,9),
                  scale_elev = scale(rep(data$occ.covs$elev,9)),
                 scale_elev2 = scale(rep(data$occ.covs$elev,9))^2,
                   site = rep(1:373,9),
                  time = rep(1:9,each = 373),
                  Int_occ = 1,
                 spatialfield = NA) %>%
  mutate(scale_time = scale(time))


elev_raster= rast(data.frame(x =  hbefElev$Easting/1000,
                                    y = hbefElev$Northing/1000,
                                    elev =  hbefElev$val))

elev_raster$scale_elev = (elev_raster$elev - mean(data$occ.covs$elev)) / sd(data$occ.covs$elev)
elev_raster$scale_elev2 = elev_raster$scale_elev^2 

all_elev = inla.group(c(Xocc$elev, values(elev_raster$elev)))

val_elev = sort(na.omit(unique(all_elev)))
Xocc$group_elev = all_elev[1:dim(Xocc)[1]]
elev_group_raster = elev_raster$elev
values(elev_group_raster) <- all_elev[-c(1:dim(Xocc)[1])]


data_list = as.list(Xocc)
data_list$Y = Y
data_list$X = Xdet



boundary = inla.nonconvex.hull(points = data$coords, convex = .3)
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


#####
h.spec <- list(rho = list(prior = 'pc.cor0', param = c(0.5, 0.1)))
spde <- inla.spde2.pcmatern(
  mesh = mesh, 
  prior.range = c(5, 0.7),
  prior.sigma = c(1, 0.5),
  constr = T) 

A_sp2 <- inla.spde.make.A(mesh = mesh, group = Xocc$time,
                         loc = cbind(Xocc$x, Xocc$y))

formula5<- inla.mdata(Y,X) ~ 
  -1   + Int_occ +             
  scale_elev +  scale_elev2 +   
  f(time, model = "ar1") +      
  f(spatialfield,                
    model=spde,                  
    A.local = A_sp2,             
    group = time,               
    control.group = list(model = 'ar1',hyper = h.spec))

model5 <- inla(formula5, #the formula
                                     data=data_list,  #the data stack
                                     family= 'occupancy',   #which family the data comes from
                                     control.fixed =  list(prec = 1, prec.intercept = 1),
                                     control.compute = list(waic = TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
                                     verbose = F,
                                     control.family = list(control.link = list(model = "logit"), 
                                                                                    link.simple = "logit",
                                                           hyper = list(beta1 = list(param = c(0,10), 
                                                                                   initial = 0),
                                                                      beta2 = list(param = c(0,10)),
                                                                      beta3 = list(param = c(0,10))
                                                           )))



summary(model5)

loocv_m1_auto <- inla.group.cv(result =  model5,num.level.sets = 3)
ULOOCV_auto  = mean(log(loocv_m1_auto$cv),na.rm=T)


####==================BTBW

data <- hbefTrends

data$y <- data$y[sp.names == "BTBW", , , ]
data$coords = data$coords/1000 #get this in km!


Y = data.frame(data$y[,1,])


Xocc = data.frame(x = rep(data$coords[,1],9),
                  y = rep(data$coords[,2],9),
                  elev = rep(data$occ.covs$elev,9),
                  scale_elev = scale(rep(data$occ.covs$elev,9)),
                 scale_elev2 = scale(rep(data$occ.covs$elev,9))^2,
                   site = rep(1:373,9),
                  time = rep(1:9,each = 373),
                  Int_occ = 1,
                 spatialfield = NA) %>%
  mutate(scale_time = scale(time))


elev_raster= rast(data.frame(x =  hbefElev$Easting/1000,
                                    y = hbefElev$Northing/1000,
                                    elev =  hbefElev$val))

elev_raster$scale_elev = (elev_raster$elev - mean(data$occ.covs$elev)) / sd(data$occ.covs$elev)
elev_raster$scale_elev2 = elev_raster$scale_elev^2 

all_elev = inla.group(c(Xocc$elev, values(elev_raster$elev)))

val_elev = sort(na.omit(unique(all_elev)))
Xocc$group_elev = all_elev[1:dim(Xocc)[1]]
elev_group_raster = elev_raster$elev
values(elev_group_raster) <- all_elev[-c(1:dim(Xocc)[1])]

data_list = as.list(Xocc)
data_list$Y = Y
data_list$X = Xdet


boundary = inla.nonconvex.hull(points = data$coords, convex = .3)
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


 
#####
h.spec <- list(rho = list(prior = 'pc.cor0', param = c(0.5, 0.1)))
spde <- inla.spde2.pcmatern(
  mesh = mesh, 
  prior.range = c(5, 0.7),
  prior.sigma = c(1, 0.5),
  constr = T) 

A_sp2 <- inla.spde.make.A(mesh = mesh, group = Xocc$time,
                         loc = cbind(Xocc$x, Xocc$y))

formula3_3 <- inla.mdata(Y,X) ~ 
  -1   + Int_occ +               
  scale_elev +  scale_elev2 +    
  f(time, model = "ar1") +       
  f(spatialfield,               
    model=spde,                 
    A.local = A_sp2,             
    group = time,               
    control.group = list(model = 'ar1',hyper = h.spec))

model6 <- inla(formula6, #the formula
                                     data=data_list,  #the data stack
                                     family= 'occupancy',   #which family the data comes from
                                     control.fixed =  list(prec = 1, prec.intercept = 1),
                                     control.compute = list(waic = TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
                                     verbose = F,
                                     control.family = list(control.link = list(model = "logit"), 
                                                                                    link.simple = "logit",
                                                           hyper = list(beta1 = list(param = c(0,10), 
                                                                                   initial = 0),
                                                                      beta2 = list(param = c(0,10)),
                                                                      beta3 = list(param = c(0,10))
                                                           )))



summary(model6)

loocv_m1_auto <- inla.group.cv(result =  model6,num.level.sets = 3)
ULOOCV_auto  = mean(log(loocv_m1_auto$cv),na.rm=T)

####===================BTNW

data <- hbefTrends
data$y <- data$y[sp.names == "BTNW", , , ]
data$coords = data$coords/1000 #get this in km!


Y = data.frame(data$y[,1,])



Xocc = data.frame(x = rep(data$coords[,1],9),
                  y = rep(data$coords[,2],9),
                  elev = rep(data$occ.covs$elev,9),
                  scale_elev = scale(rep(data$occ.covs$elev,9)),
                 scale_elev2 = scale(rep(data$occ.covs$elev,9))^2,
                   site = rep(1:373,9),
                  time = rep(1:9,each = 373),
                  Int_occ = 1,
                 spatialfield = NA) %>%
  mutate(scale_time = scale(time))


elev_raster= rast(data.frame(x =  hbefElev$Easting/1000,
                                    y = hbefElev$Northing/1000,
                                    elev =  hbefElev$val))

elev_raster$scale_elev = (elev_raster$elev - mean(data$occ.covs$elev)) / sd(data$occ.covs$elev)
elev_raster$scale_elev2 = elev_raster$scale_elev^2 

all_elev = inla.group(c(Xocc$elev, values(elev_raster$elev)))

val_elev = sort(na.omit(unique(all_elev)))
Xocc$group_elev = all_elev[1:dim(Xocc)[1]]
elev_group_raster = elev_raster$elev
values(elev_group_raster) <- all_elev[-c(1:dim(Xocc)[1])]

data_list = as.list(Xocc)
data_list$Y = Y
data_list$X = Xdet




boundary = inla.nonconvex.hull(points = data$coords, convex = .3)
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

#####
h.spec <- list(rho = list(prior = 'pc.cor0', param = c(0.5, 0.1)))
spde <- inla.spde2.pcmatern(
  mesh = mesh, 
  prior.range = c(5, 0.7),
  prior.sigma = c(1, 0.5),
  constr = T) 

A_sp2 <- inla.spde.make.A(mesh = mesh, group = Xocc$time,
                         loc = cbind(Xocc$x, Xocc$y))

formula7 <- inla.mdata(Y,X) ~ 
  -1   + Int_occ +              
  scale_elev +  scale_elev2 +    
  f(time, model = "ar1") +      
  f(spatialfield,                
    model=spde,                 
    A.local = A_sp2,            
    group = time,                
    control.group = list(model = 'ar1',hyper = h.spec))

model7 <- inla(formula7, #the formula
                                     data=data_list,  #the data stack
                                     family= 'occupancy',   #which family the data comes from
                                     control.fixed =  list(prec = 1, prec.intercept = 1),
                                     control.compute = list(waic = TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
                                     verbose = F,
                                     control.family = list(control.link = list(model = "logit"), 
                                                                                    link.simple = "logit",
                                                           hyper = list(beta1 = list(param = c(0,10), 
                                                                                   initial = 0),
                                                                      beta2 = list(param = c(0,10)),
                                                                      beta3 = list(param = c(0,10))
                                                           )))



summary(model7)

loocv_m1_auto <- inla.group.cv(result =  model7,num.level.sets = 3)
ULOOCV_auto  = mean(log(loocv_m1_auto$cv),na.rm=T)


####CAWA


data <- hbefTrends
data$coords = data$coords/1000 #get this in km!


Y = data.frame(data$y[,1,])


Xocc = data.frame(x = rep(data$coords[,1],9),
                  y = rep(data$coords[,2],9),
                  elev = rep(data$occ.covs$elev,9),
                  scale_elev = scale(rep(data$occ.covs$elev,9)),
                 scale_elev2 = scale(rep(data$occ.covs$elev,9))^2,
                   site = rep(1:373,9),
                  time = rep(1:9,each = 373),
                  Int_occ = 1,
                 spatialfield = NA) %>%
  mutate(scale_time = scale(time))


elev_raster= rast(data.frame(x =  hbefElev$Easting/1000,
                                    y = hbefElev$Northing/1000,
                                    elev =  hbefElev$val))

elev_raster$scale_elev = (elev_raster$elev - mean(data$occ.covs$elev)) / sd(data$occ.covs$elev)
elev_raster$scale_elev2 = elev_raster$scale_elev^2 

all_elev = inla.group(c(Xocc$elev, values(elev_raster$elev)))

val_elev = sort(na.omit(unique(all_elev)))
Xocc$group_elev = all_elev[1:dim(Xocc)[1]]
elev_group_raster = elev_raster$elev
values(elev_group_raster) <- all_elev[-c(1:dim(Xocc)[1])]


data_list = as.list(Xocc)
data_list$Y = Y
data_list$X = Xdet



boundary = inla.nonconvex.hull(points = data$coords, convex = .3)
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


 
#####
h.spec <- list(rho = list(prior = 'pc.cor0', param = c(0.5, 0.1)))
spde <- inla.spde2.pcmatern(
  mesh = mesh, 
  prior.range = c(5, 0.7),
  prior.sigma = c(1, 0.5),
  constr = T) 

A_sp2 <- inla.spde.make.A(mesh = mesh, group = Xocc$time,
                         loc = cbind(Xocc$x, Xocc$y))

formula8 <- inla.mdata(Y,X) ~ 
  -1   + Int_occ +               
  scale_elev +  scale_elev2 +    
  f(time, model = "ar1") +       
  f(spatialfield,                
    model=spde,                  
    A.local = A_sp2,             
    group = time,                
    control.group = list(model = 'ar1',hyper = h.spec))

model8 <- inla(formula8, #the formula
                                     data=data_list,  #the data stack
                                     family= 'occupancy',   #which family the data comes from
                                     control.fixed =  list(prec = 1, prec.intercept = 1),
                                     control.compute = list(waic = TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
                                     verbose = F,
                                     control.family = list(control.link = list(model = "logit"), 
                                                                                    link.simple = "logit",
                                                           hyper = list(beta1 = list(param = c(0,10), 
                                                                                   initial = 0),
                                                                      beta2 = list(param = c(0,10)),
                                                                      beta3 = list(param = c(0,10))
                                                           )))



summary(model8)

loocv_m1_auto <- inla.group.cv(result =  model8,num.level.sets = 3)
ULOOCV_auto  = mean(log(loocv_m1_auto$cv),na.rm=T)


####============MAWA

data <- hbefTrends

data$y <- data$y[sp.names == "MAWA", , , ]
data$coords = data$coords/1000 #get this in km!


Y = data.frame(data$y[,1,])


Xocc = data.frame(x = rep(data$coords[,1],9),
                  y = rep(data$coords[,2],9),
                  elev = rep(data$occ.covs$elev,9),
                  scale_elev = scale(rep(data$occ.covs$elev,9)),
                 scale_elev2 = scale(rep(data$occ.covs$elev,9))^2,
                   site = rep(1:373,9),
                  time = rep(1:9,each = 373),
                  Int_occ = 1,
                 spatialfield = NA) %>%
  mutate(scale_time = scale(time))


elev_raster= rast(data.frame(x =  hbefElev$Easting/1000,
                                    y = hbefElev$Northing/1000,
                                    elev =  hbefElev$val))

elev_raster$scale_elev = (elev_raster$elev - mean(data$occ.covs$elev)) / sd(data$occ.covs$elev)
elev_raster$scale_elev2 = elev_raster$scale_elev^2 

all_elev = inla.group(c(Xocc$elev, values(elev_raster$elev)))

val_elev = sort(na.omit(unique(all_elev)))
Xocc$group_elev = all_elev[1:dim(Xocc)[1]]
elev_group_raster = elev_raster$elev
values(elev_group_raster) <- all_elev[-c(1:dim(Xocc)[1])]


data_list = as.list(Xocc)
data_list$Y = Y
data_list$X = Xdet



boundary = inla.nonconvex.hull(points = data$coords, convex = .3)
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


 
#####
h.spec <- list(rho = list(prior = 'pc.cor0', param = c(0.5, 0.1)))
spde <- inla.spde2.pcmatern(
  mesh = mesh, 
  prior.range = c(5, 0.7),
  prior.sigma = c(1, 0.5),
  constr = T) 

A_sp2 <- inla.spde.make.A(mesh = mesh, group = Xocc$time,
                         loc = cbind(Xocc$x, Xocc$y))

formula9 <- inla.mdata(Y,X) ~ 
  -1   + Int_occ +              
  scale_elev +  scale_elev2 +    
  f(time, model = "ar1") +      
  f(spatialfield,               
    model=spde,                  
    A.local = A_sp2,             
    group = time,                
    control.group = list(model = 'ar1',hyper = h.spec))

model9 <- inla(formula9, #the formula
                                     data=data_list,  #the data stack
                                     family= 'occupancy',   #which family the data comes from
                                     control.fixed =  list(prec = 1, prec.intercept = 1),
                                     control.compute = list(waic = TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
                                     verbose = F,
                                     control.family = list(control.link = list(model = "logit"), 
                                                                                    link.simple = "logit",
                                                           hyper = list(beta1 = list(param = c(0,10), 
                                                                                   initial = 0),
                                                                      beta2 = list(param = c(0,10)),
                                                                      beta3 = list(param = c(0,10))
                                                           )))



summary(model9)


loocv_m1_auto <- inla.group.cv(result =  model9,num.level.sets = 3)
ULOOCV_auto  = mean(log(loocv_m1_auto$cv),na.rm=T)


####============NAWA

data <- hbefTrends

data$y <- data$y[sp.names == "NAWA", , , ]
data$coords = data$coords/1000 #get this in km!


Y = data.frame(data$y[,1,])


Xocc = data.frame(x = rep(data$coords[,1],9),
                  y = rep(data$coords[,2],9),
                  elev = rep(data$occ.covs$elev,9),
                  scale_elev = scale(rep(data$occ.covs$elev,9)),
                 scale_elev2 = scale(rep(data$occ.covs$elev,9))^2,
                   site = rep(1:373,9),
                  time = rep(1:9,each = 373),
                  Int_occ = 1,
                 spatialfield = NA) %>%
  mutate(scale_time = scale(time))


elev_raster= rast(data.frame(x =  hbefElev$Easting/1000,
                                    y = hbefElev$Northing/1000,
                                    elev =  hbefElev$val))

elev_raster$scale_elev = (elev_raster$elev - mean(data$occ.covs$elev)) / sd(data$occ.covs$elev)
elev_raster$scale_elev2 = elev_raster$scale_elev^2 

all_elev = inla.group(c(Xocc$elev, values(elev_raster$elev)))

val_elev = sort(na.omit(unique(all_elev)))
Xocc$group_elev = all_elev[1:dim(Xocc)[1]]
elev_group_raster = elev_raster$elev
values(elev_group_raster) <- all_elev[-c(1:dim(Xocc)[1])]



data_list = as.list(Xocc)
data_list$Y = Y
data_list$X = Xdet



boundary = inla.nonconvex.hull(points = data$coords, convex = .3)
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


 
#####
h.spec <- list(rho = list(prior = 'pc.cor0', param = c(0.5, 0.1)))
spde <- inla.spde2.pcmatern(
  mesh = mesh, 
  prior.range = c(5, 0.7),
  prior.sigma = c(1, 0.5),
  constr = T) 

A_sp2 <- inla.spde.make.A(mesh = mesh, group = Xocc$time,
                         loc = cbind(Xocc$x, Xocc$y))

formula10 <- inla.mdata(Y,X) ~ 
  -1   + Int_occ +               
  scale_elev +  scale_elev2 +    
  f(time, model = "ar1") +       
  f(spatialfield,               
    model=spde,                  
    A.local = A_sp2,            
    group = time,                
    control.group = list(model = 'ar1',hyper = h.spec))

model10 <- inla(formula10, #the formula
                                     data=data_list,  #the data stack
                                     family= 'occupancy',   #which family the data comes from
                                     control.fixed =  list(prec = 1, prec.intercept = 1),
                                     control.compute = list(waic = TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
                                     verbose = F,
                                     control.family = list(control.link = list(model = "logit"), 
                                                                                    link.simple = "logit",
                                                           hyper = list(beta1 = list(param = c(0,10), 
                                                                                   initial = 0),
                                                                      beta2 = list(param = c(0,10)),
                                                                      beta3 = list(param = c(0,10))
                                                           ))))



summary(model10)

loocv_m1_auto <- inla.group.cv(result =  model10,num.level.sets = 3)
ULOOCV_auto  = mean(log(loocv_m1_auto$cv),na.rm=T)


####=================OVEN


data <- hbefTrends

data$y <- data$y[sp.names == "OVEN", , , ]
data$coords = data$coords/1000 #get this in km!


Xocc = data.frame(x = rep(data$coords[,1],9),
                  y = rep(data$coords[,2],9),
                  elev = rep(data$occ.covs$elev,9),
                  scale_elev = scale(rep(data$occ.covs$elev,9)),
                 scale_elev2 = scale(rep(data$occ.covs$elev,9))^2,
                   site = rep(1:373,9),
                  time = rep(1:9,each = 373),
                  Int_occ = 1,
                 spatialfield = NA) %>%
  mutate(scale_time = scale(time))


elev_raster= rast(data.frame(x =  hbefElev$Easting/1000,
                                    y = hbefElev$Northing/1000,
                                    elev =  hbefElev$val))

elev_raster$scale_elev = (elev_raster$elev - mean(data$occ.covs$elev)) / sd(data$occ.covs$elev)
elev_raster$scale_elev2 = elev_raster$scale_elev^2 

all_elev = inla.group(c(Xocc$elev, values(elev_raster$elev)))

val_elev = sort(na.omit(unique(all_elev)))
Xocc$group_elev = all_elev[1:dim(Xocc)[1]]
elev_group_raster = elev_raster$elev
values(elev_group_raster) <- all_elev[-c(1:dim(Xocc)[1])]


data_list = as.list(Xocc)
data_list$Y = Y
data_list$X = Xdet



boundary = inla.nonconvex.hull(points = data$coords, convex = .3)
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

 
#####
h.spec <- list(rho = list(prior = 'pc.cor0', param = c(0.5, 0.1)))
spde <- inla.spde2.pcmatern(
  mesh = mesh, 
  prior.range = c(5, 0.7),
  prior.sigma = c(1, 0.5),
  constr = T) 

A_sp2 <- inla.spde.make.A(mesh = mesh, group = Xocc$time,
                         loc = cbind(Xocc$x, Xocc$y))

formula11 <- inla.mdata(Y,X) ~ 
  -1   + Int_occ +               
  scale_elev +  scale_elev2 +   
  f(time, model = "ar1") +   
  f(spatialfield,               
    model=spde,                 
    A.local = A_sp2,           
    group = time,               
    control.group = list(model = 'ar1',hyper = h.spec))

model11 <- inla(formula11, #the formula
                                     data=data_list,  #the data stack
                                     family= 'occupancy',   #which family the data comes from
                                     control.fixed =  list(prec = 1, prec.intercept = 1),
                                     control.compute = list(waic = TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
                                     verbose = F,
                                     control.family = list(control.link = list(model = "logit"), 
                                                                                    link.simple = "logit",
                                                           hyper = list(beta1 = list(param = c(0,10), 
                                                                                   initial = 0),
                                                                      beta2 = list(param = c(0,10)),
                                                                      beta3 = list(param = c(0,10))
                                                           )))



summary(model11)

loocv_m1_auto <- inla.group.cv(result =  model11,num.level.sets = 3)
ULOOCV_auto  = mean(log(loocv_m1_auto$cv),na.rm=T)


#####============REVI

data <- hbefTrends
data$y <- data$y[sp.names == 'REVI', , , ]
data$coords = data$coords/1000 #get this in km!


Y = data.frame(data$y[,1,])


Xocc = data.frame(x = rep(data$coords[,1],9),
                  y = rep(data$coords[,2],9),
                  elev = rep(data$occ.covs$elev,9),
                  scale_elev = scale(rep(data$occ.covs$elev,9)),
                 scale_elev2 = scale(rep(data$occ.covs$elev,9))^2,
                   site = rep(1:373,9),
                  time = rep(1:9,each = 373),
                  Int_occ = 1,
                 spatialfield = NA) %>%
  mutate(scale_time = scale(time))


elev_raster= rast(data.frame(x =  hbefElev$Easting/1000,
                                    y = hbefElev$Northing/1000,
                                    elev =  hbefElev$val))

elev_raster$scale_elev = (elev_raster$elev - mean(data$occ.covs$elev)) / sd(data$occ.covs$elev)
elev_raster$scale_elev2 = elev_raster$scale_elev^2 

all_elev = inla.group(c(Xocc$elev, values(elev_raster$elev)))

val_elev = sort(na.omit(unique(all_elev)))
Xocc$group_elev = all_elev[1:dim(Xocc)[1]]
elev_group_raster = elev_raster$elev
values(elev_group_raster) <- all_elev[-c(1:dim(Xocc)[1])]


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

formula12 <- inla.mdata(Y,X) ~ 
  -1   + Int_occ +               
  scale_elev +  scale_elev2 +    
  f(time, model = "ar1") +       
  f(spatialfield,              
    model=spde,                 
    A.local = A_sp2,           
    group = time,               
    control.group = list(model = 'ar1',hyper = h.spec))

model12 <- inla(formula12, #the formula
                                     data=data_list,  #the data stack
                                     family= 'occupancy',   #which family the data comes from
                                     control.fixed =  list(prec = 1, prec.intercept = 1),
                                     control.compute = list(waic = TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
                                     verbose = F,
                                     control.family = list(control.link = list(model = "logit"), 
                                                                                    link.simple = "logit",
                                                           hyper = list(beta1 = list(param = c(0,10), 
                                                                                   initial = 0),
                                                                      beta2 = list(param = c(0,10)),
                                                                      beta3 = list(param = c(0,10))
                                                           )))



summary(model12)

loocv_m1_auto <- inla.group.cv(result =  model12,num.level.sets = 3)
ULOOCV_auto  = mean(log(loocv_m1_auto$cv),na.rm=T)





