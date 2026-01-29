
 library(INLA)
library(inlabru)
library(fmesher)
library(tidyverse)
library(sf)
library(terra)
library(dplyr)

inla.Occupancy_detCov <- function(X_det){

  if(class(X_det)=="list"){
    if(length(X_det)>10){
      warning("exceeded number of detection covariates, numerical issues may occur")
    }

    if(lapply(X_det, ncol)%>%unlist()%>%unique()%>%length()>2){
      stop("inconsistent number of visits in provided detection covariates")
    }
    if(length(lapply(X_det, nrow) %>% unlist() %>% unique())>1){
      stop("inconsistent number of sites in provided detection covariates")
    }
    K<- lapply(X_det, ncol) %>% unlist() %>% max() # Max num of visits
    M<- lapply(X_det, nrow) %>% unlist() %>% unique() # Number of sites
    P <- length(X_det)

    if(lapply(X_det, ncol)%>%unlist()%>%unique()%>%length()==2 & 
       1 %in% lapply(X_det, ncol)%>%unlist()%>%unique()){
      warning(paste("At least one covariate of dimension [",M,",1] has been provided, values for this covariate will be repeated over the max numver of visits",sep=""))
      for(l in which(lapply(X_det, ncol) %>% unlist() < K)){
        X_det[[l]] <- do.call("cbind",replicate(K,X_det[[l]]))

      }
    }
    covariates <- do.call("cbind", lapply(1:K, function(i) {
      do.call("cbind", lapply(X_det, function(mat) mat[, i]))
    }))

  }

  if(is.data.frame(X_det)|is.matrix(X_det)){
    K<- ncol(X_det)
    M<- nrow(X_det)
    P <- 1
    covariates <- as.matrix(X_det)
  }

  X_mat <- matrix(NA,nrow=M,ncol=K*(P+1))
  X_mat[,seq(1,(K*(P+1)),by=(P+1))]<-1 # add Intercept at the begining of each visit-specific covariate matrix
  X_mat[, which(!(1:(K*(P+1)) %in% seq(1,(K*(P+1)),by=(P+1))))] <- covariates
  return(X_mat)

}


# read the data
space_time_data <- read.csv("Occ_data_multi.csv")
Occ_data_multi

# raster data
x_covariate <- terra::rast('xs.tif')

 


space_time_data_sf <- space_time_data %>%
  st_as_sf(coords = c('x.loc','y.loc')) %>%
  mutate(site = as.factor(cellid),
         species = as.factor(species), 
         terra::extract(x_covariate, st_coordinates(.)))


levels(space_time_data_sf$site) <- 1:nlevels(space_time_data_sf$site) 

spat_data  = space_time_data_sf %>% 
  dplyr::select(-cellid) %>%
  mutate(x.loc= st_coordinates(space_time_data_sf)[,1],
         y.loc= st_coordinates(space_time_data_sf)[,2]) %>% st_drop_geometry()



nT <- length(spat_data$time %>% unique) # number of time points
ncells <- length(values(x_covariate$x_s)) # number of cells in the whole area

# scale and centre the covariate value in the data input

spat_data <- spat_data %>%
  mutate(scale_x_s = scale(x_s)  %>% c())



mesh = fm_mesh_2d(loc.domain = st_coordinates(boundary_sf)[,1:2],
                    offset = c(-0.1, -.2),
                    max.edge = c(15, 30))

spde <- inla.spde2.pcmatern(mesh = mesh,
                              prior.range = c(100, 0.5),
                              prior.sigma = c(0.5, 0.5))
A_sp <- inla.spde.make.A(mesh = mesh,
                         loc = as.matrix(spat_data[,c("x.loc","y.loc")]),
                         group = spat_data$time)

# Detection /non-detection data
Y_mat <- spat_data %>%
  dplyr::select(num_range("y",1:3)) %>% 
  as.matrix()
# Detection covariates matrix
X_det <- spat_data %>%
  dplyr::select(num_range("g",1:3)) %>% 
  inla.Occupancy_detCov()

X_cov <- spat_data %>%
  dplyr::select(c(time,site,scale_x_s )) %>%
  mutate(Int_occ = 1,
         spatialfield = rep(NA,nrow(spat_data)))

head(X_cov)


spat_data$species_id <- as.numeric(factor(spat_data$species))



# Append to a list:
data_list_msom <- as.list(X_cov) 
data_list_msom$Y <- Y_mat
data_list_msom$X <- X_det
data_list_msom$species_id <- spat_data$species_id

data_list_msom$species_slope<-spat_data$species_id
data_list_msom$scale_x_s2<- spat_data$scale_x_s^2


# PC prior for the temporal autocorrelation
h.spec <- list(rho = list(prior = 'pc.cor0', param = c(0.5, 0.3)))


Q_species <- matrix(1, n_species, n_species)
diag(Q_species) <- 2 



 formula_species <- inla.mdata(Y, X) ~ -1 + Int_occ +
     f(species_id, model = "generic0", Cmatrix = Q_species, constr=TRUE) +
     f(species_slope, I(scale_x_s2), model="iid", constr=TRUE)+
     f(spatialfield, model = spde,
       A.local = A_sp,
       group = time,
       control.group = list(model = 'ar1', hyper = h.spec))
 
   model_msom <- inla(formula_species ,
                      data = data_list_14,
                      family = 'occupancy',
                      control.fixed = list(prec = 1/2.72, prec.intercept = 1/2.72),
                      control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
                      control.inla = list(int.strategy = "eb"),
                      verbose = FALSE,
                      control.family = list(  
                        control.link = list(model = "logit"),
                        link.simple = "logit"
                      ))
  summary(model_msom )


loocv_multi_auto<- inla.group.cv(result =  model_msom,num.level.sets = 3)
ULOOCV_auto_multi  = mean(log(loocv_multi_auto$cv),na.rm=T)


table = data.frame(  DIC = c(model_msom$dic$dic),
                    WAIC_obs = c(model_msom$waic$waic/120000),
                    mlik = c(model_msom$mlik[1,1]),
                    ULOOCV=ULOOCV_auto_multi
                    )
table 



        
# fixed
fixed_int <- model_msom$summary.fixed["Int_occ","mean"]

##====plot gamma

gamma_true  <- gamma_species
summary_gamma <- model_msom$summary.random$species_id

gamma_est   <- summary_gamma$mean-0.5
gamma_lwr   <- summary_gamma$`0.025quant`-0.5
gamma_upr   <- summary_gamma$`0.975quant`-0.5

species_id <- 1:length(gamma_true)

y_lim <- range(c(gamma_true, gamma_lwr, gamma_upr))

pdf("gama.pdf", width=6, height=6) 

plot(species_id, gamma_est,
     type="n", ylim=y_lim,
     xlab="Species", ylab="Gamma (species effect)") 

polygon(c(species_id, rev(species_id)),
        c(gamma_lwr, rev(gamma_upr)),
        col=rgb(0,0,1,0.2), border=NA)  

lines(species_id, gamma_est, col="blue", lwd=2)
points(species_id, gamma_est, col="blue", pch=16)
lines(species_id, gamma_true, col="red", lwd=2, lty=2)
points(species_id, gamma_true, col="red", pch=17)

legend("topleft",
       legend=c("Estimated gamma (CI)", "True gamma"),
       col=c("blue","red"), lwd=2, pch=c(16,17),
       lty=c(1,2))

dev.off()

####=============== Single-Species

 beta_est<-c()

 #============================= species 1
 spat_data1= spat_data[spat_data$species == 1, ]
  

 A_sp1 <- inla.spde.make.A(mesh = mesh,
                          loc = as.matrix(spat_data1[,c("x.loc","y.loc")]),
                          group = spat_data1$time)
# Detection /non-detection data
 Y_mat <- spat_data1 %>%
   dplyr::select(num_range("y",1:3)) %>% 
   as.matrix()


 # Detection covariates matrix
 X_det <- spat_data1 %>%
   dplyr::select(num_range("g",1:3)) %>% 
   inla.Occupancy_detCov()
 
 
 X_cov <- spat_data1 %>%
  dplyr::select(c(time,scale_x_s)) %>%
   mutate(Int_occ = 1,
          spatialfield = rep(NA,nrow(spat_data1)))
 

 # Append to a list:
 data_list1 <- as.list(X_cov) 
 data_list1$Y <- Y_mat
 data_list1$X <- X_det

#=======
 formula_z1 <- inla.mdata(Y, X) ~ -1 + Int_occ +
  scale_x_s^2  +
  f(spatialfield, model = spde,
     A.local = A_sp1,
     group = time,
     control.group = list(model = 'ar1', hyper = h.spec))

model_single1 <- inla(formula_z1,
                    data = data_list1,
                    family = 'occupancy',
                    control.fixed = list(prec = 1/2.72, prec.intercept = 1/2.72),
                    control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
                    control.inla = list(int.strategy = "eb"),
                    verbose = FALSE,
                    control.family = list(
                      control.link = list(model = "logit"),
                      link.simple = "logit"
                    ))

summary(model_single1)


 loocv_m1_auto_z1 <- inla.group.cv(result =  model_single1,num.level.sets = 3)
  ULOOCV_auto_z1  = mean(log(loocv_m1_auto_z1$cv),na.rm=T)




table = data.frame(  DIC = c(model_single1$dic$dic),
                    WAIC = c(model_single1$waic$waic/10000),
                    mlik = c(model_single1$mlik[1,1]),
                    LGOCV_auto = c(ULOOCV_auto_z1)
)
table 




  beta_est[1]  <- model_single1$summary.fixed["scale_x_s","mean"]


 #============================= species 2
 spat_data2= spat_data[spat_data$species == 2, ]
  

 A_sp2 <- inla.spde.make.A(mesh = mesh,
                          loc = as.matrix(spat_data2[,c("x.loc","y.loc")]),
                          group = spat_data2$time)
# Detection /non-detection data
 Y_mat <- spat_data2 %>%
   dplyr::select(num_range("y",1:3)) %>% 
   as.matrix()


 # Detection covariates matrix
 X_det <- spat_data2 %>%
   dplyr::select(num_range("g",1:3)) %>% 
   inla.Occupancy_detCov()
 
 
 X_cov <- spat_data2 %>%
  dplyr::select(c(time,scale_x_s)) %>%
   mutate(Int_occ = 1,
          spatialfield = rep(NA,nrow(spat_data2)))
 

 
 # Append to a list:
 data_list2 <- as.list(X_cov) 
 data_list2$Y <- Y_mat
 data_list2$X <- X_det


#=======
 formula_z2 <- inla.mdata(Y, X) ~ -1 + Int_occ +
  scale_x_s^2  +
  f(spatialfield, model = spde,
     A.local = A_sp2,
     group = time,
     control.group = list(model = 'ar1', hyper = h.spec))

model_single2 <- inla(formula_z2,
                    data = data_list2,
                    family = 'occupancy',
                    control.fixed = list(prec = 1/2.72, prec.intercept = 1/2.72),
                    control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
                    control.inla = list(int.strategy = "eb"),
                    verbose = FALSE,
                    control.family = list(
                      control.link = list(model = "logit"),
                      link.simple = "logit"
                    ))

summary(model_single2)

 loocv_m1_auto_z2 <- inla.group.cv(result =  model_single2,num.level.sets = 3)
  ULOOCV_auto_z2  = mean(log(loocv_m1_auto_z2$cv),na.rm=T)


table = data.frame(  DIC = c(model_single2$dic$dic),
                    WAIC = c(model_single2$waic$waic/10000),
                    mlik = c(model_single2$mlik[1,1]),
                    LGOCV_auto = c(ULOOCV_auto_z2)
)
table 


  beta_est[2]  <- model_single2$summary.fixed["scale_x_s","mean"]

################ species 3
 
 spat_data3= spat_data[spat_data$species ==3,] 
 

 A_sp3 <- inla.spde.make.A(mesh = mesh,
                          loc = as.matrix(spat_data3[,c("x.loc","y.loc")]),
                          group = spat_data3$time)
 # Detection /non-detection data
 Y_mat <- spat_data3 %>%
   dplyr::select(num_range("y",1:3)) %>% 
   as.matrix()
 # Detection covariates matrix
 X_det <- spat_data3 %>%
   dplyr::select(num_range("g",1:3)) %>% 
   inla.Occupancy_detCov()
 
 
 X_cov <- spat_data3 %>%
  dplyr::select(c(time,scale_x_s)) %>%
   mutate(Int_occ = 1,
          spatialfield = rep(NA,nrow(spat_data3)))



 # Append to a list:
 data_list3 <- as.list(X_cov) 
 data_list3$Y <- Y_mat
 data_list3$X <- X_det

formula_z3 <- inla.mdata(Y, X) ~ -1 + Int_occ +
   scale_x_s^2 +
  f(spatialfield, model = spde,
     A.local = A_sp3,
     group = time,
     control.group = list(model = 'ar1', hyper = h.spec))



 model_single3<- inla(formula_z3,
                    data = data_list3,
                    family = 'occupancy',
                   control.fixed = list(prec = 1/2.72, prec.intercept = 1/2.72),
                    control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
                    control.inla = list(int.strategy = "eb"),
                    verbose = FALSE,
                   control.family = list(
                     control.link = list(model = "logit"),
                     link.simple = "logit"
                   ))

 summary(model_single3)
 loocv_m1_auto_z3<- inla.group.cv(result =  model_single3,num.level.sets = 3)
  ULOOCV_auto_z3  = mean(log(loocv_m1_auto_z3$cv),na.rm=T)


table = data.frame(  DIC = c(model_single3$dic$dic),
                     WAIC = c(model_single3$waic$waic/10000),
                     mlik = c(model_single3$mlik[1,1]),
                     LGOCV_auto = c(ULOOCV_auto_z3)
 )
table 



  beta_est[3]  <- model_single3$summary.fixed["scale_x_s","mean"]

################# species 4
 
 spat_data4= spat_data[spat_data$species ==4,] 
 
 A_sp4 <- inla.spde.make.A(mesh = mesh,
                          loc = as.matrix(spat_data4[,c("x.loc","y.loc")]),
                          group = spat_data4$time)
 # Detection /non-detection data
 Y_mat <- spat_data4 %>%
   dplyr::select(num_range("y",1:3)) %>% 
   as.matrix()
 # Detection covariates matrix
 X_det <- spat_data4 %>%
   dplyr::select(num_range("g",1:3)) %>% 
   inla.Occupancy_detCov()
 
 
 X_cov <- spat_data4 %>%
  dplyr::select(c(time,scale_x_s)) %>%
   mutate(Int_occ = 1,
          spatialfield = rep(NA,nrow(spat_data4)))

 # Append to a list:
 data_list4 <- as.list(X_cov) 
 data_list4$Y <- Y_mat
 data_list4$X <- X_det

 
formula_z4 <- inla.mdata(Y, X) ~ -1 + Int_occ +
   scale_x_s^2 +
  f(spatialfield, model = spde,
     A.local = A_sp4,
     group = time,
     control.group = list(model = 'ar1', hyper = h.spec))



 model_single4<- inla(formula_z4,
                    data = data_list4,
                    family = 'occupancy',
                   control.fixed = list(prec = 1/2.72, prec.intercept = 1/2.72),
                    control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
                    control.inla = list(int.strategy = "eb"),
                    verbose = FALSE,
                   control.family = list(
                     control.link = list(model = "logit"),
                     link.simple = "logit"
                   ))

 summary(model_single4)

 loocv_m1_auto_z4<- inla.group.cv(result =  model_single4,num.level.sets = 3)
  ULOOCV_auto_z4  = mean(log(loocv_m1_auto_z4$cv),na.rm=T)

table = data.frame(  DIC = c(model_single4$dic$dic),
                     WAIC = c(model_single4$waic$waic/10000),
                     mlik = c(model_single4$mlik[1,1]),
                     LGOCV_auto = c(ULOOCV_auto_z4)
 )
table 



  beta_est[4]  <- model_single4$summary.fixed["scale_x_s","mean"]

################# species 5


  spat_data5 <-spat_data[spat_data$species ==5, ]

 A_sp5 <- inla.spde.make.A(mesh = mesh,
                          loc = as.matrix(spat_data5[,c("x.loc","y.loc")]),
                          group = spat_data5$time)
 # Detection /non-detection data
 Y_mat <- spat_data5 %>%
   dplyr::select(num_range("y",1:3)) %>% 
   as.matrix()
 # Detection covariates matrix
 X_det <- spat_data5 %>%
   dplyr::select(num_range("g",1:3)) %>% 
   inla.Occupancy_detCov()
 
 
 X_cov <- spat_data5 %>%
  dplyr::select(c(time,scale_x_s)) %>%
   mutate(Int_occ = 1,
          spatialfield = rep(NA,nrow(spat_data5)))


 # Append to a list:
 data_list5 <- as.list(X_cov) 
 data_list5$Y <- Y_mat
 data_list5$X <- X_det

formula_z5 <- inla.mdata(Y, X) ~ -1 + Int_occ +
   scale_x_s^2 +
  f(spatialfield, model = spde,
     A.local = A_sp5,
     group = time,
     control.group = list(model = 'ar1', hyper = h.spec))



 model_single5<- inla(formula_z5,
                    data = data_list5,
                    family = 'occupancy',
                   control.fixed = list(prec = 1/2.72, prec.intercept = 1/2.72),
                    control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
                    control.inla = list(int.strategy = "eb"),
                    verbose = FALSE,
                   control.family = list(
                     control.link = list(model = "logit"),
                     link.simple = "logit"
                   ))

 summary(model_single5)
 loocv_m1_auto_z5<- inla.group.cv(result =  model_single5,num.level.sets = 3)
  ULOOCV_auto_z5  = mean(log(loocv_m1_auto_z5$cv),na.rm=T)


table = data.frame(  DIC = c(model_single5$dic$dic),
                     WAIC = c(model_single5$waic$waic/10000),
                     mlik = c(model_single5$mlik[1,1]),
                     LGOCV_auto = c(ULOOCV_auto_z5)
 )
table 

  beta_est[5]  <- model_single5$summary.fixed["scale_x_s","mean"]

################# species 6

 spat_data6 <-spat_data[spat_data$species ==6, ]


 A_sp6 <- inla.spde.make.A(mesh = mesh,
                          loc = as.matrix(spat_data6[,c("x.loc","y.loc")]),
                          group = spat_data6$time)

 # Detection /non-detection data
 Y_mat <- spat_data6 %>%
   dplyr::select(num_range("y",1:3)) %>% 
   as.matrix()
 # Detection covariates matrix
 X_det <- spat_data6 %>%
   dplyr::select(num_range("g",1:3)) %>% 
   inla.Occupancy_detCov()
 
 
 X_cov <- spat_data6 %>%
  dplyr::select(c(time,scale_x_s)) %>%
   mutate(Int_occ = 1,
          spatialfield = rep(NA,nrow(spat_data6)))


 # Append to a list:
 data_list6 <- as.list(X_cov) 
 data_list6$Y <- Y_mat
 data_list6$X <- X_det

 
formula_z6 <- inla.mdata(Y, X) ~ -1 + Int_occ +
   scale_x_s^2 +
  f(spatialfield, model = spde,
     A.local = A_sp6,
     group = time,
     control.group = list(model = 'ar1', hyper = h.spec))



 model_single6<- inla(formula_z6,
                    data = data_list6,
                    family = 'occupancy',
                   control.fixed = list(prec = 1/2.72, prec.intercept = 1/2.72),
                    control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
                    control.inla = list(int.strategy = "eb"),
                    verbose = FALSE,
                   control.family = list(
                     control.link = list(model = "logit"),
                     link.simple = "logit"
                   ))

 summary(model_single6)

 loocv_m1_auto_z6<- inla.group.cv(result =  model_single6,num.level.sets = 3)
  ULOOCV_auto_z6  = mean(log(loocv_m1_auto_z6$cv),na.rm=T)



table = data.frame(  DIC = c(model_single6$dic$dic),
                     WAIC = c(model_single6$waic$waic/10000),
                     mlik = c(model_single6$mlik[1,1]),
                     LGOCV_auto = c(ULOOCV_auto_z6)
 )
table 



  beta_est[6]  <- model_single6$summary.fixed["scale_x_s","mean"]



################# species 7

 spat_data7 <-spat_data[spat_data$species ==7, ]

  

 A_sp7 <- inla.spde.make.A(mesh = mesh,
                          loc = as.matrix(spat_data7[,c("x.loc","y.loc")]),
                          group = spat_data7$time)
 # Detection /non-detection data
 Y_mat <- spat_data7 %>%
   dplyr::select(num_range("y",1:3)) %>% 
   as.matrix()
 # Detection covariates matrix
 X_det <- spat_data7 %>%
   dplyr::select(num_range("g",1:3)) %>% 
   inla.Occupancy_detCov()
 
 
 X_cov <- spat_data7 %>%
  dplyr::select(c(time,scale_x_s)) %>%
   mutate(Int_occ = 1,
          spatialfield = rep(NA,nrow(spat_data7)))


 # Append to a list:
 data_list7 <- as.list(X_cov) 
 data_list7$Y <- Y_mat
 data_list7$X <- X_det

formula_z7 <- inla.mdata(Y, X) ~ -1 + Int_occ +
   scale_x_s^2 +
  f(spatialfield, model = spde,
     A.local = A_sp7,
     group = time,
     control.group = list(model = 'ar1', hyper = h.spec))



 model_single7<- inla(formula_z7,
                    data = data_list7,
                    family = 'occupancy',
                   control.fixed = list(prec = 1/2.72, prec.intercept = 1/2.72),
                    control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
                    control.inla = list(int.strategy = "eb"),
                    verbose = FALSE,
                   control.family = list(
                     control.link = list(model = "logit"),
                     link.simple = "logit"
                   ))

 summary(model_single7)

 loocv_m1_auto_z7<- inla.group.cv(result =  model_single7,num.level.sets = 3)
  ULOOCV_auto_z7  = mean(log(loocv_m1_auto_z7$cv),na.rm=T)


table = data.frame(  DIC = c(model_single7$dic$dic),
                     WAIC = c(model_single7$waic$waic/10000),
                     mlik = c(model_single7$mlik[1,1]),
                     LGOCV_auto = c(ULOOCV_auto_z7)
 )
table 


  beta_est[7]  <- model_single7$summary.fixed["scale_x_s","mean"]

################# species 8


  spat_data8 <-spat_data[spat_data$species ==8, ]

  A_sp8 <- inla.spde.make.A(mesh = mesh,
                          loc = as.matrix(spat_data8[,c("x.loc","y.loc")]),
                          group = spat_data8$time)
 # Detection /non-detection data
 Y_mat <- spat_data8 %>%
   dplyr::select(num_range("y",1:3)) %>% 
   as.matrix()
 # Detection covariates matrix
 X_det <- spat_data8 %>%
   dplyr::select(num_range("g",1:3)) %>% 
   inla.Occupancy_detCov()
 
 
 X_cov <- spat_data8 %>%
  dplyr::select(c(time,scale_x_s)) %>%
   mutate(Int_occ = 1,
          spatialfield = rep(NA,nrow(spat_data8)))


 # Append to a list:
 data_list8 <- as.list(X_cov) 
 data_list8$Y <- Y_mat
 data_list8$X <- X_det

 
formula_z8 <- inla.mdata(Y, X) ~ -1 + Int_occ +
   scale_x_s^2 +
  f(spatialfield, model = spde,
     A.local = A_sp8,
     group = time,
     control.group = list(model = 'ar1', hyper = h.spec))



 model_single8<- inla(formula_z8,
                    data = data_list8,
                    family = 'occupancy',
                   control.fixed = list(prec = 1/2.72, prec.intercept = 1/2.72),
                    control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
                    control.inla = list(int.strategy = "eb"),
                    verbose = FALSE,
                   control.family = list(
                     control.link = list(model = "logit"),
                     link.simple = "logit"
                   ))

 summary(model_single8)

 loocv_m1_auto_z8<- inla.group.cv(result =  model_single8,num.level.sets = 3)
  ULOOCV_auto_z8  = mean(log(loocv_m1_auto_z8$cv),na.rm=T)

table = data.frame(  DIC = c(model_single8$dic$dic),
                     WAIC = c(model_single8$waic$waic/10000),
                     mlik = c(model_single8$mlik[1,1]),
                     LGOCV_auto = c(ULOOCV_auto_z8)
 )
table 


  beta_est[8]  <- model_single8$summary.fixed["scale_x_s","mean"]

################# species 9

  spat_data9 <-spat_data[spat_data$species ==9, ]

  A_sp9 <- inla.spde.make.A(mesh = mesh,
                          loc = as.matrix(spat_data9[,c("x.loc","y.loc")]),
                          group = spat_data9$time)
 # Detection /non-detection data
 Y_mat <- spat_data9 %>%
   dplyr::select(num_range("y",1:3)) %>% 
   as.matrix()
 # Detection covariates matrix
 X_det <- spat_data9 %>%
   dplyr::select(num_range("g",1:3)) %>% 
   inla.Occupancy_detCov()
 
 
 X_cov <- spat_data9 %>%
  dplyr::select(c(time,scale_x_s)) %>%
   mutate(Int_occ = 1,
          spatialfield = rep(NA,nrow(spat_data9)))


 # Append to a list:
 data_list9 <- as.list(X_cov) 
 data_list9$Y <- Y_mat
 data_list9$X <- X_det

formula_z9 <- inla.mdata(Y, X) ~ -1 + Int_occ +
   scale_x_s^2 +
  f(spatialfield, model = spde,
     A.local = A_sp9,
     group = time,
     control.group = list(model = 'ar1', hyper = h.spec))



 model_single9<- inla(formula_z9,
                    data = data_list9,
                    family = 'occupancy',
                   control.fixed = list(prec = 1/2.72, prec.intercept = 1/2.72),
                    control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
                    control.inla = list(int.strategy = "eb"),
                    verbose = FALSE,
                   control.family = list(
                     control.link = list(model = "logit"),
                     link.simple = "logit"
                   ))

 summary(model_single9)

 loocv_m1_auto_z9<- inla.group.cv(result =  model_single9,num.level.sets = 3)
  ULOOCV_auto_z9  = mean(log(loocv_m1_auto_z9$cv),na.rm=T)

table = data.frame(  DIC = c(model_single9$dic$dic),
                     WAIC = c(model_single9$waic$waic/10000),
                     mlik = c(model_single9$mlik[1,1]),
                     LGOCV_auto = c(ULOOCV_auto_z9)
 )

table 

  beta_est[9]  <- model_single9$summary.fixed["scale_x_s","mean"]

################# species 10
 
 spat_data10= spat_data[spat_data$species ==10,] 
 
 
 A_sp10 <- inla.spde.make.A(mesh = mesh,
                          loc = as.matrix(spat_data10[,c("x.loc","y.loc")]),
                          group = spat_data10$time)
 # Detection /non-detection data
 Y_mat <- spat_data10 %>%
   dplyr::select(num_range("y",1:3)) %>% 
   as.matrix()
 # Detection covariates matrix
 X_det <- spat_data10 %>%
   dplyr::select(num_range("g",1:3)) %>% 
   inla.Occupancy_detCov()
 
 
 X_cov <- spat_data10 %>%
  dplyr::select(c(time,scale_x_s)) %>%
   mutate(Int_occ = 1,
          spatialfield = rep(NA,nrow(spat_data10)))


 # Append to a list:
 data_list10 <- as.list(X_cov) 
 data_list10$Y <- Y_mat
 data_list10$X <- X_det

formula_z10 <- inla.mdata(Y, X) ~ -1 + Int_occ +
   scale_x_s^2+
  f(spatialfield, model = spde,
     A.local = A_sp10,
     group = time,
     control.group = list(model = 'ar1', hyper = h.spec))



 model_single10<- inla(formula_z10,
                    data = data_list10,
                    family = 'occupancy',
                   control.fixed = list(prec = 1/2.72, prec.intercept = 1/2.72),
                    control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
                    control.inla = list(int.strategy = "eb"),
                    verbose = FALSE,
                   control.family = list(
                     control.link = list(model = "logit"),
                     link.simple = "logit"
                   ))

 summary(model_single10)


 loocv_m1_auto_z10<- inla.group.cv(result =  model_single10,num.level.sets = 3)
  ULOOCV_auto_z10  <- mean(log(loocv_m1_auto_z10$cv),na.rm=T)


table = data.frame(  DIC = c(model_single10$dic$dic),
                     WAIC = c(model_single10$waic$waic/10000),
                     mlik = c(model_single10$mlik[1,1]),
                     LGOCV_auto = c(ULOOCV_auto_z10)
 )

table 


  beta_est[10]  <- model_single10$summary.fixed["scale_x_s","mean"]


################# species 11

  spat_data11 <-spat_data[spat_data$species ==11, ]

 


 A_sp11 <- inla.spde.make.A(mesh = mesh,
                          loc = as.matrix(spat_data11[,c("x.loc","y.loc")]),
                          group = spat_data11$time)
 # Detection /non-detection data
 Y_mat <- spat_data11 %>%
   dplyr::select(num_range("y",1:3)) %>% 
   as.matrix()
 # Detection covariates matrix
 X_det <- spat_data11 %>%
   dplyr::select(num_range("g",1:3)) %>% 
   inla.Occupancy_detCov()
 
 
 X_cov <- spat_data11 %>%
  dplyr::select(c(time,scale_x_s)) %>%
   mutate(Int_occ = 1,
          spatialfield = rep(NA,nrow(spat_data11)))


 # Append to a list:
 data_list11 <- as.list(X_cov) 
 data_list11$Y <- Y_mat
 data_list11$X <- X_det

 
formula_z11 <- inla.mdata(Y, X) ~ -1 + Int_occ +
   scale_x_s^2 +
  f(spatialfield, model = spde,
     A.local = A_sp11,
     group = time,
     control.group = list(model = 'ar1', hyper = h.spec))



 model_single11<- inla(formula_z11,
                    data = data_list11,
                    family = 'occupancy',
                   control.fixed = list(prec = 1/2.72, prec.intercept = 1/2.72),
                    control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
                    control.inla = list(int.strategy = "eb"),
                    verbose = FALSE,
                   control.family = list(
                     control.link = list(model = "logit"),
                     link.simple = "logit"
                   ))

 summary(model_single11)

 loocv_m1_auto_z11<- inla.group.cv(result =  model_single11,num.level.sets = 3)
  ULOOCV_auto_z11  = mean(log(loocv_m1_auto_z11$cv),na.rm=T)


table = data.frame(  DIC = c(model_single11$dic$dic),
                     WAIC = c(model_single11$waic$waic/10000),
                     mlik = c(model_single11$mlik[1,1]),
                     LGOCV_auto = c(ULOOCV_auto_z11)
 )

table 


  beta_est[11]  <- model_single11$summary.fixed["scale_x_s","mean"]

################# species 12

  spat_data12 <-spat_data[spat_data$species ==12, ]

 

 A_sp12 <- inla.spde.make.A(mesh = mesh,
                          loc = as.matrix(spat_data12[,c("x.loc","y.loc")]),
                          group = spat_data12$time)
 # Detection /non-detection data
 Y_mat <- spat_data12 %>%
   dplyr::select(num_range("y",1:3)) %>% 
   as.matrix()
 # Detection covariates matrix
 X_det <- spat_data12 %>%
   dplyr::select(num_range("g",1:3)) %>% 
   inla.Occupancy_detCov()
 
 
 X_cov <- spat_data12 %>%
  dplyr::select(c(time,scale_x_s)) %>%
   mutate(Int_occ = 1,
          spatialfield = rep(NA,nrow(spat_data12)))

 # Append to a list:
 data_list12 <- as.list(X_cov) 
 data_list12$Y <- Y_mat
 data_list12$X <- X_det

 
formula_z12 <- inla.mdata(Y, X) ~ -1 + Int_occ +
   scale_x_s^2 +
  f(spatialfield, model = spde,
     A.local = A_sp12,
     group = time,
     control.group = list(model = 'ar1', hyper = h.spec))



 model_single12<- inla(formula_z12,
                    data = data_list12,
                    family = 'occupancy',
                   control.fixed = list(prec = 1/2.72, prec.intercept = 1/2.72),
                    control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
                    control.inla = list(int.strategy = "eb"),
                    verbose = FALSE,
                   control.family = list(
                     control.link = list(model = "logit"),
                     link.simple = "logit"
                   ))

 summary(model_single12)

 loocv_m1_auto_z12<- inla.group.cv(result =  model_single12,num.level.sets = 3)
  ULOOCV_auto_z12  = mean(log(loocv_m1_auto_z12$cv),na.rm=T)


table = data.frame(  DIC = c(model_single12$dic$dic),
                     WAIC = c(model_single12$waic$waic/10000),
                     mlik = c(model_single12$mlik[1,1]),
                     LGOCV_auto = c(ULOOCV_auto_z12)
 )

table 


  beta_est[12]  <- model_single12$summary.fixed["scale_x_s","mean"]

##=====================

#MSRE_msom

beta1_msom_est <- model_msom$summary.random$species_slope$mean
beta1_msom_est 

msre_msom <- mean(
  ((beta1_msom_est - beta1) )^2
)


#MSRE_single
msre_single <- mean((beta1-beta_est)^2)


data.frame(
  species = 1:length(beta1),
  beta_true = beta1,
  beta_msom = beta1_msom_est,
  beta_single = beta_est,
  MSRE_single = ((beta_est - beta1))^2,
  MSRE_msom   = ((beta1_msom_est   - beta1))^2
)


####============================
#Predicted occupancy probability


pred_df = data.frame(x = rep(crds(x_covariate)[,1],nT),
                     y = rep(crds(x_covariate)[,2],nT),
                     x_s = rep(as.vector(scale(values(x_covariate$ x_s ))),nT),
                     time = rep(c(1:nT), each= ncells))



pred_dfnew_expanded <- pred_df[rep(1:nrow(pred_df), each = n_species), ]
pred_dfnew_expanded$species_id <- rep(1:n_species, times = nrow(pred_df))
pred_dfnew_expanded$xs_2<-(pred_dfnew_expanded$x_s)^2


pred1 = pred_dfnew_expanded %>% filter(time == 5)


A_pred <- inla.spde.make.A(
  mesh = mesh,
  loc = cbind(pred1$x, pred1$y),
  group = pred1$time
)

sample = inla.posterior.sample(1000, model_msom)

func_pred <- function(...) {
  linear_predictor <- (
     species_slope[pred1$species_id] * pred1$xs_2+
    species_id[pred1$species_id] +
    (A_pred %*% spatialfield)[,1]
  )
  linear_predictor
}


samples_pred <- inla.posterior.sample.eval(func_pred, sample)


pred1$mean <- apply(samples_pred, 1, mean)
pred1$sd <- apply(samples_pred, 1, sd)

pred1$psi <- plogis(pred1$mean)

 
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
                   )


pdf("psis.pdf", width=5, height=4)  

ggplot(pred1, aes(x = x, y = y, fill = psi)) +
  geom_tile() + 
  facet_wrap(~species_id) +
  scale_fill_scico(name = "",    direction = -1 ) +
  coord_equal() +
  theme_maps
dev.off()


pdf("sdS.pdf", width=5, height=4)  

ggplot(pred1, aes(x = x, y = y, fill = sd)) +
  geom_tile() + 
  facet_wrap(~species_id) +
  scale_fill_scico(name = ,    direction = -1 ) +
  coord_equal() +
  theme_maps
dev.off()


