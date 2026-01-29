
space_time_data <- read.csv("Occ_data_ZA_detect.csv")
 head(space_time_data)


 nT <- length(space_time_data$time %>% unique) # number of time points
 
 # raster data
 x_covariate <- terra::rast('xs_detect.tif')
 
 
 
 x_s_vector <- as.vector(x_covariate)  

 # evaluate covariate at each cell
 space_time_data <- space_time_data %>%
   mutate(x_s = x_s_vector)

 n_rows <- nrow(space_time_data) / nT  
 
 
 space_time_data_sf <- space_time_data %>%
   st_as_sf(coords = c('x.loc','y.loc')) 
 
 
 space_time_data_sf <- space_time_data_sf %>%
         mutate(site = as.factor(cellid))
 
 # (optional) if an i.i.d random effect by sites is to be used.
 levels(space_time_data_sf$site) <- 1:nlevels(space_time_data_sf$site) 
 
 spat_data  = space_time_data_sf %>% 
   dplyr::select(-cellid) %>%
   mutate(x.loc= st_coordinates(space_time_data_sf)[,1],
          y.loc= st_coordinates(space_time_data_sf)[,2]) %>% st_drop_geometry()
 
 

 
 # number of cells in the whole area
 ncells <- length( values(x_covariate))
 
 spat_data <- spat_data %>%
   mutate(scale_x_s = scale(x_s)  %>% c())
 scale_x_s<-scale(x_s)
 


pred_df = data.frame(x = rep(crds(x_covariate)[,1],nT),
                     y = rep(crds(x_covariate)[,2],nT),
                     x_s = rep(as.vector(scale(values(x_covariate$ x_s ))),nT),
                     time = rep(c(1:nT), each= ncells))



 boundary_sf = st_bbox(c(xmin = 0, xmax = 300, ymax = 0, ymin = 300)) %>%
   st_as_sfc() %>% st_as_sf()
 
 mesh = fm_mesh_2d(loc.domain = st_coordinates(boundary_sf)[,1:2],
                     offset = c(-0.1, -.2),
                     max.edge = c(15, 30))
 
 spde <- inla.spde2.pcmatern(mesh = mesh,
                              prior.range = c(100, 0.5),
                               prior.sigma = c(0.5, 0.5))

 A_sp <- inla.spde.make.A(mesh = mesh,
                          loc = as.matrix(spat_data[,c("x.loc","y.loc")]),
                        group = spat_data$time)

 Y_mat <- spat_data %>%
   dplyr::select(num_range("y",1:3)) %>% 
   as.matrix()
 # Detection covariates matrix
 X_det <- spat_data %>%
   dplyr::select(num_range("g2.",1:3),) %>% 
   inla.Occupancy_detCov()

 head(X_det)



 X_cov <- spat_data %>%
  dplyr::select(c(time,site,scale_x_s  )) %>%
  mutate(Int_occ = 1,
         spatialfield = rep(NA,nrow(spat_data)))


 # Append to a list:
 data_list <- as.list(X_cov) 
 data_list$Y <- Y_mat
 data_list$X <- X_det

 #---------------------------
 # PC prior for the temporal autocorrelation
 h.spec <- list(rho = list(prior = 'pc.cor0', param = c(0.5, 0.3)))



formula_A <- inla.mdata(Y, X) ~
  -1 + Int_occ +
  scale_x_s^2 +
  f(spatialfield,
    model = spde,
    A.local = A_sp,
    group = time,
    control.group = list(model = "ar1", hyper = h.spec))


model_A <- inla(formula_A, # the formula
                   data=data_list,  # the data list
                   family= 'occupancy',   
                   control.fixed =  list(prec = 1, prec.intercept = 1),
                   control.compute = list(dic= TRUE,waic = TRUE,config = TRUE),
                   control.inla = list(int.strategy = "eb"),
                   verbose = F,
                   control.family = list(control.link = list(model = "logit"),
                                         link.simple = "logit",
                           hyper = list(beta1 = list(param = c(0,1/3),
                                                     initial = 0),
                                        beta2 = list(param = c(0,1/3)))
                           ))


summary(model_A )




loocv_sim_autoA <- inla.group.cv(result = model_A, num.level.sets = 3)
UlogCV_simA <- mean(log(loocv_sim_autoA$cv), na.rm=TRUE)


#####===========with species B
 
Xdet_2 = cbind(
   X_det[,1:2],Occ_data_B[,1],
   X_det[,3:4], Occ_data_B[,2], 
   X_det[,5:6], Occ_data_B[,3]
)

head(Xdet_2)



 # Append to a list:
 data_list2 <- as.list(X_cov) 
 data_list2$Y <- Y_mat
 data_list2$X <- Xdet_2



formula_B <- inla.mdata(Y, X) ~
  -1 + Int_occ +
  scale_x_s^2 +
  f(spatialfield,
    model = spde,
    A.local = A_sp,
    group = time,
    control.group = list(model = "ar1", hyper = h.spec))


model_B <- inla(formula_B, # the formula
                   data=data_list2,  # the data list
                   family= 'occupancy',   
                   control.fixed =  list(prec = 1, prec.intercept = 1),
                   control.compute = list(dic= TRUE,waic = TRUE,config = TRUE),
                   control.inla = list(int.strategy = "eb"),
                   verbose = F,
                   control.family = list(control.link = list(model = "logit"),
                                         link.simple = "logit",
                           hyper = list(beta1 = list(param = c(0,1/3),
                                                     initial = 0),
                                        beta2 = list(param = c(0,1/3)))
                           ))

  
summary(model_B )




loocv_B_auto <- inla.group.cv(result =  model_xs,num.level.sets = 3)
ULOOCV_auto_B  = mean(log(loocv_B_auto$cv),na.rm=T)





