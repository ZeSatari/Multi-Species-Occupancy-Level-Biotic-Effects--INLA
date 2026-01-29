
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

### LOAD AND PREPARE THE DATA

sp.names <- dimnames(hbefTrends$y)[[1]]
n_species <- length(sp.names)


hozor<-c()
for(i in 1:n_species){

  sp <- sp.names[i]
  y_sp <- hbefTrends$y[sp, , , ]
 
  Y_sp <- data.frame(y_sp[,1,])
  for (t in 2:9) {
    Y_sp <- rbind(Y_sp, data.frame(y_sp[,t,]))
  }


total_obs <- nrow(Y_sp )
  count_ones <- sum(Y_sp[,1] == 1, na.rm = TRUE)

hozor[i]<-round(count_ones/total_obs*100, 1)
}


################
Y_list <- list()
Xocc_list <- list()
Xdet_list <- list()

for (i in 1:n_species) {
  sp <- sp.names[i]
  y_sp <- hbefTrends$y[sp, , , ]
  coords <- hbefTrends$coords / 1000
  elev <- hbefTrends$occ.covs$elev
  
  Y_sp <- data.frame(y_sp[,1,])
  for (j in 2:9) {
    Y_sp <- rbind(Y_sp, data.frame(y_sp[,j,]))
  }
  
  Xocc_sp <- data.frame(
    species = rep(sp, 373*9),
    x = rep(coords[,1], 9),
    y = rep(coords[,2], 9),
    elev = rep(elev, 9),
    scale_elev = scale(rep(elev, 9)),
    scale_elev2 = scale(rep(elev, 9))^2,
    site = rep(1:373, 9),
    time = rep(1:9, each = 373),
    Int_occ = 1
  )
  
  Y_list[[sp]] <- Y_sp
  Xocc_list[[sp]] <- Xocc_sp
}

Y_all <- do.call(rbind, Y_list)
Xocc_all <- do.call(rbind, Xocc_list)
Xocc_all$species_id <- as.numeric(factor(Xocc_all$species))
Xocc_all$spatialfield <- NA


 head(Xocc_all)


###=====With identical visit-specific detection covariates

data <- hbefTrends

data$y <- data$y[sp.names == 'REVI', , , ]
data$coords = data$coords/1000 



Y = data.frame(data$y[,1,])
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
##############



data_list_msom <- as.list(Xocc_all)
data_list_msom$Y <- Y_all
data_list_msom$X <- Xdet

boundary = inla.nonconvex.hull(points = data$coords, convex = .3)
mesh = inla.mesh.2d(boundary = boundary,
                    #   loc = cbind(data$X, data$Y),
                    max.edge = c(0.1,0.7),
                    min.angle = 20,
                    offset = c(.01, 1),
                    cutoff = 0.12,
)




h.spec <- list(rho = list(prior = 'pc.cor0', param = c(0.5, 0.1)))

spde <- inla.spde2.pcmatern(
  mesh = mesh, 
  prior.range = c(5, 0.7),
  prior.sigma = c(1, 0.5),
  constr = T) 


A_sp_msom <- inla.spde.make.A(
  mesh = mesh,
  loc = cbind(Xocc_all$x, Xocc_all$y),
  group = Xocc_all$time
)


formula_msom <- inla.mdata(Y_all, Xdet) ~ 
 -1 + Int_occ +
 scale_elev + scale_elev2 +
 f(species_id, model = "iid") +
 f(time, model = "ar1") + 
 f(spatialfield, model = spde, A.local = A_sp_msom, group = time,
 control.group = list(model = 'ar1', hyper = h.spec))

model_msom <- inla(formula_msom, 
 data = data_list_msom, 
 family = 'occupancy', 
 control.fixed = list(prec = 1, prec.intercept = 1), 
 control.compute = list(dic = TRUE,waic = TRUE, config = TRUE), 
 verbose = FALSE, 
 control.family = list( 
 control.link = list(model = "logit"), 
 link.simple = "logit"  ))

summary(model_msom )



species_effects <- model_msom$summary.random$species_id
plot_inla_effects(species_effects)



species_effects <- model_msom$summary.random$species_id

ggplot(species_effects, aes(x = ID, y = mean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = `0.025quant`, ymax = `0.975quant`), width = 0.2) +
  labs(x = "Species ID", y = "γᵢ (Species-level effect)") +
  theme_minimal()

species_effects$species <- sp.names
ggplot(species_effects, aes(x = species, y = mean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = `0.025quant`, ymax = `0.975quant`), width = 0.2) +
  labs(x = "Species", y = "γᵢ (Species-level effect)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



species_effects$species <- sp.names  

species_effects[, c("species", "mean", "sd", "0.025quant", "0.975quant")]

ggplot(species_effects, aes(x = reorder(species, sd), y = sd)) +
  geom_col(fill = "steelblue") +
  labs(x = "Species", y = "Uncertainty (SD of γᵢ)", title = "Uncertainty in Species Effects") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



#LGOC
df_sf <- Xocc_all %>% mutate(id = 1:nrow(Xocc_all)) %>% st_as_sf(coords =c("x","y"))

# create buffer of size 300 (based on estimated range) centred at each site

buffer_25 <- st_buffer(df_sf, dist = 0.75) 

# empty lists to include the indexes of the leave-out-group for each observation i
I_i <- list()

# loop though each observation and store the leave-out-group based on the buffer
for( i in 1:nrow(df_sf)){
  
  # Temporal filtering of data within a 2 years of span of  observation i
  df_sf_subset <- df_sf %>% 
    filter( between(time,left = df_sf$time[i]-2, right = df_sf$time[i]+2)) 
  # Spatial filtering of the observations that are within the buffer of the ith observation
  Buffer_i <-df_sf_subset %>% st_intersects(buffer_25[i,],sparse = FALSE) %>% # identify 
    unlist()
  
  # obtain the indexes of the leave out group
  I_i[[i]] <-  df_sf_subset[Buffer_i,] %>%  pull(id)
  
}


loocv_m1_auto <- inla.group.cv(result =  model_msom,num.level.sets = 3)
ULOOCV_auto  = mean(log(loocv_m1_auto$cv),na.rm=T)



loocv_m1_manual <- inla.group.cv(result =  model_msom,groups = I_i )
ULOOCV_manual = mean(log(loocv_m1_manual$cv),na.rm=T)


 data.frame(auto= ULOOCV_auto,
            manual=ULOOCV_manual)

table = data.frame(  DIC = c(model_msom$dic$dic),
                    WAIC = c(model_msom$waic$waic),
                    mlik = c(model_msom$mlik[1,1]),
                    LGOCV_auto = c(ULOOCV_auto)
)

###########################


elev_raster= rast(data.frame(x =  hbefElev$Easting/1000,
                                    y = hbefElev$Northing/1000,
                                    elev =  hbefElev$val))

elev_raster$scale_elev = (elev_raster$elev - mean(data$occ.covs$elev)) / sd(data$occ.covs$elev)
elev_raster$scale_elev2 = elev_raster$scale_elev^2 

all_elev = inla.group(c(Xocc_all$elev, values(elev_raster$elev)))

val_elev = sort(na.omit(unique(all_elev)))
Xocc_all$group_elev = all_elev[1:dim(Xocc_all)[1]]
elev_group_raster = elev_raster$elev
values(elev_group_raster) <- all_elev[-c(1:dim(Xocc_all)[1])]


pred_df = data.frame(x = rep(crds(elev_raster, na.rm = F)[,1],9),
                     y = rep(crds(elev_raster, na.rm = F)[,2],9),
                     scale_elev = rep(values(elev_raster$scale_elev),9),
                     scale_elev2 = rep(values(elev_raster$scale_elev2),9),
                     group_elev = rep(values(elev_group_raster),9),
                     time = rep(c(1:9), each= length(values(elev_raster$elev)))
                     ) %>%
  dplyr::filter(!is.na(scale_elev))


n_species <- length(sp.names)
pred_df_expanded <- pred_df[rep(1:nrow(pred_df), each = n_species), ]
pred_df_expanded$species_id <- rep(1:n_species, times = nrow(pred_df))


pred1 = pred_df_expanded %>% filter(time == 9)



A_pred <- inla.spde.make.A(
  mesh = mesh,
  loc = cbind(pred1$x, pred1$y),
  group = pred1$time
)

sample = inla.posterior.sample(1000, model_msom)

func_pred <- function(...) {
  linear_predictor <- (
    Int_occ +
    scale_elev * pred1$scale_elev +
    scale_elev2 * pred1$scale_elev2 +
    species_id[pred1$species_id] +
    time[pred1$time] +
    (A_pred %*% spatialfield)[,1]
  )
  linear_predictor
}


samples_pred <- inla.posterior.sample.eval(func_pred, sample)


pred1$mean <- apply(samples_pred, 1, mean)
pred1$sd <- apply(samples_pred, 1, sd)

pred1$psi <- plogis(pred1$mean)
pred1$species_name <- sp.names[pred1$species_id]


 pdf("psiR.pdf", width=8, height=4) 

ggplot(pred1, aes(x = x, y = y, fill = psi)) +
  geom_tile() + 
  facet_wrap(~species_name) +
  scale_fill_scico(name = "",    direction = -1 ) +
  coord_equal() +
  theme_maps

dev.off()


pdf("sdR.pdf", width=8, height=4)  

ggplot(pred1, aes(x = x, y = y, fill = sd)) +
  geom_tile() + 
  facet_wrap(~species_name) +
  scale_fill_scico(name = "Occupancy probability",    direction = -1 ) +
  coord_equal() +
  theme_maps
dev.off()


