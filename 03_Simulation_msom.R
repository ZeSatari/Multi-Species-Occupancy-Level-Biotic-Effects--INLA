
#===========library=============

library(INLA)
library(fmesher)-
library(tidyverse)
library(spatstat)
library(sf)
library(terra)
library(ggplot2)
library(gt)
library(dplyr)
library(viridis)
library(viridisLite)
library(scico)
library(patchwork)
library(MASS)

#==================================

set.seed(123)
#==========Data Simulation=========
#==========1 Set up================
# Define spatial domain
win <- owin(c(0,300), c(0,300))
npix <- 1000

Domain <- rast(nrows=npix, ncols=npix,
               xmax=win$xrange[2],xmin=win$xrange[1],
               ymax = win$yrange[2],ymin=win$yrange[1])


values(Domain) <- 1:ncell(Domain)
xy <- crds(Domain)

# Define regular grid
cell_size = 3
customGrid <- st_make_grid(Domain,cellsize = c(cell_size,cell_size)) %>% 
  st_cast("MULTIPOLYGON") %>%
  st_sf() %>%
  mutate(cellid = row_number())

# number of cells
ncells <- nrow(customGrid)
#==========2 Simulate data for a simple spatial occupancy model========
# Spatial boundary
boundary_sf = st_bbox(c(xmin = 0, xmax = 300, ymin= 0, ymax  = 300)) %>%
  st_as_sfc()

# Create a fine mesh
mesh_sim = fm_mesh_2d(loc.domain = st_coordinates(boundary_sf)[,1:2],
                offset = c(-0.1, -.2),
                 max.edge = c(4, 50))


# Matern model
matern_sim <- inla.spde2.pcmatern(mesh_sim,
                                  prior.range = c(100, 0.5),
                                  prior.sigma = c(0.5, 0.5))


range_spde = 100
sigma_spde = 0.5

params<- c(sigma=sigma_spde,range=range_spde)


# Precision matrix
Q1 = inla.spde.precision(matern_sim, theta = c(log(range_spde),
                                              log(sigma_spde)))


# Simulate three spatial fields
seed = 123
sim_field = inla.qsample(n = 1, Q = Q1, seed = seed)

#----------------

# Obtain the centroid of each cell
coord_grid  = st_coordinates(customGrid %>% st_centroid())
# A matrix
A_proj = inla.spde.make.A(mesh_sim, loc = coord_grid)

# Spatial components
omega_s = (A_proj %*% sim_field)[,1] # spatial random field
x_s = (A_proj %*% sim_field)[,1] # spatial environmental covariate
head(x_s)
 
# create rasters
x_rast = rast(data.frame(x = coord_grid[,1], y = coord_grid[,2],x_s))

# save raster data

writeRaster(x_rast,file='xs.tif',overwrite=TRUE)

#==============inla.Occupancy_detCov==========
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

    K<- lapply(X_det, ncol) %>% unlist() %>% max()
    M<- lapply(X_det, nrow) %>% unlist() %>% unique()
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

  X_mat <- matrix(NA, nrow = M, ncol = K * (P + 1))
  X_mat[, seq(1, (K * (P + 1)), by = (P + 1))] <- 1  # intercepts
  X_mat[, which(!(1:(K * (P + 1)) %in% seq(1, (K * (P + 1)), by = (P + 1))))] <- covariates
  return(X_mat)
}

x_s_scale=scale(x_s)

xs2<-x_s_scale^2

# -------------------------
# 3. Simulate space-time occupancy data 
# -------------------------

# load helping functions 
 source('spde-book-functions.R')

n_species <- 12
nT <- 5
K <- 3
seed <- 1890

set.seed(seed)
nsites = round(ncells*.20)

R2 <- matrix(0.5, n_species, n_species)  
R2[,4]<-R2[4,]<-0.25
R2[1,4]<-R2[4,1]<-0.4
diag(R2) <- 1  

L2 <- chol(R2) 

sigma <- 0.8  
cov_matrix2 <- (sigma^2) * R2  

set.seed(seed)
gamma_species <- mvrnorm(1, 
                        mu = rep(-0.5, n_species),  
                        Sigma = cov_matrix2)       


mu_beta0= -0.5

beta1 <- rnorm(n_species, mean=0.7, sd=1)
mu_beta1<-mean(beta1)


rho<- 0.6

 epsilon.t <- book.rspde(coord_grid, range = range_spde, seed = seed,
                           sigma = sigma_spde, n = nT, mesh = mesh_sim,
                           return.attributes = TRUE) 
  omega_st=epsilon.t
for (t in 2:nT){
  omega_st[, t] <- rho * omega_st[, t - 1] + sqrt(1 - rho^2) * epsilon.t[, t]
}


z_array <- array(NA, dim = c(ncells, nT, n_species))
psi_array <- array(NA, dim = c(ncells, nT, n_species))

for (s in 1:n_species) {
  for (t in 1:nT) {
    psi_array[,t,s] <- inla.link.logit(mu_beta0 +  beta1[s] * xs2+ 
                                               + gamma_species[s]+ 
                                               omega_st[,t], inverse = TRUE)
    z_array[,t,s] <- rbinom(ncells, size = 1, prob = psi_array[,t,s])
  }
}


alpha0 <- rnorm(n_species, -0.9, 0.2)


alpha1 <- rnorm(n_species, -1, 0.3)


if (!exists("nsites") || !exists("site_id")) {
  set.seed(seed)                
  nsites <- round(ncells * 0.20)
  site_id <- sample(1:ncells, size = nsites, replace = FALSE)
  customGrid$sample <- ifelse(customGrid$cellid %in% site_id, 1, 0)
} else {
  stopifnot(all(site_id %in% customGrid$cellid))
  if (is.null(customGrid$sample)) customGrid$sample <- ifelse(customGrid$cellid %in% site_id, 1, 0)
}

y_list <- list()
missing.p <- 0.1

Occ_data_multi_list <- list()

for (s in 1:n_species) {
  g.t <- array(runif(n = nsites * K * nT, -1, 1), dim = c(nsites, K, nT))
  y.t <- list()
  
  for (t in 1:nT) {
    y_mat <- matrix(NA, nrow = nsites, ncol = K)
    for (j in 1:K) {
      p <- inla.link.logit(alpha0[s] + alpha1[s]*g.t[,j,t], inverse = T)
      y_mat[,j] <- rbinom(nsites, 1, prob = z_array[site_id,t,s] * p)
      turnNA <- rbinom(nsites, 1, missing.p)
      y_mat[, j] <- ifelse(turnNA == 1, NA, y_mat[, j])
    }
    
    df <- data.frame(
      cellid = site_id,
      time = t,
      y1 = y_mat[,1],
      y2 = y_mat[,2],
      y3 = y_mat[,3],
      g1 = g.t[,1,t],
      g2 = g.t[,2,t],
      g3 = g.t[,3,t],
      species = s
    )
    
    y.t[[t]] <- df
  }
  
  Occ_data_multi_list[[s]] <- do.call("rbind", y.t)
}

Occ_data_multi <- do.call("rbind", Occ_data_multi_list)

centroids <- st_coordinates(st_centroid(customGrid[site_id,]))
Occ_data_multi$x.loc <- centroids[,1]
Occ_data_multi$y.loc <- centroids[,2]



write.csv(Occ_data_multi, file = "Occ_data_multi.csv", row.names = FALSE)

hozor2<-c()
for(j in 1:n_species){
data <- Occ_data_multi[Occ_data_multi$species == j, ]


total_obs <- nrow(data)
  count_ones <- sum(data[3] == 1, na.rm = TRUE)

  
  cat(paste("y1",j, ":\n"))
  cat(paste("  - تعداد 1ها:", count_ones, 
            "(", round(count_ones/total_obs*100, 1), "%)\n"))
hozor2[j]<-count_ones/total_obs*100
}


mean_hozor2 <- mean(hozor2)

range(hozor2)


