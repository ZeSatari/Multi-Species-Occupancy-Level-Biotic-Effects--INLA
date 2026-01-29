
#===========library=============

library(INLA)
library(fmesher)
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
#==================================

set.seed(121)
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
boundary_sf = st_bbox(c(xmin = 0, xmax = 300, ymax = 0, ymin = 300)) %>%
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

# Precision matrix

Q = inla.spde.precision(matern_sim, theta = c(log(range_spde),
                                              log(sigma_spde)))
# Simulate three spatial fields
seed = 121
sim_field = inla.qsample(n = 1, Q = Q, seed = seed)
sim_field2 = inla.qsample(n = 1, Q = Q, seed = 111)
#----------------

# Obtain the centroid of each cell
coord_grid  = st_coordinates(customGrid %>% st_centroid())
# A matrix
A_proj = inla.spde.make.A(mesh_sim, loc = coord_grid)

# Spatial components
omega_s = (A_proj %*% sim_field)[,1] # spatial random field
x_s = (A_proj %*% sim_field)[,1] # spatial environmental covariate

x_s2 = (A_proj %*% sim_field2)[,1] # spatial environmental covariate


# create rasters
x_rast = rast(data.frame(x = coord_grid[,1], y = coord_grid[,2],x_s))

x_rast2 = rast(data.frame(x = coord_grid[,1], y = coord_grid[,2],x_s2))

# save raster data

writeRaster(x_rast,file='xs_detect.tif',overwrite=TRUE)
writeRaster(x_rast2,file='xs2_detect.tif',overwrite=TRUE)






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



# -------------------------
# 3. Simulate space-time occupancy data 
# -------------------------

# load helping functions 
 source('spde-book-functions.R')

# Time points
nT <- 5

seed=12345
if (!exists("nsites") || !exists("site_id")) {
  set.seed(12345)                
  nsites <- round(ncells * 0.20)
  site_id <- sample(1:ncells, size = nsites, replace = FALSE)
  customGrid$sample <- ifelse(customGrid$cellid %in% site_id, 1, 0)
} else {
  stopifnot(all(site_id %in% customGrid$cellid))
  if (is.null(customGrid$sample)) customGrid$sample <- ifelse(customGrid$cellid %in% site_id, 1, 0)
}



# parameters 
params <- c(sigma = sigma_spde, range = range_spde)

# generate temporal samples for shared and species-specific fields
epsilon_shared <- book.rspde(coord_grid, range = params[2], seed = 111,
                            sigma = params[1], n = nT, mesh = mesh_sim,
                            return.attributes = TRUE)

rho <- 0.65
omega_shared <- epsilon_shared
for (t in 2:nT) {
  omega_shared[, t] <- rho * omega_shared[, t - 1] + sqrt(1 - rho^2) * epsilon_shared[, t]
}



# -------------------------
# State (occupancy) process: 
# -------------------------
beta0_A <- -0.5
beta1_A <- 0.3

beta0_B <- -0.5
beta1_B <- 0.4

nc <- nrow(coord_grid)
z2.mat <- psi2.mat <- matrix(NA, nrow = nc, ncol = nT)   # species B
z1.mat <- psi1.mat <- matrix(NA, nrow = nc, ncol = nT)   # species A

x_s2 <- scale(x_s2)
x_s <- scale(x_s)

  # species B 
for (t in 1:nT) {
  psi2.mat[, t] <- inla.link.logit(beta0_B + 
   beta1_B * x_s2^2 + 
omega_shared[, t], inverse = TRUE)
  z2.mat[, t]   <- rbinom(n = nc, size = 1, prob = psi2.mat[, t])
}


  # species A
for (t in 1:nT) {
  psi1.mat[, t] <- inla.link.logit(
    beta0_A +
    beta1_A * x_s^2 +
    omega_shared[, t],
    inverse = TRUE
  )

  z1.mat[, t] <- rbinom(nc, 1, psi1.mat[, t])
}

# -------------------------
# Observation (detection) process 
# -------------------------
K <- 3

alpha_A <- c(qlogis(0.3), 0.4)

alpha_B <- c(qlogis(0.25), 0.5)


# survey-level covariate: 
g.t <- array(runif(n = nsites * K * nT, -1, 1), dim = c(nsites, K, nT))
g.t2 <- array(runif(n = nsites * K * nT, -1, 1), dim = c(nsites, K, nT))

yA.t <- vector("list", length = nT)
yB.t <- vector("list", length = nT)

alpha3 <- 0.7

missing.p <- 0.1

for (t in 1:nT) {
  Y_A <- matrix(NA, nrow = nsites, ncol = K)
  Y_B <- matrix(NA, nrow = nsites, ncol = K)
  for (j in 1:K) {
    p_B <- inla.link.logit(alpha_B[1] + alpha_B[2] * g.t[, j, t], inverse = TRUE)  # detection prob for B

    Y_B[, j] <- rbinom(n = nsites, size = 1, prob = z2.mat[site_id, t] * p_B)

    turnNA <- rbinom(nsites, 1, missing.p)
    Y_B[, j] <- ifelse(turnNA == 1, NA, Y_B[, j])

       Y_B_eff <- ifelse(is.na(Y_B[, j]), 0, Y_B[, j])

         p_A <- inla.link.logit(
           alpha_A[1] +
           alpha_A[2] * g.t2[, j, t] +
           alpha3 * Y_B_eff,
           inverse = TRUE
              )

    Y_A[, j] <- rbinom(n = nsites, size = 1, prob = z1.mat[site_id, t] * p_A)


    turnNA <- rbinom(nsites, 1, missing.p)
    Y_A[, j] <- ifelse(turnNA == 1, NA, Y_A[, j])
  }
 yB.t[[t]] <- Y_B
  yA.t[[t]] <- Y_A
}

# -------------------------
# species B
Obs_data_B<- do.call("rbind", yB.t) %>% as.data.frame()
colnames(Obs_data_B) <- paste0("y", 1:K)

# reshape covariates g.t 
g_long <- matrix(NA, nrow = nsites * nT, ncol = K)
for (t in 1:nT) {
  idx <- ((t - 1) * nsites + 1):(t * nsites)
  g_long[idx, ] <- g.t[, , t]
}
colnames(g_long) <- paste0("g2.", 1:K)

Obs_data_B <- cbind(Obs_data_B, g_long)
Obs_data_B$cellid <- rep(site_id, times = nT)      
Obs_data_B$time   <- rep(1:nT, each = nsites)

Z_B_vec <-psi_B_vec <-  c()
for (t in 1:nT) {
  Z_B_vec <- c(Z_B_vec, z2.mat[site_id, t])
  psi_B_vec <- c(psi_B_vec, psi2.mat[site_id, t])

}


Occ_cells <- customGrid %>% st_centroid() %>% filter(sample == 1) %>% dplyr::select(cellid)
coords <- st_coordinates(Occ_cells)
Occ_cells <- Occ_cells %>% st_drop_geometry() %>% mutate(x.loc = coords[,1], y.loc = coords[,2])


Occ_data_B <- left_join(Obs_data_B, Occ_cells, by = "cellid") %>%
               arrange(cellid, time)



write.csv(Occ_data_B, file = "Occ_data_zb_detect.csv", row.names = FALSE)


# -------------------------
# species A


Obs_data_A<- do.call("rbind", yA.t) %>% as.data.frame()
colnames(Obs_data_A) <- paste0("y", 1:K)


g_long <- matrix(NA, nrow = nsites * nT, ncol = K)
for (t in 1:nT) {
  idx <- ((t - 1) * nsites + 1):(t * nsites)
  g_long[idx, ] <- g.t2[, , t]
}
colnames(g_long) <- paste0("g2.", 1:K)

Obs_data_A <- cbind(Obs_data_A, g_long)
Obs_data_A$cellid <- rep(site_id, times = nT)      
Obs_data_A$time   <- rep(1:nT, each = nsites)

Z_A_vec <-psi_A_vec <-  c()
for (t in 1:nT) {
  Z_A_vec <- c(Z_A_vec, z2.mat[site_id, t])
  psi_A_vec <- c(psi_A_vec, psi2.mat[site_id, t])

}



Occ_cells <- customGrid %>% st_centroid() %>% filter(sample == 1) %>% dplyr::select(cellid)
coords <- st_coordinates(Occ_cells)
Occ_cells <- Occ_cells %>% st_drop_geometry() %>% mutate(x.loc = coords[,1], y.loc = coords[,2])

Occ_data_A <- left_join(Obs_data_A, Occ_cells, by = "cellid") %>%
               arrange(cellid, time)



write.csv(Occ_data_A, file = "Occ_data_ZA_detect.csv", row.names = FALSE)





