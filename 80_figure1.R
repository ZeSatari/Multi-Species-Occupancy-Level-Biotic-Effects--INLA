## =====================================================================
## 80_figure1.R -- predicted occupancy and posterior uncertainty
## =====================================================================
## Reads results/real/fits.rds; does not refit anything.
##
## Produces TWO panels per species, not one. Reviewers of spatial models
## ask where the estimate is uncertain, and a posterior mean surface on
## its own does not say. The manuscript's earlier figure showed only the
## mean.
##
## Uncertainty is summarised by the posterior standard deviation of psi
## on the PROBABILITY scale, not the logit scale. A logit-scale SD is
## large wherever psi is near 0 or 1 regardless of how well determined it
## is, so it maps sampling geometry rather than knowledge.
## =====================================================================

source("R/00_setup.R")

suppressPackageStartupMessages({
  library(ggplot2)
  library(spOccupancy)
  library(patchwork)
})

data(hbefTrends)
OUTDIR <- "results/real"
FIGDIR <- "figures"
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

fits   <- readRDS(file.path(OUTDIR, "fits.rds"))
m_msom <- fits$msom

SP_NAMES <- dimnames(hbefTrends$y)[[1]]
N_SP     <- length(SP_NAMES)
N_SITE   <- dim(hbefTrends$y)[2]
N_YEAR   <- dim(hbefTrends$y)[3]
## Years are read from the data rather than assumed. If the storage order
## ever changes, the figure would otherwise silently show the wrong slice
## while the caption still said 2018.
YEARS <- dimnames(hbefTrends$y)[[3]]
stopifnot(length(YEARS) == N_YEAR, !is.null(YEARS))

TARGET_YEAR_INDEX <- N_YEAR
TARGET_YEAR       <- YEARS[TARGET_YEAR_INDEX]
message("Mapping survey year ", TARGET_YEAR,
        " (slice ", TARGET_YEAR_INDEX, " of ", N_YEAR, ")")
stopifnot(TARGET_YEAR == "2018")   # the year named in the manuscript

coords <- hbefTrends$coords / 1000
elev_s <- as.numeric(scale(hbefTrends$occ.covs$elev))

## ---------------------------------------------------------------------
## Posterior draws of psi at each site, for one species and one year.
##
## inla.posterior.sample() returns the joint posterior, so the intercept,
## the species effect, the elevation terms and the field are drawn
## together and their dependence is preserved. Combining marginal
## summaries instead would understate or overstate the uncertainty
## depending on their correlation.
##
## Components are located BY NAME. The order in the latent vector is not
## contractual, and an earlier version of this pipeline assumed it was,
## which silently attached species effects to the wrong species.
## ---------------------------------------------------------------------
N_DRAWS <- 1000

message("Drawing ", N_DRAWS, " joint posterior samples...")
samp <- INLA::inla.posterior.sample(n = N_DRAWS, result = m_msom)
nm   <- rownames(samp[[1]]$latent)

idx_of <- function(pattern, expected = 1L) {
  i <- grep(pattern, nm)
  if (length(i) != expected) {
    stop(sprintf("idx_of('%s'): matched %d, expected %d.\n  first names: %s",
                 pattern, length(i), expected,
                 paste(head(nm, 15), collapse = " | ")))
  }
  i
}

i_int  <- idx_of("^Int_occ:")
i_ele  <- idx_of("^elev_s:")
i_ele2 <- idx_of("^elev_s2:")
i_sp   <- vapply(seq_len(N_SP),
                 \(k) idx_of(sprintf("^species_id:%d$", k)), integer(1))

## Field nodes for the target year. The field is stored as
## n.spde blocks of one season each, so season t occupies
## ((t-1)*n.spde + 1):(t*n.spde).
i_field_all <- grep("^spatialfield:", nm)
n_spde <- length(i_field_all) / N_YEAR
stopifnot(n_spde == floor(n_spde))
i_field_t <- i_field_all[((TARGET_YEAR_INDEX - 1) * n_spde + 1):
                         (TARGET_YEAR_INDEX * n_spde)]

## Projector from mesh nodes to site locations, rebuilt to match the fit.
bnd  <- fmesher::fm_nonconvex_hull(coords, convex = 0.3)
mesh <- fmesher::fm_mesh_2d(boundary = bnd, max.edge = c(0.35, 1.2),
                            min.angle = 21, offset = c(0.05, 1.5),
                            cutoff = 0.35)
stopifnot(mesh$n == n_spde)
A_site <- INLA::inla.spde.make.A(mesh, loc = as.matrix(coords))

psi_draws <- function(sp_index) {
  out <- matrix(NA_real_, N_SITE, N_DRAWS)
  for (d in seq_len(N_DRAWS)) {
    v   <- samp[[d]]$latent
    eta <- v[i_int] + v[i_sp[sp_index]] +
           v[i_ele] * elev_s + v[i_ele2] * elev_s^2 +
           as.numeric(A_site %*% v[i_field_t])
    out[, d] <- plogis(eta)
  }
  out
}

## ---------------------------------------------------------------------
## Assemble.
## ---------------------------------------------------------------------
naive_occ <- vapply(SP_NAMES, function(sp) {
  y  <- hbefTrends$y[sp, , , ]
  Y  <- do.call(rbind, lapply(seq_len(N_YEAR), \(t) y[, t, ]))
  ok <- rowSums(!is.na(Y)) > 0
  mean(rowSums(Y[ok, , drop = FALSE] == 1, na.rm = TRUE) > 0)
}, numeric(1))

## Order panels by prevalence so the sparse species, where uncertainty is
## largest and the hierarchical model matters most, are adjacent.
sp_order <- SP_NAMES[order(naive_occ)]

dat <- do.call(rbind, lapply(seq_len(N_SP), function(k) {
  d <- psi_draws(k)
  data.frame(species = SP_NAMES[k],
             x = coords[, 1], y = coords[, 2],
             mean = rowMeans(d),
             sd   = apply(d, 1, sd))
}))
dat$species <- factor(dat$species, levels = sp_order)

saveRDS(dat, file.path(OUTDIR, "figure1_data.rds"))

message(sprintf("psi range %.3f-%.3f | posterior SD range %.3f-%.3f",
                min(dat$mean), max(dat$mean), min(dat$sd), max(dat$sd)))

## ---------------------------------------------------------------------
## Panels. Viridis for the mean (perceptually uniform, colour-blind safe)
## and a single-hue ramp for the SD, so the two are not confused at a
## glance. A shared scale across species is essential: per-panel scales
## would make a rare species look as though it were as widely distributed
## as a common one.
## ---------------------------------------------------------------------
base_theme <- theme_minimal(base_size = 9) +
  theme(panel.grid = element_blank(),
        axis.text  = element_text(size = 6),
        axis.title = element_blank(),
        strip.text = element_text(size = 8, face = "bold"),
        legend.position = "bottom",
        legend.key.width = unit(1.4, "cm"),
        legend.key.height = unit(0.3, "cm"))

p_mean <- ggplot(dat, aes(x, y, colour = mean)) +
  geom_point(size = 0.8) +
  scale_colour_viridis_c(name = "Occupancy probability",
                         limits = c(0, 1), option = "D") +
  coord_equal() +
  facet_wrap(~ species, ncol = 4) +
  labs(subtitle = "(a) Posterior mean") +
  base_theme

p_sd <- ggplot(dat, aes(x, y, colour = sd)) +
  geom_point(size = 0.8) +
  scale_colour_distiller(name = "Posterior standard deviation",
                         palette = "YlOrRd", direction = 1,
                         limits = c(0, max(dat$sd))) +
  coord_equal() +
  facet_wrap(~ species, ncol = 4) +
  labs(subtitle = "(b) Posterior standard deviation") +
  base_theme

fig <- p_mean / p_sd

ggsave(file.path(FIGDIR, "figure1.pdf"), fig, width = 9, height = 11)
ggsave(file.path(FIGDIR, "figure1.png"), fig, width = 9, height = 11, dpi = 300)
message("Wrote ", file.path(FIGDIR, "figure1.pdf"))

## ---------------------------------------------------------------------
## The relationship worth stating in the caption: uncertainty is highest
## for the species with fewest detections, which is where the
## hierarchical structure is doing the most work.
## ---------------------------------------------------------------------
summ <- aggregate(cbind(mean, sd) ~ species, dat, mean)
summ$naive_occ <- naive_occ[as.character(summ$species)]
summ <- summ[order(summ$naive_occ), ]
print(summ, row.names = FALSE, digits = 3)
cat(sprintf("\ncorrelation of mean posterior SD with naive occupancy: %.3f\n",
            cor(summ$sd, summ$naive_occ)))
write.csv(summ, file.path(OUTDIR, "figure1_summary.csv"), row.names = FALSE)
