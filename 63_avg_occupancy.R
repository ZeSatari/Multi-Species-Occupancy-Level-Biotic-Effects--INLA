## =====================================================================
## 63_avg_occupancy.R -- is the inadmissibility claim actually supported?
## =====================================================================
## The manuscript compares naive occupancy against
##      inv.logit(beta0 + gamma_i),
## the occupancy probability at the covariate and field means. That is
## not the quantity naive occupancy estimates. By Jensen's inequality
##      inv.logit(E[eta]) != E[inv.logit(eta)],
## and the two can differ substantially when the linear predictor varies
## across sites, which is exactly what a spatial field is for.
##
## The comparable quantity is the average occupancy probability over the
## site-years actually surveyed,
##      psi_bar_i = (1/ST) sum_{s,t} E[psi_ist | y],
## and, better still, the posterior predictive probability that a species
## is recorded at least once in K visits,
##      P(recorded)_i = (1/ST) sum_{s,t} E[ psi_ist * (1-(1-p)^K) | y ],
## which is what naive occupancy is a realisation of.
##
## This script computes both, for the hierarchical model and for the
## independent fits, and compares them with the observed frequency.
##
## The claim in the manuscript stands only if P(recorded) under the
## independent fits falls materially below the observed frequency. If it
## does not, the claim must be withdrawn.
## =====================================================================

source("R/00_setup.R")
source("R/10_metrics.R")
source("R/30_fit_occupancy.R")

suppressPackageStartupMessages(library(spOccupancy))
data(hbefTrends)

OUTDIR <- "results/real"
fits <- readRDS(file.path(OUTDIR, "fits.rds"))
m_msom <- fits$msom
ssom   <- fits$ssom

SP_NAMES <- dimnames(hbefTrends$y)[[1]]
N_SP   <- length(SP_NAMES)
N_SITE <- dim(hbefTrends$y)[2]
N_YEAR <- dim(hbefTrends$y)[3]
K      <- dim(hbefTrends$y)[4]
N_ROW  <- N_SITE * N_YEAR
N_DRAWS <- 1000

stack_species <- function(sp) {
  y <- hbefTrends$y[sp, , , ]
  do.call(rbind, lapply(seq_len(N_YEAR), \(t) y[, t, ]))
}
stack_det_cov <- function(nm) {
  z <- hbefTrends$det.covs[[nm]]
  M <- do.call(rbind, lapply(seq_len(N_YEAR), \(t) z[, t, ]))
  (M - mean(M, na.rm = TRUE)) / sd(M, na.rm = TRUE)
}

coords <- hbefTrends$coords / 1000
elev_s <- as.numeric(scale(hbefTrends$occ.covs$elev))
elev_rep  <- rep(elev_s,   times = N_YEAR)
time_rep  <- rep(seq_len(N_YEAR), each = N_SITE)

## Observed frequency: proportion of site-years with at least one
## detection, over site-years that were actually surveyed.
observed <- vapply(SP_NAMES, function(sp) {
  Y <- stack_species(sp)
  ok <- rowSums(!is.na(Y)) > 0
  mean(rowSums(Y[ok, , drop = FALSE] == 1, na.rm = TRUE) > 0)
}, numeric(1))

## ---------------------------------------------------------------------
## Rebuild the projector used in the fit.
## ---------------------------------------------------------------------
bnd  <- fmesher::fm_nonconvex_hull(coords, convex = 0.3)
mesh <- fmesher::fm_mesh_2d(boundary = bnd, max.edge = c(0.35, 1.2),
                            min.angle = 21, offset = c(0.05, 1.5),
                            cutoff = 0.35)
A_site <- INLA::inla.spde.make.A(mesh, loc = as.matrix(coords))

## ---------------------------------------------------------------------
## Joint posterior draws, components located BY NAME.
## ---------------------------------------------------------------------
draw_psi <- function(fit, species_index = NULL) {
  s  <- INLA::inla.posterior.sample(n = N_DRAWS, result = fit)
  nm <- rownames(s[[1]]$latent)

  idx1 <- function(pat) {
    i <- grep(pat, nm)
    if (length(i) != 1L) stop("pattern '", pat, "' matched ", length(i))
    i
  }
  i_int  <- idx1("^Int_occ:")
  i_ele  <- idx1("^elev_s:")
  i_ele2 <- idx1("^elev_s2:")
  i_sp   <- if (is.null(species_index)) NA_integer_ else
              idx1(sprintf("^species_id:%d$", species_index))

  i_field <- grep("^spatialfield:", nm)
  n_spde  <- length(i_field) / N_YEAR
  stopifnot(n_spde == mesh$n)

  out <- matrix(NA_real_, N_ROW, N_DRAWS)
  for (d in seq_len(N_DRAWS)) {
    v <- s[[d]]$latent
    fixed <- v[i_int] + v[i_ele] * elev_rep + v[i_ele2] * elev_rep^2
    if (!is.na(i_sp)) fixed <- fixed + v[i_sp]
    ## field value at each site-year: pick the slice matching its year
    fld <- numeric(N_ROW)
    for (tt in seq_len(N_YEAR)) {
      sl <- i_field[((tt - 1) * n_spde + 1):(tt * n_spde)]
      fld[time_rep == tt] <- as.numeric(A_site %*% v[sl])
    }
    out[, d] <- plogis(fixed + fld)
  }
  out
}

## Detection probability, averaged over the draws of the detection
## hyperparameters and over the observed survey covariates.
p_detect_mean <- function(fit, Xday, Xtod) {
  dr <- det_draws(fit, n = N_DRAWS)          # alpha0, alpha1, ...
  ## columns: intercept, day, tod (biotic term absent from these fits)
  pm <- matrix(NA_real_, N_ROW, K)
  for (j in seq_len(K)) {
    eta <- outer(rep(1, N_ROW), dr[, 1]) +
           outer(Xday[, j], dr[, 2]) +
           outer(Xtod[, j], dr[, 3])
    pm[, j] <- rowMeans(plogis(eta))
  }
  pm
}

Xday <- stack_det_cov("day"); Xtod <- stack_det_cov("tod")
Xday[is.na(Xday)] <- 0;       Xtod[is.na(Xtod)] <- 0

## ---------------------------------------------------------------------
## Compute for every species, under both models.
## ---------------------------------------------------------------------
res <- list()
for (i in seq_len(N_SP)) {
  sp <- SP_NAMES[i]
  message("species ", sp)
  Y  <- stack_species(sp)
  surveyed <- rowSums(!is.na(Y)) > 0

  psi_h <- draw_psi(m_msom, species_index = i)
  psi_s <- draw_psi(ssom[[sp]], species_index = NULL)

  p_h <- p_detect_mean(m_msom,   Xday, Xtod)
  p_s <- p_detect_mean(ssom[[sp]], Xday, Xtod)

  ## P(at least one detection in the visits actually made)
  nvis <- rowSums(!is.na(Y))
  det_any <- function(pm) 1 - apply(1 - pm, 1, function(z) prod(z[seq_len(K)]))

  rec_h <- rowMeans(psi_h) * det_any(p_h)
  rec_s <- rowMeans(psi_s) * det_any(p_s)

  res[[i]] <- data.frame(
    species        = sp,
    observed       = observed[[sp]],
    psi_bar_msom   = mean(rowMeans(psi_h)[surveyed]),
    psi_bar_ssom   = mean(rowMeans(psi_s)[surveyed]),
    psi_mean_msom  = plogis(m_msom$summary.fixed["Int_occ", "mean"] +
                            m_msom$summary.random$species_id$mean[i]),
    psi_mean_ssom  = plogis(ssom[[sp]]$summary.fixed["Int_occ", "mean"]),
    prec_msom      = mean(rec_h[surveyed]),
    prec_ssom      = mean(rec_s[surveyed])
  )
}
out <- do.call(rbind, res)
out <- out[order(out$observed), ]

cat("\n", strrep("=", 92), "\n", sep = "")
cat("psi_mean_* : occupancy at the covariate and field means (what the paper now reports)\n")
cat("psi_bar_*  : average posterior occupancy over surveyed site-years (Jensen-correct)\n")
cat("prec_*     : posterior predictive P(recorded at least once), comparable to 'observed'\n\n")
print(out, row.names = FALSE, digits = 3)

cat("\n--- does the independent fit fall below the observed frequency? ---\n")
out$viol_mean <- out$psi_mean_ssom < out$observed
out$viol_bar  <- out$psi_bar_ssom  < out$observed
out$viol_prec <- out$prec_ssom     < out$observed
cat(sprintf("  using psi at means      : %d of %d species\n", sum(out$viol_mean), nrow(out)))
cat(sprintf("  using averaged psi      : %d of %d species\n", sum(out$viol_bar),  nrow(out)))
cat(sprintf("  using P(recorded)       : %d of %d species\n", sum(out$viol_prec), nrow(out)))
cat("\n  same three checks for the hierarchical model:\n")
cat(sprintf("  using psi at means      : %d\n", sum(out$psi_mean_msom < out$observed)))
cat(sprintf("  using averaged psi      : %d\n", sum(out$psi_bar_msom  < out$observed)))
cat(sprintf("  using P(recorded)       : %d\n", sum(out$prec_msom     < out$observed)))

cat("\nThe manuscript's claim requires the P(recorded) row to show a clear\n")
cat("shortfall for the independent fits and none for the hierarchical model.\n")
cat("If it does not, the claim must be withdrawn or restated.\n")

write.csv(out, file.path(OUTDIR, "avg_occupancy.csv"), row.names = FALSE)
