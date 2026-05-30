# Compare alternative random-effects structures for the sugarcane MA.
# All models share the SAME fixed effects (PCY ~ delta_T_C) and V = 0,
# so REML-based logLik/AIC/BIC ARE comparable across them.

suppressPackageStartupMessages(library(metafor))

dat <- read.csv("ma_data_clean.csv", stringsAsFactors = FALSE)
dat <- dat[dat$include == "yes", ]
dat$PCY_pct   <- as.numeric(dat$PCY_pct)
dat$delta_T_C <- as.numeric(dat$delta_T_C)
dat <- dat[!is.na(dat$PCY_pct) & !is.na(dat$delta_T_C), ]
dat$obs <- seq_len(nrow(dat))
dat$es_id <- factor(dat$obs)   # observation-level residual (obs is globally unique)
dat$vi  <- 0

cat(sprintf("n = %d effects, %d studies, %d sites, %d crop models\n\n",
            nrow(dat), length(unique(dat$study_id)),
            length(unique(dat$site)), length(unique(dat$crop_model))))

fit <- function(label, randstr) {
  m <- try(rma.mv(PCY_pct, vi, mods = ~ delta_T_C, random = randstr,
                  data = dat, method = "REML"), silent = TRUE)
  if (inherits(m, "try-error")) {
    cat("------------------------------------------------------------\n")
    cat(label, "\n  FAILED:", attr(m, "condition")$message, "\n\n")
    return(invisible(NULL))
  }
  b   <- coef(m)["delta_T_C"]
  ci  <- c(m$ci.lb[2], m$ci.ub[2])
  vc  <- paste(sprintf("%s=%.1f", m$s.names, m$sigma2), collapse = ", ")
  cat("------------------------------------------------------------\n")
  cat(label, "\n")
  cat(sprintf("  delta_T slope = %+.2f%%/C  (95%% CI %+.2f to %+.2f, QM p = %.4g)\n",
              b, ci[1], ci[2], m$QMp))
  cat(sprintf("  variance components: %s\n", vc))
  cat(sprintf("  logLik = %.2f | AIC = %.2f | BIC = %.2f\n",
              logLik(m), AIC(m), BIC(m)))
  invisible(m)
}

# (1) site only  -- NB: with V=0 there is no observation-level residual term
fit("(1) site only:                random = ~ 1 | site",
    ~ 1 | site)

# (2) site + observation residual  (study_id:obs == per-effect residual here)
fit("(2) site + obs residual:      list(~ 1 | site, ~ 1 | es_id)",
    list(~ 1 | site, ~ 1 | es_id))

# (3) site + nested study/obs
fit("(3) site + study/obs:         list(~ 1 | site, ~ 1 | study_id/obs)",
    list(~ 1 | site, ~ 1 | study_id/obs))

# (4) lecture-canonical 3-level: effects nested in studies (no site, no crop)
fit("(4) study/obs (3-level, lecture): random = ~ 1 | study_id/obs",
    ~ 1 | study_id/obs)

# (5) current model: study/obs + crossed crop_model
fit("(5) study/obs + crop_model (current): list(~ 1 | study_id/obs, ~ 1 | crop_model)",
    list(~ 1 | study_id/obs, ~ 1 | crop_model))

cat("------------------------------------------------------------\n")
cat("Note: lower AIC/BIC = better fit (comparable: same fixed effects, REML).\n")

cat("\n\n############ APPLES-TO-APPLES: crop_model held in ALL models ############\n")
cat("(only the clustering term varies; crop_model crossed in every model)\n\n")

# A) site + crop_model
fit("(A) site + crop_model:                 list(~1|site, ~1|crop_model)",
    list(~ 1 | site, ~ 1 | crop_model))

# B) site + obs residual + crop_model
fit("(B) site + obs + crop_model:           list(~1|site, ~1|es_id, ~1|crop_model)",
    list(~ 1 | site, ~ 1 | es_id, ~ 1 | crop_model))

# C) site + study/obs + crop_model
fit("(C) site + study/obs + crop_model:     list(~1|site, ~1|study_id/obs, ~1|crop_model)",
    list(~ 1 | site, ~ 1 | study_id/obs, ~ 1 | crop_model))

# D) study/obs + crop_model  (= model 5, the reference)
fit("(D) study/obs + crop_model (= model 5): list(~1|study_id/obs, ~1|crop_model)",
    list(~ 1 | study_id/obs, ~ 1 | crop_model))

cat("------------------------------------------------------------\n")
cat("Note: lower AIC/BIC = better fit (comparable: same fixed effects, REML).\n")
