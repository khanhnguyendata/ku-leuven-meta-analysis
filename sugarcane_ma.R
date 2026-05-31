library(metafor)

dat <- read.csv("ma_data_clean.csv", stringsAsFactors = FALSE)

dat <- dat[dat$include == "yes", ]
dat$PCY_pct   <- as.numeric(dat$PCY_pct)
dat$delta_T_C <- as.numeric(dat$delta_T_C)
dat$study_id <- as.factor(dat$study_id)
dat <- dat[!is.na(dat$PCY_pct) & !is.na(dat$delta_T_C), ]
dat$sim <- factor(seq_len(nrow(dat)))
dat$vi  <- 0

m_sim <- rma.mv(yi = PCY_pct, V = vi, mods = ~ delta_T_C,
                 random = list(~ 1 | sim),
                 data = dat, method = "REML")
summary(m_sim)

m_li <- rma.mv(yi = PCY_pct, V = vi, mods = ~ delta_T_C,
                random = list(~ 1 | sim, ~ 1 | site, ~ 1 | crop_model),
                data = dat, method = "REML")
summary(m_li)
anova(m_sim, m_li)

m_full <- rma.mv(yi = PCY_pct, V = vi, mods = ~ delta_T_C,
                random = list(~ 1 | study_id/sim, ~ 1 | site, ~ 1 | crop_model),
                data = dat, method = "REML")
summary(m_full)
anova(m_li, m_full)

m_final <- rma.mv(yi = PCY_pct, V = vi, mods = ~ delta_T_C,
                random = list(~ 1 | study_id/sim, ~ 1 | crop_model),
                data = dat, method = "REML")
summary(m_final)
anova(m_final, m_full)

# --- Check: should crop_model vary the SLOPE (random slope) or be a fixed effect? ---
# Only 3 crop models (DSSAT-CANEGRO=39 sim/7 studies, APSIM-Sugar=4 sim/1 study,
# EPIC=2 sim/1 study with a single delta_T). crop_model is thus nearly nested in
# study_id, so any slope difference is confounded with study/design.
m_rslope <- tryCatch(
  rma.mv(yi = PCY_pct, V = vi, mods = ~ delta_T_C,
         random = list(~ 1 | study_id/sim, ~ delta_T_C | crop_model),
         struct = "GEN", data = dat, method = "REML"),
  error = function(e) e)
if (inherits(m_rslope, "error")) {
  cat("m_rslope problem:", conditionMessage(m_rslope), "\n")
} else {
  print(summary(m_rslope))
  print(anova(m_final, m_rslope))
}

# crop_model (and its interaction with delta_T) as FIXED effects instead.
m_fixint <- tryCatch(
  rma.mv(yi = PCY_pct, V = vi, mods = ~ delta_T_C * crop_model,
         random = list(~ 1 | study_id/sim), data = dat, method = "REML", test = "t"),
  error = function(e) e)
if (inherits(m_fixint, "error"))
  cat("m_fixint problem:", conditionMessage(m_fixint), "\n") else summary(m_fixint)

# Conclusion of the crop_model slope check:
#  - m_rslope "fits" and beats m_final by LRT (p~0.0003), BUT the intercept-slope
#    correlation hits the boundary rho = 1.00 -> singular / over-parameterized fit
#    (expected with only 3 outer levels). Not a trustworthy/generalizable slope variance.
#  - m_fixint: the delta_T x crop_model interaction is significant, but R drops EPIC's
#    slope (single delta_T) and the APSIM slope (-5.4 %/C) comes entirely from one study
#    (XA3Z3WRK), so the interaction is really a single-study/design contrast, not a
#    generalizable crop-model effect.
#  - The steep APSIM vs flat DSSAT difference is real WITHIN XA3Z3WRK (same site,
#    Piracicaba), but rests on that one study -> report descriptively, not as a model term.
# DECISION: keep crop_model as a random INTERCEPT (m_final); do NOT use a random slope
# or a fixed delta_T x crop_model interaction. Pooled slope ~= the DSSAT-CANEGRO slope.

# --- Check: should the delta_T slope vary by STUDY (random slope)? ---
# 9 studies, 8 with >=2 distinct delta_T (only AADRMRIW is single-valued), so better
# informed than the crop_model slope. NOTE: ~ delta_T_C | study_id already supplies the
# study random intercept, so use ~ 1 | sim (not study_id/sim) to avoid double-counting it.
# Watch for a boundary rho (singular fit) and check the LRT / effect on the slope CI.
m_sslope <- tryCatch(
  rma.mv(yi = PCY_pct, V = vi, mods = ~ delta_T_C,
         random = list(~ delta_T_C | study_id, ~ 1 | sim, ~ 1 | crop_model),
         struct = "GEN", data = dat, method = "REML"),
  error = function(e) e)
if (inherits(m_sslope, "error")) {
  cat("m_sslope problem:", conditionMessage(m_sslope), "\n")
} else {
  print(summary(m_sslope))
  print(anova(m_final, m_sslope))
}

# Robust (cluster-robust / RVE) SEs on the random-intercept model as an alternative
# safeguard against within-study dependence without estimating a slope variance.
robust(m_final, cluster = dat$study_id)

# Conclusion of the study-slope check:
#  - m_sslope: random slope by study is NOT supported -- intercept-slope rho = -1.00
#    (singular fit), LRT vs m_final p = 0.28, and AIC worse (348.7 vs 347.3). So, like
#    crop_model, the slope should stay fixed; keep m_final's random-intercept structure.
#  - robust(): cluster-robust SEs (9 study clusters, df = 7) give delta_T_C = -1.29,
#    SE = 0.69, p = 0.105, CI [-2.92, 0.35] -- i.e. NOT significant, vs the model-based
#    p = 0.005. With only 9 studies the pooled temperature effect is suggestive but not
#    robustly significant; report both and temper the conclusion accordingly.

data.frame(model  = c("sim", "sim + site + crop", "study/sim + site + crop", "study/sim + crop"),
           AIC    = c(AIC(m_sim), AIC(m_li), AIC(m_full), AIC(m_final)),
           BIC    = c(BIC(m_sim), BIC(m_li), BIC(m_full), BIC(m_final)),
           logLik = c(logLik(m_sim), logLik(m_li), logLik(m_full), logLik(m_final)),
           row.names = NULL)
confint(m_final, sigma2 = 1)

# --- Normality check: simulation-level (observation) residuals of m_final ---
# Each row is one simulation, so rstandard() returns one standardized residual
# per simulation. No qqnorm() method exists for rma.mv, so build the QQ plot by
# hand and back it up with a Shapiro-Wilk test.
res_z <- rstandard(m_final)$z
fitted_vals <- fitted(m_final)
png("ma_qqresid.png", width = 7, height = 6, units = "in", res = 300)
qqnorm(res_z, pch = 19, main = "Normal Q-Q: m_final standardized residuals")
qqline(res_z, col = "red", lwd = 2)
dev.off()
print(shapiro.test(res_z))

# Residuals vs fitted: check linearity of the delta_T_C relationship and
# homoscedasticity of the simulation-level residuals.
png("ma_residfit.png", width = 7, height = 6, units = "in", res = 300)
plot(fitted_vals, res_z, pch = 19,
     xlab = "Fitted PCY (%)", ylab = "Standardized residual",
     main = "Residuals vs fitted: m_final")
abline(h = 0, col = "grey60", lty = 2)
abline(h = c(-2, 2), col = "grey80", lty = 3)
lines(lowess(fitted_vals, res_z), col = "red", lwd = 2)
dev.off()

# --- Influence diagnostics: spot outlying / influential simulations ---
# No influence() method exists for rma.mv, so compute the pieces directly.
# cooks.distance() refits the model once per simulation (deletion), so it is
# the slow part; rstudent() gives deleted (externally studentized) residuals.
res_stud <- rstudent(m_final)$z
cookd    <- cooks.distance(m_final)
hatv     <- hatvalues(m_final)
idx      <- seq_along(res_stud)
# Common cut-offs: |studentized| > ~1.96 (or Bonferroni), Cook's D > 4/n.
cook_cut <- 4 / length(cookd)
flag     <- which(abs(res_stud) > qnorm(0.975) | cookd > cook_cut)

# Highlight all Marin (2015) / XA3Z3WRK simulations in red; label only the most
# influential one (largest Cook's distance).
out_i  <- which.max(cookd)
out_lab <- sprintf("Marin (2015) | %s | %g \u00b0C", dat$crop_model[out_i], dat$delta_T_C[out_i])
pt_col  <- ifelse(dat$study_id == "XA3Z3WRK", "red", "black")
png("ma_influence.png", width = 8, height = 8, units = "in", res = 300)
op <- par(mfrow = c(3, 1), mar = c(4, 4.5, 1.5, 1))
plot(idx, res_stud, type = "o", pch = 19, col = pt_col, xlab = "", ylab = "Studentized residual")
abline(h = c(-1, 1) * qnorm(0.975), col = "grey70", lty = 3); abline(h = 0, col = "grey85")
plot(idx, cookd, type = "o", pch = 19, col = pt_col, xlab = "", ylab = "Cook's distance")
abline(h = cook_cut, col = "grey70", lty = 3)
text(out_i, cookd[out_i], labels = out_lab, pos = 2, cex = 0.7, col = "red")
plot(idx, hatv, type = "o", pch = 19, col = pt_col, xlab = "Simulation (row order)", ylab = "Leverage (hat)")
abline(h = 2 * mean(hatv), col = "grey70", lty = 3)
par(op)
dev.off()

# Print the flagged simulations for the report.
infl_tab <- data.frame(row = idx, study_id = dat$study_id, crop_model = dat$crop_model,
                       delta_T_C = dat$delta_T_C, PCY_pct = dat$PCY_pct,
                       stud_res = round(res_stud, 2), cooks_d = round(cookd, 3),
                       hat = round(hatv, 3))
cat("\nFlagged influential/outlying simulations (|stud.res|>1.96 or Cook's D >",
    round(cook_cut, 3), "):\n")
print(infl_tab[flag, ], row.names = FALSE)

# --- Sensitivity: refit m_final WITHOUT the influential study XA3Z3WRK (Marin 2015) ---
# This is the study that drives the influence diagnostics (the steep APSIM slope).
# Compare the pooled intercept/slope to gauge how much it depends on that one study.
dat_nox <- droplevels(dat[dat$study_id != "XA3Z3WRK", ])
m_nox <- rma.mv(yi = PCY_pct, V = vi, mods = ~ delta_T_C,
                random = list(~ 1 | study_id/sim, ~ 1 | crop_model),
                data = dat_nox, method = "REML")
cat("\n--- m_final vs m_nox (XA3Z3WRK removed) ---\n")
comp <- data.frame(
  term     = c("intercept", "delta_T_C"),
  full     = round(coef(m_final), 3),
  full_se  = round(m_final$se, 3),
  full_p   = round(m_final$pval, 4),
  nox      = round(coef(m_nox), 3),
  nox_se   = round(m_nox$se, 3),
  nox_p    = round(m_nox$pval, 4))
print(comp, row.names = FALSE)
cat(sprintf("n: full = %d sims / %d studies; nox = %d sims / %d studies\n",
            m_final$k, length(unique(dat$study_id)),
            m_nox$k, length(unique(dat_nox$study_id))))
print(summary(m_nox))


pdf("ma_regplot.pdf", width = 7, height = 6)
regplot(m_final, mod = "delta_T_C", psize = 1, xlab = "Temperature increase (°C)", ylab = "PCY (%)")
dev.off()

# --- Sensitivity regplot: full fit vs fit without XA3Z3WRK ---
# Highlight the XA3Z3WRK simulations and overlay both pooled regression lines so
# the reader can see how much that one study pulls the slope.
is_x <- dat$study_id == "XA3Z3WRK"
pt_col <- ifelse(is_x, "red", "grey50")
# Slope CI half-widths (95%) for the +/- legend labels.
mar_full <- (m_final$ci.ub[2] - m_final$ci.lb[2]) / 2
mar_nox  <- (m_nox$ci.ub[2]   - m_nox$ci.lb[2])   / 2
png("ma_regplot_nox.png", width = 7, height = 6, units = "in", res = 300)
plot(dat$delta_T_C, dat$PCY_pct, pch = 19, col = pt_col,
     xlab = expression("Temperature increase (" * degree * "C)"), ylab = "PCY (%)")
abline(h = 0, col = "grey80", lty = 3)
abline(a = coef(m_final)[1], b = coef(m_final)[2], col = "black", lwd = 3)
abline(a = coef(m_nox)[1],   b = coef(m_nox)[2],   col = "red",   lwd = 3)
legend("bottomleft", bg = "white", cex = 0.8,
       legend = c("Marin 2015 simulations", "other simulations",
                  sprintf("full fit (slope %.2f \u00b1 %.2f)", coef(m_final)[2], mar_full),
                  sprintf("without Marin 2015 (slope %.2f \u00b1 %.2f)", coef(m_nox)[2], mar_nox)),
       pch = c(19, 19, NA, NA), lwd = c(NA, NA, 3, 3),
       col = c("red", "grey50", "black", "red"))
dev.off()

studies <- sort(unique(dat$study_id))
study_labels <- c(
  "65VYUJZ8" = "Verma 2023",
  "7PCFLLSK" = "Guhan 2024",
  "AADRMRIW" = "Sanders 2019",
  "FWLHHVEY" = "Jones 2015",
  "HP54EZKH" = "Jones & Singels 2018",
  "PCXRP7XW" = "Marin 2013",
  "Q6N2E27T" = "Ramachandran 2017",
  "XA3Z3WRK" = "Marin 2015",
  "Z5PF5XNQ" = "Sonkar 2020"
)
cols <- setNames(palette.colors(length(studies), "Tableau 10"), studies)
pdf("ma_studies.pdf", width = 8, height = 6)
plot(dat$delta_T_C, dat$PCY_pct, pch = 19, col = cols[dat$study_id],
     xlab = "Temperature increase (°C)", ylab = "PCY (%)")
abline(h = 0, col = "grey80", lty = 3)
for (s in studies) {
  ds <- dat[dat$study_id == s, ]
  if (length(unique(ds$delta_T_C)) >= 2) abline(lm(PCY_pct ~ delta_T_C, ds), col = cols[s], lwd = 1.5)
}
abline(a = coef(m_final)[1], b = coef(m_final)[2], col = "black", lwd = 3)
legend("bottomleft", legend = study_labels[studies], col = cols[studies],
       pch = 19, lwd = 1.5, cex = 0.7, ncol = 2, bg = "white")
dev.off()

# --- Exploratory scatter (no regression lines) ---
# Show the raw simulations vs temperature change, coloured by study, and how many
# fall above vs below 0 % yield change. No overall or study-specific fit lines.
n_total <- nrow(dat)
n_below <- sum(dat$PCY_pct < 0)
n_above <- sum(dat$PCY_pct > 0)
n_zero  <- sum(dat$PCY_pct == 0)
pct_below <- 100 * n_below / n_total
pct_above <- 100 * n_above / n_total
cat(sprintf("\nSimulations: n=%d | below 0: %d (%.1f%%) | above 0: %d (%.1f%%) | at 0: %d\n",
            n_total, n_below, pct_below, n_above, pct_above, n_zero))

# Simple (unweighted) mean PCY across simulations and its 95% CI.
pcy_tt   <- t.test(dat$PCY_pct)
pcy_mean <- unname(pcy_tt$estimate)
pcy_ci   <- pcy_tt$conf.int
cat(sprintf("Mean PCY = %.2f%% (95%% CI %.2f to %.2f)\n", pcy_mean, pcy_ci[1], pcy_ci[2]))

pdf("ma_explore.pdf", width = 8, height = 6)
plot(dat$delta_T_C, dat$PCY_pct, pch = 19, col = cols[dat$study_id],
     xlab = expression("Temperature increase (" * degree * "C)"), ylab = "PCY (%)")
usr <- par("usr")
rect(usr[1], pcy_ci[1], usr[2], pcy_ci[2], col = rgb(0, 0, 0, 0.08), border = NA)
abline(h = pcy_mean, col = "black", lwd = 2)
abline(h = pcy_ci, col = "black", lwd = 1, lty = 2)
abline(h = 0, col = "grey60", lty = 3, lwd = 1.5)
legend("bottomleft", legend = study_labels[studies], col = cols[studies],
       pch = 19, cex = 0.7, ncol = 2, bg = "white")
legend("topright",
       legend = c(sprintf("Mean PCY = %.1f%%", pcy_mean),
                  sprintf("95%% CI [%.1f, %.1f]", pcy_ci[1], pcy_ci[2])),
       lty = c(1, 2), lwd = c(2, 1), col = "black", cex = 0.7, bg = "white")
dev.off()

crop_models <- sort(unique(dat$crop_model))
cm_cols <- setNames(palette.colors(length(crop_models), "Tableau 10"), crop_models)
pdf("ma_cropmodels.pdf", width = 8, height = 6)
plot(dat$delta_T_C, dat$PCY_pct, pch = 19, col = cm_cols[dat$crop_model],
     xlab = "Temperature increase (°C)", ylab = "PCY (%)")
abline(h = 0, col = "grey80", lty = 3)
for (cm in crop_models) {
  dc <- dat[dat$crop_model == cm, ]
  if (length(unique(dc$delta_T_C)) >= 2) abline(lm(PCY_pct ~ delta_T_C, dc), col = cm_cols[cm], lwd = 1.5)
}
abline(a = coef(m_final)[1], b = coef(m_final)[2], col = "black", lwd = 3)
legend("bottomleft", legend = crop_models, col = cm_cols, pch = 19, lwd = 1.5, cex = 0.7, bg = "white")
dev.off()

for (cm in crop_models) {
  dc <- dat[dat$crop_model == cm, ]
  if (length(unique(dc$delta_T_C)) >= 2) {
    cf <- coef(lm(PCY_pct ~ delta_T_C, dc))
    cat(sprintf("%-22s n=%2d  intercept=%7.2f  slope=%7.3f\n", cm, nrow(dc), cf[1], cf[2]))
  } else {
    cat(sprintf("%-22s n=%2d  (single delta_T = %.1f, no slope)\n", cm, nrow(dc), unique(dc$delta_T_C)))
  }
}
