library(metafor)

dat <- read.csv("ma_data_clean.csv", stringsAsFactors = FALSE)

dat <- dat[dat$include == "yes", ]
dat$PCY_pct   <- as.numeric(dat$PCY_pct)
dat$delta_T_C <- as.numeric(dat$delta_T_C)
dat$study_id <- as.factor(dat$study_id)
dat <- dat[!is.na(dat$PCY_pct) & !is.na(dat$delta_T_C), ]
dat$obs <- factor(seq_len(nrow(dat)))
dat$vi  <- 0

m_base <- rma.mv(yi = PCY_pct, V = vi, mods = ~ delta_T_C,
                 random = list(~ 1 | obs, ~ 1 | crop_model),
                 data = dat, method = "REML")
summary(m_base)

m_site <- rma.mv(yi = PCY_pct, V = vi, mods = ~ delta_T_C,
                 random = list(~ 1 | obs, ~ 1 | site, ~ 1 | crop_model),
                 data = dat, method = "REML")
summary(m_site)
anova(m_base, m_site)

m_ext <- rma.mv(yi = PCY_pct, V = vi, mods = ~ delta_T_C,
                random = list(~ 1 | study_id/obs, ~ 1 | crop_model),
                data = dat, method = "REML")
summary(m_ext)
anova(m_base, m_ext)

data.frame(model  = c("obs + crop", "obs + site + crop", "study/obs + crop"),
           AIC    = c(AIC(m_base), AIC(m_site), AIC(m_ext)),
           BIC    = c(BIC(m_base), BIC(m_site), BIC(m_ext)),
           logLik = c(logLik(m_base), logLik(m_site), logLik(m_ext)),
           row.names = NULL)
confint(m_ext, sigma2 = 1)


pdf("ma_regplot.pdf", width = 7, height = 6)
regplot(m_ext, mod = "delta_T_C", psize = 1, xlab = "Temperature increase (°C)", ylab = "PCY (%)")
dev.off()

studies <- sort(unique(dat$study_id))
cols <- setNames(palette.colors(length(studies), "Tableau 10"), studies)
pdf("ma_studies.pdf", width = 8, height = 6)
plot(dat$delta_T_C, dat$PCY_pct, pch = 19, col = cols[dat$study_id],
     xlab = "Temperature increase (°C)", ylab = "PCY (%)")
abline(h = 0, col = "grey80", lty = 3)
for (s in studies) {
  ds <- dat[dat$study_id == s, ]
  if (length(unique(ds$delta_T_C)) >= 2) abline(lm(PCY_pct ~ delta_T_C, ds), col = cols[s], lwd = 1.5)
}
abline(a = coef(m_ext)[1], b = coef(m_ext)[2], col = "black", lwd = 3)
legend("bottomleft", legend = studies, col = cols, pch = 19, lwd = 1.5, cex = 0.7, ncol = 2, bg = "white")
dev.off()
