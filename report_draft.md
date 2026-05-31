# Statistical Analysis, Results, and Discussion

---

## Statistical Analysis

The outcome of interest is the percent change in sugarcane yield (PCY, %) relative to each study's baseline period, and the predictor is the temperature change applied in each simulation (ΔT, °C). After excluding observations with missing PCY or ΔT, the final dataset comprised 45 simulation observations drawn from 9 published studies.

**Sampling variance.** The primary studies are deterministic crop-model simulations, each reporting a single projected yield change per site–scenario combination without a standard error or replicated measurements. An observation-level sampling variance cannot therefore be computed from the available data. Sampling variance was accordingly set to zero (V = 0), which treats each simulation as a point estimate and gives every simulation equal weight; all heterogeneity is consequently absorbed by the random-effects variance components described below.

**Model specification.** A multilevel meta-regression was fitted using `metafor::rma.mv` with REML estimation. Temperature change (ΔT) entered as the sole fixed moderator. Four random-effects structures were compared in a sequential ladder:

1. **m_obs** (`~ 1 | obs`): a single-level model treating all 45 simulations as independent. This serves as a baseline reference that ignores the known clustering of multiple estimates within each study; it is not a candidate inferential model.

2. **m_li** (`~ 1 | obs` + `~ 1 | site` + `~ 1 | crop_model`): a direct replication of the random-effects structure used by Li et al. (2021) for their cotton meta-analysis, which included study site and crop model as crossed random intercepts. Li's two additional random effects — GCM/climate model and climate scenario — could not be replicated here because approximately half of the sugarcane studies are one-at-a-time sensitivity analyses that apply prescribed temperature increments without using a GCM or an emission scenario.

3. **m_full** (`~ 1 | study_id/obs` + `~ 1 | site` + `~ 1 | crop_model`): extends m_li by adding a hierarchical study-level random intercept. Multiple simulation outputs drawn from the same publication are not independent — they share authors, study design, site selection, and baseline climatology — yet this within-study dependence was not explicitly modelled by Li et al. (2021). The nested term `study_id/obs` partitions variance into a between-study component and a within-study residual.

4. **m_final** (`~ 1 | study_id/obs` + `~ 1 | crop_model`): obtained by dropping `site` from m_full after its variance component estimated to effectively zero in every model that included it alongside `study_id`, indicating complete redundancy once study clustering was accounted for.

**Slope heterogeneity checks.** Allowing the temperature slope to vary across crop models (random slope or fixed ΔT × crop_model interaction) and across studies (random slope by study) was also explored. Both were rejected: in each case the estimated intercept–slope correlation reached the boundary of the parameter space (ρ = ±1.00), indicating a singular, over-parameterised fit. The practical cause is sparse data — APSIM-Sugar and EPIC each appear in a single study (EPIC at only one ΔT value) — so the slope cannot be reliably distinguished from study and design effects. The temperature slope was therefore kept fixed across all groups, as in Li et al. (2021).

Models were compared using nested likelihood-ratio tests (LRT) where applicable, and by AIC and BIC across the full ladder.

---

## Results

**Model selection.** Table 1 summarises fit statistics for the four models. Adding site and crop-model random intercepts to the null (m_obs → m_li) significantly improved fit (LRT p = 0.015, ΔAIC = −4.4). Introducing the hierarchical study-level intercept (m_li → m_full) yielded a further significant improvement (LRT p = 0.006, ΔAIC = −5.5), confirming that within-study non-independence is a meaningful source of heterogeneity that the Li-style structure does not capture. Dropping site from m_full did not alter the fit at all (LRT p = 1.00), and the parsimonious m_final (`study_id/obs` + `crop_model`) achieved the lowest AIC and BIC of all four models.

| Model | Random structure | AIC | BIC | log-lik |
|---|---|---|---|---|
| m_obs | obs | 359.2 | 364.4 | −176.6 |
| m_li | obs + site + crop_model | 354.8 | 363.6 | −172.4 |
| m_full | study/obs + site + crop_model | 349.3 | 359.8 | −168.6 |
| **m_final** | **study/obs + crop_model** | **347.3** | **356.1** | **−168.6** |

*Table 1. Model fit statistics. m_final (bold) is the preferred model.*

**Variance components.** In m_final, the estimated standard deviations are 8.1 PCY percentage points at the study level, 10.4 at the within-study residual level, and 14.2 across crop models, reflecting substantial heterogeneity at all levels. The large crop-model component mirrors the well-documented sensitivity of sugarcane yield projections to the choice of simulation model (Figure 3).

**Effect of temperature on sugarcane yield.** m_final estimates a pooled slope of **−1.29% PCY per +1 °C** (95% CI: −2.19 to −0.38, p = 0.005), indicating that each additional degree of warming reduces simulated sugarcane yield by approximately 1.3 percentage points. Notably, accounting for study-level dependence sharpened this estimate considerably: in m_li, which ignores within-study clustering, the slope was −1.07%/°C and statistically marginal (p = 0.09). The improvement demonstrates that the extra variance soaked up by the study intercept in m_final reduces noise around the fixed-effect estimate.

Figure 1 shows the fitted m_final regression line overlaid on all 45 observations. The negative relationship is consistent across the ΔT range covered by the corpus (−3 °C to +9 °C), though with wide scatter. Figure 2 displays individual OLS regression lines per study, illustrating the substantial between-study heterogeneity in both intercept and slope that the random-effects structure is designed to accommodate. Figure 3 presents per-crop-model lines: APSIM-Sugar shows a markedly steeper decline (~−5.4%/°C) compared to the near-flat DSSAT-CANEGRO response (~−0.3%/°C per OLS), though the former rests on a single study at one site (Piracicaba, Brazil) and should not be interpreted as a generalised crop-model effect.

---

## Model validation & sensitivity analysis

A normal quantile–quantile plot of the simulation-level standardized residuals from m_final (Figure 4) shows the points falling close to the reference line across the full range, with only minor deviation near the centre and well-behaved tails. A Shapiro–Wilk test was consistent with normality (W = 0.98, p = 0.70), so the model's residual assumption is well supported.

Per-simulation influence measures (standardized deleted residuals, Cook's distance, and leverage; Figure 5) identified simulation 34 — an APSIM-Sugar run at ΔT = −3 °C from Marin et al. (2015) (study XA3Z3WRK), which projected a +17% yield change — as a clear outlier. Its Cook's distance (≈1.6) is an order of magnitude larger than that of any other simulation, combining high leverage with a large residual. The remaining high-influence points also belong to this study, reflecting that APSIM-Sugar appears only within Marin et al. (2015) and produces the steepest temperature response in the corpus.

Because this study dominated the diagnostics, m_final was refitted without it. The pooled temperature slope became markedly less negative, from −1.29%/°C (95% CI: −2.19 to −0.38, p = 0.005) to −0.78%/°C (95% CI: −1.66 to 0.10, p = 0.082) — no longer significant at the 5% level, with the interval now spanning zero (Figure 6). The significance of the pooled effect thus leans substantially on a single study, echoing the cluster-robust analysis (p = 0.11). The direction is unchanged, but this must be weighed against the loss of information: removal reduces the dataset from 45 to 37 simulations and from 9 to 8 studies, and eliminates the entire APSIM-Sugar level. The full-data estimate is therefore retained as primary.

---

## Discussion

The pooled estimate of −1.3%/°C is broadly consistent with Li et al. (2021), who found a −1.64%/°C response for cotton using the same meta-regression approach. The milder sensitivity found here is ecologically plausible: sugarcane is a C4 crop with a higher temperature optimum than the C3 cotton modelled by Li et al., and C4 species are generally less responsive to moderate warming (Sage & Kubien, 2007). Notably, Li et al.'s restricted-dataset estimate (−7.79%/°C) is far steeper, but that model conditioned on complete covariate records — a restriction that would reduce the present dataset to an unusably small subset.

A methodological improvement over Li et al. (2021) is the explicit publication-level random intercept (`study_id`), which their specification omitted. The m_li → m_full LRT (p = 0.006) confirms that within-study non-independence is a meaningful source of heterogeneity; accounting for it sharpened the temperature slope from −1.07%/°C (p = 0.09) to −1.29%/°C (p = 0.005).

**Limitations.** The analysis is restricted to temperature change as the sole moderator because it is the only variable reported consistently across all nine studies. The additional drivers modelled by Li et al. — precipitation change, CO₂ concentration, and adaptation measures — are available in only a subset of the sugarcane papers and could not be included without an unacceptable loss of observations. Future work should systematically extract these variables from the existing studies where they are reported, and prioritise new sugarcane simulation studies that provide the full set of climate-driver information, enabling a richer multi-moderator model and more nuanced conclusions about the drivers of sugarcane yield change. The nine included studies are also geographically concentrated in a small number of regions (South Asia, Brazil, and southern Africa), so the pooled estimate may not generalise to other major sugarcane-producing areas. Finally, potential publication bias cannot be ruled out: researchers may be more likely to publish simulation results that show large yield changes, and because crop-model simulations do not produce standard errors or replicated measurements, it is hard to quantify this bias using tools like a funnel plot.
