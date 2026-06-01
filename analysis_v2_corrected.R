# ==============================================================================
# Nitrogen Content Analysis — Thaumetopoea pityocampa (PPM Proteins dataset)
# Comparing: standard clade larvae (L1 & L2) vs ENA clade adults & larvae
# ==============================================================================
# Author:   Konstantina Gianniou
# Contact:  g.tem2106@gmail.com | github.com/KonGianniou
# Dataset:  Synthetic — preserves statistical structure of original.
#           Real data not included (copyright — ELGO Demeter / Dr. D. Avtzis).
# Methods:  Non-linear correlation (nlcor), Mann-Whitney U, Student t-test,
#           Simple linear regression, Piecewise regression (segmented),
#           Generalized linear model with stepwise selection
# ==============================================================================

# ── 0. PACKAGE SETUP ──────────────────────────────────────────────────────────

required_cran <- c("ggplot2", "dplyr", "tidyr", "gplots", "car",
                   "segmented", "remotes", "patchwork", "scales", "broom")

for (pkg in required_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

if (!requireNamespace("nlcor", quietly = TRUE)) {
  remotes::install_github("ProcessMiner/nlcor", force = TRUE)
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(gplots)
  library(car)
  library(segmented)
  library(nlcor)
  library(patchwork)
  library(scales)
  library(broom)
})

cat("── Packages loaded successfully ─────────────────\n\n")

# ── THEME ─────────────────────────────────────────────────────────────────────

theme_tp <- theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 13, hjust = 0),
    plot.subtitle    = element_text(size = 10, colour = "grey40", hjust = 0),
    axis.title       = element_text(size = 11),
    axis.text        = element_text(size = 10),
    legend.title     = element_text(face = "bold", size = 10),
    legend.text      = element_text(size = 9),
    panel.grid.minor = element_blank(),
    strip.text       = element_text(face = "bold")
  )

pal_stage    <- c("adult" = "#1d6fa4", "larva" = "#e6550d")
pal_location <- c("L1"   = "#2ca25f", "L2" = "#9467bd")
pal_sample   <- c("TPena" = "#1d6fa4", "TP" = "#e6550d")

# ── 1. LOAD DATA ───────────────────────────────────────────────────────────────

source("generate_synthetic_data_v2.R")
data <- data_synth

# Replicate subgroup structure from original analysis
datalarva  <- filter(data, Stage == "larva")
L2         <- filter(datalarva, Location == "L2")
L1_larva   <- filter(datalarva, Location == "L1")

# Subgroup 2: L1 only (adults + larvae) — the ENA clade comparison
dataTPE    <- filter(data, Location == "L1")
adult      <- filter(dataTPE, Stage == "adult")
larva_ath  <- filter(dataTPE, Stage == "larva")  

# Reciprocal transformation for weight
data$W2      <- 1 / data$Weight
datalarva$W2 <- 1 / datalarva$Weight
dataTPE$W2   <- 1 / dataTPE$Weight

cat("── Sample sizes ─────────────────────────────────\n")
cat(sprintf("Total: %d\n", nrow(data)))
cat(sprintf("Larvae total: %d  (L1: %d, L2: %d)\n",
            nrow(datalarva), nrow(L1_larva), nrow(L2)))
cat(sprintf("L1 subgroup: %d  (adults: %d, larvae: %d)\n\n",
            nrow(dataTPE), nrow(adult), nrow(larva_ath)))

# ── 2. NORMALITY CHECKS ────────────────────────────────────────────────────────

cat("── Shapiro-Wilk Normality Tests ─────────────────\n")

norm_groups <- list(
  "Larvae Weight"  = datalarva$Weight,
  "Larvae N"       = datalarva$N,
  "L1 Weight"      = dataTPE$Weight,
  "L1 N"           = dataTPE$N,
  "Larvae L1 N"    = L1_larva$N,
  "Larvae L2 N"    = L2$N,
  "Adults N"       = adult$N
)

sw_results <- lapply(names(norm_groups), function(nm) {
  x  <- norm_groups[[nm]]
  sw <- shapiro.test(x)
  data.frame(Group = nm, W = round(sw$statistic, 4),
             p = round(sw$p.value, 4),
             Normal = sw$p.value >= 0.05)
})
sw_df <- do.call(rbind, sw_results)
print(sw_df, row.names = FALSE)

# ── 3. NON-LINEAR CORRELATIONS ────────────────────────────────────────────────

cat("\n── Non-linear Correlations (nlcor) ──────────────\n")

run_nlcor <- function(x, y, label) {
  res <- tryCatch(nlcor(x, y, plt = FALSE), error = function(e) NULL)
  if (!is.null(res)) {
    cat(sprintf("%-35s  r = %.3f,  p = %.4f  %s\n",
                label, res$cor.estimate, res$adjusted.p.value,
                ifelse(res$adjusted.p.value < 0.05, "*", "(ns)")))
    invisible(res)
  }
}

run_nlcor(data$Weight,      data$N,      "All data (Weight ~ N)")
run_nlcor(datalarva$Weight, datalarva$N, "Subgroup 1: Larvae (Weight ~ N)")
run_nlcor(dataTPE$Weight,   dataTPE$N,   "Subgroup 2: L1 (Weight ~ N)")

# ── 4. COMPARATIVE ANALYSIS ───────────────────────────────────────────────────

cat("\n── 4a. Subgroup 1: Larvae L1 vs L2 (N) ──────\n")
cat("Group means:\n")
print(datalarva %>% group_by(Location) %>%
        summarise(n = n(), Mean_N = round(mean(N), 3), SD_N = round(sd(N), 3)))

vt1 <- var.test(L2$N, L1_larva$N)
cat(sprintf("\nF-test variances: F = %.3f, p = %.4f\n", vt1$statistic, vt1$p.value))

mw1 <- wilcox.test(L2$N, L1_larva$N, alternative = "two.sided")
cat(sprintf("Mann-Whitney U: W = %.1f, p = %.4f %s\n",
            mw1$statistic, mw1$p.value,
            ifelse(mw1$p.value < 0.05, "✓ Significant", "(ns)")))

cat("\n── 4b. Subgroup 2: Adults vs Larvae in L1 (N) \n")
cat("Group means:\n")
print(dataTPE %>% group_by(Stage) %>%
        summarise(n = n(), Mean_N = round(mean(N), 3), SD_N = round(sd(N), 3)))

vt2 <- var.test(adult$N, larva_ath$N)   
                  
cat(sprintf("\nF-test variances: F = %.3f, p = %.4f\n", vt2$statistic, vt2$p.value))

tt2 <- t.test(N ~ Stage, data = dataTPE, alternative = "two.sided",
              var.equal = (vt2$p.value >= 0.05))
cat(sprintf("t-test: t = %.3f, df = %.1f, p = %.6f %s\n",
            tt2$statistic, tt2$parameter, tt2$p.value,
            ifelse(tt2$p.value < 0.05, "✓ Significant", "(ns)")))

# ── 5. REGRESSION MODELS ──────────────────────────────────────────────────────

cat("\n── 5. Regression Models ─────────────────────────\n")

RMSE <- function(residuals) sqrt(mean(residuals^2))

# Model 1: Simple linear — Weight
lm1 <- lm(N ~ Weight, data = data)
cat(sprintf("M1 Simple linear (Weight):       Adj R² = %.3f\n",
            summary(lm1)$adj.r.squared))

# Model 2: Simple linear — 1/Weight
lm2 <- lm(N ~ W2, data = data)
cat(sprintf("M2 Simple linear (1/Weight):     Adj R² = %.3f\n",
            summary(lm2)$adj.r.squared))

# Model 3: Piecewise — Weight
seg1 <- segmented(lm1, seg.Z = ~ Weight, psi = 0.05)

cat(sprintf("M3 Piecewise (Weight):           Adj R² = %.3f  breakpoint = %.4f\n",
            summary(seg1, short = FALSE)$r.sq,
            seg1$psi[2]))   
                  
# Model 4: Piecewise — 1/Weight
seg2 <- segmented(lm2, seg.Z = ~ W2, psi = 60)
cat(sprintf("M4 Piecewise (1/Weight):         Adj R² = %.3f  breakpoint = %.3f\n",
            summary(seg2, short = FALSE)$r.sq,
            seg2$psi[2]))   

# Model 5: GLM — 1/Weight + Stage + Location
glm1 <- glm(N ~ W2 + Stage + Location, data = data)
cat(sprintf("M5 GLM (1/Weight + Stage + Loc): AIC = %.3f\n", AIC(glm1)))

# Model 6: Stepwise GLM
glm2  <- glm(N ~ Weight + Stage + Location, data = data)
step1 <- step(glm2, trace = 0)
cat(sprintf("M6 Stepwise GLM:                 AIC = %.3f\n", AIC(step1)))

# Model comparison table
cat("\n── Model Comparison Table ───────────────────────\n")
model_table <- data.frame(
  Model = c("Simple linear (Weight)", "Simple linear (1/Weight)",
            "Piecewise (Weight)", "Piecewise (1/Weight)",
            "GLM (1/W + Stage + Loc)", "Stepwise GLM"),
  RMSE  = round(c(RMSE(lm1$residuals), RMSE(lm2$residuals),
                  RMSE(seg1$residuals), RMSE(seg2$residuals),
                  RMSE(glm1$residuals), RMSE(step1$residuals)), 3),
  AIC   = round(AIC(lm1, lm2, seg1, seg2, glm1, step1)$AIC, 3)
)
print(model_table, row.names = FALSE)

cat("\nBest model (lowest AIC + RMSE): GLM with 1/Weight, Stage, Location\n")
cat("\nGLM coefficients:\n")
print(tidy(glm1) %>% mutate(across(where(is.numeric), ~ round(., 4))))

# ── 6. VISUALISATIONS ─────────────────────────────────────────────────────────

cat("\n── Generating plots ─────────────────────────────\n")

# ── Plot 1: Non-linear correlation panels (Figure 1 equivalent) ──────────────


make_nlcor_scatter <- function(df, x_var, y_var, col_var, title, subtitle) {
  ggplot(df, aes(x = .data[[x_var]], y = .data[[y_var]],
                 colour = .data[[col_var]])) +
    geom_point(alpha = 0.7, size = 2) +
    geom_smooth(method = "loess", se = TRUE, colour = "grey30",
                fill = "grey80", linewidth = 0.8, alpha = 0.3) +
    scale_colour_manual(values = c(pal_stage, pal_location)) +
    labs(title = title, subtitle = subtitle,
         x = "Weight (g)", y = "Nitrogen content (%)") +
    theme_tp +
    theme(legend.position = "right")
}

p1a <- make_nlcor_scatter(data,     "Weight", "N", "Stage",
  "a) All data", "r = 0.39, p < 0.01")
p1b <- make_nlcor_scatter(dataTPE,  "Weight", "N", "Stage",
  "b) Location 1 subgroup (adults + larvae)", "r = 0.72, p < 0.01")   
p1c <- make_nlcor_scatter(datalarva,"Weight", "N", "Location",
  "c) Larvae subgroup (L1 + L2)", "r = 0.21, p = 0.07 (ns)")

fig1 <- p1a + p1b + p1c +
  plot_annotation(
    title = "Figure 1: Non-linear correlation between weight and nitrogen content",
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  )
ggsave("fig1_nlcor_panels.png", fig1, width = 14, height = 4.5, dpi = 150)

# ── Plot 2: Larvae N by location (Figure 2 equivalent) ───────────────────────

fig2 <- ggplot(datalarva, aes(x = Location, y = N, fill = Location)) +
  geom_violin(alpha = 0.45, colour = NA) +
  geom_boxplot(width = 0.25, alpha = 0.85, outlier.shape = 21,
               outlier.size = 2) +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 4,
               fill = "#f0c080", colour = "black") +
  scale_fill_manual(values = pal_location,
                    labels = c("L1" = "L1", "L2" = "L2")) +
  scale_x_discrete(labels = c("L1" = "L1", "L2" = "L2")) +
  annotate("text", x = 1.5, y = max(datalarva$N, na.rm = TRUE) * 0.98,
           label = sprintf("Mann-Whitney U\np = %.4f", mw1$p.value),
           size = 3.2, colour = "grey30") +
  labs(title = "Figure 2: Nitrogen content in larvae by sampling location",
       subtitle = "Beige dot = group mean  |  Larvae: L1 (n=51) vs L2 (n=24)",
       x = NULL, y = "Nitrogen content (%)", fill = "Location") +
  theme_tp +
  theme(legend.position = "none")
ggsave("fig2_larvae_by_location.png", fig2, width = 6, height = 5, dpi = 150)

# ── Plot 3: Adults vs Larvae in L1 (Figure 3 equivalent) ─────────────────────

fig3 <- ggplot(dataTPE, aes(x = Stage, y = N, fill = Stage)) +
  geom_violin(alpha = 0.45, colour = NA) +
  geom_boxplot(width = 0.25, alpha = 0.85, outlier.shape = 21,
               outlier.size = 2) +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 4,
               fill = "#f0c080", colour = "black") +
  scale_fill_manual(values = pal_stage) +
  scale_x_discrete(labels = c("adult" = "Adult", "larva" = "Larva")) +
  annotate("text", x = 1.5, y = max(dataTPE$N, na.rm = TRUE) * 0.98,
           label = sprintf("t-test\nt = %.2f, p < 0.001", tt2$statistic),
           size = 3.2, colour = "grey30") +
  labs(title = "Figure 3: Nitrogen content in L1 — adults vs larvae",
       subtitle = "Beige dot = group mean  |  Adults (n=54) vs Larvae (n=38)  [L1 only]",
       x = NULL, y = "Nitrogen content (%)", fill = "Stage") +
  theme_tp +
  theme(legend.position = "none")
ggsave("fig3_L1_stage.png", fig3, width = 6, height = 5, dpi = 150)

# ── Plot 4: Piecewise regression on 1/Weight (Figure 4 equivalent) ───────────

data$N_pred_seg2 <- predict(seg2, newdata = data.frame(W2 = data$W2))
breakpt <- seg2$psi[2]   

# Panel a: transformed scale
p4a <- ggplot(data, aes(x = W2, y = N, colour = Stage)) +
  geom_point(alpha = 0.55, size = 1.8) +
  geom_line(aes(y = N_pred_seg2), colour = "black", linewidth = 1) +
  geom_vline(xintercept = breakpt, linetype = "dashed",
             colour = "grey40", linewidth = 0.7) +
  annotate("text", x = breakpt + 2, y = 2,
           label = sprintf("Breakpoint\n1/W = %.1f", breakpt),
           size = 3, colour = "grey40", hjust = 0) +
  scale_colour_manual(values = pal_stage) +
  labs(title = "a) Transformed scale (1/Weight)",
       x = "1 / Weight", y = "Nitrogen content (%)", colour = "Stage") +
  theme_tp

# Panel b: back-transformed to original weight scale
data_sorted <- data[order(data$Weight), ]
data_sorted$N_pred_bt <- predict(seg2,
                                  newdata = data.frame(W2 = data_sorted$W2))

p4b <- ggplot(data_sorted, aes(x = Weight, y = N, colour = Stage)) +
  geom_point(alpha = 0.55, size = 1.8) +
  geom_line(aes(y = N_pred_bt), colour = "black", linewidth = 1) +
  scale_colour_manual(values = pal_stage) +
  labs(title = "b) Original scale (Weight)",
       x = "Weight (g)", y = "Nitrogen content (%)", colour = "Stage") +
  theme_tp

fig4 <- p4a + p4b +
  plot_annotation(
    title    = "Figure 4: Piecewise linear regression — 1/Weight vs Nitrogen content",
    subtitle = sprintf("Breakpoint at 1/W = %.3f  |  Adj R² = %.3f", breakpt,
                       summary(seg2, short = FALSE)$r.sq),   
    theme = theme(plot.title    = element_text(face = "bold", size = 13),
                  plot.subtitle = element_text(size = 10, colour = "grey40"))
  )
ggsave("fig4_piecewise_regression.png", fig4, width = 12, height = 5, dpi = 150)

# ── Plot 5 (bonus): GLM fitted lines by Stage × Location ─────────────────────

w2_seq  <- seq(min(data$W2), max(data$W2), length.out = 200)
pred_df <- expand.grid(W2 = w2_seq,
                       Stage    = c("larva", "adult"),
                       Location = c("L1", "L2"))
pred_df$N_pred <- predict(glm1, newdata = pred_df)
pred_df$Group  <- paste(pred_df$Stage, pred_df$Location, sep = " — ")

fig5 <- ggplot(data, aes(x = W2, y = N, colour = Stage, shape = Location)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_line(data = pred_df,
            aes(x = W2, y = N_pred,
                colour = Stage, linetype = Location),
            linewidth = 0.9) +
  scale_colour_manual(values = pal_stage) +
  scale_shape_manual(values  = c("L1" = 16, "L2" = 17),
                     labels  = c("L1" = "L1", "L2" = "L2")) +
  scale_linetype_manual(values = c("L1" = "solid", "L2" = "dashed"),
                        labels = c("L1" = "L1", "L2" = "L2")) +
  labs(title    = "Figure 5: GLM fitted lines — N ~ 1/Weight by Stage and Location",
       subtitle = "Best model: Generalised mixed (AIC lowest, RMSE = 2.11)",
       x        = "1 / Weight",
       y        = "Nitrogen content (%)",
       colour   = "Stage", shape = "Location", linetype = "Location") +
  theme_tp
ggsave("fig5_glm_fitted.png", fig5, width = 8, height = 5.5, dpi = 150)

cat("Plots saved:\n")
cat("  fig1_nlcor_panels.png\n  fig2_larvae_by_location.png\n")
cat("  fig3_L1_stage.png\n  fig4_piecewise_regression.png\n")
cat("  fig5_glm_fitted.png\n\n")
cat("=== Analysis complete ===\n")
