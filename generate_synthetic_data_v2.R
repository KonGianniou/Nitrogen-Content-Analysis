# ==============================================================================
# Synthetic Data Generator — Thaumetopoea pityocampa Nitrogen Study
# Full dataset: PPM Proteins analysis
# ==============================================================================
# Calibrated to the exact statistical structure of the original study.
# Real measurements are not included (copyright — ELGO Demeter / Dr. D. Avtzis).
# The synthetic data preserves: sample sizes per group, group means & SDs,
# weight ranges, CP~N correlation structure (~perfect: tau=0.999),
# and the bimodal weight distribution of the ENA adult group.
# ==============================================================================

set.seed(2022)

generate_tp_data <- function() {

  rows <- list()

  # ── Helper: simulate a group ──────────────────────────────────────────────
  add_group <- function(sample_label, stage, location,
                        n, w_mean, w_sd, w_min, w_max,
                        N_mean, N_sd, code_prefix, has_CP = TRUE) {

    # Weight: truncated normal
    W <- pmax(w_min, pmin(w_max,
          rnorm(n, w_mean, w_sd)))

    # N: based on observed mean/SD, with small individual noise
    N_base <- rnorm(n, N_mean, N_sd)

    # CP = N * 6.25 (standard nitrogen-to-protein conversion) + tiny noise
    CP <- if (has_CP) N_base * 6.25 + rnorm(n, 0, 0.05) else rep(NA_real_, n)

    for (i in seq_len(n)) {
      rows[[length(rows) + 1]] <<- data.frame(
        Sample   = sample_label,
        Stage    = stage,
        Location = location,
        Weight   = round(W[i], 4),
        CP       = if (has_CP) round(CP[i], 3) else NA_real_,
        N        = round(N_base[i], 3),
        Code     = paste0(code_prefix, i),
        stringsAsFactors = FALSE
      )
    }
  }

  # ── GROUP DEFINITIONS (calibrated to original summary stats) ──────────────
  #
  # Tp standard larvae — L1
  # n=13, weight range ~0.05–0.20, N ~9.8 (but two outliers ~1.6 & 2.4 visible)
  # We model the main cluster (n=11) and two low-N outliers separately
  add_group("T. pityocampa pure", "larva", "L1",
            n=11, w_mean=0.118, w_sd=0.050, w_min=0.050, w_max=0.210,
            N_mean=11.25, N_sd=0.90, code_prefix="Tplarva")

  # Two low-N outlier larvae (L1) — likely contaminated/moulting samples
  rows[[length(rows) + 1]] <- data.frame(
    Sample="T. pityocampa pure", Stage="larva", Location="L1",
    Weight=0.1827, CP=15.061, N=2.410, Code="Tplarva_outL1",
    stringsAsFactors=FALSE)
  rows[[length(rows) + 1]] <- data.frame(
    Sample="T. pityocampa pure", Stage="larva", Location="L1",
    Weight=0.0549, CP=10.001, N=1.600, Code="Tplarva_outL2",
    stringsAsFactors=FALSE)

  # Tp standard larvae — L2
  # n=19+5=24 (early batches + κάμπιες), weight ~0.011–0.020, N ~10.26, SD=2.46
  add_group("T. pityocampa", "larva", "L2",
            n=5, w_mean=0.052, w_sd=0.010, w_min=0.035, w_max=0.070,
            N_mean=10.06, N_sd=0.50, code_prefix="Tplarva_LA")

  add_group("T. pityocampa (kampia)", "larva", "L2",
            n=19, w_mean=0.016, w_sd=0.003, w_min=0.011, w_max=0.022,
            N_mean=10.38, N_sd=1.80, code_prefix="Tplarva_LB")

  # TPena larvae — L1 (ENA clade)
  # n=38, weight ~0.011–0.017, N mean=11.15, SD=1.31
  # No CP values recorded (NA in original)
  add_group("TPena", "larva", "L1",
            n=38, w_mean=0.0138, w_sd=0.0015, w_min=0.0110, w_max=0.0170,
            N_mean=11.15, N_sd=1.31, code_prefix="TpENAlarva", has_CP=FALSE)

  # TPena adults — L1 (ENA clade)
  # n=54, weight ~0.010–0.013, N mean=20.06, SD=0.95
  # Very tight weight range, very high N — the key biological signal
  add_group("TPena", "adult", "L1",
            n=8, w_mean=0.0645, w_sd=0.006, w_min=0.054, w_max=0.072,
            N_mean=9.53, N_sd=0.65, code_prefix="TpENAadult_std")

  add_group("TPena", "adult", "L1",
            n=46, w_mean=0.01125, w_sd=0.0008, w_min=0.0100, w_max=0.0132,
            N_mean=20.06, N_sd=0.95, code_prefix="TpENAadult")

  do.call(rbind, rows)
}

data_synth <- generate_tp_data()

# Standardise location labels
data_synth$Location[data_synth$Location == "L1"] <- "L1"
data_synth$Location[data_synth$Location == "L2"] <- "L2"

# Add Sample column used in original script
data_synth$Sample_group <- ifelse(
  grepl("TPena", data_synth$Code) & data_synth$Stage == "adult", "TPena",
  ifelse(grepl("TPena", data_synth$Code) & data_synth$Stage == "larva", "TPena",
         "TP"))

cat("=== Synthetic Data Verification ===\n\n")

cat("N by Stage x Location:\n")
print(aggregate(N ~ Stage + Location, data_synth,
  FUN = function(x) c(n=length(x), mean=round(mean(x),2), sd=round(sd(x),2))))

cat("\nWeight by Stage x Location:\n")
print(aggregate(Weight ~ Stage + Location, data_synth,
  FUN = function(x) c(n=length(x), mean=round(mean(x),4), sd=round(sd(x),4))))

write.csv(data_synth, "Data_TP_synthetic.csv", row.names = FALSE)
cat("\nSaved: Data_TP_synthetic.csv  |  Total rows:", nrow(data_synth), "\n")
