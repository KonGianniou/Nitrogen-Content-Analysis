# Nitrogen Content Analysis — *Thaumetopoea pityocampa*

**Author:** Konstantina Gianniou  
**Affiliation:** ELGO Demeter – Forest Research Institute (collaboration with Dr. D. Avtzis)  
**Status:** Manuscript submitted for peer-reviewed publication

> ⚠️ \\\*\\\*Public Repository Disclaimer\\\*\\\*
>
> This repository is a \\\*\\\*public-facing, adapted version\\\*\\\* of the analysis accompanying a manuscript submitted for peer-reviewed publication.
>
> Two deliberate modifications have been made relative to the submitted work:
>
> 1. \\\*\\\*Synthetic data.\\\*\\\* The original measurements are the property of ELGO Demeter / Forest Research Institute and cannot be shared publicly. The dataset included here (`generate\\\_synthetic\\\_data.R`) is synthetic: it preserves the statistical structure of the real data (sample sizes, group means, standard deviations, distributional shape) but contains no actual measurements. All results produced by running this repository will therefore differ numerically from those reported in the manuscript.
>
> 2. \\\*\\\*Modified analysis code.\\\*\\\* The R script (`analysis.R`) has been rewritten relative to the version used in the submitted manuscript in order to avoid any copyright concerns. The statistical approach, choice of tests, and interpretation are identical; only the specific implementation details have been changed.
>
> If you are a reviewer or editor and require access to the original data and code, please contact the author directly.

\---

## Overview

This repository contains the full statistical pipeline for a comparative entomological study examining nitrogen content in larvae and adults of *Thaumetopoea pityocampa* (pine processionary moth) and the ENA clade, across two Greek sampling locations (L1 and L2).

The study asks:

* Does nitrogen content differ significantly between developmental stages (larva vs. adult)?
* Does sampling location (L1 vs. L2) affect nitrogen content?
* Is there a non-linear relationship between specimen weight and nitrogen content?

\---

## Repository Structure

```
tpityocampa-nitrogen/
├── generate\\\_synthetic\\\_data.R   # Generates synthetic dataset (see Data note)
├── analysis.R                  # Full statistical analysis pipeline
├── plot\\\_N\\\_by\\\_stage.png         # Output: violin-boxplot by stage
├── plot\\\_N\\\_stage\\\_x\\\_location.png # Output: violin-boxplot by stage × location
├── plot\\\_N\\\_by\\\_location.png      # Output: violin-boxplot by location
└── README.md
```

\---

## Data Note

> \\\*\\\*Real data are not included in this repository\\\*\\\* as they are the property of ELGO Demeter / Forest Research Institute. The file `generate\\\_synthetic\\\_data.R` produces a synthetic dataset that preserves the statistical structure of the original (sample sizes, group means, standard deviations, distributional shape) without exposing any actual measurements. All analysis code runs on this synthetic data and is fully reproducible.

\---

## Methods

|Step|Method|R package|
|-|-|-|
|Normality assessment|Shapiro-Wilk test|`stats`|
|Non-linear correlation (Weight \~ N)|nlcor|`nlcor`|
|Stage comparison (larvae vs adults)|Mann-Whitney U|`stats`|
|Location comparison|Mann-Whitney U|`stats`|
|Adults Ath vs Thess (small n)|Permutation test|`coin`|
|Visualisation|Violin + boxplot|`ggplot2`|

Non-parametric tests were chosen throughout due to the non-normal distribution of nitrogen values (confirmed by Shapiro-Wilk, p < 0.05 for most subgroups). For the adult L2 subgroup (n = 18), a permutation test was used instead of Mann-Whitney U due to small sample size.

\---

## Key Findings (from original data)

* Adults showed significantly higher mean nitrogen content than larvae (p < 0.01)
* Samples from L1 exhibited significantly higher mean N content than L2
* A strong non-linear correlation between weight and N content was found for adults in L1 (r = 0.96, p < 0.01)
* No significant correlation was found for larvae subgroups

\---

## How to Run

```r
# 1. Install dependencies (handled automatically in analysis.R)
# 2. Generate synthetic data
source("generate\\\_synthetic\\\_data.R")
# 3. Run full analysis
source("analysis.R")
```

Requires R ≥ 4.0. The `nlcor` package installs from GitHub automatically.

\---

## Contact

Konstantina Gianniou  
g.tem2106@gmail.com | LinkedIn | GitHub

