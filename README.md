# qPCR Analysis Pipeline (Hypoxia/PTSD study)

Processing and statistical analysis of qPCR data (HIF-1α, HIF-2α, HIF-3α, PACAP, PAI-1) in hippocampus and mPFC of rats (PTSD model + intermittent hypoxic training).

## Structure

```
project/
├── preprocess.py             # Ct → ΔCt → Grubbs → ΔΔCt → [C]rel
├── molecular_analysis_v2.py  # group statistics + boxplot figures
├── correlation_analysis.py   # Spearman correlations, gene × behaviour
└── data/
    ├── PCR_Results_Corrected.xlsx
    ├── OF_DF_ALL.xlsx
    ├── EP_DF_ALL.xlsx
    └── DL_DF_ALL.xlsx
```

## Installation

```bash
pip install numpy pandas scipy matplotlib openpyxl
```

## Running

`preprocess.py` is imported automatically by the other two scripts — no need to run it separately.

```bash
python molecular_analysis.py    # group statistics + figures
python correlation_analysis.py     # gene–behaviour correlations
```

Results are saved to `results/<timestamp>/`.

## What each script does

**preprocess.py**
Raw Ct → paired with Actin by animal order (within group+treat+tissue) → ΔCt → outlier removal (iterative Grubbs' test, α=0.05) → ΔΔCt (calibrated to Control+none, per tissue) → `[C]rel = 2^(−ΔΔCt)` (Livak & Schmittgen, 2001).

**molecular_analysis_v2.py**
One-way ANOVA/Kruskal-Wallis + pairwise post-hoc (t-test/Mann-Whitney), Benjamini-Hochberg correction (25 comparisons per gene×tissue). Generates boxplots (HIF-1/2/3, PACAP, PAI-1) and `Art2_qPCR_stats_v2.csv`.

**correlation_analysis.py**
Spearman correlations of group means (n=7) of [C]rel vs behavioural parameters (EPM, OF, DLB). Outputs a figure of significant correlations and `Art2_correlations_stats.csv`.

## Notes

- Cortex tissue is labeled **mPFC** throughout.
- "Undetermined" Ct values are excluded as non-amplified.
- Each run creates a new `results/<timestamp>/` folder — previous results are not overwritten.
