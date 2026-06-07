# Mode-Dependent Effects of Chemical, Intermittent, and Hyperbaric Hypoxia on HIF-1α, HIF-2α, HIF-3α, PACAP, and PAI mRNA Expression in the Single Prolonged Stress Model of PTSD

Statistical analysis pipeline for a rat PTSD model study comparing three
hypoxic interventions — cobalt chloride (CoCl₂), intermittent normobaric
hypoxia (IHT), and hypobaric chamber exposure — on anxiety-like behaviour
and gene expression in the medial prefrontal cortex and hippocampus.

## Authors

Denys Porkhalo1, Denys Pashevin1, Yana Naumenko2, Roman Koval3, Mariia Kozlovska4,5, Mariia Kamkina6, Yan Tytarenko7, Alla Portnychenko4,5, Victor Dosenko1 
1. Department of General and Molecular Pathophysiology, Bogomoletz Institute of Physiology of NAS of Ukraine, Kyiv, Ukraine
2. Department of Biophysics of Sensory Signalling, Bogomoletz Institute of Physiology of NAS of Ukraine, Kyiv
3. Department of Cytology, Bogomoletz Institute of Physiology of NAS of Ukraine, Kyiv
4. Department of Hypoxia, Bogomoletz Institute of Physiology of NAS of Ukraine, Kyiv
5. International Centre for Astronomical and Medico-Ecological Research, NAS of Ukraine, Kyiv, Ukraine
6. Kyiv Medical University, Kyiv, Ukraine
7. Bogomolets National Medical University, Kyiv, Ukraine


## Study Overview

Adult male rats (n = 65; 7 experimental groups) were subjected to the
Single Prolonged Stress (SPS) PTSD protocol followed by one of three
hypoxic interventions. Anxiety-like behaviour was evaluated using the
Elevated Plus Maze (EPM), Open Field Test (OFT), and Dark-Light Box (DLB).
Gene expression of HIF-1α, HIF-2α, HIF-3α, PACAP, and PAI-1 was
quantified in the medial prefrontal cortex and hippocampus by RT-qPCR
(Applied Biosystems 7500 Fast, 2^(−ΔΔCt) method, β-actin reference gene).

**Experimental groups:**

| Group | n |
|---|---|
| Control (no intervention) | 18 |
| PTSD (no intervention) | 18 |
| Control + CoCl₂ | 12 |
| PTSD + CoCl₂ | 11 |
| Control + IHT | 12 |
| PTSD + IHT | 10 |
| Control + Barochamber | 12 |
| PTSD + Barochamber | 12 |

## Setup

```bash
pip install numpy scipy matplotlib openpyxl
```

## Data Files

Place in `data/` folder:

**Molecular data (qPCR):**
- `PTSD_Hypoxia_Cortex.xlsx` — gene expression, medial prefrontal cortex
- `PTSD_Hypoxia_Hipo.xlsx` — gene expression, hippocampus

Each file contains 5 sheets (HIF1, HIF2, HIF3, PACAP, PAI1).
Groups are arranged in columns (6 columns per group, [C]rel at offset+4):
Control (0), PTSD (6), Control+CoCl₂ (12), PTSD+CoCl₂ (18),
Control+IHT (24), PTSD+IHT (30), Control+Bar (36), PTSD+Bar (42).

**Behavioural data:**
- `OF_DF_ALL.xlsx` — Open Field Test
- `EP_DF_ALL.xlsx` — Elevated Plus Maze
- `DL_DF_ALL.xlsx` — Dark-Light Box

Columns: `animal №`, `group`, `treatment`, followed by measured parameters.
Groups: `control` / `ptsd` / `cont`. Treatments: `non`, `stage_1` (CoCl₂),
`stage_2` (IHT), `stage_3` (barochamber).

## Usage

```bash
# Molecular analysis (qPCR)
python molecular.py
python molecular.py --no-outliers

# Behavioural analysis
python behavioural.py
python behavioural.py --test EPM
python behavioural.py --no-outliers

# Spearman correlations (selected gene-behaviour pairs)
python correlation.py
```

## Statistical Pipeline

1. **Outlier removal** — IQR 1.5× per group
2. **Normality** — Shapiro-Wilk test per group
3. **Overall test** — one-way ANOVA (all normal) or Kruskal-Wallis
4. **Post-hoc** — Student's t-test (after ANOVA) or Mann-Whitney U (after Kruskal)
5. **Multiple comparisons** — Benjamini-Hochberg FDR correction per gene
   (molecular data only; 12 pairwise comparisons per gene)

Behavioural comparisons between Control and PTSD groups use raw p-values.
Correlations are exploratory Spearman rank analyses (selected pairs, uncorrected).
