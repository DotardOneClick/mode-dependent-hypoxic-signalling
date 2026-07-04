"""
Preprocessing pipeline — Article 2, raw amplification data → [C]rel
Steps:
    1. Load raw long-format Ct data (PCR_Results_Corrected.xlsx)
    2. Pair each target-gene Ct with Actin (reference) Ct by animal order
       within (group, treat, tissue) → delta Ct = Ct(target) - Ct(Actin)
    3. Remove delta-Ct outliers with iterative Grubbs' test (alpha=0.05)
    4. delta-delta Ct = delta Ct(sample) - mean(delta Ct) of Control+none,
       calibrated separately per tissue (cortex/mPFC and hippocampus have
       different baseline expression)
    5. [C]rel = 2^(-delta-delta Ct)   (Livak & Schmittgen, 2001)
Outputs:
    data/Art2_per_animal_Crel.csv   — long format, one row per animal/gene
    data/Art2_grubbs_outliers.csv   — log of removed points
"""
import numpy as np
import pandas as pd
import openpyxl
from scipy import stats

RAW_FILE = 'data/PCR_Results_Corrected.xlsx'

TARGET_GENES = ['HIF-1', 'HIF-2', 'HIF-3', 'PACAP', 'PAI']
REF_GENE = 'Actin'

GROUP_MAP = {
    ('control', 'none'):   'Control',
    ('ptsd',    'none'):   'PTSD',
    ('control', 'cobalt'): 'Control+CoCl₂',
    ('ptsd',    'cobalt'): 'PTSD+CoCl₂',
    ('control', 'iht'):    'Control+IHT',
    ('ptsd',    'iht'):    'PTSD+IHT',
    ('control', 'bar'):    'Control+Bar',
    ('ptsd',    'bar'):    'PTSD+Bar',
}
TISSUE_MAP = {'hippo': 'Hippocampus', 'mpfc': 'mPFC'}


def load_raw(filepath):
    wb = openpyxl.load_workbook(filepath, data_only=True)
    ws = wb['Аркуш1']
    rows = list(ws.iter_rows(values_only=True))
    header = rows[0]
    idx = {h: i for i, h in enumerate(header)}
    data = []
    for r in rows[1:]:
        if r[idx['Group']] == 'ntc' or r[idx['Group']] is None:
            continue
        ct = r[idx['Cт']]
        if not isinstance(ct, (int, float)):
            continue  
        data.append({
            'group':  r[idx['Group']],
            'treat':  r[idx['Treat']],
            'tissue': r[idx['Sample Name']],
            'target': r[idx['Target Name']],
            'ct':     ct,
        })
    return data


def build_blocks(data):
    """blocks[(group,treat,tissue,target)] = [ct1, ct2, ...] in original row order"""
    blocks = {}
    for rec in data:
        key = (rec['group'], rec['treat'], rec['tissue'], rec['target'])
        blocks.setdefault(key, []).append(rec['ct'])
    return blocks


def grubbs_remove(vals, alpha=0.05):
    """Iterative two-sided Grubbs' test. Returns (kept_vals, removed_vals)."""
    vals = list(vals)
    removed = []
    while len(vals) >= 3:
        mean = np.mean(vals)
        sd = np.std(vals, ddof=1)
        if sd == 0:
            break
        abs_dev = [abs(v - mean) for v in vals]
        i_max = int(np.argmax(abs_dev))
        G = abs_dev[i_max] / sd
        n = len(vals)
        t_crit = stats.t.ppf(1 - alpha / (2 * n), n - 2)
        G_crit = ((n - 1) / np.sqrt(n)) * np.sqrt(t_crit ** 2 / (n - 2 + t_crit ** 2))
        if G > G_crit:
            removed.append(vals.pop(i_max))
        else:
            break
    return vals, removed


def run():
    data = load_raw(RAW_FILE)
    blocks = build_blocks(data)

    # 1-2. delta Ct per animal, paired positionally within (group,treat,tissue)
    delta_rows = []
    for (group, treat, tissue), _ in {(k[0], k[1], k[2]): None for k in blocks}.items():
        actin_vals = blocks.get((group, treat, tissue, REF_GENE), [])
        for gene in TARGET_GENES:
            target_vals = blocks.get((group, treat, tissue, gene), [])
            n = min(len(target_vals), len(actin_vals))
            for i in range(n):
                delta_rows.append({
                    'group': group, 'treat': treat, 'tissue': tissue, 'gene': gene,
                    'animal_idx': i + 1,
                    'ct_target': target_vals[i], 'ct_actin': actin_vals[i],
                    'delta_ct': target_vals[i] - actin_vals[i],
                })
    ddf = pd.DataFrame(delta_rows)

    # 3. Grubbs outlier removal per (group,treat,tissue,gene) on delta_ct
    outlier_log = []
    keep_mask = pd.Series(True, index=ddf.index)
    for (group, treat, tissue, gene), sub in ddf.groupby(['group', 'treat', 'tissue', 'gene']):
        kept_vals, removed_vals = grubbs_remove(sub['delta_ct'].tolist())
        if removed_vals:
            for rv in removed_vals:
                ridx = sub.index[sub['delta_ct'] == rv][0]
                keep_mask[ridx] = False
                outlier_log.append({
                    'group': group, 'treat': treat, 'tissue': tissue, 'gene': gene,
                    'removed_delta_ct': round(rv, 4), 'n_before': len(sub),
                })
    outliers_df = pd.DataFrame(outlier_log)
    ddf_clean = ddf[keep_mask].copy()

    # 4. delta-delta Ct, calibrated to Control+none, per tissue+gene
    crel_rows = []
    for (tissue, gene), sub in ddf_clean.groupby(['tissue', 'gene']):
        ctrl_mean = sub.loc[(sub['group'] == 'control') & (sub['treat'] == 'none'), 'delta_ct'].mean()
        for _, r in sub.iterrows():
            ddct = r['delta_ct'] - ctrl_mean
            crel = 2 ** (-ddct)
            crel_rows.append({
                'group_label': GROUP_MAP[(r['group'], r['treat'])],
                'tissue_label': TISSUE_MAP[r['tissue']],
                'gene': gene,
                'animal_idx': r['animal_idx'],
                'delta_ct': round(r['delta_ct'], 4),
                'delta_delta_ct': round(ddct, 4),
                'Crel': round(crel, 6),
            })
    crel_df = pd.DataFrame(crel_rows)

    crel_df.to_csv('data/Art2_per_animal_Crel.csv', index=False)
    outliers_df.to_csv('data/Art2_grubbs_outliers.csv', index=False)

    print(f"Total animal-gene datapoints (paired): {len(ddf)}")
    print(f"Grubbs outliers removed: {len(outliers_df)}")
    print(f"Final per-animal Crel rows: {len(crel_df)}")
    print("\nOutliers removed:")
    if len(outliers_df):
        print(outliers_df.to_string(index=False))
    return crel_df, outliers_df


if __name__ == '__main__':
    run()
