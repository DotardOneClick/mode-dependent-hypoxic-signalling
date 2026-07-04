"""
qPCR Statistical (v2, recalculated from raw Ct)
Pipeline:
    1. Raw Ct (long format) -> delta Ct vs Actin, paired by animal order
       within group+treat+tissue
    2. Grubbs' test (alpha=0.05, iterative) removes delta-Ct outliers
    3. delta-delta Ct vs Control+none (calibrated per tissue+gene)
    4. [C]rel = 2^(-delta-delta Ct)                (Livak & Schmittgen, 2001)
    5. Group comparisons: one-way ANOVA/Kruskal-Wallis + t-test/Mann-Whitney
       post-hoc, same pairwise scheme as Article 2 original analysis
Figures (cortex tissue labeled "mPFC" everywhere):
    Art2_HIF123_comparison.png/pdf   - 2 rows (mPFC/Hippocampus) x 3 cols (HIF-1/2/3)
    Art2_PACAP_comparison.png/pdf    - 2 rows x 1 col
    Art2_PAI_comparison.png/pdf      - 2 rows x 1 col
Usage:
    python molecular_art2_v2.py
"""
import os
from datetime import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

from preprocess import run as run_preprocess

# ── Config ──────────────────────────────────────────────────────────────────
TISSUES     = ['mPFC', 'Hippocampus']
GENE_LABELS = {'HIF-1': 'HIF-1α', 'HIF-2': 'HIF-2α', 'HIF-3': 'HIF-3α',
               'PACAP': 'PACAP', 'PAI': 'PAI-1'}

REF_GROUPS = ['Control', 'PTSD']
BOX_GROUPS = ['Control+CoCl₂', 'PTSD+CoCl₂',
              'Control+IHT',   'PTSD+IHT',
              'Control+Bar',   'PTSD+Bar']
BOX_POSITIONS = [0, 1, 3, 4, 6, 7]
THERAPY_COLORS = {
    'Control+CoCl₂': '#64B5F6', 'PTSD+CoCl₂': '#1565C0',
    'Control+IHT':   '#81C784', 'PTSD+IHT':   '#2E7D32',
    'Control+Bar':   '#FFD54F', 'PTSD+Bar':   '#E65100',
}
CTRL_LINE_COLOR = '#388E3C'
PTSD_LINE_COLOR = '#7B1FA2'
COMPARISON_PAIRS = [
    ('Control+CoCl₂', 'Control',       'Control+CoCl₂ vs Control'),
    ('PTSD+CoCl₂',    'PTSD',          'PTSD+CoCl₂ vs PTSD'),
    ('PTSD+CoCl₂',    'Control',       'PTSD+CoCl₂ vs Control'),
    ('PTSD+CoCl₂',    'Control+CoCl₂', 'PTSD+CoCl₂ vs Control+CoCl₂'),
    ('Control+IHT',   'Control',       'Control+IHT vs Control'),
    ('PTSD+IHT',      'PTSD',          'PTSD+IHT vs PTSD'),
    ('PTSD+IHT',      'Control',       'PTSD+IHT vs Control'),
    ('PTSD+IHT',      'Control+IHT',   'PTSD+IHT vs Control+IHT'),
    ('Control+Bar',   'Control',       'Control+Bar vs Control'),
    ('PTSD+Bar',      'PTSD',          'PTSD+Bar vs PTSD'),
    ('PTSD+Bar',      'Control',       'PTSD+Bar vs Control'),
    ('PTSD+Bar',      'Control+Bar',   'PTSD+Bar vs Control+Bar'),
    ('Control+IHT',   'Control+CoCl₂', 'Control+IHT vs Control+CoCl₂'),
    ('Control+Bar',   'Control+CoCl₂', 'Control+Bar vs Control+CoCl₂'),
    ('Control+Bar',   'Control+IHT',   'Control+Bar vs Control+IHT'),
    ('PTSD+IHT',      'PTSD+CoCl₂',    'PTSD+IHT vs PTSD+CoCl₂'),
    ('PTSD+Bar',      'PTSD+CoCl₂',    'PTSD+Bar vs PTSD+CoCl₂'),
    ('PTSD+Bar',      'PTSD+IHT',      'PTSD+Bar vs PTSD+IHT'),
    ('PTSD+IHT',      'Control+CoCl₂', 'PTSD+IHT vs Control+CoCl₂'),
    ('PTSD+Bar',      'Control+CoCl₂', 'PTSD+Bar vs Control+CoCl₂'),
    ('Control+IHT',   'PTSD+CoCl₂',    'Control+IHT vs PTSD+CoCl₂'),
    ('Control+Bar',   'PTSD+CoCl₂',    'Control+Bar vs PTSD+CoCl₂'),
    ('PTSD+Bar',      'Control+IHT',   'PTSD+Bar vs Control+IHT'),
    ('Control+Bar',   'PTSD+IHT',      'Control+Bar vs PTSD+IHT'),
    ('PTSD',          'Control',       'PTSD vs Control'),
]
plt.rcParams.update({'font.family': 'DejaVu Sans', 'pdf.fonttype': 42, 'ps.fonttype': 42})

# ── Statistics (unchanged pairwise scheme) ────────────────────────────────────
def shapiro_wilk(v):
    if len(v) < 3: return False, None
    _, p = stats.shapiro(v)
    return p >= 0.05, round(p, 4)

def run_overall(groups_dict):
    vals = list(groups_dict.values())
    norm = all(shapiro_wilk(v)[0] for v in vals if len(v) >= 3)
    if norm:
        stat, p = stats.f_oneway(*vals)
        return {'test': 'One-way ANOVA', 'statistic': round(stat, 4), 'p': round(p, 6), 'all_normal': True}
    stat, p = stats.kruskal(*vals)
    return {'test': 'Kruskal-Wallis', 'statistic': round(stat, 4), 'p': round(p, 6), 'all_normal': False}

def run_posthoc(v1, v2, use_ttest):
    if len(v1) < 3 or len(v2) < 3: return None, 'n/a'
    if use_ttest:
        _, p = stats.ttest_ind(v1, v2)
        return round(p, 6), 't-test'
    _, p = stats.mannwhitneyu(v1, v2, alternative='two-sided')
    return round(p, 6), 'Mann-Whitney U'

def sig_label(p):
    if p is None: return ''
    if p < 0.001: return '***'
    if p < 0.01:  return '**'
    if p < 0.05:  return '*'
    return 'ns'

def benjamini_hochberg(pvals):
    """BH-FDR adjusted p-values. Input/output preserve order; None passed through."""
    idx = [i for i, p in enumerate(pvals) if p is not None]
    ps = [pvals[i] for i in idx]
    m = len(ps)
    adj = [None] * len(pvals)
    if m == 0:
        return adj
    order = sorted(range(m), key=lambda k: ps[k])
    ranked = [None] * m
    prev = 1.0
    for rank in range(m - 1, -1, -1):
        k = order[rank]
        val = ps[k] * m / (rank + 1)
        prev = min(prev, val)
        ranked[k] = min(prev, 1.0)
    for j, i in enumerate(idx):
        adj[i] = ranked[j]
    return adj

def _get_x(group_name):
    if group_name in BOX_GROUPS:
        return BOX_POSITIONS[BOX_GROUPS.index(group_name)]
    return None

# ── Plot (single gene panel) ──────────────────────────────────────────────────
def plot_gene(ax, ref_data, box_data, gene_label, posthoc, show_ylabel=False,
              show_xlabels=True, panel_title=None):
    x_min, x_max = -0.7, 7.7
    ctrl_vals = ref_data.get('Control', [])
    ptsd_vals = ref_data.get('PTSD', [])
    ctrl_median = np.median(ctrl_vals) if ctrl_vals else None
    ptsd_median = np.median(ptsd_vals) if ptsd_vals else None

    if ctrl_median is not None:
        ax.hlines(ctrl_median, x_min, x_max, colors=CTRL_LINE_COLOR, linestyles='--', linewidth=2.0, zorder=6)
    if ptsd_median is not None:
        ax.hlines(ptsd_median, x_min, x_max, colors=PTSD_LINE_COLOR, linestyles=':', linewidth=2.0, zorder=6)

    all_vals = list(ctrl_vals) + list(ptsd_vals)
    for pos, grp in zip(BOX_POSITIONS, BOX_GROUPS):
        vals = box_data.get(grp, [])
        if not vals: continue
        all_vals += vals
        color = THERAPY_COLORS[grp]
        bp = ax.boxplot([vals], positions=[pos], widths=0.55, patch_artist=True, notch=False,
                         medianprops=dict(color='#212121', linewidth=2.0),
                         whiskerprops=dict(color='#424242', linewidth=1.2),
                         capprops=dict(color='#424242', linewidth=1.2),
                         flierprops=dict(marker='', markersize=0),
                         boxprops=dict(linewidth=1.2), zorder=3)
        bp['boxes'][0].set_facecolor(color)
        bp['boxes'][0].set_alpha(0.75)
        bp['boxes'][0].set_edgecolor('#212121')
        np.random.seed(42 + pos)
        jitter = np.random.uniform(-0.15, 0.15, len(vals))
        ax.scatter(pos + jitter, vals, color='#212121', s=22, alpha=0.65, linewidths=0, zorder=5)

    if not all_vals: return
    data_max, data_min = max(all_vals), min(all_vals)
    all_pos = [v for v in all_vals if v > 0]
    use_log = all_pos and (max(all_pos) / max(min(all_pos), 0.001)) > 8
    top = data_max * 1.60
    bottom = -top * 0.06 if data_min >= 0 else data_min * 1.1
    ax.set_xlim(x_min, x_max)
    if use_log:
        ax.set_yscale('log')
        ax.set_ylim(bottom=min(all_pos) * 0.5)
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda val, _: f'{val:.2g}'))
    else:
        ax.set_ylim(bottom, top)

    # significance brackets
    sig_pairs = [(g1, g2, r) for (g1, g2), r in posthoc.items() if r.get('sig', '') not in ('', 'ns')]
    if use_log:
        base_y = ax.get_ylim()[1]
        step = base_y * 0.28
    else:
        base_y = top * 0.62
        step = top * 0.14
    level_tops = {}
    for g1, g2, result in sig_pairs:
        sl = result.get('sig', '')
        p1, p2 = _get_x(g1), _get_x(g2)
        if p1 is None or p2 is None: continue
        lo, hi = min(p1, p2), max(p1, p2)
        y = base_y
        for col in [lo, hi]:
            if col in level_tops:
                y = max(y, level_tops[col] + step)
        for col in [lo, hi]:
            level_tops[col] = y + step * 0.6
        needed = y + step * 1.4
        if needed > ax.get_ylim()[1]:
            if use_log: ax.set_ylim(top=needed)
            else: ax.set_ylim(bottom, needed)
        h = step * 0.22
        ax.plot([p1, p1, p2, p2], [y, y+h, y+h, y], lw=1.0, color='#212121', clip_on=False, zorder=8)
        ax.text((p1+p2)/2, y + h*1.3, sl, ha='center', va='bottom', fontsize=12,
                color='#212121', fontweight='bold', clip_on=False, zorder=9)

    ax.set_xticks(BOX_POSITIONS)
    if show_xlabels:
        ax.set_xticklabels([g.replace('+', '\n+') for g in BOX_GROUPS], fontsize=11,
                            fontweight='bold', color='#212121', linespacing=1.2)
    else:
        ax.set_xticklabels([])
    ax.set_title(panel_title if panel_title else gene_label,
                 fontsize=15, fontweight='bold', color='#212121', pad=8)
    if show_ylabel:
        ax.set_ylabel('[C]rel', fontsize=14, fontweight='bold', color='#424242', labelpad=10)
    if show_xlabels:
        ax.set_xlabel('Group + treatment', fontsize=14, fontweight='bold', color='#424242', labelpad=8)
    ax.set_facecolor('#EBEBEB')
    ax.yaxis.grid(True, color='white', linewidth=1.0, zorder=0)
    ax.xaxis.grid(False)
    ax.set_axisbelow(True)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_color('#BDBDBD')
    ax.spines['bottom'].set_linewidth(0.8)
    ax.tick_params(axis='x', length=0, pad=8)
    ax.tick_params(axis='y', length=0, labelsize=12)


def tissue_row_label(ax, text):
    ax.annotate(text, xy=(-0.30, 0.5), xycoords='axes fraction', rotation=90,
                ha='center', va='center', fontsize=12, fontweight='bold', color='#212121')


def make_multi_figure(genes, crel_data, stats_data, title, out_prefix):
    single_gene = len(genes) == 1

    if single_gene:
        # tissues side by side: 1 row x len(TISSUES) cols
        gene = genes[0]
        n_rows, n_cols = 1, len(TISSUES)
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(6.6*n_cols, 6.6), squeeze=False)
        fig.patch.set_facecolor('white')
        for ci, tissue in enumerate(TISSUES):
            data = crel_data[tissue][gene]
            ref_dict = {g: data.get(g, []) for g in REF_GROUPS}
            box_dict = {g: data.get(g, []) for g in BOX_GROUPS}
            posthoc = stats_data[tissue][gene]['posthoc']
            plot_gene(axes[0][ci], ref_dict, box_dict, GENE_LABELS[gene], posthoc,
                      show_ylabel=(ci == 0), show_xlabels=True,
                      panel_title=f'{GENE_LABELS[gene]} — {tissue}')
    else:
        # gene grid: rows = tissue, cols = gene
        n_rows, n_cols = len(TISSUES), len(genes)
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(5.6*n_cols, 5.8*n_rows), squeeze=False)
        fig.patch.set_facecolor('white')
        for ri, tissue in enumerate(TISSUES):
            for ci, gene in enumerate(genes):
                data = crel_data[tissue][gene]
                ref_dict = {g: data.get(g, []) for g in REF_GROUPS}
                box_dict = {g: data.get(g, []) for g in BOX_GROUPS}
                posthoc = stats_data[tissue][gene]['posthoc']
                plot_gene(axes[ri][ci], ref_dict, box_dict, GENE_LABELS[gene], posthoc,
                          show_ylabel=True, show_xlabels=True,
                          panel_title=f'{GENE_LABELS[gene]} — {tissue}')

    legend_elements = [
        mlines.Line2D([0], [0], color=CTRL_LINE_COLOR, linestyle='--', linewidth=2, label='Control (median)'),
        mlines.Line2D([0], [0], color=PTSD_LINE_COLOR, linestyle=':', linewidth=2, label='PTSD (median)'),
    ] + [mpatches.Patch(facecolor=THERAPY_COLORS[g], edgecolor='#212121', alpha=0.75, label=g) for g in BOX_GROUPS]
    legend_y = -0.10 if single_gene else -0.06
    fig.legend(handles=legend_elements, loc='lower center', ncol=4, fontsize=12, frameon=True,
               fancybox=False, edgecolor='#E0E0E0', bbox_to_anchor=(0.5, legend_y))

    fig.suptitle(title, fontsize=18, fontweight='bold', color='#212121', y=1.02)
    fig.text(0.5, (-0.17 if single_gene else -0.11),
              '[C]rel = 2^(-ΔΔCt) vs Control(none), Δ actin-normalized. '
              'ΔCt outliers removed by Grubbs\' test. '
              'Dashed/dotted lines = Control/PTSD baseline medians. '
              'Significance from Benjamini-Hochberg FDR-adjusted p (per gene×tissue): '
              '* p<0.05, ** p<0.01, *** p<0.001.',
              ha='center', fontsize=10, color='#9E9E9E', style='italic')

    if single_gene:
        plt.tight_layout(rect=[0.01, 0.02, 0.99, 0.99])
        plt.subplots_adjust(wspace=0.22)
    else:
        plt.tight_layout(rect=[0.01, 0.02, 0.99, 0.99])
        plt.subplots_adjust(wspace=0.28, hspace=0.4)
    fig.savefig(f'{out_prefix}.png', dpi=300, bbox_inches='tight', facecolor='white', pad_inches=0.1)
    fig.savefig(f'{out_prefix}.pdf', bbox_inches='tight', facecolor='white', pad_inches=0.1)
    plt.close(fig)
    print(f"  Figure -> {out_prefix}.png / .pdf")


# ── Main ────────────────────────────────────────────────────────────────────
def run():
    ts = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    out_dir = os.path.join('results', ts)
    os.makedirs(out_dir, exist_ok=True)

    print("="*60)
    print("qPCR Analysis — Article 2 (v2, recalculated)")
    print(f"Output -> {out_dir}/")
    print("="*60)

    crel_df, outliers_df = run_preprocess()
    outliers_df.to_csv(os.path.join(out_dir, 'Art2_grubbs_outliers.csv'), index=False)
    crel_df.to_csv(os.path.join(out_dir, 'Art2_per_animal_Crel.csv'), index=False)

    genes = ['HIF-1', 'HIF-2', 'HIF-3', 'PACAP', 'PAI']
    all_groups = REF_GROUPS + BOX_GROUPS

    # nested: crel_data[tissue][gene][group_label] = [values]
    crel_data = {t: {g: {} for g in genes} for t in TISSUES}
    stats_data = {t: {g: {} for g in genes} for t in TISSUES}
    csv_rows = []

    for tissue in TISSUES:
        for gene in genes:
            sub = crel_df[(crel_df.tissue_label == tissue) & (crel_df.gene == gene)]
            all_data = {g: sub.loc[sub.group_label == g, 'Crel'].tolist() for g in all_groups}
            crel_data[tissue][gene] = all_data

            non_empty = {g: v for g, v in all_data.items() if len(v) >= 3}
            if len(non_empty) < 2:
                stats_data[tissue][gene] = {'posthoc': {}}
                continue
            overall = run_overall(non_empty)

            # 1st pass: raw post-hoc p per pair
            posthoc = {}
            pair_order = []
            raw_ps = []
            for g1, g2, pair_label in COMPARISON_PAIRS:
                v1, v2 = all_data.get(g1, []), all_data.get(g2, [])
                if len(v1) < 3 or len(v2) < 3: continue
                p, test = run_posthoc(v1, v2, overall['all_normal'])
                posthoc[(g1, g2)] = {'p': p, 'test': test, 'label': pair_label}
                pair_order.append((g1, g2))
                raw_ps.append(p)

            # BH-FDR across this gene x tissue family (25 pairs)
            p_adj = benjamini_hochberg(raw_ps)
            for (g1, g2), pa in zip(pair_order, p_adj):
                posthoc[(g1, g2)]['p_adj'] = pa
                sl = sig_label(pa)
                posthoc[(g1, g2)]['sig'] = sl
                if sl not in ('ns', ''):
                    v1, v2 = all_data.get(g1, []), all_data.get(g2, [])
                    d = '↑' if np.mean(v1) > np.mean(v2) else '↓'
                    print(f"  {tissue:<12} {GENE_LABELS[gene]:<8} {posthoc[(g1,g2)]['label']:<35} "
                          f"{sl} p={posthoc[(g1,g2)]['p']:.4f} p_adj={pa:.4f} {d}")

            stats_data[tissue][gene] = {'posthoc': posthoc, 'overall': overall}

            for g1, g2, pair_label in COMPARISON_PAIRS:
                v1, v2 = all_data.get(g1, []), all_data.get(g2, [])
                if not v1 or not v2: continue
                ph = posthoc.get((g1, g2), {})
                m1, m2 = np.mean(v1), np.mean(v2)
                s1, s2 = np.std(v1, ddof=1), np.std(v2, ddof=1)
                csv_rows.append({
                    'tissue': tissue, 'gene': GENE_LABELS[gene], 'comparison': pair_label,
                    'group1': g1, 'group2': g2,
                    'overall_test': overall['test'], 'overall_p': overall['p'],
                    'posthoc_test': ph.get('test', ''), 'posthoc_p': ph.get('p', ''),
                    'posthoc_p_adj_BH': (round(ph['p_adj'], 6) if ph.get('p_adj') is not None else ''),
                    'significance': ph.get('sig', ''),
                    'direction': ('↑' if m1 > m2 else '↓') if ph.get('sig', '') not in ('ns', '') else '—',
                    'mean_g1_Crel': round(m1, 4), 'sd_g1': round(s1, 4),
                    'sem_g1': round(s1/np.sqrt(len(v1)), 4), 'n_g1': len(v1),
                    'mean_g2_Crel': round(m2, 4), 'sd_g2': round(s2, 4),
                    'sem_g2': round(s2/np.sqrt(len(v2)), 4), 'n_g2': len(v2),
                })

    pd.DataFrame(csv_rows).to_csv(os.path.join(out_dir, 'Art2_qPCR_stats_v2.csv'), index=False)
    print(f"\n  CSV    -> Art2_qPCR_stats_v2.csv")

    make_multi_figure(['HIF-1', 'HIF-2', 'HIF-3'], crel_data, stats_data,
                       'HIF-1/2/3 expression',
                       os.path.join(out_dir, 'Art2_HIF123_comparison'))
    make_multi_figure(['PACAP'], crel_data, stats_data,
                       'PACAP expression',
                       os.path.join(out_dir, 'Art2_PACAP_comparison'))
    make_multi_figure(['PAI'], crel_data, stats_data,
                       'PAI-1 expression',
                       os.path.join(out_dir, 'Art2_PAI_comparison'))

    print(f"\nDone! -> {out_dir}/")
    return out_dir


if __name__ == '__main__':
    run()
