"""
qPCR Statistical Analysis
Control (non) and PTSD (non) shown as reference median lines.
6 therapy groups shown as boxplots.
Statistical pipeline:
    1. Outlier removal — IQR 1.5×
    2. Normality — Shapiro-Wilk per group
    3. Overall — one-way ANOVA (all normal) or Kruskal-Wallis
    4. Post-hoc — t-test (after ANOVA) or Mann-Whitney U (after Kruskal)
    5. Multiple comparisons — Benjamini-Hochberg FDR correction (per gene)
    ALL 28 pairwise comparisons (all groups vs all groups)
Usage:
    python molecular.py
    python molecular.py --no-outliers
Data files:
    PTSD_Hypoxia_Cortex.xlsx
    PTSD_Hypoxia_Hipo.xlsx
"""
import os, csv, argparse
from itertools import combinations
from datetime import datetime
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from scipy import stats
import openpyxl
import warnings
warnings.filterwarnings('ignore')

# ── Config ──────────────────────────────────────────────────────────────────────
DATA_FILES = {
    'Cortex':      'PTSD_Hypoxia_Cortex.xlsx',
    'Hippocampus': 'PTSD_Hypoxia_Hipo.xlsx',
}
SHEETS      = ['HIF1', 'HIF2', 'HIF3', 'PACAP', 'PAI1']
GENE_LABELS = {
    'HIF1': 'HIF-1α', 'HIF2': 'HIF-2α', 'HIF3': 'HIF-3α',
    'PACAP': 'PACAP',  'PAI1': 'PAI-1',
}
GROUP_OFFSETS = {
    'Control':       0,
    'PTSD':          6,
    'Control+CoCl₂': 12,
    'PTSD+CoCl₂':    18,
    'Control+IHT':   24,
    'PTSD+IHT':      30,
    'Control+Bar':   36,
    'PTSD+Bar':      42,
}
DATA_ROW = 6
ALL_GROUPS = ['Control', 'PTSD',
              'Control+CoCl₂', 'PTSD+CoCl₂',
              'Control+IHT',   'PTSD+IHT',
              'Control+Bar',   'PTSD+Bar']
REF_GROUPS = ['Control', 'PTSD']
BOX_GROUPS = ['Control+CoCl₂', 'PTSD+CoCl₂',
              'Control+IHT',   'PTSD+IHT',
              'Control+Bar',   'PTSD+Bar']
BOX_POSITIONS = [0, 1, 3, 4, 6, 7]

ALL_POSITIONS = {
    'Control':       -2,
    'PTSD':          -1,
    'Control+CoCl₂': 0,
    'PTSD+CoCl₂':    1,
    'Control+IHT':   3,
    'PTSD+IHT':      4,
    'Control+Bar':   6,
    'PTSD+Bar':      7,
}

THERAPY_COLORS = {
    'Control+CoCl₂': '#64B5F6',
    'PTSD+CoCl₂':    '#1565C0',
    'Control+IHT':   '#81C784',
    'PTSD+IHT':      '#2E7D32',
    'Control+Bar':   '#FFD54F',
    'PTSD+Bar':      '#E65100',
}
CTRL_LINE_COLOR = '#388E3C'
PTSD_LINE_COLOR = '#7B1FA2'

# ALL 28 pairwise comparisons
COMPARISON_PAIRS = [
    (g1, g2, f'{g2} vs {g1}')
    for g1, g2 in combinations(ALL_GROUPS, 2)
]

plt.rcParams.update({
    'font.family':  'DejaVu Sans',
    'pdf.fonttype': 42,
    'ps.fonttype':  42,
})

# ── Data loading ────────────────────────────────────────────────────────────────
def load_vals(filepath, sheet, group):
    wb   = openpyxl.load_workbook(filepath, data_only=True)
    ws   = wb[sheet]
    rows = list(ws.iter_rows(values_only=True))
    wb.close()
    off       = GROUP_OFFSETS[group]
    actin_col = off + 1
    gene_col  = off + 2
    vals = []
    for r in rows[DATA_ROW:]:
        if len(r) <= gene_col:
            continue
        ac, gc = r[actin_col], r[gene_col]
        if (ac is not None and gc is not None
                and isinstance(ac, (int, float))
                and isinstance(gc, (int, float))):
            vals.append(2 ** (float(ac) - float(gc)) * 100)
    return vals

def remove_iqr(vals, factor=1.5):
    if len(vals) < 4: return vals
    q1, q3 = np.percentile(vals, 25), np.percentile(vals, 75)
    iqr = q3 - q1
    return [v for v in vals if q1 - factor*iqr <= v <= q3 + factor*iqr]

# ── Statistics ──────────────────────────────────────────────────────────────────
def shapiro_wilk(v):
    if len(v) < 3: return False, None
    _, p = stats.shapiro(v)
    return p >= 0.05, round(p, 4)

def run_overall(groups_dict):
    vals = list(groups_dict.values())
    norm = all(shapiro_wilk(v)[0] for v in vals if len(v) >= 3)
    if norm:
        stat, p = stats.f_oneway(*vals)
        return {'test': 'One-way ANOVA', 'statistic': round(stat, 4),
                'p': round(p, 6), 'all_normal': True}
    stat, p = stats.kruskal(*vals)
    return {'test': 'Kruskal-Wallis', 'statistic': round(stat, 4),
            'p': round(p, 6), 'all_normal': False}

def run_posthoc(v1, v2, use_ttest):
    if len(v1) < 3 or len(v2) < 3: return None, 'n/a'
    if use_ttest:
        _, p = stats.ttest_ind(v1, v2)
        return round(p, 6), 't-test'
    _, p = stats.mannwhitneyu(v1, v2, alternative='two-sided')
    return round(p, 6), 'Mann-Whitney U'

def bh_correction(p_values, alpha=0.05):
    valid   = [(k, p) for k, p in p_values if p is not None]
    n       = len(valid)
    if n == 0:
        return {}
    sorted_pairs = sorted(valid, key=lambda x: x[1])
    adj_list = [None] * n
    prev_adj = 1.0
    for rank in range(n - 1, -1, -1):
        key, p_raw = sorted_pairs[rank]
        adj = min(prev_adj, p_raw * n / (rank + 1))
        adj = min(adj, 1.0)
        adj_list[rank] = (key, round(adj, 6))
        prev_adj = adj
    return {key: adj for key, adj in adj_list}

def sig_label(p):
    if p is None: return ''
    if p < 0.001: return '***'
    if p < 0.01:  return '**'
    if p < 0.05:  return '*'
    return 'ns'

def _get_x(group_name):
    return ALL_POSITIONS.get(group_name, None)

# ── Plot ─────────────────────────────────────────────────────────────────────────
def plot_gene(ax, ref_data, box_data, gene_label, posthoc):
    x_min, x_max = -0.7, 7.7

    # ── Reference lines ────────────────────────────────────────────────────────
    ctrl_vals   = ref_data.get('Control', [])
    ptsd_vals   = ref_data.get('PTSD', [])
    ctrl_median = np.median(ctrl_vals) if ctrl_vals else None
    ptsd_median = np.median(ptsd_vals) if ptsd_vals else None

    if ctrl_median is not None:
        ax.hlines(ctrl_median, x_min, x_max,
                  colors=CTRL_LINE_COLOR, linestyles='--',
                  linewidth=2.0, zorder=6, label='Control median')
    if ptsd_median is not None:
        ax.hlines(ptsd_median, x_min, x_max,
                  colors=PTSD_LINE_COLOR, linestyles=':',
                  linewidth=2.0, zorder=6, label='PTSD median')

    # ── Boxplots ───────────────────────────────────────────────────────────────
    all_vals = list(ctrl_vals) + list(ptsd_vals)
    for pos, grp in zip(BOX_POSITIONS, BOX_GROUPS):
        vals = box_data.get(grp, [])
        if not vals: continue
        all_vals += vals
        color = THERAPY_COLORS[grp]

        bp = ax.boxplot(
            [vals], positions=[pos], widths=0.55,
            patch_artist=True, notch=False,
            medianprops=dict(color='#212121', linewidth=2.0),
            whiskerprops=dict(color='#424242', linewidth=1.2),
            capprops=dict(color='#424242', linewidth=1.2),
            flierprops=dict(marker='', markersize=0),
            boxprops=dict(linewidth=1.2),
            zorder=3,
        )
        bp['boxes'][0].set_facecolor(color)
        bp['boxes'][0].set_alpha(0.75)
        bp['boxes'][0].set_edgecolor('#212121')

        np.random.seed(42 + pos)
        jitter = np.random.uniform(-0.15, 0.15, len(vals))
        ax.scatter(pos + jitter, vals, color='#212121',
                   s=22, alpha=0.65, linewidths=0, zorder=5)

    # ── Y limits (initial) ─────────────────────────────────────────────────────
    if not all_vals: return
    data_max  = max(all_vals)
    data_min  = min(all_vals)
    all_pos_v = [v for v in all_vals if v > 0]
    use_log   = (all_pos_v and
                 (max(all_pos_v) / max(min(all_pos_v), 0.001)) > 8)

    ax.set_xlim(x_min, x_max)
    if use_log:
        ax.set_yscale('log')
        ax.set_ylim(bottom=min(all_pos_v) * 0.5)
        ax.yaxis.set_major_formatter(
            plt.FuncFormatter(lambda val, _: f'{val:.2g}')
        )
    else:
        top    = data_max * 1.55
        bottom = -top * 0.04 if data_min >= 0 else data_min * 1.1
        ax.set_ylim(bottom, top)

    # ── Significance brackets ─────────────────────────────────────────────────
    sig_pairs = [
        (g1, g2, r) for (g1, g2), r in posthoc.items()
        if r.get('sig_bh', '') not in ('ns', '')
        and g1 in BOX_GROUPS and g2 in BOX_GROUPS
    ]
    sig_pairs.sort(key=lambda item: abs(
        ALL_POSITIONS.get(item[0], 0) - ALL_POSITIONS.get(item[1], 0)
    ))

    y_lo, y_hi = ax.get_ylim()
    data_span  = y_hi - y_lo if not use_log else None

    top_at = {}

    for g1, g2, result in sig_pairs:
        sl = result.get('sig_bh', '')
        p1 = ALL_POSITIONS.get(g1)
        p2 = ALL_POSITIONS.get(g2)
        if p1 is None or p2 is None: continue
        lo_x, hi_x = min(p1, p2), max(p1, p2)

        if use_log:
            cur_top_y = ax.get_ylim()[1]
            base = cur_top_y * 0.85
            step = cur_top_y * 0.12
        else:
            cur_y_lo, cur_y_hi = ax.get_ylim()
            cur_span = cur_y_hi - cur_y_lo
            base = cur_y_lo + cur_span * 0.82
            step = cur_span * 0.055

        y_start = base
        for x_col in range(int(lo_x), int(hi_x) + 1):
            if x_col in top_at:
                y_start = max(y_start, top_at[x_col] + step * 0.3)

        bracket_y = y_start
        h = step * 0.25

        needed = bracket_y + h + step * 0.8
        if not use_log:
            cur_y_lo2, cur_y_hi2 = ax.get_ylim()
            if needed > cur_y_hi2:
                ax.set_ylim(cur_y_lo2, needed * 1.05)
        else:
            _, cur_top2 = ax.get_ylim()
            if needed > cur_top2:
                ax.set_ylim(top=needed * 1.1)

        ax.plot([p1, p1, p2, p2],
                [bracket_y, bracket_y + h, bracket_y + h, bracket_y],
                lw=1.1, color='#212121', clip_on=False, zorder=8)
        ax.text((p1 + p2) / 2, bracket_y + h * 1.15, sl,
                ha='center', va='bottom', fontsize=11,
                color='#212121', fontweight='bold',
                clip_on=False, zorder=9)

        for x_col in range(int(lo_x), int(hi_x) + 1):
            top_at[x_col] = max(top_at.get(x_col, 0),
                                bracket_y + h + step * 0.4)

    # ── Formatting ────────────────────────────────────────────────────────────
    ax.set_xticks(BOX_POSITIONS)
    ax.set_xticklabels(
        [g.replace('+', '\n+') for g in BOX_GROUPS],
        fontsize=10, fontweight='bold', color='#212121', linespacing=1.2
    )
    ax.set_title(gene_label, fontsize=13, fontweight='bold',
                 color='#212121', pad=10)
    ax.set_ylabel('[C]rel', fontsize=12, fontweight='bold',
                  color='#424242', labelpad=12)

    ax.set_facecolor('#EBEBEB')
    ax.yaxis.grid(True, color='white', linewidth=1.0, zorder=0)
    ax.xaxis.grid(False)
    ax.set_axisbelow(True)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_color('#BDBDBD')
    ax.spines['bottom'].set_linewidth(0.8)
    ax.tick_params(axis='x', length=0, pad=10)
    ax.tick_params(axis='y', length=0, labelsize=10)


def _draw_panel(sheets_subset, ref_dict, box_dict, posthoc_dict,
                letter_start, out_prefix, panel_title):
    LETTERS = list('ABCDEFGHIJKLMNOPQRSTUVWXYZ')
    n = len(sheets_subset)
    if n == 0:
        return None, letter_start

    if n <= 3:
        ncols, nrows = n, 1
    else:
        ncols = -(-n // 2)
        nrows = 2

    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(ncols * 5.0, nrows * 6.8),
        squeeze=False
    )
    fig.patch.set_facecolor('white')

    letter_idx = letter_start
    for gi, sheet in enumerate(sheets_subset):
        row, col = divmod(gi, ncols)
        ax = axes[row][col]
        plot_gene(ax,
                  ref_dict.get(sheet, {}),
                  box_dict.get(sheet, {}),
                  GENE_LABELS[sheet],
                  posthoc_dict.get(sheet, {}))

        # ── Panel letter in upper-left corner ──────────────────────────────
        letter = LETTERS[letter_idx] if letter_idx < len(LETTERS) else '?'
        ax.text(-0.08, 1.06, letter,
                transform=ax.transAxes,
                fontsize=16, fontweight='bold', color='#111111',
                va='top', ha='left', clip_on=False,
                fontfamily='DejaVu Sans')
        letter_idx += 1

    for empty in range(n, nrows * ncols):
        row, col = divmod(empty, ncols)
        axes[row][col].set_visible(False)

    fig.suptitle(panel_title,
                 fontsize=13, fontweight='bold', color='#212121', y=1.02)

    plt.tight_layout(rect=[0, 0, 1, 1])
    plt.subplots_adjust(wspace=0.38, hspace=0.52)
    fig.savefig(f'{out_prefix}.png', dpi=300,
                bbox_inches='tight', facecolor='white')
    fig.savefig(f'{out_prefix}.pdf', bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'  Panel "{panel_title}" ({n} genes) → {out_prefix}.png/.pdf')
    return fig, letter_idx


def make_figure(tissue, ref_dict, box_dict, posthoc_dict, output_prefix):
    HIF_SHEETS   = ['HIF1', 'HIF2', 'HIF3']
    OTHER_SHEETS = ['PACAP', 'PAI1']

    def has_sig(sheet):
        return any(r.get('sig_bh', '') not in ('ns', '')
                   for r in posthoc_dict.get(sheet, {}).values())

    hif_draw   = [s for s in HIF_SHEETS   if has_sig(s)]
    other_draw = [s for s in OTHER_SHEETS if has_sig(s)]

    if not hif_draw and not other_draw:
        print(f'  [{tissue}] No significant genes — figures skipped.')
        return

    letter_idx = 0

    if hif_draw:
        _, letter_idx = _draw_panel(
            hif_draw, ref_dict, box_dict, posthoc_dict,
            letter_start=letter_idx,
            out_prefix=f'{output_prefix}_HIF',
            panel_title=f'HIF isoforms — {tissue}',
        )

    if other_draw:
        _draw_panel(
            other_draw, ref_dict, box_dict, posthoc_dict,
            letter_start=letter_idx,
            out_prefix=f'{output_prefix}_PACAP_PAI',
            panel_title=f'PACAP & PAI-1 — {tissue}',
        )


# ── CSV ─────────────────────────────────────────────────────────────────────────
def save_csv(rows, filepath):
    if not rows: return
    fields = list(rows[0].keys())
    with open(filepath, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f'  CSV    → {filepath}')


# ── Main ─────────────────────────────────────────────────────────────────────────
def run(remove_outliers=True):
    ts      = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    out_dir = os.path.join('results', ts)
    os.makedirs(out_dir, exist_ok=True)

    print(f"\n{'='*60}")
    print(f'qPCR Analysis — All-vs-All (28 pairs per gene)')
    print(f"Outlier removal: {'ON (IQR 1.5×)' if remove_outliers else 'OFF'}")
    print(f'Multiple comparisons: Benjamini-Hochberg FDR (per gene)')
    print(f'Output → {out_dir}/')
    print('='*60)

    all_csv = []

    for tissue, filepath in DATA_FILES.items():
        print(f"\n{'─'*60}")
        print(f'  {tissue}')
        print(f"{'─'*60}")

        ref_dict     = {}
        box_dict     = {}
        posthoc_dict = {}

        for sheet in SHEETS:
            gene = GENE_LABELS[sheet]

            all_data = {}
            for g in ALL_GROUPS:
                raw  = load_vals(filepath, sheet, g)
                vals = remove_iqr(raw) if remove_outliers else raw
                all_data[g] = vals

            ref_dict[sheet] = {g: all_data[g] for g in REF_GROUPS}
            box_dict[sheet] = {g: all_data[g] for g in BOX_GROUPS}

            non_empty = {g: v for g, v in all_data.items() if len(v) >= 3}
            if len(non_empty) < 2: continue
            overall = run_overall(non_empty)

            # ── Step 1: compute all 28 raw p-values ───────────────────────────
            raw_posthoc = {}
            for g1, g2, pair_label in COMPARISON_PAIRS:
                v1, v2 = all_data.get(g1, []), all_data.get(g2, [])
                if len(v1) < 3 or len(v2) < 3: continue
                p, test = run_posthoc(v1, v2, overall['all_normal'])
                raw_posthoc[(g1, g2)] = {
                    'p': p, 'test': test, 'label': pair_label
                }

            # ── Step 2: BH correction across all 28 pairs for this gene ───────
            p_adj = bh_correction(
                [((g1, g2), v['p']) for (g1, g2), v in raw_posthoc.items()]
            )

            # ── Step 3: build final posthoc with corrected sig ────────────────
            posthoc = {}
            for (g1, g2), v in raw_posthoc.items():
                p_bh  = p_adj.get((g1, g2))
                sl_bh = sig_label(p_bh)
                posthoc[(g1, g2)] = {
                    'p':      v['p'],
                    'p_bh':   p_bh,
                    'sig':    sig_label(v['p']),
                    'sig_bh': sl_bh,
                    'test':   v['test'],
                    'label':  v['label'],
                }
                if sl_bh not in ('ns', ''):
                    d = '↑' if np.mean(all_data[g2]) > np.mean(all_data[g1]) else '↓'
                    print(f'  {gene:<8} {v["label"]:<40} '
                          f'p={v["p"]:.4f} p_BH={p_bh:.4f} {sl_bh} {d}')

            posthoc_dict[sheet] = posthoc

            # ── CSV ───────────────────────────────────────────────────────────
            for g1, g2, pair_label in COMPARISON_PAIRS:
                v1, v2 = all_data.get(g1, []), all_data.get(g2, [])
                if not v1 or not v2: continue
                ph = posthoc.get((g1, g2), {})
                m1, m2 = np.mean(v1), np.mean(v2)
                s1, s2 = np.std(v1, ddof=1), np.std(v2, ddof=1)
                all_csv.append({
                    'tissue':           tissue,
                    'gene':             gene,
                    'comparison':       pair_label,
                    'group1':           g1,
                    'group2':           g2,
                    'overall_test':     overall['test'],
                    'overall_p':        overall['p'],
                    'posthoc_test':     ph.get('test', ''),
                    'posthoc_p_raw':    ph.get('p', ''),
                    'posthoc_p_BH':     ph.get('p_bh', ''),
                    'sig_raw':          ph.get('sig', ''),
                    'sig_BH':           ph.get('sig_bh', ''),
                    'direction':        ('↑' if m2 > m1 else '↓') if ph.get('sig_bh', '') not in ('ns', '') else '—',
                    'mean_g1':          round(m1, 4),
                    'sd_g1':            round(s1, 4),
                    'sem_g1':           round(s1 / np.sqrt(len(v1)), 4),
                    'n_g1':             len(v1),
                    'mean_g2':          round(m2, 4),
                    'sd_g2':            round(s2, 4),
                    'sem_g2':           round(s2 / np.sqrt(len(v2)), 4),
                    'n_g2':             len(v2),
                    'outliers_rm_g1':   len(load_vals(filepath, sheet, g1)) - len(v1),
                    'outliers_rm_g2':   len(load_vals(filepath, sheet, g2)) - len(v2),
                })

        make_figure(tissue, ref_dict, box_dict, posthoc_dict,
                    os.path.join(out_dir, f'qPCR_{tissue}'))

    save_csv(all_csv, os.path.join(out_dir, 'qPCR_all_vs_all_stats.csv'))
    print(f'\nDone! → {out_dir}/')
    return out_dir


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--no-outliers', action='store_true')
    args = parser.parse_args()
    run(remove_outliers=not args.no_outliers)