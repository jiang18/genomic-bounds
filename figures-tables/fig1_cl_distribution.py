"""
fig1_cl_distribution.py
-----------------------
Figure 1: Distribution of per-SNP GRAMMAR-Gamma coefficients c_l across
three livestock chip datasets, for both full-GRM and LOCO mixed models.

Layout: 2 x 3 grid (double-column width).
  Row 1 (full-GRM c_l):  Chinese Holstein | Karacabey Merino | German Holstein
  Row 2 (LOCO  c_l):     Chinese Holstein | Karacabey Merino | (empty)

Each panel shows kernel density estimates of c_l at h2 = 0.1, 0.3, 0.5,
together with a per-h2 CV annotation in the upper-left corner.  The narrow,
concentrated distributions confirm that the scalar approximation
c_l ≈ c-bar is accurate for standard chip panels.

Data files:
    Per-SNP coefficient TSV files from the GRAMMAR-Gamma analysis.

Usage:
    python fig1_cl_distribution.py \
        --coeff_dir ../grammar-gamma-coeffs/analysis_outputs/grammar_gamma/results \
        --out figures/fig1_cl_distribution.pdf
"""

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

import style
from style import DOUBLE_COL, DPI, clean_axes, panel_label

# Layout: rows indexed 0=full-GRM, 1=LOCO; columns indexed 0,1,2.
# Each entry is (dataset key, label, file basename, has_loco)
DATASETS = {
    'chinese_holstein': {
        'name':     'Chinese Holstein',
        'subtitle': r'($N = 2{,}510$ and $M = 42{,}775$)',
        'full':     'chinese_holstein_coeffs.tsv',
        'loco':     'chinese_holstein_loco_coeffs.tsv',
    },
    'karacabey_merino': {
        'name':     'Karacabey Merino',
        'subtitle': r'($N = 734$ and $M = 35{,}774$)',
        'full':     'karacabey_merino_coeffs.tsv',
        'loco':     'karacabey_merino_loco_coeffs.tsv',
    },
    'german_holstein': {
        'name':     'German Holstein',
        'subtitle': r'($N = 5{,}024$ and $M = 42{,}216$)',
        'full':     'german_holstein_coeffs.tsv',
        'loco':     None,    # no chromosome map for LOCO
    },
}

# Position grid: row, col -> (dataset key, mode 'full'/'loco')
GRID = {
    (0, 0): ('chinese_holstein', 'full'),
    (0, 1): ('karacabey_merino', 'full'),
    (0, 2): ('german_holstein',  'full'),
    (1, 0): ('chinese_holstein', 'loco'),
    (1, 1): ('karacabey_merino', 'loco'),
    (1, 2): None,
}

PANEL_LABELS = {
    (0, 0): '(A)', (0, 1): '(B)', (0, 2): '(C)',
    (1, 0): '(D)', (1, 1): '(E)',
}

H2_VALUES = [0.1, 0.3, 0.5]

# One colour per h2 (same palette as S(N) and other figures)
H2_COLORS = {
    0.1: '#CC6677',   # rose
    0.3: '#44AA99',   # teal
    0.5: '#332288',   # indigo
}


def plot_panel(ax, df, mode):
    """
    Plot KDEs for the three h2 values on *ax*, using c_l columns from *df*.

    mode : 'full' uses c_h2_<h2>; 'loco' uses c_loco_h2_<h2>.

    Returns dict {h2: cv} for the legend-cum-annotation.
    """
    cv_by_h2 = {}
    prefix = 'c_h2_' if mode == 'full' else 'c_loco_h2_'

    for h2 in H2_VALUES:
        col = f'{prefix}{h2}'
        vals = df[col].dropna().values
        color = H2_COLORS[h2]
        cv = float(np.std(vals) / np.mean(vals))
        cv_by_h2[h2] = cv

        kde = gaussian_kde(vals, bw_method=0.15)
        x_grid = np.linspace(
            max(0, np.percentile(vals, 0.5)),
            np.percentile(vals, 99.5),
            300,
        )
        density = kde(x_grid)

        ax.fill_between(x_grid, density, alpha=0.25, color=color, zorder=2)
        ax.plot(x_grid, density, color=color, lw=1.0, zorder=3)
        ax.axvline(np.mean(vals), color=color, ls='--', lw=0.7,
                   alpha=0.7, zorder=2)

    return cv_by_h2


def add_topleft_legend(ax, cv_by_h2):
    """Stacked legend-style annotation in the top-left of the panel.

    Layout:
        h2     CV
        ■ 0.1   0.08
        ■ 0.3   0.13
        ■ 0.5   0.16
    """
    base_x, base_y = 0.04, 0.96
    line_dy = 0.10
    annot_fs = 6.5

    # Header row
    ax.text(base_x, base_y, r'$h^2$', transform=ax.transAxes,
            fontsize=annot_fs, va='top', ha='left', color='0.2')
    ax.text(base_x + 0.21, base_y, r'CV', transform=ax.transAxes,
            fontsize=annot_fs, va='top', ha='left', color='0.2')

    for i, h2 in enumerate(H2_VALUES, start=1):
        y = base_y - i * line_dy
        # Coloured square swatch
        ax.text(base_x, y, '■', transform=ax.transAxes,
                fontsize=annot_fs, va='top', ha='left',
                color=H2_COLORS[h2])
        # h2 value
        ax.text(base_x + 0.045, y, rf'{h2}', transform=ax.transAxes,
                fontsize=annot_fs, va='top', ha='left', color='0.15')
        # CV value (no "CV=" prefix)
        ax.text(base_x + 0.21, y, rf'{cv_by_h2[h2]:.2f}',
                transform=ax.transAxes,
                fontsize=annot_fs, va='top', ha='left', color='0.15')


def make_figure(coeff_dir, out_path):
    fig_width  = DOUBLE_COL
    fig_height = DOUBLE_COL * 0.62           # taller for 2-row grid
    fig, axes  = plt.subplots(2, 3, figsize=(fig_width, fig_height),
                              layout='constrained')

    # Cache loaded data frames
    loaded = {}
    for ds_key, ds in DATASETS.items():
        loaded[ds_key] = {}
        full_path = os.path.join(coeff_dir, ds['full'])
        loaded[ds_key]['full'] = pd.read_csv(full_path, sep='\t')
        if ds['loco']:
            loco_path = os.path.join(coeff_dir, ds['loco'])
            loaded[ds_key]['loco'] = pd.read_csv(loco_path, sep='\t')

    for (row, col), entry in GRID.items():
        ax = axes[row, col]

        if entry is None:
            ax.set_visible(False)
            continue

        ds_key, mode = entry
        df = loaded[ds_key][mode]

        cv_by_h2 = plot_panel(ax, df, mode)
        add_topleft_legend(ax, cv_by_h2)

        # Two-line title: name+mode on line 1, sample sizes on line 2.
        # Capitalised at title start (matches journal convention).
        mode_tag = 'Full-GRM' if mode == 'full' else 'LOCO'
        title = (rf'{mode_tag}: {DATASETS[ds_key]["name"]}'
                 + '\n'
                 + DATASETS[ds_key]['subtitle'])
        ax.set_title(title, fontsize=7, pad=4)

        # Axis labels
        if row == 1 or (row == 0 and col == 2):
            ax.set_xlabel(r'Per-SNP coefficient $c_\ell$')
        if col == 0:
            ax.set_ylabel('Density')

        ax.set_xlim(left=0)
        ax.set_ylim(bottom=0)

        panel_label(ax, PANEL_LABELS[(row, col)], x=-0.15, y=1.12)
        clean_axes(ax)

    os.makedirs(os.path.dirname(out_path) or '.', exist_ok=True)
    fig.savefig(out_path, dpi=DPI)
    print(f'Saved: {out_path}')
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Figure 1: Distribution of per-SNP GRAMMAR-Gamma '
                    'coefficients in three livestock chip datasets '
                    '(full-GRM and LOCO).')
    parser.add_argument(
        '--coeff_dir',
        default='../grammar-gamma-coeffs/analysis_outputs/grammar_gamma/results',
        help='Directory containing per-SNP coefficient TSV files')
    parser.add_argument(
        '--out', default='figures/fig1_cl_distribution.pdf',
        help='Output file path')
    args = parser.parse_args()
    make_figure(args.coeff_dir, args.out)
