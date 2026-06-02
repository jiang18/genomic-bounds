"""
fig4_cbar_divergence.py
-----------------------
Figure 4 (main): Divergence of full-GRM and LOCO average GRAMMAR-Gamma
coefficients with sample size.

Layout: 3 panels side by side (double-column width), one per Ne (20, 50, 100).

Each panel plots the average GRAMMAR-Gamma coefficient (cbar) against sample
size N on a log scale.  For each heritability value (h2 = 0.1, 0.3, 0.5):
  - cbar_full  is shown as a dashed line (approaches 0 as N grows),
  - cbar_loco  is shown as a solid line with markers (climbs toward an
    asymptotic upper bound as N grows),
  - loco_upper_bound is shown as a thin dotted horizontal line in each h2 colour.

A light grey horizontal line at y = 1 marks the full-GRM / LOCO boundary.

Data files are pre-computed CSV tables in ../msprime/:
    ../msprime/cbar_Ne20.csv, ../msprime/cbar_Ne50.csv, ../msprime/cbar_Ne100.csv

Usage:
    python fig4_cbar_divergence.py --out figures/fig4_cbar_divergence.pdf
"""

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.lines as mlines

import style
from style import DOUBLE_COL, DPI, clean_axes, panel_label

# ── Configuration ──────────────────────────────────────────────────────────
NE_VALUES    = [20, 50, 100]
H2_VALUES    = [0.1, 0.3, 0.5]
PANEL_LABELS = ['(A)', '(B)', '(C)']

# One colour per h2 value (Paul Tol muted set, same hues as Ne palette)
H2_COLORS = {
    0.1: '#CC6677',   # rose
    0.3: '#44AA99',   # teal
    0.5: '#332288',   # indigo
}


def load_data(ne, script_dir):
    """Load cbar CSV for a given Ne.  Returns a pandas DataFrame."""
    path = os.path.join(script_dir, '..', 'msprime', f'cbar_Ne{ne}.csv')
    return pd.read_csv(path)


def make_figure(out_path):
    script_dir = os.path.dirname(os.path.abspath(__file__))

    fig_width  = DOUBLE_COL
    fig_height = DOUBLE_COL * 0.45                # raised for legend room
    fig, axes  = plt.subplots(1, 3, figsize=(fig_width, fig_height),
                              layout='constrained')

    for col_idx, ne in enumerate(NE_VALUES):
        ax = axes[col_idx]
        df = load_data(ne, script_dir)

        # Reference line at y = 1
        ax.axhline(1.0, color='#DDDDDD', lw=0.6, zorder=0)

        for h2 in H2_VALUES:
            sub = df[np.isclose(df['h2'], h2)].sort_values('N')
            N_arr = sub['N'].values
            color = H2_COLORS[h2]

            # LOCO upper bound (constant across N for each h2) -- dotted, heavy
            bound = sub['loco_upper_bound'].iloc[0]
            ax.axhline(bound, color=color, lw=1.2, ls=':',
                       alpha=0.7, zorder=1)

            # Full-GRM cbar -- dashed, no markers (original style)
            ax.plot(N_arr, sub['cbar_full'].values,
                    color=color, ls='--', lw=1.0, zorder=2)

            # LOCO cbar — solid with markers
            ax.plot(N_arr, sub['cbar_loco'].values,
                    color=color, ls='-', lw=1.2,
                    marker='o', markersize=3, zorder=3,
                    label=rf'$h^2={h2}$' if col_idx == 0 else None)

        # ── Axes formatting ────────────────────────────────────────────
        ax.set_xscale('log')
        all_N = df['N'].unique()
        ax.set_xlim(min(all_N) * 0.7, max(all_N) * 1.5)
        ax.xaxis.set_major_formatter(ticker.FuncFormatter(
            lambda x, _: f'{int(x/1000)}k' if x >= 1000 else f'{int(x)}'))
        ax.xaxis.set_minor_formatter(ticker.NullFormatter())

        ax.set_ylim(0, 2.6)
        ax.set_xlabel(r'Sample size $N$')
        if col_idx == 0:
            ax.set_ylabel('Average coefficient')

        ax.set_title(rf'$N_e = {ne}$', pad=3)
        panel_label(ax, PANEL_LABELS[col_idx], x=-0.18, y=1.08)
        clean_axes(ax)

    # ── Legend (first panel only) ──────────────────────────────────────────
    ax0 = axes[0]

    # Collect h2-colour handles already placed by the LOCO plot calls
    h2_handles, h2_labels = ax0.get_legend_handles_labels()

    # Style legend entries
    loco_handle = mlines.Line2D([], [], color='grey', ls='-', lw=1.0,
                                marker='o', markersize=2.5, label='LOCO')
    full_handle = mlines.Line2D([], [], color='grey', ls='--', lw=1.0,
                                label='Full-GRM')
    bound_handle = mlines.Line2D([], [], color='grey', ls=':', lw=1.2,
                                 alpha=0.7, label='LOCO bound')

    all_handles = h2_handles + [loco_handle, full_handle, bound_handle]
    all_labels  = h2_labels  + ['LOCO', 'Full-GRM', 'LOCO bound']

    ax0.legend(handles=all_handles, labels=all_labels,
               loc='upper left', fontsize=6,
               handlelength=1.6, labelspacing=0.3,
               ncol=2, columnspacing=0.8)

    # ── Save ──────────────────────────────────────────────────────────────
    fig.savefig(out_path, dpi=DPI)
    print(f'Saved: {out_path}')
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Figure 4: Divergence of full-GRM and LOCO '
                    'GRAMMAR-Gamma coefficients with sample size.')
    parser.add_argument('--out', default='figures/fig4_cbar_divergence.pdf',
                        help='Output file path (default: figures/fig4_cbar_divergence.pdf)')
    args = parser.parse_args()
    make_figure(args.out)
