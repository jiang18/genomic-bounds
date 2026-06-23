"""
fig5_gblup_reliability.py
-------------------------
Figure 5: Empirical vs. theoretical in-sample GBLUP reliability.

Single-panel scatter comparing the empirical in-sample reliability (squared
correlation between the GBLUP prediction and the true genetic value, averaged
over 100 replicates) with the predicted R^2_in = 1 - lambda*S(N)/N (Equation 17)
across all 36 (Ne, h2, N) scenarios.

Visual encoding:
  - Colour    = Ne (20, 50, 100)
  - Marker    = h2 (circle 0.1, triangle 0.3, square 0.5)
  - y errorbar = SD of the empirical reliability across the 100 replicates

Data source:
    gblup_reliability.csv  (columns: Ne, h2, N, nrep, R2_random_mean,
    R2_full_mean, R2_full_sd, R2_theory, ratio_full_over_theory), produced by
    msprime/compute_gblup_reliability.py --nrep 100.

Usage:
    python fig5_gblup_reliability.py --csv ../msprime/gblup_reliability.csv \
                                     --out figures/fig5_gblup_reliability.pdf
"""

import argparse
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

import style
from style import COLORS, DOUBLE_COL, DPI, clean_axes

NE_ORDER = [20, 50, 100]
NE_KEYS = {20: 'Ne20', 50: 'Ne50', 100: 'Ne100'}
H2_MARKERS = {0.1: 'o', 0.3: '^', 0.5: 's'}   # circle, triangle, square


def main():
    ap = argparse.ArgumentParser(
        description="Figure 5: empirical vs. theoretical GBLUP reliability")
    ap.add_argument('--csv', default='../msprime/gblup_reliability.csv',
                    help='Summary CSV from compute_gblup_reliability.py')
    ap.add_argument('--out', default='figures/fig5_gblup_reliability.pdf')
    args = ap.parse_args()

    df = pd.read_csv(args.csv)
    df['Ne_int'] = df['Ne'].str.replace('Ne', '', regex=False).astype(int)
    print(f"Loaded {len(df)} scenarios")

    ratio = df['R2_full_mean'] / df['R2_theory']
    ratio_mean, ratio_min, ratio_max = ratio.mean(), ratio.min(), ratio.max()

    fig_w = 0.55 * DOUBLE_COL
    fig, ax = plt.subplots(1, 1, figsize=(fig_w, fig_w))

    for ne in NE_ORDER:
        color = COLORS[NE_KEYS[ne]]
        for h2, marker in H2_MARKERS.items():
            sub = df[(df['Ne_int'] == ne) & (np.isclose(df['h2'], h2))]
            if sub.empty:
                continue
            ax.errorbar(sub['R2_theory'], sub['R2_full_mean'],
                        yerr=sub['R2_full_sd'],
                        fmt=marker, ms=4.5, mfc=color, mec='white', mew=0.3,
                        ecolor=color, elinewidth=0.6, capsize=1.5,
                        alpha=0.85, linestyle='None', zorder=3)

    # Identity line
    lo = min(df['R2_theory'].min(), df['R2_full_mean'].min()) - 0.03
    hi = max(df['R2_theory'].max(), df['R2_full_mean'].max()) + 0.03
    ax.plot([lo, hi], [lo, hi], ls='--', lw=0.8, color='grey', zorder=1)
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_xlabel(r'Predicted $\bar{R}^2_{\mathrm{in}} = 1 - \lambda S(N)/N$ (Eq. 17)')
    ax.set_ylabel(r'Empirical $\bar{R}^2_{\mathrm{in}}$')

    stats_text = (f"Mean(ratio) = {ratio_mean:.2f}\n"
                  r"Bars: $\pm$1 SD (100 reps)" "\n"
                  rf"{len(df)} scenarios = $3\,N_e \times 3\,h^2 \times 4\,N$")
    ax.text(0.03, 0.97, stats_text, transform=ax.transAxes,
            va='top', ha='left', fontsize=6.5,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                      edgecolor='0.7', alpha=0.9))

    ne_handles = [
        mlines.Line2D([], [], color=COLORS[NE_KEYS[ne]], marker='o',
                      linestyle='None', markersize=4.5, label=rf'$N_e = {ne}$')
        for ne in NE_ORDER
    ]
    h2_handles = [
        mlines.Line2D([], [], color='0.35', marker=H2_MARKERS[h2],
                      linestyle='None', markersize=4.5, label=rf'$h^2 = {h2}$')
        for h2 in sorted(H2_MARKERS)
    ]
    ax.legend(handles=ne_handles + h2_handles, loc='lower right',
              fontsize=6.5, markerscale=1.0, handletextpad=0.3,
              labelspacing=0.35, ncol=2, columnspacing=0.8)

    clean_axes(ax)
    ax.set_aspect('equal', adjustable='box')

    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    fig.savefig(args.out, dpi=DPI)
    plt.close(fig)
    print(f"Saved {args.out}  ({len(df)} scenarios, mean ratio = {ratio_mean:.3f})")


if __name__ == '__main__':
    main()
