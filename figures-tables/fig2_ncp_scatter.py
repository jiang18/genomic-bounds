"""
fig2_ncp_scatter.py
--------------------
Figure 2: Theoretical vs. empirical non-centrality parameter.

Single-panel scatter plot comparing per-SNP theoretical NCP (using the
eigenvalue-specific coefficient c_l) to the empirical NCP (mean chi-squared
minus 1) across all focal-SNP simulation scenarios.

Visual encoding:
  - Colour  = Ne (20, 50, 100)
  - Marker  = h2 (circle 0.1, triangle 0.3, square 0.5)
Each Ne contributes 120 points: 3 h2 x 2 beta2 x 4 N x 5 focal SNPs
(one SNP per MAF bin per chromosome), for 360 total.

Data source:
    plink_out/Ne{ne}/h2_{h2}__beta2_{beta2}/N_{N}/summary.csv

Usage:
    python fig2_ncp_scatter.py --plink_dir ../plink_out \
                               --out figures/fig2_ncp_scatter.pdf
"""

import argparse
import os
import re

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

import style
from style import COLORS, DOUBLE_COL, DPI, clean_axes

# ---------------------------------------------------------------------------
# Regex patterns for parsing the directory hierarchy
# ---------------------------------------------------------------------------
RE_NE    = re.compile(r'Ne(\d+)')
RE_SCEN  = re.compile(r'h2_([\d.]+)__beta2_([\d.]+)')
RE_N     = re.compile(r'N_(\d+)')

# ---------------------------------------------------------------------------
# Visual-encoding maps
# ---------------------------------------------------------------------------
NE_ORDER = [20, 50, 100]
NE_KEYS  = {20: 'Ne20', 50: 'Ne50', 100: 'Ne100'}

H2_MARKERS = {0.1: 'o', 0.3: '^', 0.5: 's'}   # circle, triangle, square


def collect_summaries(plink_dir):
    """
    Walk *plink_dir* and load every ``summary.csv``, annotating each row
    with Ne, h2, beta2, and N parsed from the enclosing directory names.

    Returns
    -------
    pd.DataFrame
        Combined data with added columns Ne, h2, beta2, N.
    """
    frames = []

    for ne_dir in sorted(os.listdir(plink_dir)):
        m_ne = RE_NE.match(ne_dir)
        if not m_ne:
            continue
        ne = int(m_ne.group(1))
        ne_path = os.path.join(plink_dir, ne_dir)
        if not os.path.isdir(ne_path):
            continue

        for scenario_dir in sorted(os.listdir(ne_path)):
            m_sc = RE_SCEN.match(scenario_dir)
            if not m_sc:
                continue
            h2    = float(m_sc.group(1))
            beta2 = float(m_sc.group(2))
            scenario_path = os.path.join(ne_path, scenario_dir)
            if not os.path.isdir(scenario_path):
                continue

            for n_dir in sorted(os.listdir(scenario_path)):
                m_n = RE_N.match(n_dir)
                if not m_n:
                    continue
                N = int(m_n.group(1))
                csv_path = os.path.join(scenario_path, n_dir, 'summary.csv')
                if not os.path.isfile(csv_path):
                    continue

                df = pd.read_csv(csv_path)
                df['Ne']    = ne
                df['h2']    = h2
                df['beta2'] = beta2
                df['N']     = N
                frames.append(df)

    if not frames:
        raise FileNotFoundError(
            f"No summary.csv files found under {plink_dir}")

    return pd.concat(frames, ignore_index=True)


def main():
    ap = argparse.ArgumentParser(
        description="Figure 2: Theoretical vs. empirical NCP scatter")
    ap.add_argument('--plink_dir', required=True,
                    help='Path to plink_out/ directory')
    ap.add_argument('--out', default='figures/fig2_ncp_scatter.pdf',
                    help='Output file path')
    args = ap.parse_args()

    # ------------------------------------------------------------------
    # Load and clean data
    # ------------------------------------------------------------------
    df = collect_summaries(args.plink_dir)
    df = df.dropna(subset=['ncp_theoretical_cl', 'ncp_empirical', 'ncp_ratio'])
    df = df[df['ncp_theoretical_cl'] > 0].copy()
    print(f"Loaded {len(df)} focal-SNP observations")

    # ------------------------------------------------------------------
    # Summary statistics from ncp_ratio = empirical / theoretical_cl
    # ------------------------------------------------------------------
    ratio_mean = df['ncp_ratio'].mean()
    ratio_sd   = df['ncp_ratio'].std()
    n_total    = len(df)

    # ------------------------------------------------------------------
    # Figure: single square panel
    # ------------------------------------------------------------------
    fig_w = 0.55 * DOUBLE_COL                       # ~3.86 in
    fig, ax = plt.subplots(1, 1, figsize=(fig_w, fig_w))

    # Plot each Ne x h2 combination separately to get distinct markers
    for ne in NE_ORDER:
        color = COLORS[NE_KEYS[ne]]
        for h2, marker in H2_MARKERS.items():
            sub = df[(df['Ne'] == ne) & (np.isclose(df['h2'], h2))]
            if sub.empty:
                continue
            ax.scatter(sub['ncp_theoretical_cl'], sub['ncp_empirical'],
                       s=14, alpha=0.6, linewidths=0.3, edgecolors='white',
                       color=color, marker=marker, zorder=3)

    # Identity line (grey dashed, from 0 to max)
    xy_max = max(df['ncp_theoretical_cl'].max(),
                 df['ncp_empirical'].max()) * 1.05
    ax.plot([0, xy_max], [0, xy_max], ls='--', lw=0.8,
            color='grey', zorder=1)

    # Axis limits and labels
    ax.set_xlim(0, xy_max)
    ax.set_ylim(0, xy_max)
    ax.set_xlabel(r'Theoretical NCP (per-SNP $c_\ell$, Eq. 1)')
    ax.set_ylabel(r'Empirical NCP (mean $\chi^2 - 1$)')

    # Stats annotation --- upper left, white rounded box
    stats_text = (f"Mean(ratio) = {ratio_mean:.2f}\n"
                  f"SD(ratio) = {ratio_sd:.2f}\n"
                  rf"{n_total} points = $3\,N_e \times 3\,h^2 \times "
                  rf"2\,\beta^2 \times 4\,N \times 5$ SNPs")
    ax.text(0.03, 0.97, stats_text, transform=ax.transAxes,
            va='top', ha='left', fontsize=6.5,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                      edgecolor='0.7', alpha=0.9))

    # ------------------------------------------------------------------
    # Two-part legend: Ne (colour) + h2 (marker shape)
    # ------------------------------------------------------------------
    # Ne colour handles
    ne_handles = [
        mlines.Line2D([], [], color=COLORS[NE_KEYS[ne]], marker='o',
                       linestyle='None', markersize=4.5,
                       label=rf'$N_e = {ne}$')
        for ne in NE_ORDER
    ]
    # h2 marker-shape handles (grey to emphasise shape, not colour)
    h2_handles = [
        mlines.Line2D([], [], color='0.35', marker=H2_MARKERS[h2],
                       linestyle='None', markersize=4.5,
                       label=rf'$h^2 = {h2}$')
        for h2 in sorted(H2_MARKERS)
    ]

    all_handles = ne_handles + h2_handles
    ax.legend(handles=all_handles, loc='lower right',
              fontsize=6.5, markerscale=1.0,
              handletextpad=0.3, labelspacing=0.35,
              ncol=2, columnspacing=0.8)

    clean_axes(ax)
    ax.set_aspect('equal', adjustable='box')

    # ------------------------------------------------------------------
    # Save
    # ------------------------------------------------------------------
    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    fig.savefig(args.out, dpi=DPI)
    plt.close(fig)
    print(f"Saved {args.out}  ({n_total} points, "
          f"mean ratio = {ratio_mean:.3f})")


if __name__ == '__main__':
    main()
