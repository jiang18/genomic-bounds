"""
fig_supp_sb_sa_ratio.py
-----------------------
Supplementary Figure S1: S_B(N) / S_A(N) ratio across sample sizes.

Layout: 3 panels side by side (double-column width), one per Ne (20, 50, 100).
Each panel has three lines, one per heritability h2 = 0.1, 0.3, 0.5.
Y-axis is shared across the three panels so that the between-Ne contrast
is immediate.

Result: S_B/S_A is small for Ne=50 and Ne=100 (< ~0.01), and rises with N
for Ne=20 to ~0.11 at h2=0.5, N=100,000 --- still small enough that the
practical ceiling is set by the top Me eigenmodes, but visible at the
smallest Ne and largest h2.

Note: S_B/S_A is identically zero whenever the simulated N <= Me, since
the LD-matrix eigenvalue array has length min(N, M) and contains no
'Group B' eigenvalues below the Me cutoff in that case.  This affects
the leftmost points of the Ne=50 and Ne=100 panels.

Usage:
    python fig_supp_sb_sa_ratio.py --eigen_dir ../eigen_out \
        --out figures/fig_supp_sb_sa_ratio.pdf
"""

import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

import style
from style import (DOUBLE_COL, DPI, NE_VALUES, me_theory,
                   compute_S_decomposed, panel_label, clean_axes, LS)

H2_VALUES    = [0.1, 0.3, 0.5]
PANEL_LABELS = ['(A)', '(B)', '(C)']

# One colour per h2 (matches Fig 3, Fig 4)
H2_COLORS = {
    0.1: '#CC6677',   # rose
    0.3: '#44AA99',   # teal
    0.5: '#332288',   # indigo
}


def discover_files(eigen_dir):
    """Scan eigen_dir for eigenvalue .npy files, grouped by Ne and N."""
    data = {}
    for ne in NE_VALUES:
        ne_dir = os.path.join(eigen_dir, f'Ne{ne}')
        if not os.path.isdir(ne_dir):
            continue
        data[ne] = {}
        for fname in sorted(os.listdir(ne_dir)):
            if fname.startswith('eigenvalues_N') and fname.endswith('.npy'):
                n_str = fname.replace('eigenvalues_N', '').replace('.npy', '')
                try:
                    data[ne][int(n_str)] = os.path.join(ne_dir, fname)
                except ValueError:
                    pass
    return data


def infer_M(files_by_ne):
    """Infer total SNP count M per Ne from the eigenvalue array at the largest N."""
    M_by_ne = {}
    for ne, files in files_by_ne.items():
        largest_N = max(files.keys())
        eigs = np.load(files[largest_N])
        M_by_ne[ne] = len(eigs)
    return M_by_ne


def make_figure(eigen_dir, out_path):
    files_by_ne = discover_files(eigen_dir)
    M_by_ne     = infer_M(files_by_ne)

    fig_width  = DOUBLE_COL
    fig_height = DOUBLE_COL * 0.40

    fig, axes = plt.subplots(1, 3, figsize=(fig_width, fig_height),
                             layout='constrained', sharey=True)

    for col_idx, ne in enumerate(NE_VALUES):
        ax = axes[col_idx]

        if ne not in files_by_ne:
            ax.set_visible(False)
            continue

        M    = M_by_ne[ne]
        Me = min(me_theory(ne), M)

        n_values = sorted(files_by_ne[ne].keys())
        N_arr    = np.array(n_values, dtype=float)

        for h2 in H2_VALUES:
            SB_vals, SA_vals = [], []
            for n_val in n_values:
                ev = np.load(files_by_ne[ne][n_val])
                ev = np.maximum(ev, 0.0)
                if len(ev) < M:
                    ev = np.concatenate([ev, np.zeros(M - len(ev))])
                Mc = min(Me, len(ev))
                _, SA, SB, _, _ = compute_S_decomposed(ev, n_val, h2, M, Mc)
                SA_vals.append(SA)
                SB_vals.append(SB)

            SA_arr = np.array(SA_vals)
            SB_arr = np.array(SB_vals)
            ratio  = SB_arr / np.where(SA_arr > 0, SA_arr, np.nan)

            ax.plot(N_arr, ratio,
                    color=H2_COLORS[h2], ls=LS['solid'], lw=1.2,
                    marker='o', markersize=2.5,
                    label=rf'$h^2={h2}$', zorder=3)

        # X-axis: log scale
        ax.set_xscale('log')
        ax.set_xlim(N_arr.min() * 0.8, N_arr.max() * 1.3)
        ax.xaxis.set_major_formatter(ticker.FuncFormatter(
            lambda x, _: f'{int(x/1000)}k' if x >= 1000 else f'{int(x)}'))

        # Y-axis is shared across panels but show labels on every panel
        # so readers can eyeball values without tracing back to panel A.
        ax.set_ylim(bottom=0)
        ax.tick_params(axis='y', labelleft=True)

        ax.set_title(rf'$N_e = {ne}$  ($M_e = {Me:,}$)', pad=3)
        ax.set_xlabel(r'Sample size $N$')
        if col_idx == 0:
            ax.set_ylabel(r'$S_B(N)\,/\,S_A(N)$')

        panel_label(ax, PANEL_LABELS[col_idx], x=-0.18, y=1.08)
        clean_axes(ax)

        # Legend in first panel only
        if col_idx == 0:
            ax.legend(loc='upper left', fontsize=6,
                      handlelength=1.5, labelspacing=0.3)

    fig.savefig(out_path, dpi=DPI)
    print(f'Saved: {out_path}')
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Supplementary Figure S1: S_B/S_A ratio.')
    parser.add_argument('--eigen_dir', required=True,
                        help='Path to eigen_out directory')
    parser.add_argument('--out', default='figures/fig_supp_sb_sa_ratio.pdf',
                        help='Output file path')
    args = parser.parse_args()
    make_figure(args.eigen_dir, args.out)
