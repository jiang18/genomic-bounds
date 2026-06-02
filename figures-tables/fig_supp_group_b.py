"""
fig_supp_group_b.py
-------------------
Supplementary Figure S2: S_B(N) vs U_B(N) (Jensen tightness).

Jensen's inequality is tight for Group B (U_B is close to S_B), validating
U_B as a proxy for S_B in the ratio bound derivation in the main text.

Layout: 3 panels side by side (one per Ne = 20, 50, 100), with a shared
y-axis to make the between-Ne comparison direct.  Each panel overlays
three h2 values (0.1, 0.3, 0.5) with the same colour palette as
Fig 3 / Fig 4.  S_B is plotted with solid lines and U_B with dashed
lines of the same colour.

Note: S_B and U_B are *both* identically zero whenever the simulated
N <= Me.  This is a structural artifact of the sample LD matrix
having rank min(N, M): positions beyond Me in the eigenvalue array
are zero-padded, so rho_B = 0 and the Jensen formula returns
U_B = N * 0 * M_B / denom = 0.  The population-level U_B is positive
in this regime; we just cannot estimate it from samples with N <= Me.
This affects the leftmost points of the Ne=50 and Ne=100 panels.

Usage:
    python fig_supp_group_b.py --eigen_dir ../eigen_out \
        --out figures/fig_supp_group_b.pdf
"""

import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.lines as mlines

import style
from style import (LS, DOUBLE_COL, DPI, NE_VALUES,
                   me_theory, compute_S_decomposed,
                   panel_label, clean_axes)

H2_VALUES    = [0.1, 0.3, 0.5]
PANEL_LABELS = ['(A)', '(B)', '(C)']

# Same h2 palette as Fig 3, Fig 4
H2_COLORS = {
    0.1: '#CC6677',   # rose
    0.3: '#44AA99',   # teal
    0.5: '#332288',   # indigo
}


def discover_files(eigen_dir):
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
            SB_vals, UB_vals = [], []
            for n_val in n_values:
                ev = np.load(files_by_ne[ne][n_val])
                ev = np.maximum(ev, 0.0)
                if len(ev) < M:
                    ev = np.concatenate([ev, np.zeros(M - len(ev))])
                Mc = min(Me, len(ev))
                _, _, SB, _, UB = compute_S_decomposed(ev, n_val, h2, M, Mc)
                SB_vals.append(SB)
                UB_vals.append(UB)

            SB_arr = np.array(SB_vals)
            UB_arr = np.array(UB_vals)

            color = H2_COLORS[h2]
            ax.plot(N_arr, SB_arr, color=color, ls='-', lw=1.2, zorder=3)
            ax.plot(N_arr, UB_arr, color=color, ls='--', lw=1.0, zorder=2)

        # X-axis
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
            ax.set_ylabel(r'$S_B(N)$ and $U_B(N)$')

        panel_label(ax, PANEL_LABELS[col_idx], x=-0.20, y=1.10)
        clean_axes(ax)

        # Two-part legend in the first panel only:
        # h2 colours + line-style key (solid = S_B, dashed = U_B)
        if col_idx == 0:
            h2_handles = [
                mlines.Line2D([], [], color=H2_COLORS[h2], ls='-', lw=1.2,
                              label=rf'$h^2={h2}$')
                for h2 in H2_VALUES
            ]
            style_handles = [
                mlines.Line2D([], [], color='0.35', ls='-', lw=1.2,
                              label=r'$S_B$'),
                mlines.Line2D([], [], color='0.35', ls='--', lw=1.0,
                              label=r'$U_B$'),
            ]
            ax.legend(handles=h2_handles + style_handles,
                      loc='upper left', fontsize=6,
                      handlelength=1.5, labelspacing=0.3,
                      ncol=1)

    fig.savefig(out_path, dpi=DPI)
    print(f'Saved: {out_path}')
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Supplementary Figure S2: S_B vs U_B (Jensen tightness).')
    parser.add_argument('--eigen_dir', required=True)
    parser.add_argument('--out', default='figures/fig_supp_group_b.pdf')
    args = parser.parse_args()
    make_figure(args.eigen_dir, args.out)
