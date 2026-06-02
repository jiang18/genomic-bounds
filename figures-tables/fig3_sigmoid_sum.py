"""
fig3_sigmoid_sum.py
--------------------
Figure 3: Saturation of the sigmoid sum S(N) at Me (top row), and
the (N, h^2) exchange in S(N) (bottom row).

Layout: 2 x 3 grid (double-column width).
  Row 1 (S vs N at fixed h^2):   one panel per Ne (20, 50, 100),
                                  curves for h^2 = 0.1, 0.3, 0.5.
  Row 2 (S vs h^2 at fixed N):   one panel per Ne (20, 50, 100),
                                  curves for N = 10,000; 20,000; 50,000
                                  (typical DRP-cohort scale in livestock).
                                  h^2 grid is restricted to [0.05, 0.95]
                                  to cover realistic phenotype heritabilities
                                  and DRP reliabilities.

Each panel marks the theoretical ceiling Me = 4 * Ne * L as a
horizontal dashed line. A light grey horizontal reference line at
Me/2 in the bottom row aids reading the (N, h^2) trade-off:
points where the three N-curves cross this line are (N, h^2) pairs
with equal S(N).

The bottom-row curves use the N-specific eigenvalue spectrum at each
fixed N (the same per-N spectra as the top row).

Usage:
    python fig3_sigmoid_sum.py --eigen_dir ../eigen_out --out figures/fig3_sigmoid_sum.pdf
"""

import argparse
import os
import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

import style
from style import (DOUBLE_COL, DPI, NE_VALUES, me_theory,
                   compute_S_decomposed, compute_rho, panel_label, clean_axes, LS)

H2_VALUES        = [0.1, 0.3, 0.5]
N_FIXED_VALUES   = [10_000, 20_000, 50_000]   # bottom-row curves; typical DRP-cohort scale
PANEL_LABELS_TOP = ['(A)', '(B)', '(C)']
PANEL_LABELS_BOT = ['(D)', '(E)', '(F)']

# One colour per h2 (top row)
H2_COLORS = {
    0.1: '#CC6677',   # rose
    0.3: '#44AA99',   # teal
    0.5: '#332288',   # indigo
}

# One colour per N (bottom row); reuse the same hue family for visual continuity
N_COLORS = {
    10_000:  '#CC6677',   # rose
    20_000:  '#44AA99',   # teal
    50_000:  '#332288',   # indigo
}


def discover_files(eigen_dir):
    """Scan eigen_dir/Ne{ne}/ for eigenvalue .npy files and return nested dict."""
    data = {}
    for ne in NE_VALUES:
        ne_dir = os.path.join(eigen_dir, f'Ne{ne}')
        if not os.path.isdir(ne_dir):
            continue
        data[ne] = {}
        for fname in sorted(os.listdir(ne_dir)):
            m = re.match(r'eigenvalues_N(\d+)\.npy$', fname)
            if m:
                n_val = int(m.group(1))
                data[ne][n_val] = os.path.join(ne_dir, fname)
    return data


def infer_M(files_by_ne):
    """For each Ne, M = length of eigenvalue array at the largest available N."""
    M_by_ne = {}
    for ne, files in files_by_ne.items():
        largest_N = max(files.keys())
        eigs = np.load(files[largest_N])
        M_by_ne[ne] = len(eigs)
    return M_by_ne


def nice_yticks(Me):
    """
    Return [0, rounded_half, Me] where rounded_half is a round number
    close to Me / 2 for clean y-axis labelling.
    """
    half = Me / 2
    for step in [5000, 2000, 1000, 500, 200, 100]:
        candidate = round(half / step) * step
        if 0 < candidate < Me:
            return [0, candidate, Me]
    return [0, Me]


def _fmt_N(x, _pos):
    """Format sample size as e.g. '5k', '10k', '200k'."""
    if x >= 1000:
        return f'{int(x / 1000)}k'
    return f'{int(x)}'


def _load_eigs(path, M):
    """Load eigenvalue array, clip negatives, zero-pad to length M."""
    ev = np.load(path)
    ev = np.maximum(ev, 0.0)
    if len(ev) < M:
        ev = np.concatenate([ev, np.zeros(M - len(ev))])
    return ev


def plot_top_panel(ax, ne, M, Me, files_by_ne):
    """Top-row: S(N) vs N for three h^2 values, using per-N eigenvalue spectra."""
    n_values = sorted(files_by_ne[ne].keys())
    N_arr    = np.array(n_values, dtype=float)

    x_lo = N_arr.min() * 0.8
    x_hi = N_arr.max() * 1.3

    for h2 in H2_VALUES:
        S_vals = []
        for n_val in n_values:
            ev = _load_eigs(files_by_ne[ne][n_val], M)
            Mc = min(Me, len(ev))
            S, _, _, _, _ = compute_S_decomposed(ev, n_val, h2, M, Mc)
            S_vals.append(S)
        ax.plot(N_arr, np.array(S_vals),
                color=H2_COLORS[h2], ls=LS['solid'], lw=1.3,
                label=rf'$h^2={h2}$', zorder=3)

    ax.axhline(Me, color='#444444', ls=LS['dashed'], lw=0.9, zorder=2)
    ax.set_xscale('log')
    ax.set_xlim(x_lo, x_hi)
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(_fmt_N))

    ymax = Me * 1.12
    ax.set_ylim(0, ymax)
    ticks = nice_yticks(Me)
    ax.set_yticks(ticks)
    ax.set_yticklabels([str(t) for t in ticks])

    ax.set_xlabel(r'Sample size $N$')


def plot_bottom_panel(ax, ne, M, Me, files_by_ne):
    """Bottom-row: S(N) vs h^2 for three fixed N values, using per-N eigenvalue spectra."""
    h2_grid = np.linspace(0.05, 0.95, 200)

    # Only plot N values that are available in this Ne's eigenvalue files
    available_N = sorted(files_by_ne[ne].keys())
    n_used      = [n for n in N_FIXED_VALUES if n in available_N]

    for n_val in n_used:
        ev = _load_eigs(files_by_ne[ne][n_val], M)
        Mc = min(Me, len(ev))
        S_vals = []
        for h2 in h2_grid:
            S, _, _, _, _ = compute_S_decomposed(ev, n_val, h2, M, Mc)
            S_vals.append(S)
        ax.plot(h2_grid, np.array(S_vals),
                color=N_COLORS[n_val], ls=LS['solid'], lw=1.3,
                label=rf'$N={n_val:,}$', zorder=3)

    # Reference lines: Me ceiling, and a half-Me guide
    ax.axhline(Me, color='#444444', ls=LS['dashed'], lw=0.9, zorder=2)
    ax.axhline(Me / 2, color='#888888', ls=':', lw=0.7, zorder=1)

    ax.set_xlim(0.05, 0.95)
    ax.set_xticks([0.1, 0.3, 0.5, 0.7, 0.9])
    ax.set_xticklabels(['0.1', '0.3', '0.5', '0.7', '0.9'])

    ymax = Me * 1.12
    ax.set_ylim(0, ymax)
    ticks = nice_yticks(Me)
    ax.set_yticks(ticks)
    ax.set_yticklabels([str(t) for t in ticks])

    ax.set_xlabel(r'Heritability $h^2$')


def make_figure(eigen_dir, out_path):
    files_by_ne = discover_files(eigen_dir)
    M_by_ne     = infer_M(files_by_ne)

    fig_width  = DOUBLE_COL
    fig_height = DOUBLE_COL * 0.72

    fig, axes = plt.subplots(2, 3, figsize=(fig_width, fig_height),
                             layout='constrained')

    for col_idx, ne in enumerate(NE_VALUES):
        ax_top = axes[0, col_idx]
        ax_bot = axes[1, col_idx]

        if ne not in files_by_ne:
            ax_top.set_visible(False)
            ax_bot.set_visible(False)
            continue

        M    = M_by_ne[ne]
        Me = min(me_theory(ne), M)

        plot_top_panel(ax_top, ne, M, Me, files_by_ne)
        plot_bottom_panel(ax_bot, ne, M, Me, files_by_ne)

        ax_top.set_title(rf'$N_e = {ne}$  ($M_e={Me:,}$)', pad=3)
        if col_idx == 0:
            ax_top.set_ylabel(r'Sigmoid sum $S(N)$')
            ax_bot.set_ylabel(r'Sigmoid sum $S(N)$')

        panel_label(ax_top, PANEL_LABELS_TOP[col_idx], x=-0.18, y=1.08)
        panel_label(ax_bot, PANEL_LABELS_BOT[col_idx], x=-0.18, y=1.08)
        clean_axes(ax_top)
        clean_axes(ax_bot)

        # Legends in first column only
        if col_idx == 0:
            ax_top.legend(loc='upper left', fontsize=6,
                          handlelength=1.4, labelspacing=0.3,
                          bbox_to_anchor=(0.02, 0.85))
            ax_bot.legend(loc='upper left', fontsize=6,
                          handlelength=1.4, labelspacing=0.3,
                          bbox_to_anchor=(0.02, 0.85))

    fig.savefig(out_path, dpi=DPI)
    print(f'Saved: {out_path}')
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Figure 3: S(N) saturation at Me (top) and (N, h^2) exchange (bottom).')
    parser.add_argument('--eigen_dir', required=True,
                        help='Path to eigen_out directory')
    parser.add_argument('--out', default='figures/fig3_sigmoid_sum.pdf',
                        help='Output file path (default: figures/fig3_sigmoid_sum.pdf)')
    args = parser.parse_args()
    make_figure(args.eigen_dir, args.out)
