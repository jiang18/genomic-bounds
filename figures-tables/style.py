"""
style.py
--------
Shared matplotlib style settings and helper constants for all figures.
Targets Genetics (GSA/OUP) figure conventions:
  - Sans-serif font (Arial with Helvetica fallback)
  - Single-column width: 8.7 cm = 3.43 in
  - Double-column width: 17.8 cm = 7.01 in
  - Ticks inward, minor ticks visible
  - No top/right spines
  - Legend without frame
  - 300 DPI output
  - Colorblind-safe palette (Paul Tol's muted set)
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------------------
# Figure dimensions (inches)
# ---------------------------------------------------------------------------
SINGLE_COL = 3.43   # 8.7 cm
DOUBLE_COL = 7.01   # 17.8 cm
# Heights are set per figure; a common aspect ratio is ~0.75 x width per panel

# ---------------------------------------------------------------------------
# Colorblind-safe palette (Paul Tol's muted, 6 colours used here)
# ---------------------------------------------------------------------------
# Usage: COLORS['Ne20'], COLORS['Ne50'], COLORS['Ne100']
#        COLORS['S'], COLORS['SA'], COLORS['SB']
COLORS = {
    'Ne20':  '#CC6677',   # rose
    'Ne50':  '#44AA99',   # teal
    'Ne100': '#332288',   # indigo
    # component colours for sigmoid sum decomposition
    'S':     '#222222',   # near-black (total)
    'SA':    '#0077BB',   # blue (Group A)
    'SB':    '#EE7733',   # orange (Group B)
    'UB':    '#EE7733',   # same hue as SB, dashed
    'bound': '#BBBBBB',   # grey (theoretical bound)
}

# Line styles
LS = {
    'solid':       '-',
    'dashed':      '--',
    'dotted':      ':',
    'dashdotted':  '-.',
}

DPI = 300

# ---------------------------------------------------------------------------
# rcParams — apply once at import
# ---------------------------------------------------------------------------
RC = {
    # Font
    'font.family':          'sans-serif',
    'font.sans-serif':      ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size':            7,
    'axes.titlesize':       7,
    'axes.labelsize':       7,
    'xtick.labelsize':      7,
    'ytick.labelsize':      7,
    'legend.fontsize':      7,
    # Ticks
    'xtick.direction':      'in',
    'ytick.direction':      'in',
    'xtick.minor.visible':  True,
    'ytick.minor.visible':  True,
    'xtick.major.size':     3.0,
    'xtick.minor.size':     1.5,
    'ytick.major.size':     3.0,
    'ytick.minor.size':     1.5,
    'xtick.major.width':    0.6,
    'ytick.major.width':    0.6,
    # Axes
    'axes.linewidth':       0.6,
    'axes.spines.top':      False,
    'axes.spines.right':    False,
    # Legend
    'legend.frameon':       False,
    'legend.handlelength':  1.5,
    'legend.handletextpad': 0.4,
    # Lines
    'lines.linewidth':      1.2,
    # Math text — use dejavusans (always available, clean sans-serif)
    'mathtext.fontset':     'dejavusans',
    # Font embedding for PDF/PS — Type 42 = TrueType embedded (required for journals)
    'pdf.fonttype':         42,
    'ps.fonttype':          42,
    # Output
    'savefig.dpi':          DPI,
    'savefig.bbox':         'tight',
    'savefig.pad_inches':   0.02,
}

mpl.rcParams.update(RC)


# ---------------------------------------------------------------------------
# Helper: clean spines (belt-and-suspenders given rcParams)
# ---------------------------------------------------------------------------
def clean_axes(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


# ---------------------------------------------------------------------------
# Panel label helper  (A), (B), (C) ...
# ---------------------------------------------------------------------------
def panel_label(ax, label, x=-0.18, y=1.05, fontsize=8, fontweight='bold'):
    ax.text(x, y, label, transform=ax.transAxes,
            fontsize=fontsize, fontweight=fontweight,
            va='top', ha='left')


# ---------------------------------------------------------------------------
# Population parameters
# ---------------------------------------------------------------------------
# Genome length in Morgans (30 chromosomes × ~1 Morgan each, typical bovine)
L_MORGANS = 25.0

NE_VALUES = [20, 50, 100]

def me_theory(Ne):
    """Theoretical Me = 4 * Ne * L."""
    return int(4 * Ne * L_MORGANS)


# ---------------------------------------------------------------------------
# Sigmoid sum components
# ---------------------------------------------------------------------------
def compute_S(eigenvalues, N, h2, M):
    """
    Compute S(N), S_A(N), S_B(N) for a given eigenvalue array.

    Parameters
    ----------
    eigenvalues : np.ndarray
        Sorted descending eigenvalues of the LD matrix R = Z'Z/N.
        Length = min(N_sim, M) where N_sim is the sample size used
        in the eigenvalue computation. Here we reweight for target N.
    N : int or float
        Target sample size for the sigmoid sum.
    h2 : float
        Heritability.
    M : int
        Total number of SNPs (so that alpha = M*(1-h2)/h2).

    Returns
    -------
    S, S_A, S_B : floats
    """
    alpha = M * (1.0 - h2) / h2
    lam = eigenvalues  # already LD eigenvalues
    # Clip negatives from numerical noise
    lam = np.maximum(lam, 0.0)
    sigm = N * lam / (N * lam + alpha)
    S = sigm.sum()
    return S


def compute_S_decomposed(eigenvalues, N, h2, M, Me):
    """
    Return S, S_A, S_B, U_A, U_B for a given N.

    Me : int
        Number of Group A eigenvalues (top Me).
    """
    alpha  = M * (1.0 - h2) / h2
    lam    = np.maximum(eigenvalues, 0.0)
    sigm   = N * lam / (N * lam + alpha)

    S      = sigm.sum()
    S_A    = sigm[:Me].sum()
    S_B    = sigm[Me:].sum()

    # Jensen upper bounds
    lam_lambda = (1.0 - h2) / h2   # = alpha / M
    rho    = lam[:Me].sum() / M
    rho_B  = lam[Me:].sum() / M
    MB     = M - Me

    # U_A
    denom_A = N * rho + Me * lam_lambda
    U_A = (N * rho * Me / denom_A) if denom_A > 0 else Me

    # U_B
    if MB > 0:
        denom_B = N * rho_B + MB * lam_lambda
        U_B = (N * rho_B * MB / denom_B) if denom_B > 0 else 0.0
    else:
        U_B = 0.0

    return S, S_A, S_B, U_A, U_B


def compute_rho(eigenvalues, Me):
    """Fraction of total eigenvalue mass in top Me eigenvalues."""
    total = eigenvalues.sum()
    if total == 0:
        return np.nan
    return eigenvalues[:Me].sum() / total


def compute_NA_star(Me, h2, rho):
    """N_A* = Me * lambda / rho,  lambda = (1-h2)/h2."""
    lam = (1.0 - h2) / h2
    return Me * lam / rho


def compute_NB_star(M, Me, h2, rho):
    """N_B* = M_B * lambda / (1-rho)."""
    lam  = (1.0 - h2) / h2
    MB   = M - Me
    rho_B = 1.0 - rho
    if rho_B <= 0:
        return np.inf
    return MB * lam / rho_B
