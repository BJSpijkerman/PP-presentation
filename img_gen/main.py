import matplotlib.pyplot as plt
import numpy as np

# -----------------------------
# LaTeX + scientific formatting
# -----------------------------
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern"],

    # --- PDF font embedding ---
    "pdf.fonttype": 42,   # Embed TrueType fonts (recommended)
    "ps.fonttype": 42,

    # Avoid Type 3 fonts (journals often reject these)
    "text.latex.preamble": r"\usepackage{amsmath}",
})

# -----------------------------
# Parameters
# -----------------------------
N_modes = 3
L = 100e-6
max_wavelength = 420e-9
frequency_in = 25e12

k_max = 2 * np.pi / max_wavelength
omega_c_in = frequency_in / (2.99e8)

# -----------------------------
# Data
# -----------------------------
k_x = np.linspace(0, k_max/100, 10**5)
vacuum_omega_c = k_x

omega_c = []
mode_lines = []

for i in range(N_modes):
    mode_number = i + 1
    omega_c_mode = np.sqrt(k_x**2 + (mode_number * np.pi / L)**2)
    mode_line = mode_number * np.pi / L

    omega_c.append(omega_c_mode)
    mode_lines.append(mode_line)

# -----------------------------
# Plot
# -----------------------------
fig, ax = plt.subplots(figsize=(7, 7))

# Vacuum dispersion (reference)
ax.plot(
    k_x,
    vacuum_omega_c,
    linestyle="--",
    linewidth=2.5,
    color="black",
    label=r"Vacuum: $\omega/c = k_x$"
)

# Mode dispersions
for i in range(N_modes):
    ax.plot(
        k_x,
        omega_c[i],
        linewidth=2,
        label=rf"Mode $m={i+1}$"
    )

for i in range(2):
    ax.hlines(
        mode_lines[i],
        0,
        k_max/100,
        linestyle="--",
        color="black",
        alpha=0.2
    )

# -----------------------------
# Excited modes
# -----------------------------
idx_max_mode = abs(np.array(mode_lines) - omega_c_in).argmin()

k_x_excited = []
for i in range(idx_max_mode):
    idx_excited = abs(omega_c[i] - omega_c_in).argmin()
    kx_val = k_x[idx_excited]

    k_x_excited.append(kx_val)

    ax.scatter(
        kx_val,
        omega_c_in,
        zorder=5
    )

# Input frequency line
ax.hlines(
    omega_c_in,
    0,
    k_max/100,
    linestyle="--",
    linewidth=2,
    label=r"Input $\omega_{\mathrm{in}}/c$"
)

ax.vlines(
    k_x_excited,
    0,
    omega_c_in,
    linestyle="--",
    alpha=0.5
)

# -----------------------------
# Axes formatting
# -----------------------------
ax.set_xlim(0, k_max/100)
ax.set_ylim(0, (N_modes + 1)*np.pi/L)

# Y ticks: mode cutoffs
ax.set_yticks(mode_lines)
ax.set_yticklabels([
    rf"$\frac{{{i+1}\pi}}{{L}}$" for i in range(N_modes)
])

# Remove x ticks for cleaner look
ax.set_xticks([])

ax.set_xlabel(r"$k_x$", fontsize=18)
ax.set_ylabel(r"$\omega / c$", fontsize=18)

# Subtle grid
#ax.grid(alpha=0.15)

# Legend
ax.legend(frameon=False)

# No propagation label
ax.text(
    0.7, 0.1,
    r"No propagation",
    transform=ax.transAxes,
    ha="left",
    va="top"
)

# Single mode propagation label
ax.text(
    0.7, 0.4,
    r"Single mode propagation",
    transform=ax.transAxes,
    ha="left",
    va="top"
)

# Multi mode propagation label
ax.text(
    0.7, 0.6,
    r"Multi mode propagation",
    transform=ax.transAxes,
    ha="left",
    va="top"
)

# Layout
fig.tight_layout()
plt.show()

# Save as publication-quality PDF
fig.savefig(
    "dispersion_diagram.pdf",
    format="pdf",
    bbox_inches="tight",
    pad_inches=0.02
)

plt.show()
