import numpy as np
import matplotlib.pyplot as plt
import os
import scienceplots

v = 3.4 # um/s

plt.style.use('science')
plt.rcParams.update({"font.size":28})          # specify font size here

output_dir = "../../fig/fig-final"

# Load simulation data
data = {}
with open("results_vs_D_T.txt", "r") as file:
    lines = file.readlines()
    T_current = None
    D_vals, v_means, v_stds = [], [], []
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith("n_bar ="):
            n_bar_fixed = float(line.split("=")[1])
        elif line.startswith("T ="):
            if T_current is not None:
                data[T_current] = (np.array(D_vals), np.array(v_means), np.array(v_stds))
            T_current = int(float(line.split("=")[1].strip(" :\n")))
            D_vals, v_means, v_stds = [], [], []
        else:
            parts = line.split()
            if len(parts) == 3:
                D_vals.append(float(parts[0]))
                v_means.append(float(parts[1]))
                v_stds.append(float(parts[2]))
    if T_current is not None:
        data[T_current] = (np.array(D_vals), np.array(v_means), np.array(v_stds))

# Constants for vertical lines (Stokes' Law)
eta = 0.001  # Pa.s
kBT = 4.1e-21  # J
radius_min = 30e-9      # meters
radius_max = 150e-9     # meters

def stokes_D(a):
    D_m2_per_s = kBT / (6 * np.pi * eta * a)
    D_um2_per_min = D_m2_per_s * 1e12 * 60  # Convert to μm²/min
    return D_um2_per_min

D_min = stokes_D(radius_max)
D_max = stokes_D(radius_min)

# Plot
plt.figure(figsize=(12, 8))
colors = plt.cm.viridis(np.linspace(0, 1, len(data)))
for idx, (T, (Ds, means, errs)) in enumerate(sorted(data.items())):
    plt.errorbar(Ds, means*60, yerr=errs*60, fmt='o-', label=f'T = {round(T/60)} hours', capsize=4, color=colors[idx])

# Add vertical lines for D_min and D_max
plt.axvline(D_min, color='gray', linestyle='--', linewidth=1.5, label='Experimental range')
plt.axvline(D_max, color='gray', linestyle='--', linewidth=1.5)

# Shaded region (optional)
plt.axvspan(D_min, D_max, color='gray', alpha=0.1)

plt.xlabel(r'Exosome Diffusivity $D$ ($\mu$m$^2$/min)')
plt.ylabel(r'Average Follower Velocity $\bar{v}$ ($\mu$m/hr)')
# plt.title(rf'Follower velocity vs diffusion constant $(\bar{{n}} = {int(n_bar_fixed)})$')
plt.legend(frameon=False, fontsize=18)
plt.ylim(0,v)
plt.tight_layout()
# plt.savefig(f'../fig/fig-final/figure4b.pdf', dpi=300)
plt.savefig(os.path.join(output_dir, "figure4b.pdf"), dpi=300)
plt.show()