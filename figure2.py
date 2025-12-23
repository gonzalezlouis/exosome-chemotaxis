import numpy as np
import matplotlib.pyplot as plt
import os
import scienceplots

plt.style.use('science')
plt.rcParams.update({"font.size": 26})

# Fixed parameters
H_values = [1, 3]
T_values = [0.5, 5, 60]  # Decay times (min)
data_dir = "../../dat/"
output_dir = "../../fig/fig-final"
n_bar_label = r'$\bar{n}$ (mean molecules per exosome)'
vel_label = r'Average Follower Velocity $\bar{v}$ ($\mu$m/hr)'

# --- Step 1: Find global min and max velocities ---
vmin, vmax = np.inf, -np.inf
for H in H_values:
    for T_decay in T_values:
        filename = os.path.join(data_dir, f"results_2d_H{H}_T{int(T_decay)}.txt")
        with open(filename, "r") as f:
            for line in f:
                line = line.strip()
                if any(line.startswith(prefix) for prefix in ["K =", "H =", "tau =", "nu ="]):
                    continue
                parts = line.split()
                if len(parts) == 2:
                    v = float(parts[1]) * 60.0  # Convert to μm/hr
                    vmin = min(vmin, v)
                    vmax = max(vmax, v)

# Optional: Add padding
padding = 0.05 * (vmax - vmin)
vmin -= padding
vmax += padding

# --- Step 2: Plot using fixed y-limits ---
for H in H_values:
    fig, ax = plt.subplots(figsize=(12, 8))

    for T_decay in T_values:
        filename = os.path.join(data_dir, f"results_2d_H{H}_T{int(T_decay)}.txt")

        n_vals, v_vals = [], []
        with open(filename, "r") as f:
            for line in f:
                line = line.strip()
                if any(line.startswith(prefix) for prefix in ["K =", "H =", "tau =", "nu ="]):
                    continue
                parts = line.split()
                if len(parts) == 2:
                    n_vals.append(float(parts[0]))
                    v_vals.append(float(parts[1]) * 60.0)

        n_vals = np.array(n_vals)
        v_vals = np.array(v_vals)

        plt.plot(n_vals, v_vals, marker='o', linestyle='None', label=f'$\\mathcal{{T}}={T_decay}$ min')

    plt.xlabel(n_bar_label)
    plt.ylabel(vel_label)
    # ax.set_yscale('log')
    plt.title(f'Hill coefficient $H={H}$')
    plt.ylim(vmin, vmax)
    plt.legend(fontsize=18,markerscale=3)
    plt.tight_layout()

    fig_suffix = 'a' if H == 1 else 'b'
    plt.savefig(os.path.join(output_dir, f"figure2{fig_suffix}.pdf"), dpi=300)

plt.show()

