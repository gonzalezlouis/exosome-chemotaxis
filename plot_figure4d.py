import numpy as np
import matplotlib.pyplot as plt
import scienceplots

v_max = 3.4 # um/hr (input as per cpp but label says um/hr usually? Check units)
# In cpp: v = 3.4/60 (um/min). So v_max in um/min is 0.056. 
# But the plot Y-axis is usually scaled to um/hr or normalized.
# The previous plot multiplied mean by 60 -> um/hr. 
# 3.4/60 * 60 = 3.4 um/hr.

plt.style.use('science')
plt.rcParams.update({"font.size":28})     # specify font size here

data = {}
current_label = None

with open("results_walls.txt", "r") as file:
    lines = file.readlines()
    D_vals, v_means, v_stds = [], [], []
    
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith("n_bar"):
            continue
        
        # Check for new block
        if line.startswith("Wall Dist"):
            # Save previous block
            if current_label is not None:
                data[current_label] = (np.array(D_vals), np.array(v_means), np.array(v_stds))
                D_vals, v_means, v_stds = [], [], []
            
            # Parse new label
            # Format: "Wall Dist = X cells:" or "Wall Dist = Infinity:"
            label_part = line.split("=")[1].strip(" :\n")
            current_label = label_part
            
        else:
            # Data lines
            parts = line.split()
            if len(parts) == 3:
                D_vals.append(float(parts[0]))
                v_means.append(float(parts[1]))
                v_stds.append(float(parts[2]))

    # Save last block
    if current_label is not None:
        data[current_label] = (np.array(D_vals), np.array(v_means), np.array(v_stds))

# --- Plotting ---
plt.figure(figsize=(12, 8))

# Define colors manually or via map
colors = ['black', 'firebrick', 'orange', 'green'] # Inf, 100, 10, 5
labels_map = {
    "Infinity": "Unbounded",
    "100 cells": r"$L = 200a$",
    "10 cells": r"$L = 20a$",
    "5 cells":  r"$L = 10a$"
}

# Sort keys to ensure consistent order: Infinity first, then descending size
sorted_keys = ["Infinity", "100 cells", "10 cells", "5 cells"]

for idx, key in enumerate(sorted_keys):
    if key not in data: continue
    
    Ds, means, errs = data[key]
    
    # Convert velocity to um/hr (multiply by 60)
    plt.errorbar(Ds, means * 60, yerr=errs * 60, 
                 fmt='o-', label=labels_map[key], 
                 color=colors[idx], capsize=4, markersize=5, linewidth=1.5)

# Vertical Lines (Stokes Law)
eta = 0.001
kBT = 4.1e-21
def stokes_D(a): # radius in meters
    D_m2_s = kBT / (6 * np.pi * eta * a)
    return D_m2_s * 1e12 * 60 # um^2/min

D_min = stokes_D(150e-9)
D_max = stokes_D(30e-9)

# Add vertical lines for D_min and D_max
plt.axvline(D_min, color='gray', linestyle='--', linewidth=1.5, label='Experimental range')
plt.axvline(D_max, color='gray', linestyle='--', linewidth=1.5)
plt.axvspan(D_min, D_max, color='gray', alpha=0.1)

plt.xlabel(r'Exosome Diffusivity $D$ ($\mu$m$^2$/min)')
plt.ylabel(r'Average Follower Velocity $\bar{v}$ ($\mu$m/hr)')
plt.legend(frameon=False, fontsize=18)
# plt.title(r'Effect of Confinement ($T=1000$ min)')
# plt.xlim(0, 1000)

plt.tight_layout()
plt.savefig("figure4d.pdf", dpi=300)
plt.show()