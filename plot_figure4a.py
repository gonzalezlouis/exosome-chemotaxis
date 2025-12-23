import matplotlib.pyplot as plt
import numpy as np
import scienceplots
import os

plt.style.use('science')
plt.rcParams.update({"font.size": 18})

# Parameters
D_values = [0.0, 50.0, 250.0, 500.0, 1000.0]
colors = ['gray', 'orange', 'red', 'blue', 'green', 'purple']  # Match D_values length
labels = [rf'$D$ = {int(D)}' for D in D_values]
cell_radius = 10

# Plot setup
fig, ax = plt.subplots(figsize=(10, 6))
ax.set_aspect('equal')

# Plot each follower trajectory
for D, color, label in zip(D_values, colors, labels):
    filename = f"follower_trajectory_{int(D)}.txt"
    if not os.path.exists(filename):
        print(f"Warning: File {filename} not found. Skipping.")
        continue

    data = np.loadtxt(filename)
    if data.ndim < 2 or data.shape[1] < 2:
        print(f"Warning: File {filename} does not contain valid 2D trajectory data.")
        continue

    x, y = data[:, 0], data[:, 1]
    x_shifted = x - x[0]
    y_shifted = y - y[0]
    ax.plot(x_shifted * 60, y_shifted * 60, color=color, label=label, linewidth=2)

# Draw initial follower cell as a filled circle
origin_circle = plt.Circle((0, 0), cell_radius, edgecolor='black',
                           facecolor='skyblue', linewidth=1.5, zorder=10)
ax.add_patch(origin_circle)

# Remove frame and ticks for aesthetics
for spine in ax.spines.values():
    spine.set_visible(False)
ax.set_facecolor('white')
fig.patch.set_facecolor('white')

# Axis labels and legend
ax.set_xlabel(r'$X$ ($\mu$m)')
ax.set_ylabel(r'$Y$ ($\mu$m)')
ax.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), frameon=False)

# Save figure
plt.savefig("../../fig/fig-final/figure4a.pdf", dpi=300, bbox_inches='tight')
plt.show()
