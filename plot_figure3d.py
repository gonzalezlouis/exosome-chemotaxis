import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scienceplots

plt.style.use('science')
plt.rcParams.update({"font.size": 24})


# Parameters
v0 = 3.4 / 60       # µm/min, source speed
K = 25             # threshold molecules
tau_values = [1, 5, 25, 50]  # memory time constants (min)

# Define the colormap
cmap = cm.get_cmap("rainbow", len(tau_values))
colors = [cmap(i) for i in range(len(tau_values))]


# File pattern (adjust if needed)
filename_template = "results_2d_tau_{}.txt"

# Function to load data from file
def load_simulation_data(filename):
    with open(filename, "r") as file:
        lines = file.readlines()

    n_vals = []
    v_vals = []
    for line in lines:
        line = line.strip()
        if not line or line.startswith(("K =", "H =", "tau")):
            continue
        parts = line.split()
        if len(parts) == 2:
            n_vals.append(float(parts[0]))
            v_vals.append(float(parts[1]))
    return np.array(n_vals), np.array(v_vals)

plt.figure(figsize=(12, 8))

for idx,tau in enumerate(tau_values):
    filename = filename_template.format(int(tau))
    try:
        n_bar, velocities = load_simulation_data(filename)
    except FileNotFoundError:
        print(f"Warning: File {filename} not found, skipping.")
        continue

    # Normalize x-axis: n_bar / K
    x = n_bar / K

    # Normalize y-axis: v / v0
    y = velocities / v0

    plt.plot(x, y, marker='o', label=f"$\\tau = {tau}$ min", color=colors[idx])

plt.xlabel(r"Scaled mean cargo size $x = \bar{n} / K$")
plt.ylabel(r"Scaled velocity $\bar{v} / v_0$")
plt.ylim(0,1)
plt.legend(fontsize = 20, loc='best')
plt.tight_layout()
plt.savefig("../../fig/fig-final/figure3d.pdf", dpi=300)
plt.show()
