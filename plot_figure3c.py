import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scienceplots

plt.style.use('science')
plt.rcParams.update({"font.size": 24})

# Parameters
K = 25
H = 3
v0 = 3.4 / 60.0
nu_values = [60, 100, 250, 600]

# Define the colormap
cmap = cm.get_cmap("rainbow", len(nu_values))
colors = [cmap(i) for i in range(len(nu_values))]


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

for idx,nu in enumerate(nu_values):
    fname = f"results_2d_nu{int(nu)}.txt"
    try:
        n_bar, velocities = load_simulation_data(fname)
    except FileNotFoundError:
        print(f"Warning: File {fname} not found, skipping.")
        continue

    x = n_bar / K 
    y = velocities / v0
    plt.plot(x, y, 'o', linestyle='-', label=rf"$\nu={nu}$",color=colors[idx])



plt.xlabel(r"Scaled mean cargo size $x = \bar{n} / K$")
plt.ylabel(r"Scaled velocity $\bar{v} / v_0$")
plt.ylim(0,1)
plt.tight_layout()
plt.legend(loc='lower center', bbox_to_anchor=(0.25, 0.0),  fontsize=24, frameon=False)
plt.savefig("../../fig/fig-final/figure3c.pdf", dpi=300)
plt.show()

