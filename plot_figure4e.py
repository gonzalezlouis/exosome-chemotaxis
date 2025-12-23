import numpy as np
import matplotlib.pyplot as plt
import os
import scienceplots

plt.style.use('science')
plt.rcParams.update({"font.size": 28})

# Parameters
v0 = 3.4 / 60       # µm/min, source speed
K = 25             # threshold molecules
output_dir = "../../fig/fig-final"

# Load simulation data from the updated protocol comparison file
nbar_values = []
data_A = []  # Constant D
data_B = []  # Stokes-Einstein D

with open("results_protocol_comparison.txt", "r") as file:
    lines = file.readlines()
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith("n_bar ="):
            nbar = float(line.split("=")[1].strip(" :\n"))
            nbar_values.append(nbar)

            # Next two lines contain ConstantD and StokesEinsteinD results
            line_A = lines[i + 1].strip()
            line_B = lines[i + 2].strip()

            parts_A = line_A.split()
            parts_B = line_B.split()
            assert parts_A[0] == "ConstantD"
            assert parts_B[0] == "StokesEinsteinD"

            data_A.append((float(parts_A[1]), float(parts_A[2])))
            data_B.append((float(parts_B[1]), float(parts_B[2])))

            i += 4  # skip blank line too
        else:
            i += 1

nbar_values = np.array(nbar_values)
data_A = np.array(data_A)
data_B = np.array(data_B)

means_A, errs_A = data_A[:, 0], data_A[:, 1]
means_B, errs_B = data_B[:, 0], data_B[:, 1]

# Plot
plt.figure(figsize=(12, 8))
plt.errorbar(nbar_values, means_A * 60, yerr=errs_A * 60, fmt='o-', label=r'Constant $D$', capsize=4)
plt.errorbar(nbar_values, means_B * 60, yerr=errs_B * 60, fmt='s--', label=r'Stokes-Einstein $D(\bar{n})$', capsize=4)

plt.xlabel(r'$\bar{n}$ (mean molecules per exosome)')
plt.ylabel(r'Average Follower Velocity $\bar{v}$ ($\mu$m/hr)')
plt.legend(frameon=False)
plt.ylim(0,v0 * 60)
# plt.savefig('../fig/fig-final/figure4c.pdf', dpi=300)
plt.savefig(os.path.join(output_dir, "figure4c.pdf"), dpi=300)
plt.show()

