# import numpy as np
# import matplotlib.pyplot as plt
# import os
# import scienceplots

# v = 3.4 # um/s reference

# try:
#     plt.style.use('science')
# except:
#     pass # Fallback if scienceplots isn't installed
# plt.rcParams.update({"font.size":16})

# output_dir = "."
# if not os.path.exists(output_dir):
#     os.makedirs(output_dir)

# # Load simulation data
# data = {}
# l0_list = []
# D_star_list = []
# max_v_list = []

# filename = "results_vs_D_l0.txt" # New filename

# if not os.path.exists(filename):
#     print(f"Error: {filename} not found. Run the C++ simulation first.")
#     exit()

# with open(filename, "r") as file:
#     lines = file.readlines()
#     l0_current = None
#     D_vals, v_means, v_stds = [], [], []
    
#     for line in lines:
#         line = line.strip()
#         if not line:
#             continue
#         if line.startswith("n_bar") or line.startswith("T_fixed"):
#             continue
#         elif line.startswith("l0_factor ="):
#             if l0_current is not None:
#                 # Store previous block
#                 data[l0_current] = (np.array(D_vals), np.array(v_means), np.array(v_stds))
            
#             l0_current = float(line.split("=")[1].strip(" :"))
#             D_vals, v_means, v_stds = [], [], []
#         else:
#             parts = line.split()
#             if len(parts) == 3:
#                 D_vals.append(float(parts[0]))
#                 v_means.append(float(parts[1]))
#                 v_stds.append(float(parts[2]))
    
#     # Store last block
#     if l0_current is not None:
#         data[l0_current] = (np.array(D_vals), np.array(v_means), np.array(v_stds))




# # Analyze Optimal D*
# print("Analyzing Optimal D*...")
# sorted_l0 = sorted(data.keys())
# for l0_factor in sorted_l0:
#     Ds, means, errs = data[l0_factor]
    
#     idx_max = np.argmax(means)
#     D_star = Ds[idx_max]
#     v_max = means[idx_max]
    
#     # --- CRITICAL FIX FOR UNITS ---
#     # Convert the factor (e.g., 8) to actual microns (80.0)
#     real_distance_um = l0_factor * 10.0 
    
#     l0_list.append(real_distance_um)      # Append ACTUAL distance
#     D_star_list.append(D_star)
#     max_v_list.append(v_max)
#     print(f"Dist={real_distance_um} um: Max V={v_max:.4f} at D={D_star}")

# # Plot 1: v vs D for different l0
# plt.figure(figsize=(10, 6))
# colors = plt.cm.viridis(np.linspace(0, 1, len(data)))
# for idx, l0 in enumerate(sorted_l0):
#     Ds, means, errs = data[l0]
#     plt.plot(Ds, means*60, 'o-', label=f'$l_0 = {int(l0)}$ sep', color=colors[idx], markersize=4)
#     # plt.fill_between(Ds, (means-errs)*60, (means+errs)*60, color=colors[idx], alpha=0.1)

# plt.xlabel(r'Diffusivity $D$ ($\mu$m$^2$/min)')
# plt.ylabel(r'Avg Velocity $\bar{v}$ ($\mu$m/hr)')
# plt.legend(frameon=False, title="Initial Separation")
# plt.title("Velocity vs Diffusivity for varying separation")
# plt.tight_layout()
# plt.savefig(os.path.join(output_dir, "v_vs_D_l0.pdf"), dpi=300)
# plt.show()

# # Plot 2: Optimal D* vs l0
# plt.figure(figsize=(8, 6))
# plt.plot(l0_list, D_star_list, 's--', color='black', markersize=8, label='Simulation')

# # Linear Fit Check
# slope, intercept = np.polyfit(l0_list, D_star_list, 1)
# fit_line = [slope * x + intercept for x in l0_list]
# plt.plot(l0_list, fit_line, 'r-', alpha=0.6, label=f'Linear Fit (slope={slope:.2f})')

# # Update Label to reflect reality
# plt.xlabel(r'Initial Separation $l_0$ ($\mu$m)')  # Units are now correct
# plt.ylabel(r'Optimal Diffusivity $D^*_{ex}$ ($\mu$m$^2$/min)')
# plt.title(r'Optimal Diffusivity vs Separation')
# plt.legend()
# plt.grid(True, linestyle='--', alpha=0.5)
# plt.tight_layout()
# plt.savefig(os.path.join(output_dir, "D_star_vs_l0.pdf"), dpi=300)
# plt.show()


import numpy as np
import matplotlib.pyplot as plt
import os
import scienceplots

v = 3.4/60 # um/min reference
a = 10 # cell radius, in um

try:
    plt.style.use('science')
except:
    pass
plt.rcParams.update({"font.size":28})

output_dir = "."
filename = "results_repeats.txt"

if not os.path.exists(filename):
    print(f"Error: {filename} not found.")
    exit()

# Data Structure: data[l0_factor][repeat_index] = (D_array, v_array)
data = {}

# --- 1. Parse Data with Repeats ---
with open(filename, "r") as file:
    lines = file.readlines()
    l0_current = None
    rep_current = None
    D_vals, v_means = [], []
    
    for line in lines:
        line = line.strip()
        if not line: continue
        if line.startswith("n_bar") or line.startswith("T_fixed"): continue
        
        # New separation block
        elif line.startswith("l0_factor ="):
            l0_current = float(line.split("=")[1].strip(" :"))
            if l0_current not in data:
                data[l0_current] = {}
        
        # New repeat block
        elif line.startswith("repeat ="):
            # Save previous repeat if exists
            if rep_current is not None and len(D_vals) > 0:
                data[l0_current][rep_current] = (np.array(D_vals), np.array(v_means))
            
            rep_current = int(line.split("=")[1].strip())
            D_vals, v_means = [], []
            
        # Data lines
        else:
            parts = line.split()
            if len(parts) >= 2:
                D_vals.append(float(parts[0]))
                v_means.append(float(parts[1]))

    # Save final block
    if l0_current is not None and rep_current is not None:
         data[l0_current][rep_current] = (np.array(D_vals), np.array(v_means))

# --- 2. Calculate Statistics for D* ---
print("Analyzing Statistics...")
l0_list_um = []
D_star_mean = []
D_star_err = []

sorted_l0 = sorted(data.keys())

for l0_fac in sorted_l0:
    repeats = data[l0_fac]
    best_Ds = []
    
    # For each repeat, find the optimal D
    for r_idx, (Ds, vs) in repeats.items():
        idx_max = np.argmax(vs)
        best_Ds.append(Ds[idx_max])
    
    # Calculate Mean and StdErr of the optimal D
    best_Ds = np.array(best_Ds)
    mu = np.mean(best_Ds)
    se = np.std(best_Ds, ddof=1) / np.sqrt(len(best_Ds))
    
    real_dist = l0_fac * 10.0 # Convert factor to microns
    l0_list_um.append(real_dist)
    D_star_mean.append(mu)
    D_star_err.append(se)
    
    print(f"Dist {real_dist} um: D* = {mu:.2f} +/- {se:.2f} (N={len(best_Ds)})")

# --- 3. Plot with Error Bars ---
plt.figure(figsize=(12, 8))

plt.errorbar(l0_list_um, D_star_mean, yerr=D_star_err, fmt='s-', color='black', 
             capsize=5, label='Simulation (Mean $\pm$ SE)')

# # Linear Fit (on mean values)
# l0_arr = np.array(l0_list_um)
# D_arr = np.array(D_star_mean)
# slope, intercept = np.polyfit(l0_arr, D_arr, 1)
# fit_line = slope * l0_arr + intercept

# plt.plot(l0_arr, fit_line, 'r--', alpha=0.6, label=f'Linear Fit (slope={slope:.2f})')

# Fitting analytics

l0_arr = np.linspace(min(l0_list_um), max(l0_list_um), 100)
D_star = v/(4*a) * (l0_arr + a)**2

plt.plot(l0_arr, D_star, 'r--', lw = 4, label=f'Analytic Prediction')


plt.xlabel(r'Initial Separation $l_0$ ($\mu$m)') # Explicit units
plt.ylabel(r'Optimal Diffusivity $D^*$ ($\mu$m$^2$/min)')
# plt.title(r'Optimal Diffusivity vs Separation ($T=24$h)')
plt.legend(frameon=False, fontsize=18)
# plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()
plt.savefig("D_star_vs_l0_with_errors.pdf", dpi=300)
plt.show()