import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# ---- GLOBAL STYLE (publication-like) ----
plt.rcParams.update({
    "font.size": 16,
    "axes.titlesize": 22,
    "axes.labelsize": 22,
    "xtick.labelsize": 22,
    "ytick.labelsize": 22
})

# ---- YOUR DATA LOADING (unchanged) ----
time_steps = []
temp_rows = []
time_list = []

with open("state.dat") as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        if parts[0] == "TIME":
            DX = float(parts[3][1:])


with open("temp_diag.dat") as f:
    for line in f:
        line = line.strip()
        if not line:
            continue

        parts = line.split()

        # skip headers safely
        if parts[0] in ["STATE", "TIME"]:
            continue

        try:
            t = float(parts[0]) * 1e9  # seconds → ns
            temps = list(map(float, parts[1:]))

            time_list.append(t)
            temp_rows.append(temps)

        except:
            continue
  
n_depth = len(temp_rows[0])  # number of spatial points

depth_list = np.arange(n_depth) * DX * 1e7  # nm

data = np.array(temp_rows).T

t_min = np.min(data)
t_max = np.max(data)

cmap = LinearSegmentedColormap.from_list("blue_red", ["blue", "red"])

# ---- FIGURE SETUP (KEY PART) ----
fig, ax = plt.subplots(
    figsize=(10, 6),
    constrained_layout=True
)

# lock visual proportions
ax.set_box_aspect(0.6)  # matches 10:6 ratio

times = np.array(time_list, dtype=float)
depth = np.array(depth_list, dtype=float)

# ---- MAIN IMAGE ----
im = ax.imshow(
    data,
    aspect='auto',
    origin='upper',
    cmap=cmap,
    extent=[times[0], times[-1], depth[-1], depth[0]]
)

# ---- LABELS ----
ax.set_xlabel("Time (ns)")
ax.set_ylabel("Depth (nm)")
ax.set_title("Temperature Evolution Heatmap")

# ---- COLORBAR (STABLE, LEFT SIDE) ----
cbar = fig.colorbar(
    im,
    ax=ax,
    shrink=0.9
)
cbar.set_label("Temperature")

# ---- ANNOTATION (FIGURE-LOCKED, NO DRIFT) ----
fig.text(
    0.82, 0.88,
    f"Hot (red): {t_max:.2f}\nCold (blue): {t_min:.2f}",
    ha='left',
    va='top',
    fontsize=12,
    bbox=dict(facecolor='white', alpha=0.8, edgecolor='none')
)

# ---- SHOW ----
plt.show()