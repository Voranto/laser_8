import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, ListedColormap

# ---- GLOBAL STYLE ----

# Set to -1 if you want to use the usual range (from minimum value found to maximum value found)
DEPTH_MIN = -1
DEPTH_MAX = -1
TIME_MIN = -1
TIME_MAX = -1

MANUAL_LUT_TEMP_MIN = False
MANUAL_LUT_TEMP_MAX = False
LUT_MIN = -1
LUT_MAX = -1
LUT_FILENAME = "fire.lut"

plt.rcParams.update({
    "font.size": 16,
    "axes.titlesize": 22,
    "axes.labelsize": 22,
    "xtick.labelsize": 22,
    "ytick.labelsize": 22
})


# ---- PRINTING INFORMATION ----
print("The imported LUT file that's being used is: ", LUT_FILENAME)
if MANUAL_LUT_TEMP_MIN:
    print("The minimum conductivity for the LUT has been manually set, and it's value is: ", LUT_MIN)
else:
    print("The minimum conductivity for the LUT will be the minimum value found in the data")

if MANUAL_LUT_TEMP_MAX:
    print("The maximum conductivity for the LUT has been manually set, and it's value is: ", LUT_MAX)
else:
    print("The maximum conductivity for the LUT will be the maximum value found in the data")

print("The xtick label size is: ", plt.rcParams.get("xtick.labelsize"))
print("The ytick label size is: ", plt.rcParams.get("ytick.labelsize"))


def load_imagej_lut(path):
    with open(path, "rb") as f:
        data = np.frombuffer(f.read(), dtype=np.uint8)
    if len(data) < 768:
        raise ValueError("Not a valid ImageJ LUT (too small)")
    data = data[:768].reshape((3, 256)).T  # (256, 3)
    return data / 255.0

lut_rgb = load_imagej_lut(LUT_FILENAME)
LUT_CMAP = ListedColormap(lut_rgb)

# ---- DATA LOADING ----
temp_rows = []
time_list = [0]

with open("state.dat") as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        if parts[0] == "TIME":
            DX = float(parts[3][1:])

with open("conductivity_diag.dat") as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        if parts[0] in ["STATE", "TIME"]:
            continue
        try:
            t = float(parts[0]) * 1e9  # seconds → ns
            temps = list(map(float, parts[1:]))
            time_list.append(t)
            temp_rows.append(temps)
        except:
            continue

depth_list = []
with open("depth_vector.dat") as f:
    for line in f:
        if not line:
            continue
        depth_list.append(float(line.strip()))
n_depth = len(depth_list)

# ---- PROCESS DATA ----
data = np.array(temp_rows).T          # shape: (n_depth, n_time)

times = np.array(time_list, dtype=float)
depth = np.array(depth_list, dtype=float)

# ---- APPLY CUSTOM RANGES ----
time_min_plot = TIME_MIN if TIME_MIN != -1 else times.min()
time_max_plot = TIME_MAX if TIME_MAX != -1 else times.max()
depth_min_plot = DEPTH_MIN if DEPTH_MIN != -1 else depth[:-12].min()
depth_max_plot = DEPTH_MAX if DEPTH_MAX != -1 else depth[:-12].max()
print("The range for the time is: (",TIME_MIN if TIME_MIN != -1 else times[:-12].min() ,",",TIME_MAX if TIME_MAX != -1 else times[:-12].max(),")" )
print("The range for the depth is: (",DEPTH_MIN if DEPTH_MIN != -1 else depth[:-12].min() ,",",DEPTH_MAX if DEPTH_MAX != -1 else depth[:-12].max(),")" )
print("Please note this only affects the linear subplot")


vmin = LUT_MIN if MANUAL_LUT_TEMP_MIN else np.min(data)
vmax = LUT_MAX if MANUAL_LUT_TEMP_MAX else np.max(data)

# ---- SPLIT INTO REGULAR + EXPONENTIAL REGIONS ----
N_EXP = 12
N_REG = n_depth - N_EXP

depth_reg = depth[:N_REG]
depth_exp = depth[N_REG:]
data_reg = data[:N_REG, :]    # shape (N_REG, n_time)
data_exp = data[N_REG:, :]    # shape (N_EXP, n_time)

# ---- CELL EDGE HELPERS ----
def cell_edges_linear(centers):
    """Arithmetic midpoints → linear-scale cell edges."""
    c = np.asarray(centers, dtype=float)
    mids  = 0.5 * (c[:-1] + c[1:])
    left  = c[0]  - (c[1]  - c[0])  / 2
    right = c[-1] + (c[-1] - c[-2]) / 2
    return np.concatenate([[left], mids, [right]])

def cell_edges_log(centers):
    """Geometric midpoints → log-scale cell edges (visually even on log axis)."""
    c = np.asarray(centers, dtype=float)
    mids  = np.sqrt(c[:-1] * c[1:])
    left  = c[0]  ** 2 / c[1]
    right = c[-1] ** 2 / c[-2]
    return np.concatenate([[left], mids, [right]])

depth_reg_edges = cell_edges_linear(depth_reg)
depth_exp_edges = cell_edges_log(depth_exp)

# ---- FIGURE ----
fig, (ax1, ax2) = plt.subplots(
    1, 2,
    figsize=(15, 7),
    gridspec_kw={"width_ratios": [3, 1]},
    constrained_layout=True,
)

# -- Left: regular region, linear depth axis --
im1 = ax1.pcolormesh(
    times, depth_reg_edges, data_reg,
    cmap=LUT_CMAP, vmin=vmin, vmax=vmax,
    shading="flat",
)
ax1.invert_yaxis()
ax1.set_xlabel("Time (ns)")
ax1.set_ylabel("Depth (nm)")
ax1.set_title("Conductivity Evolution")

# Dashed boundary line at the bottom of the regular region
ax1.axhline(depth_reg[-1], color="black", lw=1.2, ls="--", alpha=0.6)

# -- Right: exponential region, log depth axis --
im2 = ax2.pcolormesh(
    times, depth_exp_edges, data_exp,
    cmap=LUT_CMAP, vmin=vmin, vmax=vmax,
    shading="flat",
)
ax2.set_yscale("log")
ax2.invert_yaxis()
ax2.set_xlabel("Time (ns)")
ax2.set_title("Extrapolated\n(blurry nodes)")

# Tick labels on the right to avoid overlap with ax1
ax2.yaxis.set_label_position("right")
ax2.yaxis.tick_right()
ax2.set_ylabel("Depth (nm)", labelpad=10)

# Visual cue for lower-fidelity region
ax2.set_facecolor("#f0f0f0")
ax2.text(
    0.5, 0.97, "×2ⁿ spacing",
    transform=ax2.transAxes,
    ha="center", va="top",
    fontsize=11, color="gray", style="italic",
)

ax1.set_xlim(time_min_plot, time_max_plot)
ax1.set_ylim(depth_max_plot, depth_min_plot) 


# ---- SHARED COLORBAR (anchored to both axes) ----
cbar = fig.colorbar(im1, ax=[ax1, ax2], shrink=0.85, pad=0.12)
cbar.set_label("Thermal Conductivity")

plt.show()