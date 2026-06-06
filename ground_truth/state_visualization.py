import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap, BoundaryNorm

# ---- STYLE ----
# Set to -1 if you want to use the usual range (from minimum value found to maximum value found)
DEPTH_MIN = -1
DEPTH_MAX = -1
TIME_MIN = -1
TIME_MAX = -1

plt.rcParams.update({
    "font.size": 16,
    "axes.titlesize": 22,
    "axes.labelsize": 22,
    "xtick.labelsize": 22,
    "ytick.labelsize": 22
})

state_info = {
    1: ("CRYSTAL", "blue"),
    2: ("LARGE POLY", "orange"),
    3: ("FINE POLY", "yellow"),
    4: ("AMORPHOUS", "purple"),
    5: ("MUSHY", "brown"),
    8: ("LIQUID", "red"),
    9: ("SUPERCOOLED", "cyan"),
}

#------------------------------PRINTING INFORMATION---------------------------
print("The following colors have been assigned: ")
for name,color in state_info.values():
    print(name, "has been assigned", color)

print("The xtick label size is: ", plt.rcParams.get("xtick.labelsize"))
print("The ytick label size is: ", plt.rcParams.get("ytick.labelsize"))

letter_to_state = {
    'C': 1,
    'P': 2,
    'F': 3,
    'A': 4,
    'M': 5,
    'L': 8,
    'S': 9,
}

# ---- PARSING ----
time_steps = []
time_list = [0]

with open("state.dat") as f:
    for line in f:
        line = line.strip()
        if not line:
            continue

        parts = line.split()

        if parts[0] == "TIME":
            DX = float(parts[3][1:])

        elif all(len(p) == 1 and p.isalpha() for p in parts):
            time_steps.append(parts)

        elif len(parts) > 1 and parts[1].isalpha():
            time_steps.append(parts[1:])
            time_list.append(float(parts[0])*1e9)
time_steps.insert(0,time_steps[0])

depth_list = []  # nm
with open("depth_vector.dat") as f:
    for line in f:
        if not line:
            continue
        depth_list.append(float(line.strip()))
n_depth = len(depth_list)

# ---- PROCESS DATA ----
data = np.array(time_steps)
numeric = np.vectorize(letter_to_state.get)(data).astype(float)
numeric = numeric.T

values = sorted(state_info.keys())
value_to_index = {v: i for i, v in enumerate(values)}
indexed = np.vectorize(lambda x: value_to_index.get(x, -1))(numeric)

colors = [state_info[v][1] for v in values]
cmap = ListedColormap(colors)

norm = BoundaryNorm(np.arange(-0.5, len(values) + 0.5, 1), cmap.N)


times = np.array(time_list, dtype=float)
depth = np.array(depth_list, dtype=float)

time_min_plot = TIME_MIN if TIME_MIN != -1 else times[:-12].min()
time_max_plot = TIME_MAX if TIME_MAX != -1 else times[:-12].max()
depth_min_plot = DEPTH_MIN if DEPTH_MIN != -1 else depth[:-12].min()
depth_max_plot = DEPTH_MAX if DEPTH_MAX != -1 else depth[:-12].max()

print("The range for the time is: (",TIME_MIN if TIME_MIN != -1 else times[:-12].min() ,",",TIME_MAX if TIME_MAX != -1 else times[:-12].max(),")" )
print("The range for the depth is: (",DEPTH_MIN if DEPTH_MIN != -1 else depth[:-12].min() ,",",DEPTH_MAX if DEPTH_MAX != -1 else depth[:-12].max(),")" )
print("Please note this only affects the linear subplot")

# ---- FIGURE (KEY FIX FOR CUT-OFF ISSUES) ----
# ---- SPLIT DATA ----
N_EXP = 12
N_REG = n_depth - N_EXP  # regular linear nodes

depth_reg = depth[:N_REG]
depth_exp = depth[N_REG:]
data_reg = indexed[:N_REG, : -2]   # shape (N_REG, n_time)
data_exp = indexed[N_REG:, :-2]   # shape (N_EXP, n_time)
# ---- CELL EDGE HELPERS ----
def cell_edges_linear(centers):
    """Arithmetic midpoints → linear-scale cell edges."""
    c = np.asarray(centers, dtype=float)
    mids = 0.5 * (c[:-1] + c[1:])
    left  = c[0]  - (c[1]  - c[0])  / 2
    right = c[-1] + (c[-1] - c[-2]) / 2
    return np.concatenate([[left], mids, [right]])

def cell_edges_log(centers):
    """Geometric midpoints → log-scale cell edges (visually even on log axis)."""
    c = np.asarray(centers, dtype=float)
    mids  = np.sqrt(c[:-1] * c[1:])          # geometric mean between neighbours
    left  = c[0]  ** 2 / c[1]                # extrapolate left  in log space
    right = c[-1] ** 2 / c[-2]               # extrapolate right in log space
    return np.concatenate([[left], mids, [right]])

depth_reg_edges = cell_edges_linear(depth_reg)
depth_exp_edges = cell_edges_log(depth_exp)

# ---- FIGURE ----
fig, (ax1, ax2) = plt.subplots(
    1, 2,
    figsize=(15, 7),
    gridspec_kw={"width_ratios": [3, 1]},   # exp region gets less horizontal space
    constrained_layout=True,
)

# -- Left: regular region, linear depth axis --
im1 = ax1.pcolormesh(times, depth_reg_edges, data_reg, cmap=cmap, norm=norm)
ax1.invert_yaxis()                           # surface (0 nm) at top
ax1.set_xlabel("Time (ns)")
ax1.set_ylabel("Depth (nm)")
ax1.set_title("State Evolution Ground truth")

# Mark the boundary with a dashed line at the bottom
ax1.axhline(depth_reg[-1], color="black", lw=1.2, ls="--", alpha=0.6)

# -- Right: extrapolated exponential region, log depth axis --
print(times.shape,depth_exp_edges.shape,data_exp.shape)

im2 = ax2.pcolormesh(times, depth_exp_edges, data_exp, cmap=cmap, norm=norm)
ax2.set_yscale("log")
ax2.invert_yaxis()
ax2.set_xlabel("Time (ns)")
ax2.set_title("Extrapolated\n(blurry nodes)")

# Tick labels on the right side so they don't overlap ax1
ax2.yaxis.set_label_position("right")
ax2.yaxis.tick_right()
ax2.set_ylabel("Depth (nm)", labelpad=10)

# Shade the entire right panel to signal lower physical fidelity
ax2.set_facecolor("#f0f0f0")
ax2.text(
    0.5, 0.97, "×2ⁿ spacing",
    transform=ax2.transAxes,
    ha="center", va="top",
    fontsize=11, color="gray", style="italic",
)

# ---- LEGEND (attached to ax1) ----
patches = [
    mpatches.Patch(color=state_info[v][1], label=state_info[v][0])
    for v in values
]
ax1.legend(
    handles=patches,
    loc="upper left",
    bbox_to_anchor=(1.01, 1.0),
    frameon=True,
    fontsize=13,
)

ax1.set_xlim(time_min_plot, time_max_plot)
ax1.set_ylim(depth_max_plot, depth_min_plot)  # inverted
plt.show()