import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap, BoundaryNorm

# ---- STYLE ----
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
time_list = []

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

n_depth = len(time_steps[0])  # number of spatial points

depth_list = np.arange(n_depth) * DX * 1e7  # nm

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

# ---- FIGURE (KEY FIX FOR CUT-OFF ISSUES) ----
fig, ax = plt.subplots(figsize=(10, 6), constrained_layout=True)

im = ax.imshow(
    indexed,
    cmap=cmap,
    norm=norm,
    aspect='auto',
    origin='upper',
    extent=[times[0], times[-1], depth[-1], depth[0]]
)

# ---- LABELS ----
ax.set_xlabel("Time (ns)")
ax.set_ylabel("Depth (nm)")
ax.set_title("State Evolution")

# ---- LEGEND (FIXED + NEVER CUT OFF) ----
patches = [
    mpatches.Patch(color=state_info[v][1], label=state_info[v][0])
    for v in values
]

legend = ax.legend(
    handles=patches,
    loc='center left',
    bbox_to_anchor=(1.02, 0.5),
    frameon=True
)

# ---- EXTRA SAFETY: prevent clipping on tight layouts ----
plt.subplots_adjust(right=0.80)

plt.show()