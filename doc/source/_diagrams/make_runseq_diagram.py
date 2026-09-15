"""Generate a nested time-loop diagram for the CMEPS run sequence.

Run this script from anywhere; it writes ``run_sequence_loops.png`` next to the
documentation source (the parent directory of this ``_diagrams`` folder).
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

fig, ax = plt.subplots(figsize=(11, 7))
ax.set_xlim(0, 110)
ax.set_ylim(0, 72)
ax.axis("off")

# (x, y, w, h, edge, fill, header, body)
boxes = [
    (2,  2, 106, 68, "#1f3b73", "#eef3fb",
     "@86400   (1 day)",
     "run GLC  (land-ice)"),
    (7,  6,  96, 54, "#2a6f97", "#e7f0f7",
     "@10800   (3 hours)",
     "run ROF  (river)"),
    (12, 10, 86, 40, "#2d8f6f", "#e8f5ef",
     "@3600   (1 hour)",
     "run OCN  (ocean)"),
    (17, 14, 76, 26, "#b5651d", "#fbf1e6",
     "@1800   (30 min - fast loop)",
     "MED: atm/ocn fluxes, ocean albedo,\naccumulate ocean forcing\nrun ICE, LND, ATM  (each: prep -> run -> post)"),
]

for (x, y, w, h, ec, fc, header, body) in boxes:
    ax.add_patch(FancyBboxPatch((x, y), w, h,
                 boxstyle="round,pad=0,rounding_size=1.2",
                 linewidth=2, edgecolor=ec, facecolor=fc))
    ax.text(x + 2.2, y + h - 2.4, header, fontsize=11, fontweight="bold",
            color=ec, ha="left", va="top")
    ax.text(x + 2.2, y + h - 6.0, body, fontsize=9.5, color="#222",
            ha="left", va="top")

# innermost doubled loop (ocean coupling interval), nested inside the fast loop
xi, yi, wi, hi = 21, 15, 70, 10.0
ax.add_patch(FancyBboxPatch((xi, yi), wi, hi,
             boxstyle="round,pad=0,rounding_size=1.0",
             linewidth=2, edgecolor="#7b2cbf", facecolor="#f3ebfb",
             linestyle="--"))
ax.text(xi + 2, yi + hi - 2.2, "@@3600   (ocean coupling interval)",
        fontsize=10, fontweight="bold", color="#7b2cbf", ha="left", va="top")
ax.text(xi + 2, yi + hi - 5.2,
        "fires once per hour (every 2nd pass of the 30-min loop):\n"
        "average accumulated ocean forcing  ->  send to OCN",
        fontsize=9.5, color="#222", ha="left", va="top")

ax.set_title("CMEPS run sequence: nested coupling-interval loops",
             fontsize=13, fontweight="bold", pad=10)

# Write the PNG next to the docs source (parent of this _diagrams directory),
# regardless of the current working directory.
_here = os.path.dirname(os.path.abspath(__file__))
out = os.path.join(os.path.dirname(_here), "run_sequence_loops.png")
fig.savefig(out, dpi=150, bbox_inches="tight")
print("wrote", out)
