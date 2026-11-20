"""
Static map preview (preview.png) of the Heatwave Box Atlas for a GitHub README: draws the Robinson map, boxes and
labels from the built ../github_pages/heatwave-box-atlas/index.html.

Usage:
    python make_preview.py
"""
import json
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch

SITE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "github_pages", "heatwave-box-atlas")
s = open(os.path.join(SITE, "index.html")).read()


def mp(d):
    v, c = [], []
    for ring in d.split("Z"):
        if not ring.strip():
            continue
        pts = [tuple(map(float, p.split())) for p in ring.strip().lstrip("M").split(" L")]
        v += pts + [pts[0]]
        c += [Path.MOVETO] + [Path.LINETO] * (len(pts) - 1) + [Path.CLOSEPOLY]
    return Path(v, c)


fig, ax = plt.subplots(figsize=(14, 3.2))
fig.patch.set_facecolor("#f3f5f7")
ax.add_patch(PathPatch(mp(re.search(r'class="ocean" d="([^"]+)"', s).group(1)), fc="#e6eef5", ec="#8b98a5", lw=.6))
for g in re.findall(r'<path d="([^"]+)"/>', re.search(r'id="grid">(.*?)</g>', s).group(1)):
    xy = [tuple(map(float, p.split())) for p in g.lstrip("M").split(" L")]
    ax.plot(*zip(*xy), color="#8b98a5", lw=.4, ls=(0, (2, 2)), alpha=.6)
ax.add_patch(PathPatch(mp(re.search(r'class="land" d="([^"]+)"', s).group(1)), fc="#d5dbe0", ec="#8b98a5", lw=.3))
cols = {("S", "W"): "#fbd38d", ("S", "E"): "#e3901b", ("C", "W"): "#a8dcb4", ("C", "E"): "#2f9a52",
        ("N", "W"): "#b3cdec", ("N", "E"): "#2f68ad"}
for d in json.loads(re.search(r'const DATA = (\[.*?\]);', s, re.S).group(1)):
    ax.add_patch(PathPatch(mp(d["path"]), fc=cols[(d["band"], d["side"])], ec="white", lw=.4))
    ax.text(d["lx"], d["ly"], d["code"], fontsize=8, ha="center", va="center", family="monospace", color="#17212b",
            bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=.7))
vb = list(map(float, re.search(r'viewBox="([^"]+)" role="group"', s).group(1).split()))
ax.set_xlim(vb[0], vb[0] + vb[2]); ax.set_ylim(vb[1] + vb[3] - 4, vb[1]); ax.set_aspect("equal"); ax.axis("off")
fig.savefig(os.path.join(SITE, "preview.png"), dpi=150, bbox_inches="tight", pad_inches=0.05,
            facecolor=fig.get_facecolor())
