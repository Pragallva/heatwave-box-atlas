"""
Probability density of daily LWA (averaged over 5-10 km), JJA, MPI-GE members 1-5: all years (1850-2100) pooled,
plus historical (1850-2014) and SSP3-7.0 (2015-2100) separately. Area-weighted (cos lat). Vertical lines: 65 m/s
(blocking threshold A_c) and the 88th percentile of the all-years distribution.

Default: all grid points 38-85N (cache from make_LWA_JJA_38-85N_5-10km_hist.py) -> LWA_JJA_38-85N_5-10km_pdf.png
--regions: each of the 18 continent lat/lon boxes (cache from make_LWA_JJA_5-10km_hist_regions.py)
           -> LWA_JJA_5-10km_pdf_<box>.png and LWA_JJA_5-10km_pdf_regions_stats.json (used by the dashboard);
           all boxes share one x-range

Other latitude band / layer (all longitudes; cache from make_LWA_JJA_38-85N_5-10km_hist.py with the same numbers)
           -> LWA_JJA_<LAT0>-<LAT1>N_<Z0>-<Z1>km_pdf.png

Usage:
    python plot_LWA_JJA_38-85N_5-10km_pdf.py
    python plot_LWA_JJA_38-85N_5-10km_pdf.py --regions
    python plot_LWA_JJA_38-85N_5-10km_pdf.py --regions --nofoot   # dashboard figures without the footnote
    python plot_LWA_JJA_38-85N_5-10km_pdf.py 35 80 5 12        # 35-80N, 5-12 km
"""
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
AC = 65.0
NOFOOT = "--nofoot" in sys.argv   # no explanatory text below the plot; output name gets _nofoot
PQ = 88          # percentile of the all-years distribution marked by a line
BANDS = {"NA": [(38, 48), (48, 58), (58, 85)], "EU": [(38, 48), (48, 58), (58, 72)], "AS": [(38, 48), (48, 58), (58, 80)]}
SPLIT = {"NA": ("< 100°W", "≥ 100°W"), "EU": ("< 15°E", "≥ 15°E"), "AS": ("0–100°E", "≥ 100°E")}
CONT = {"NA": "North America", "EU": "Europe", "AS": "Asia"}

plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Nimbus Sans", "DejaVu Sans"],
                     "font.size": 15, "axes.spines.top": False, "axes.spines.right": False})


def box_name(code):
    cc, side, b = code[:2], code[3], code[4]
    lo, hi = BANDS[cc]["SCN".index(b)]
    return f"{CONT[cc]} {'West' if side == 'W' else 'East'} ({SPLIT[cc][side == 'E']}), {lo}–{hi}°N"


def percentile(h, edges, q):
    """q-th percentile from a histogram, linear interpolation within the bin."""
    cdf = np.cumsum(h) / h.sum()
    return np.interp(q / 100, np.concatenate([[0], cdf[:-1]]), edges[:-1])


def plot_pdf(edges, W, scen, title, where, out, xmax=160, zr="5–10"):
    centres, width = 0.5 * (edges[:-2] + edges[1:-1]), np.diff(edges[:-1])
    fig, ax = plt.subplots(figsize=(10, 6), constrained_layout=True)
    st = {}
    for sel, key, lab, col, fill in [(np.ones(len(scen), bool), "all", "All years", "0.25", True),
                                     (scen == "historical", "hist", "Historical", "#1E3A8A", False),
                                     (scen == "ssp370", "ssp", "SSP3-7.0", "#B45309", False)]:
        h = W[sel].sum(axis=0)
        tot = h.sum()
        pdf = h[:-1] / (tot * width)                              # last bin = overflow > 300 m/s
        st[key] = round(100 * h[edges[:-1] >= AC].sum() / tot, 2)
        st[key + "_mean"] = round((h[:-1] * centres).sum() / h[:-1].sum(), 2)
        st[key + "_median"] = round(percentile(h, edges, 50), 2)
        if fill:
            ax.fill_between(centres, pdf, color=col, alpha=0.18, lw=0)
            ax.plot(centres, pdf, color=col, lw=2.0, label=lab)
            pq = percentile(h, edges, PQ)
            st["p88"], st["pct_of_65"] = round(pq, 2), round(100 - st["all"], 2)
        else:
            ax.plot(centres, pdf, color=col, lw=1.8, ls="--", label=lab)
    ax.set_ylim(bottom=0)
    ytop = ax.get_ylim()[1]
    ax.axvline(AC, color="#C2189A", lw=2.2, zorder=5)
    left = pq < AC     # 65 label on the side away from the percentile line
    ax.text(AC + (2 if left else -2), ytop * 0.97, f"{AC:.0f} m s$^{{-1}}$", color="#C2189A",
            fontweight="bold", va="top", ha="left" if left else "right")
    ax.axvline(pq, color="#0F766E", lw=2.2, ls=(0, (5, 2)), zorder=5)
    ax.text(max(pq, AC) + 2, ytop * 0.62, f"{PQ}th percentile,\nall years: {pq:.1f} m s$^{{-1}}$",
            color="#0F766E", fontweight="bold", va="top", ha="left")   # in the empty space right of both lines
    ax.set_xlim(0, xmax)
    ax.set_xlabel(f"LWA, {zr} km mean (m s$^{{-1}}$)")
    ax.set_ylabel("Probability density (per m s$^{-1}$)")
    ax.set_title(title, loc="left", fontweight="bold")
    ax.legend(frameon=False, loc="upper right", fontsize=13)
    ax.grid(axis="y", color="0.9", lw=0.7); ax.set_axisbelow(True)
    if not NOFOOT:
        fig.text(0.01, -0.02, f"Daily LWA averaged over {zr} km height at {where}, June–August, historical + SSP3-7.0, "
                 "members 1–5; area-weighted (cos latitude); 0.5 m s$^{-1}$ bins.",
                 ha="left", va="top", fontsize=11, color="0.3", wrap=True)
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"{os.path.basename(out)}: {PQ}th pct {pq:.1f}; 65 m/s = {st['pct_of_65']:.1f}th pct; "
          f"% >= 65: all {st['all']:.1f}, hist {st['hist']:.1f}, SSP {st['ssp']:.1f}")
    return st


if "--regions" in sys.argv:
    c = np.load(os.path.join(HERE, "LWA_JJA_5-10km_hist_regions_cache.npz"))
    edges, W, scen, codes, npts = c["edges"], c["wcounts"], c["scenario"], list(c["codes"]), c["npts"]
    # common x-range: 99.5th percentile of the widest box, rounded up to 20 m/s
    xmax = 20 * np.ceil(max(percentile(W[:, k].sum(axis=0), edges, 99.5) for k in range(len(codes))) / 20)
    stats = {}
    for k, code in enumerate(codes):
        stats[code] = plot_pdf(edges, W[:, k], scen, f"{code}: {box_name(code)} — daily LWA, JJA",
                               f"the {npts[k]} grid points (1.875°) whose centre lies in {code}",
                               os.path.join(HERE, f"LWA_JJA_5-10km_pdf_{code}{'_nofoot' if NOFOOT else ''}.png"), xmax)
        stats[code]["npts"] = int(npts[k])
    json.dump(stats, open(os.path.join(HERE, "LWA_JJA_5-10km_pdf_regions_stats.json"), "w"), indent=1, default=lambda v: round(float(v), 2))
else:
    nums = [int(a) for a in sys.argv[1:] if a.isdigit()]
    lat0, lat1, z0, z1 = nums if len(nums) == 4 else (38, 85, 5, 10)
    tag = f"{lat0}-{lat1}N_{z0}-{z1}km"
    c = np.load(os.path.join(HERE, f"LWA_JJA_{tag}_hist_cache.npz"))
    plot_pdf(c["edges"], c["wcounts"], c["scenario"], f"Daily LWA at grid points {lat0}–{lat1}°N, JJA, MPI-GE members 1–5",
             f"every grid point {lat0}–{lat1}°N (all longitudes)", os.path.join(HERE, f"LWA_JJA_{tag}_pdf.png"),
             zr=f"{z0}–{z1}")
