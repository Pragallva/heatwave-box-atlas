"""
Heatwave days per decade and their regime fractions for two regions of the lat/lon-band analysis
(see ../2026-10-04/region_map_latlonbands.png):
  EU-WS : Europe (countries west of 45E), West of 15E, southern band 38-48N
  EU-WC : Europe (countries west of 45E), West of 15E, central band 48-58N
  NA-WS : North America (Greenland excl.), West of 100W, southern band 38-48N
  NA-EN : North America (Greenland excl.), East of 100W, northern band 58-85N
  SCAN  : Scandinavia, mainland Norway + Sweden + Finland + Denmark (as in regime_shift_per_decade_WesternEurope_SWEurope_Scandinavia)
  GL    : Greenland, centroid >= 60N (as in the Arctic/Greenland regime-shift figure)
A heatwave belongs to the region if its centroid lies inside (Natural Earth 1:110m countries).

Heatwave days of an event = its persistence (max number of T > T_P90 days at a grid point of the heatwave
core within +/-10 days of the block peak). Per decade: sum over events, divided by the number of members (5);
error bars = range across the 5 members.
Top: total heatwave days per decade per member. Bottom: fraction of heatwave days in each regime.

Usage:
    python plot_heatwave_days_per_decade_EU-WS_NA-WS.py              # all regions
    python plot_heatwave_days_per_decade_EU-WS_NA-WS.py NA-EN GL     # selected regions
    python plot_heatwave_days_per_decade_EU-WS_NA-WS.py EU-WS --violin   # top panel as violins of the 5 member values
                                                                         # (-> heatwave_days_per_decade_<code>_violin.png)
    python plot_heatwave_days_per_decade_EU-WS_NA-WS.py EU-WS --sum      # top panel = sum over all 5 members
                                                                         # (-> heatwave_days_per_decade_<code>_sum5members.png)
    python plot_heatwave_days_per_decade_EU-WS_NA-WS.py EU-WS --sum --t2m  # panel (c) with 2 m air temperature
                                                                         # (-> ..._sum5members_t2m.png)
    python plot_heatwave_days_per_decade_EU-WS_NA-WS.py --latlon --sum --ts --png   # all 18 lat/lon boxes, PNG only
    python plot_heatwave_days_per_decade_EU-WS_NA-WS.py EU-WS --sum --ts --lwa --png  # + panel (d) LWA, all heatwaves
    python plot_heatwave_days_per_decade_EU-WS_NA-WS.py --latlon --sum --ts --lwa --png --nofoot  # dashboard figures, no footnote
    python plot_heatwave_days_per_decade_EU-WS_NA-WS.py EU-WS --sum --ts   # panel (c) with skin temperature
                                                                         # (-> ..._sum5members_ts.png)
"""
import os
import textwrap
import sys

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cartopy.io.shapereader as shpreader
from shapely.geometry import Point, box
from shapely.ops import unary_union
from shapely.prepared import prep

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = ("/projects/NS11084K/MPI-LE/daily/post_processing/blocking_info_JJA/mid_to_upper_level/"
       "JJA_NH_Ac=65_new_heatwave_centroid/pandas_track_data_summer_5-12km_summer_heatwaves_more_area_"
       "smooth_baseline_with_land_ocean_mask_and_dry_soil_Bo_ratio_0.5_threshold.csv")
CSV = os.environ.get("HW_TRACK_CSV", CSV)   # override, e.g. a local copy (.csv or .csv.gz)
REGIMES = [("energy_limited_regime", "Energy-limited", "#2563EB"),
           ("transitional_regime", "Transitional", "#16A34A"),
           ("moisture_limited_regime", "Moisture-limited", "#DC2626"),
           ("unknown", "Atm-to-sfc heating", "#9CA3AF")]
DECADES = np.arange(1850, 2100, 10)
# per-event heatwave temperature: ta mean over sigma = 1.0-0.6, +/-10 deg box, days -3..+3 (index = CSV row)
T2M = "--t2m" in sys.argv     # 2 m air temperature (tas) instead of ta(sigma 1.0-0.6) in panel (c)
TS = "--ts" in sys.argv       # skin temperature derived from rlus/rlds (make_heatwave_event_skinT_cache.py)
if TS:
    TCACHE, TCOL = os.path.join(HERE, "heatwave_event_ts_days-3to3_cache.csv"), "ts"
else:
    TCACHE, TCOL = os.path.join(HERE, "..", "2026-10-04", "heatwave_event_tas_days-3to3_cache.csv" if T2M else
                                "heatwave_event_ta_sigma1.0-0.6_unweighted_days-3to3_cache.csv"), "tas"
WIN, MIN_EV = 10, 1
X0 = 2015          # hinge of the joined trend: centre of the 2010s, the last historical (transition) decade


def tight_ylim(ax, lo, hi, frac=0.08):
    """y-limits around the running-mean envelopes (lo/hi = min/max of mean -/+ s.e.) with a small margin."""
    pad = frac * (hi - lo) if hi > lo else 1.0
    ax.set_ylim(lo - pad, hi + pad)


def joined_trend(x, y):
    """Continuous piecewise-linear least-squares fit with a hinge at X0:
    y = a + b_h (x - X0) for x <= X0 and a + b_s (x - X0) for x > X0, so both lines share the value a at X0.
    Returns a, (b_h, p_h), (b_s, p_s) with two-sided t-test p-values of the slopes (dof = n - 3)."""
    d = x - X0
    A = np.column_stack([np.ones_like(d), np.minimum(d, 0), np.maximum(d, 0)])
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    dof = len(y) - 3
    cov = np.sum((y - A @ coef) ** 2) / dof * np.linalg.inv(A.T @ A)
    pv = [2 * stats.t.sf(abs(coef[i] / np.sqrt(cov[i, i])), dof) for i in (1, 2)]
    return coef[0], (coef[1], pv[0]), (coef[2], pv[1])
YEARS = np.arange(1850, 2101)


def window_mean(y, v):
    """Mean and standard error of all values within +/-WIN years (>= MIN_EV values)."""
    m, se = np.full(YEARS.size, np.nan), np.full(YEARS.size, np.nan)
    for i, yr in enumerate(YEARS):
        sel = np.abs(y - yr) <= WIN
        if sel.sum() >= MIN_EV:
            m[i] = v[sel].mean(); se[i] = v[sel].std(ddof=1) / np.sqrt(sel.sum())
    return m, se

recs = list(shpreader.Reader(shpreader.natural_earth("110m", "cultural", "admin_0_countries")).records())
europe = unary_union([r.geometry for r in recs if r.attributes["CONTINENT"] == "Europe"]).intersection(box(-25, 34, 45, 72))
north_am = unary_union([r.geometry for r in recs if r.attributes["CONTINENT"] == "North America"
                        and r.attributes["NAME"] != "Greenland"])
REGIONS = [("EU-WS", "Europe West (< 15°E), 38–48°N", europe.intersection(box(-25, 38, 15, 48))),
           ("EU-WC", "Europe West (< 15°E), 48–58°N", europe.intersection(box(-25, 48, 15, 58))),
           ("NA-WS", "North America West (< 100°W), 38–48°N", north_am.intersection(box(-180, 38, -100, 48))),
           ("NA-EN", "North America East (≥ 100°W), 58–85°N", north_am.intersection(box(-100, 58, -20, 85))),
           ("GL", "Greenland (≥ 60°N)", [r.geometry for r in recs if r.attributes["NAME"] == "Greenland"][0]
            .intersection(box(-80, 60, 0, 90))),
           ("SCAN", "Scandinavia (Norway, Sweden, Finland, Denmark)",
            unary_union([r.geometry for r in recs if r.attributes["NAME"] in ("Norway", "Sweden", "Finland", "Denmark")])
            .intersection(box(3, 54, 33, 72)))]
# all 18 lat/lon-band boxes of region_map_latlonbands.png (same definitions as the regime-shift figures):
# <continent>-<W|E><S|C|N>; Asia West = 0-100E (Asia polygon incl. Russia east of 45E), Asia East = east of 100E
asia = unary_union([r.geometry for r in recs if r.attributes["CONTINENT"] == "Asia"]
                   + [[r.geometry for r in recs if r.attributes["NAME"] == "Russia"][0].intersection(box(45, 30, 180, 85))])
for cc, cname, geom, bands, sname, west, east in [
        ("NA", "North America", north_am, [(38, 48), (48, 58), (58, 85)], "100°W", [(-180, -100)], [(-100, -20)]),
        ("EU", "Europe", europe, [(38, 48), (48, 58), (58, 72)], "15°E", [(-25, 15)], [(15, 45)]),
        ("AS", "Asia", asia, [(38, 48), (48, 58), (58, 80)], "100°E", [(0, 100)], [(100, 180), (-180, 0)])]:
    for side, lons, rel in [("W", west, "<"), ("E", east, "≥")]:
        for (lo, hi), b, bname in zip(bands, "SCN", ("southern", "central", "northern")):
            code_ = f"{cc}-{side}{b}"
            if code_ in [r[0] for r in REGIONS]:
                continue
            poly = unary_union([geom.intersection(box(l0, lo, l1, hi)) for l0, l1 in lons])
            REGIONS.append((code_, f"{cname} {'West' if side == 'W' else 'East'} ({rel} {sname}), {lo}–{hi}°N", poly))
LATLON = [f"{c}-{s}{b}" for c in ("NA", "EU", "AS") for s in "WE" for b in "SCN"]
if "--latlon" in sys.argv:                   # all 18 continent lat/lon boxes
    REGIONS = sorted([r for r in REGIONS if r[0] in LATLON], key=lambda r: LATLON.index(r[0]))
PNG_ONLY = "--png" in sys.argv
LWA = "--lwa" in sys.argv     # add panel (d): LWA of all heatwaves (make_heatwave_event_LWA_cache.py)
LWA_ZTOP = 10   # km; LWA averaged over 5-LWA_ZTOP km (make_heatwave_event_LWA_cache.py LWA_ZTOP)
LCACHE = os.path.join(HERE, f"heatwave_event_LWA_5-{LWA_ZTOP}km_days-3to3_cache.csv")
VIOLIN = "--violin" in sys.argv
SUM = "--sum" in sys.argv
NOFOOT = "--nofoot" in sys.argv   # no explanatory text below the panels; output name gets _nofoot
ARGS = [a for a in sys.argv[1:] if not a.startswith("--")]
if ARGS:
    REGIONS = [r for r in REGIONS if r[0] in ARGS]

df = pd.read_csv(CSV, usecols=["member", "peak_block_date", "HeatWave", "energy_state", "persistence",
                               "HW_x_loc", "HW_y_loc"])
df = df[df["HeatWave"] == True].copy()
df["decade"] = pd.to_datetime(df["peak_block_date"]).dt.year // 10 * 10
df["lon"] = ((df["HW_x_loc"] + 180) % 360) - 180
pts = [Point(x, y) for x, y in zip(df["lon"], df["HW_y_loc"])]
members = sorted(df["member"].unique())
tcache = pd.read_csv(TCACHE).dropna(subset=[TCOL]).set_index("index")
df["T"] = tcache[TCOL].reindex(df.index).values
if LWA:
    df["LWA"] = pd.read_csv(LCACHE).set_index("index")["LWA"].reindex(df.index).values
df["year"] = pd.to_datetime(df["peak_block_date"]).dt.year

plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Nimbus Sans", "DejaVu Sans"],
                     "font.size": 17, "axes.spines.top": False, "axes.spines.right": False})
for code, name, poly in REGIONS:
    pg = prep(poly)
    sub = df[[pg.contains(p) for p in pts]]
    per_mem = (sub.groupby(["decade", "member"])["persistence"].sum().unstack(fill_value=0)
               .reindex(index=DECADES, columns=members, fill_value=0))
    total = per_mem.mean(axis=1)
    days = (sub.groupby(["decade", "energy_state"])["persistence"].sum().unstack(fill_value=0)
            .reindex(index=DECADES, columns=[k for k, _, _ in REGIMES], fill_value=0))
    frac = days.div(days.sum(axis=1).replace(0, np.nan), axis=0) * 100

    fig, axs = plt.subplots(4 if LWA else 3, 1, figsize=(14, 19.5 if LWA else 15), sharex=True, constrained_layout=True)
    a0, a1, a2 = axs[:3]
    alast = axs[-1]
    if SUM:
        a0.bar(DECADES + 5, per_mem.sum(axis=1), width=8.6, color="0.3", alpha=0.85, edgecolor="white", linewidth=0.4)
    elif VIOLIN:
        ok = per_mem.std(axis=1).values > 0          # a KDE needs spread; constant decades get only dots
        vp = a0.violinplot([per_mem.loc[d].values for d in DECADES[ok]], positions=DECADES[ok] + 5, widths=8.6,
                           showextrema=False)
        for b in vp["bodies"]:
            b.set_facecolor("0.55"); b.set_edgecolor("0.25"); b.set_linewidth(0.6); b.set_alpha(0.6)
        rng = np.random.default_rng(0)
        for d in DECADES:
            v = per_mem.loc[d].values
            a0.scatter(d + 5 + rng.uniform(-1.2, 1.2, v.size), v, s=12, color="0.1", alpha=0.7,
                       edgecolors="none", zorder=3)
        a0.plot(DECADES + 5, total, color="k", lw=0.8, alpha=0.5, zorder=3)
        a0.scatter(DECADES + 5, total, marker="D", s=30, color="white", edgecolors="k", linewidths=1.0, zorder=4,
                   label="Member mean")
        a0.legend(loc="upper left", frameon=False, fontsize=15)
    else:
        a0.bar(DECADES + 5, total, width=8.6, color="0.3", alpha=0.85, edgecolor="white", linewidth=0.4)
        a0.errorbar(DECADES + 5, total, yerr=[total - per_mem.min(axis=1), per_mem.max(axis=1) - total],
                    fmt="none", ecolor="0.15", elinewidth=0.9, capsize=2.5)
    yval = per_mem.sum(axis=1) if SUM else total
    # joined linear trends of the decadal values (x = decade centre): historical and SSP3-7.0 segments meet at X0
    trend_txt = []
    a, *fits = joined_trend(yval.index.values + 5.0, yval.values.astype(float))
    for (b, p), xx, col, lab in zip(fits, [(1850, X0), (X0, 2100)], ["#1E3A8A", "#B45309"], ["Historical", "SSP3-7.0"]):
        xx = np.array(xx)
        a0.plot(xx, a + b * (xx - X0), color=col, lw=2.6, zorder=5,
                label=f"{lab} trend: {10 * b:+.1f} days per decade (p = {p:.2f})")
        trend_txt.append(f"{lab} {10 * b:+.1f}/dec (p={p:.2f})")
    top = (per_mem.sum(axis=1).max() if SUM else per_mem.max(axis=1).max()) * 1.06
    a0.set_ylim(0, top)
    # legend and period labels sit above the axes so the bars can fill the panel
    a0.legend(loc="lower left", bbox_to_anchor=(0, 1.0), frameon=False, fontsize=15, handlelength=1.8,
              borderaxespad=0.2)
    a0.set_ylabel("Heatwave days per decade\n(sum of 5 members)" if SUM else "Heatwave days per decade\nper member")
    a0.set_title(f"(a) {code}: {name} — heatwave days per decade\n(n = {len(sub)} heatwaves, "
                 f"{sub['persistence'].sum():,.0f} heatwave days)", loc="left", fontsize=18, fontweight="bold", pad=62)

    bottom = np.zeros(len(DECADES))
    for key, label, color in REGIMES:
        v = frac[key].fillna(0).values
        a1.bar(DECADES + 5, v, width=8.6, bottom=bottom, color=color, alpha=0.9, edgecolor="white",
               linewidth=0.4, label=label)
        bottom += v
    a1.set_ylim(0, 100)
    a1.set_ylabel("Fraction of heatwave days (%)")
    a1.set_title(f"(b) {code} — fraction of heatwave days by regime", loc="left", fontsize=18, fontweight="bold")
    # (c) heatwave temperature by regime (running means only; y-range fitted to them)
    env = []
    for key, label, color in REGIMES:
        g = sub[(sub["energy_state"] == key) & sub["T"].notna()]
        if len(g) == 0:
            continue
        m, se = window_mean(g["year"].values, g["T"].values)
        se = np.nan_to_num(se)
        a2.fill_between(YEARS, m - se, m + se, color=color, alpha=0.18, lw=0)
        a2.plot(YEARS, m, color=color, lw=2.2, zorder=4)
        env += [np.nanmin(m - se), np.nanmax(m + se)]
    if env:
        tight_ylim(a2, min(env), max(env))
    a2.set_ylabel("Heatwave surface\ntemperature (°C)" if TS else "Heatwave 2 m air\ntemperature (°C)" if T2M else "Heatwave temperature (°C)\nσ = 1.0–0.6 mean")
    a2.set_title(f"(c) {code} — heatwave {'surface ' if TS else '2 m air ' if T2M else ''}temperature by regime", loc="left", fontsize=18, fontweight="bold")
    if LWA:   # (d) LWA of all heatwaves in the box (no regime split), same running mean
        a3 = axs[3]
        g = sub[sub["LWA"].notna()]
        LWA_COLOR = "#C2189A"   # magenta
        m, se = window_mean(g["year"].values, g["LWA"].values)
        a3.fill_between(YEARS, m - se, m + se, color=LWA_COLOR, alpha=0.18, lw=0)
        a3.plot(YEARS, m, color=LWA_COLOR, lw=2.4, zorder=4)
        se = np.nan_to_num(se)
        tight_ylim(a3, np.nanmin(m - se), np.nanmax(m + se))
        print(f"   LWA running mean {np.nanmin(m):.1f}-{np.nanmax(m):.1f}")
        a3.set_ylabel(f"LWA, 5–{LWA_ZTOP} km mean\n(m s$^{{-1}}$)")
        a3.set_title(f"(d) {code} — local wave activity, all heatwaves (n = {len(g)})", loc="left", fontsize=18,
                     fontweight="bold")
    alast.set_xlabel("Year")
    # regime legend between panels (b) and (c)
    a1.legend(loc="upper center", bbox_to_anchor=(0.5, -0.025), ncol=4, frameon=False, fontsize=16,
              handlelength=1.6, borderaxespad=0)
    for ax in axs:
        ax.axvline(2015, color="k", lw=1.8, ls="--", zorder=6)
        ax.grid(axis="y", color="0.9", lw=0.7); ax.set_axisbelow(True)
    a1.set_xlim(1850, 2100)
    for ax in [a0, *axs[2:]]:
        ax.axvspan(2015, 2100, color="#FEF3E2", lw=0, zorder=0)
    tr = a0.get_xaxis_transform()
    a0.text(2012, 1.015, "← Historical", ha="right", va="bottom", fontsize=17, fontweight="bold", color="#1E3A8A",
            transform=tr)
    a0.text(2018, 1.015, "SSP3-7.0 →", ha="left", va="bottom", fontsize=17, fontweight="bold", color="#B45309",
            transform=tr)
    foot = (
             "Heatwave days = sum of event persistence (days with T > T$_{P90}$ in the heatwave core, ±10 days around the block "
             "peak) of heatwaves whose centroid lies in the region. Dashed line: 2015 (start of SSP3-7.0).\n"
             + ("(a) Total over all 5 members. " if SUM else
                "(a) Violins: the 5 members' decadal totals (dots: members; diamonds: member mean). " if VIOLIN else
                "(a) Mean over 5 members, error bars = member range. ")
             + "Trend lines: joined (continuous piecewise-linear) least-squares fit to the decadal values, with the historical "
             "and SSP3-7.0 segments meeting at the 2010s (transition decade, 2015).\n"
             "(b) 5 members pooled.\n(c) " + ("surface (skin) temperature Ts = [(rlus − (1−ε) rlds)/(εσ)]$^{1/4}$, ε = 0.996, from daily rlus/rlds" if TS else
                                       "2 m air temperature (tas)" if T2M else "ta averaged over model levels σ = 1.0–0.6")
             + ", ±10° box around the heatwave centroid, "
             f"days −3…+3; lines: mean of a regime's heatwaves within ±{WIN} yr (≥ {MIN_EV} event{'s' if MIN_EV > 1 else ''}), shading ±1 s.e.; "
             "y-axis fitted to the running means."
             + (f"\n(d) Local wave activity (LWA) averaged over 5–{LWA_ZTOP} km, the LWA composite box around the event and days −3…+3; "
                f"mean of all heatwaves in the box (all regimes) within ±{WIN} yr, shading ±1 s.e.; y-axis fitted to the running mean." if LWA else ""))
    foot = "\n".join(textwrap.fill(line, 135) for line in foot.split("\n"))
    if not NOFOOT:
        fig.text(0.005, -0.01, foot,
                 ha="left", va="top", fontsize=14.5, color="0.3")

    out = os.path.join(HERE, f"heatwave_days_per_decade_{code}{'_sum5members' if SUM else '_violin' if VIOLIN else ''}{'_ts' if TS else '_t2m' if T2M else ''}{'_lwa' if LWA else ''}{'_nofoot' if NOFOOT else ''}.png")
    fig.savefig(out, dpi=200, bbox_inches="tight")
    if not PNG_ONLY:
        fig.savefig(out.replace(".png", ".pdf"), bbox_inches="tight")
    plt.close(fig)
    h, s = days.loc[1850:2010].sum(), days.loc[2020:2090].sum()
    print(f"{code}: n={len(sub)}  days/decade/member hist {total.loc[1850:2010].mean():.1f}, "
          f"SSP 2020-2090 {total.loc[2020:2090].mean():.1f}")
    for k, lab, _ in REGIMES:
        print(f"   {lab:20s} hist {100 * h[k] / h.sum():5.1f}%   SSP {100 * s[k] / s.sum():5.1f}%")
    print("   trends:", "; ".join(trend_txt))
    print("Saved", out)
