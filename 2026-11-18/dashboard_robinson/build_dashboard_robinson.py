"""
Robinson-projection version of ../../2026-10-07/dashboard/build_dashboard.py (same boxes, statistics and figures;
only the map is reprojected), simplified: the sidebar shows heatwave counts only and the methods text lives in details.html.
Output: the static site ../github_pages/heatwave-box-atlas/ (index.html, details.html, figs/), also published as the
claude.ai artifact. Figures are the footnote-free versions:
    ../../2026-10-07: python plot_heatwave_days_per_decade_EU-WS_NA-WS.py --latlon --sum --ts --lwa --png --nofoot
    ../../2026-10-08: python plot_LWA_JJA_38-85N_5-10km_pdf.py --regions --nofoot

Build the interactive heatwave-box dashboard (index.html + figs/*.png).

Map: Natural Earth 1:110m land, 30-88N, equirectangular; the 18 lat/lon boxes are the land-clipped
region polygons used in plot_heatwave_days_per_decade_EU-WS_NA-WS.py (--latlon). Clicking a box shows
heatwave_days_per_decade_<box>_sum5members_ts_lwa.png, or (LWA distribution view)
../../2026-10-08/LWA_JJA_5-10km_pdf_<box>.png. Side stats: heatwave count, the share of heatwave days per regime,
historical (1850s-2010s) vs SSP3-7.0 (2020s-2090s), and LWA 88th percentile / share >= 65 m/s
(../../2026-10-08/LWA_JJA_5-10km_pdf_regions_stats.json).

Usage:
    python build_dashboard.py
"""
import json
import os
import shutil

import numpy as np
import pandas as pd
import cartopy.io.shapereader as shpreader
import shapely
from pyproj import Transformer
from shapely.geometry import Point, box
from shapely.ops import unary_union
from shapely.prepared import prep

HERE = os.path.dirname(os.path.abspath(__file__))
SITE = os.path.join(HERE, "..", "github_pages", "heatwave-box-atlas")
FIGDIR = os.path.join(HERE, "..", "..", "2026-10-07")
CSV = ("/projects/NS11084K/MPI-LE/daily/post_processing/blocking_info_JJA/mid_to_upper_level/"
       "JJA_NH_Ac=65_new_heatwave_centroid/pandas_track_data_summer_5-12km_summer_heatwaves_more_area_"
       "smooth_baseline_with_land_ocean_mask_and_dry_soil_Bo_ratio_0.5_threshold.csv")
CSV = os.environ.get("HW_TRACK_CSV", CSV)   # override, e.g. a local copy (.csv or .csv.gz)
LWADIR = os.path.join(os.path.dirname(FIGDIR), "2026-10-08")
LWASTATS = json.load(open(os.path.join(LWADIR, "LWA_JJA_5-10km_pdf_regions_stats.json")))
LAT0, LAT1 = 30, 88
REGIMES = ["energy_limited_regime", "transitional_regime", "moisture_limited_regime", "unknown"]

recs = list(shpreader.Reader(shpreader.natural_earth("110m", "cultural", "admin_0_countries")).records())
G = {r.attributes["NAME"]: r.geometry for r in recs}
europe = unary_union([r.geometry for r in recs if r.attributes["CONTINENT"] == "Europe"]).intersection(box(-25, 34, 45, 72))
north_am = unary_union([r.geometry for r in recs if r.attributes["CONTINENT"] == "North America"
                        and r.attributes["NAME"] != "Greenland"])
asia = unary_union([r.geometry for r in recs if r.attributes["CONTINENT"] == "Asia"]
                   + [G["Russia"].intersection(box(45, 30, 180, 85))])
land = unary_union([r.geometry for r in recs]).intersection(box(-180, LAT0, 180, LAT1))

BOXES = []
for cc, cname, geom, bands, sname, west, east in [
        ("NA", "North America", north_am, [(38, 48), (48, 58), (58, 85)], "100°W", [(-180, -100)], [(-100, -20)]),
        ("EU", "Europe", europe, [(38, 48), (48, 58), (58, 72)], "15°E", [(-25, 15)], [(15, 45)]),
        ("AS", "Asia", asia, [(38, 48), (48, 58), (58, 80)], "100°E", [(0, 100)], [(100, 180), (-180, 0)])]:
    for side, lons, rel in [("W", west, "west of"), ("E", east, "east of")]:
        for (lo, hi), b, bname in zip(bands, "SCN", ("Southern", "Central", "Northern")):
            poly = unary_union([geom.intersection(box(l0, lo, l1, hi)) for l0, l1 in lons])
            BOXES.append(dict(code=f"{cc}-{side}{b}", continent=cname, band=b,
                              label=f"{bname} {cname}, {'West' if side == 'W' else 'East'}",
                              where=f"{lo}–{hi}°N, {rel} {sname}", side=side, geom=poly))


ROB = Transformer.from_crs("EPSG:4326", "ESRI:54030", always_xy=True)   # Robinson, central meridian 0
SCALE = 1e5                                                              # map units = 100 km


def proj(lon, lat):
    x, y = ROB.transform(np.asarray(lon, float), np.asarray(lat, float))
    return x / SCALE, -y / SCALE                                         # SVG y points down


def path(geom, nd=2, seg=0.5):
    """SVG path for a (multi)polygon, densified to `seg` degrees and projected to Robinson."""
    out = []
    geom = shapely.segmentize(geom, seg)
    for p in getattr(geom, "geoms", [geom]):
        if p.is_empty or p.geom_type != "Polygon":
            continue
        for ring in [p.exterior, *p.interiors]:
            xy = np.asarray(ring.coords)
            X, Y = proj(xy[:, 0], xy[:, 1])
            out.append("M" + " L".join(f"{x:.{nd}f} {y:.{nd}f}" for x, y in zip(X, Y)) + "Z")
    return "".join(out)


def line(lons, lats, nd=2):
    X, Y = proj(lons, lats)
    return "M" + " L".join(f"{x:.{nd}f} {y:.{nd}f}" for x, y in zip(X, Y))


# ---- statistics per box -------------------------------------------------------------------------
df = pd.read_csv(CSV, usecols=["member", "peak_block_date", "HeatWave", "energy_state", "persistence",
                               "HW_x_loc", "HW_y_loc"])
df = df[df["HeatWave"] == True].copy()
df["decade"] = pd.to_datetime(df["peak_block_date"]).dt.year // 10 * 10
pts = [Point(((x + 180) % 360) - 180, y) for x, y in zip(df["HW_x_loc"], df["HW_y_loc"])]
data = []
for b in BOXES:
    pg = prep(b["geom"])
    sub = df[[pg.contains(p) for p in pts]]
    shares = {}
    for per, (d0, d1) in [("hist", (1850, 2010)), ("ssp", (2020, 2090))]:
        s = sub[sub["decade"].between(d0, d1)].groupby("energy_state")["persistence"].sum().reindex(REGIMES, fill_value=0)
        shares[per] = [round(100 * v / s.sum(), 1) if s.sum() else 0 for v in s.values]
    days = sub.groupby("decade")["persistence"].sum()
    src = os.path.join(FIGDIR, f"heatwave_days_per_decade_{b['code']}_sum5members_ts_lwa.png")
    ls = LWASTATS[b["code"]]
    c = b["geom"].representative_point() if b["geom"].geom_type == "Polygon" else \
        max(b["geom"].geoms, key=lambda g: g.area).representative_point()
    data.append(dict(code=b["code"], label=b["label"], where=b["where"], continent=b["continent"],
                     band=b["band"], side=b["side"], n=int(len(sub)), days=int(sub["persistence"].sum()),
                     hist=shares["hist"], ssp=shares["ssp"],
                     lwa=dict(p88=ls["p88"], hist=ls["hist"], ssp=ls["ssp"], all=ls["all"]),
                     dpd_hist=round(days.reindex(range(1850, 2020, 10), fill_value=0).mean(), 0),
                     dpd_ssp=round(days.reindex(range(2020, 2100, 10), fill_value=0).mean(), 0),
                     path=path(b["geom"]), lx=round(float(proj(c.x, c.y)[0]), 2), ly=round(float(proj(c.x, c.y)[1]), 2)))
    print(b["code"], len(sub), flush=True)

land_path = path(land.simplify(0.15), nd=1)
ocean = path(box(-180, LAT0, 180, LAT1), nd=2)
xmin = float(proj(-180, LAT0)[0]); xmax = -xmin
ymin, ymax = float(proj(0, LAT1)[1]), float(proj(0, LAT0)[1])
viewbox = f"{xmin - 1:.1f} {ymin - 1:.1f} {xmax - xmin + 2:.1f} {ymax - ymin + 4.5:.1f}"
# graticule: meridians every 30 deg (curved in Robinson), parallels at the band edges (straight)
grid, lats = [], np.linspace(LAT0, LAT1, 60)
for lon in range(-150, 180, 30):
    grid.append(f'<path d="{line(np.full(lats.size, lon), lats)}"/>')
    tx, ty = proj(lon, LAT0)
    lab = "0°" if lon == 0 else f"{abs(lon)}°{'W' if lon < 0 else 'E'}"
    grid.append(f'<text x="{float(tx):.2f}" y="{float(ty) + 2.4:.2f}" text-anchor="middle">{lab}</text>')
for lat in (38, 48, 58, 72):
    grid.append(f'<path d="{line(np.linspace(-180, 180, 5), np.full(5, lat))}"/>')
    tx, ty = proj(-180, lat)
    grid.append(f'<text x="{float(tx) + 0.8:.2f}" y="{float(ty) - 0.4:.2f}">{lat}°N</text>')
tpl = open(os.path.join(HERE, "template.html")).read()
html = (tpl.replace("__LAND__", land_path).replace("__OCEAN__", ocean).replace("__GRID__", "".join(grid))
        .replace("__VIEWBOX__", viewbox).replace("__DATA__", json.dumps(data, ensure_ascii=False)))
os.makedirs(os.path.join(SITE, "figs"), exist_ok=True)
open(os.path.join(SITE, "index.html"), "w").write(html)
shutil.copy(os.path.join(HERE, "details.html"), os.path.join(SITE, "details.html"))
for b in BOXES:
    shutil.copy(os.path.join(FIGDIR, f"heatwave_days_per_decade_{b['code']}_sum5members_ts_lwa_nofoot.png"),
                os.path.join(SITE, "figs", f"{b['code']}.png"))
    shutil.copy(os.path.join(LWADIR, f"LWA_JJA_5-10km_pdf_{b['code']}_nofoot.png"),
                os.path.join(SITE, "figs", f"{b['code']}_lwa_pdf.png"))
print("viewBox", viewbox)
print("index.html", os.path.getsize(os.path.join(SITE, "index.html")) // 1024, "kB")
