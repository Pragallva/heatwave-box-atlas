"""
Histogram cache of daily LWA (averaged over height 5-10 km), JJA, MPI-GE members 1-5, historical (1850-2014) +
SSP3-7.0 (2015-2100), for the 18 continent lat/lon boxes of the heatwave dashboard (land-clipped Natural Earth 1:110m
polygons, same definitions as 2026-10-07/dashboard/build_dashboard.py). A grid point (1.875 deg) belongs to a box
if its centre lies inside the polygon. One histogram per file and box, 0.5 m/s bins 0-300 m/s (+ overflow),
counted per grid point-day and with cos(lat) area weights.
Output: LWA_JJA_5-10km_hist_regions_cache.npz (counts/wcounts: file x box x bin; codes; npts per box).

Usage:
    python make_LWA_JJA_5-10km_hist_regions.py [nproc]
"""
import glob
import os
import re
import sys
from multiprocessing import Pool

import numpy as np
import xarray as xr
import cartopy.io.shapereader as shpreader
from shapely.geometry import box
from shapely.ops import unary_union
from shapely import vectorized

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = "/projects/NS11084K/MPI-LE/daily"
EDGES = np.append(np.arange(0, 300.01, 0.5), np.inf)
Z0, Z1 = 5000, 10000

recs = list(shpreader.Reader(shpreader.natural_earth("110m", "cultural", "admin_0_countries")).records())
G = {r.attributes["NAME"]: r.geometry for r in recs}
europe = unary_union([r.geometry for r in recs if r.attributes["CONTINENT"] == "Europe"]).intersection(box(-25, 34, 45, 72))
north_am = unary_union([r.geometry for r in recs if r.attributes["CONTINENT"] == "North America"
                        and r.attributes["NAME"] != "Greenland"])
asia = unary_union([r.geometry for r in recs if r.attributes["CONTINENT"] == "Asia"]
                   + [G["Russia"].intersection(box(45, 30, 180, 85))])
CODES, POLYS = [], []
for cc, geom, bands, west, east in [
        ("NA", north_am, [(38, 48), (48, 58), (58, 85)], [(-180, -100)], [(-100, -20)]),
        ("EU", europe, [(38, 48), (48, 58), (58, 72)], [(-25, 15)], [(15, 45)]),
        ("AS", asia, [(38, 48), (48, 58), (58, 80)], [(0, 100)], [(100, 180), (-180, 0)])]:
    for side, lons in [("W", west), ("E", east)]:
        for (lo, hi), b in zip(bands, "SCN"):
            CODES.append(f"{cc}-{side}{b}")
            POLYS.append(unary_union([geom.intersection(box(l0, lo, l1, hi)) for l0, l1 in lons]))

f0 = sorted(glob.glob(f"{BASE}/member01/historical/LWA/month06/*.nc"))[0]
with xr.open_dataset(f0) as ds:
    LAT = ds["lat"].sel(lat=slice(30, 90)).values
    LON = ds["lon"].values
lon2, lat2 = np.meshgrid(((LON + 180) % 360) - 180, LAT)
MASKS = np.array([vectorized.contains(p, lon2, lat2) for p in POLYS])        # box x lat x lon
WLAT = np.cos(np.deg2rad(LAT))


def one(f):
    with xr.open_dataset(f) as ds:
        lwa = ds["LWA"].sel(lat=slice(30, 90), height=slice(Z0, Z1)).mean("height").values  # time, lat, lon
    h = np.zeros((len(CODES), EDGES.size - 1)); hw = np.zeros_like(h)
    for k, m in enumerate(MASKS):
        v = lwa[:, m]                                                       # time x npts
        w = np.broadcast_to(WLAT[np.nonzero(m)[0]][None, :], v.shape)
        ok = np.isfinite(v)
        h[k] = np.histogram(v[ok], EDGES)[0]
        hw[k] = np.histogram(v[ok], EDGES, weights=w[ok])[0]
    return h, hw


if __name__ == "__main__":
    print(dict(zip(CODES, MASKS.sum(axis=(1, 2)))), flush=True)
    files, meta = [], []
    for m in range(1, 6):
        for scen in ("historical", "ssp370"):
            for mon in (6, 7, 8):
                for f in sorted(glob.glob(f"{BASE}/member{m:02d}/{scen}/LWA/month{mon:02d}/LWA_dailyZ_*.nc")):
                    files.append(f)
                    meta.append((m, scen, mon, int(re.search(r"_(\d{4})\d{4}-\d{8}\.nc$", f).group(1))))
    print(len(files), "files", flush=True)
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 32) as p:
        res = []
        for i, r in enumerate(p.imap(one, files, chunksize=4)):
            res.append(r)
            if i % 200 == 0:
                print(i, flush=True)
    np.savez_compressed(os.path.join(HERE, "LWA_JJA_5-10km_hist_regions_cache.npz"), edges=EDGES, codes=CODES,
                        npts=MASKS.sum(axis=(1, 2)), counts=np.array([r[0] for r in res]).astype(np.int64),
                        wcounts=np.array([r[1] for r in res]).astype(np.float32),
                        member=[x[0] for x in meta], scenario=[x[1] for x in meta], month=[x[2] for x in meta],
                        year=[x[3] for x in meta])
    print("done")
