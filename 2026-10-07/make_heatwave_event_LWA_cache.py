"""
Per-heatwave-event LWA cache: LWA_events from composite_data/LWA averaged over height 5-ZTOP km (default 5-10 km;
5-12 km = blocking-detection layer), the composite box (lat +/-4.5, lon +/-10.5 deg around the event centre) and days -3..+3 around the block peak.
Output: heatwave_event_LWA_5-<ZTOP>km_days-3to3_cache.csv with columns index (= regime-CSV row), domain, scenario,
member, LWA.

Usage:
    python make_heatwave_event_LWA_cache.py          # 5-10 km
    python make_heatwave_event_LWA_cache.py 12       # 5-12 km
"""
import glob
import os
import sys

import pandas as pd
import xarray as xr

HERE = os.path.dirname(os.path.abspath(__file__))
ZTOP = int(sys.argv[1]) if len(sys.argv) > 1 else 10   # km
SRC = ("/projects/NS11084K/MPI-LE/daily/post_processing/blocking_info_JJA/mid_to_upper_level/"
       "JJA_NH_Ac=65_new_heatwave_centroid/composite_data/LWA/")
rows = []
for f in sorted(glob.glob(SRC + "*_DAY*.nc")):
    domain, scen, mem = os.path.basename(f).split("_")[:3]
    with xr.open_dataset(f, chunks={"event": 50}) as ds:
        lwa = (ds["LWA_events"].sel(time=slice(-3, 3), height=slice(5000, ZTOP * 1000))
               .mean(dim=("time", "height", "lat", "lon"), skipna=True).compute().values)
        idx = ds["index_events"].values
    rows.append(pd.DataFrame(dict(index=idx, domain=domain, scenario=scen, member=mem, LWA=lwa)))
    print(os.path.basename(f), len(idx), flush=True)
out = pd.concat(rows, ignore_index=True)
print("events:", len(out), "unique index:", out["index"].nunique())
out.to_csv(os.path.join(HERE, f"heatwave_event_LWA_5-{ZTOP}km_days-3to3_cache.csv"), index=False)
