"""
Per-heatwave-event radiative skin temperature, derived from the stored surface-flux composites
(composite_data/SHF/all_members_more_area_smooth_baseline; no daily ts output is available on NIRD):

    Ts = [ (rlus - (1 - eps) * rlds) / (eps * sigma) ]^(1/4),   eps = 0.996 (ECHAM6 surface emissivity)

computed for each grid point and day, then averaged over the +/-10 deg box around the heatwave centroid and
days -3..+3. Output: heatwave_event_ts_days-3to3_cache.csv (index = regime-CSV row, ts and tas in deg C).

Usage:
    python make_heatwave_event_skinT_cache.py
"""
import glob
import os

import pandas as pd
import xarray as xr

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = ("/projects/NS11084K/MPI-LE/daily/post_processing/blocking_info_JJA/mid_to_upper_level/"
       "JJA_NH_Ac=65_new_heatwave_centroid/composite_data/SHF/all_members_more_area_smooth_baseline/")
EPS, SIGMA = 0.996, 5.670374419e-8
rows = []
for f in sorted(glob.glob(SRC + "*.nc")):
    domain, scen, mem = os.path.basename(f).split("_")[:3]
    with xr.open_dataset(f) as ds:
        sel = dict(time=slice(-3, 3))
        ts = ((ds["rlus_events"].sel(**sel) - (1 - EPS) * ds["rlds_events"].sel(**sel)) / (EPS * SIGMA)) ** 0.25
        rows.append(pd.DataFrame(dict(index=ds["index_events"].values, domain=domain, scenario=scen, member=mem,
                                      year=pd.to_datetime(ds["peak_block_date_events"].values).year,
                                      ts=ts.mean(("time", "lat", "lon")).values - 273.15,
                                      tas=ds["tas_events"].sel(**sel).mean(("time", "lat", "lon")).values - 273.15)))
    print(os.path.basename(f), flush=True)
out = pd.concat(rows, ignore_index=True)
out.to_csv(os.path.join(HERE, "heatwave_event_ts_days-3to3_cache.csv"), index=False)
print(out.groupby("domain")[["ts", "tas"]].describe().T)
print("land ts - tas:", (out.ts - out.tas)[out.domain == "land"].describe())
