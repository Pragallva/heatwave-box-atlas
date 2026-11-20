# Heatwave Box Atlas: files to recreate it

Everything needed to rebuild the Heatwave Box Atlas web page (`2026-11-18/github_pages/heatwave-box-atlas/`)
from the heatwave track table, plus the scripts that made the intermediate caches from the raw model output on NIRD.
Keep this folder layout: the scripts find each other's outputs through relative paths (`2026-10-07/`, `2026-10-08/`,
`2026-11-18/`).

Model: MPI-ESM1.2-LR (MPI Grand Ensemble), members 1–5, historical (1850–2014) + SSP3-7.0 (2015–2100), JJA.

## Contents

| Path | What it is |
|---|---|
| `inputs/heatwave_tracks.csv.gz` | Heatwave/blocking track table (gzipped copy of `pandas_track_data_summer_5-12km_summer_heatwaves_more_area_smooth_baseline_with_land_ocean_mask_and_dry_soil_Bo_ratio_0.5_threshold.csv` from `/projects/NS11084K/MPI-LE/daily/post_processing/blocking_info_JJA/mid_to_upper_level/JJA_NH_Ac=65_new_heatwave_centroid/`) |
| `inputs/natural_earth/` | Natural Earth 1:110m admin-0 countries (used for the land-clipped boxes and the map) |
| `2026-10-07/make_heatwave_event_skinT_cache.py` | Per-heatwave skin temperature → `heatwave_event_ts_days-3to3_cache.csv` (included; needs NIRD composites to rerun) |
| `2026-10-07/make_heatwave_event_LWA_cache.py` | Per-heatwave LWA, 5–10 km → `heatwave_event_LWA_5-10km_days-3to3_cache.csv` (included; needs NIRD composites to rerun) |
| `2026-10-07/plot_heatwave_days_per_decade_EU-WS_NA-WS.py` | The 4-panel figure per box |
| `2026-10-08/make_LWA_JJA_5-10km_hist_regions.py` | Daily LWA histograms per box → `LWA_JJA_5-10km_hist_regions_cache.npz` (included; needs NIRD daily LWA to rerun) |
| `2026-10-08/plot_LWA_JJA_38-85N_5-10km_pdf.py` | LWA distribution figure per box + `LWA_JJA_5-10km_pdf_regions_stats.json` |
| `2026-11-18/dashboard_robinson/build_dashboard_robinson.py` | Builds the web page (Robinson map, sidebar numbers) and copies figures into the site |
| `2026-11-18/dashboard_robinson/template.html`, `details.html` | Page template and the methods page |
| `2026-11-18/dashboard_robinson/make_preview.py` | `preview.png` map image for a GitHub README |
| `2026-11-18/github_pages/heatwave-box-atlas/` | The finished site: `index.html`, `details.html`, `figs/`, `preview.png`, `README.md` |

## Data not included in this repository

Two files are too large for the repository and are listed in `.gitignore`:

| File | Size | Needed for |
|---|---|---|
| `inputs/heatwave_tracks.csv.gz` | 46 MB | steps 1 and 3 below |
| `2026-10-08/LWA_JJA_5-10km_hist_regions_cache.npz` | 62 MB | step 2 below |

Put them at those paths (from NIRD, or from the author on request) before rebuilding. The finished site in
`2026-11-18/github_pages/heatwave-box-atlas/` is in the repository and is published by GitHub Pages
from the `gh-pages` branch: https://pragallva.github.io/heatwave-box-atlas/

After rebuilding, commit to `main` and update the live site with

```bash
git subtree push --prefix 2026-11-18/github_pages/heatwave-box-atlas origin gh-pages
```

## Python environment

Python 3.10 with numpy 1.26, pandas 2.3, scipy 1.14, matplotlib 3.8, shapely 2.0, cartopy 0.22, pyproj 3.7,
xarray 2024.1 (xarray only for the cache scripts). For example:

```bash
conda create -n hwatlas -c conda-forge python=3.10 numpy pandas scipy matplotlib shapely cartopy pyproj xarray
```

Also tested with Python 3.11, numpy 2.3, pandas 2.3, matplotlib 3.10, shapely 2.1, cartopy 0.25 (2026-11-20): the
numbers and HTML come out identical; figures differ only in font rendering.

Cartopy downloads the Natural Earth file on first use. Offline, copy `inputs/natural_earth/*` to
`~/.local/share/cartopy/shapefiles/natural_earth/cultural/`.

## Rebuild the site from the included inputs and caches

Run from this folder:

```bash
export HW_TRACK_CSV=$PWD/inputs/heatwave_tracks.csv.gz

# 1. 4-panel figures, 18 boxes, without footnotes -> 2026-10-07/heatwave_days_per_decade_<box>_sum5members_ts_lwa_nofoot.png
(cd 2026-10-07 && python plot_heatwave_days_per_decade_EU-WS_NA-WS.py --latlon --sum --ts --lwa --png --nofoot)

# 2. LWA distribution figures + stats -> 2026-10-08/LWA_JJA_5-10km_pdf_<box>_nofoot.png, LWA_JJA_5-10km_pdf_regions_stats.json
(cd 2026-10-08 && python plot_LWA_JJA_38-85N_5-10km_pdf.py --regions --nofoot)

# 3. Web page -> 2026-11-18/github_pages/heatwave-box-atlas/ (index.html, details.html, figs/)
(cd 2026-11-18/dashboard_robinson && python build_dashboard_robinson.py && python make_preview.py)
```

Open `2026-11-18/github_pages/heatwave-box-atlas/index.html` in a browser. The site is static; that folder is what goes
to GitHub Pages.

## Rebuild the caches from raw model output (NIRD only)

These read tens of GB of model output under `/projects/NS11084K/MPI-LE/daily/` and are only needed if the
heatwave definition or the averaging changes:

```bash
(cd 2026-10-07 && python make_heatwave_event_skinT_cache.py && python make_heatwave_event_LWA_cache.py)   # composite_data/SHF, composite_data/LWA
(cd 2026-10-08 && python make_LWA_JJA_5-10km_hist_regions.py 32)                                           # daily LWA, 3765 files, 32 processes
```

## Definitions

See `2026-11-18/github_pages/heatwave-box-atlas/details.html` (also linked from the atlas page).
