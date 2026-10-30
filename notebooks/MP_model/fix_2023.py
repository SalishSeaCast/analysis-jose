#!/usr/bin/env python
"""Convert noleap time axes to gregorian in the 2023 biol_T files, in place.
    python fix_biol_2023.py            # dry run + check
    python fix_biol_2023.py --apply    # convert + check
"""
import os, sys
from datetime import datetime, timedelta
import numpy as np
import netCDF4 as nc
import cftime
import xarray as xr

def get_Biol_path(date):
    # Make biol path
    path = '/data/jvalenti/results/'
    datestr = [date.strftime(fmt) for fmt in ('%d%b%y', '%Y%m%d-%Y%m%d')]
    path = os.path.join(path, datestr[0].lower(), f'SHEM_1h_{datestr[1]}_biol_T.nc')
    return path  

NOLEAP = {"noleap", "365_day"}
apply = "--apply" in sys.argv

def biol_files_2023():
    fs, d = [], datetime(2023, 1, 1)
    while d.year == 2023:
        fs.append(get_Biol_path(d))
        d += timedelta(days=1)
    return fs

def reencode(var, units, calendar):
    """Decode with noleap, re-encode the same dates as gregorian."""
    vals = np.ma.filled(var[:], np.nan).astype("float64")
    ok = np.isfinite(vals)
    old = cftime.num2date(vals[ok], units, calendar)
    new = [cftime.datetime(d.year, d.month, d.day, d.hour, d.minute, d.second,
                           d.microsecond, calendar="gregorian") for d in old]
    out = vals.copy()
    out[ok] = cftime.date2num(new, units, "gregorian")
    if not np.issubdtype(var.dtype, np.floating):
        out = np.round(out)
    return out.astype(var.dtype)

def fix(path_):
    done = []
    with nc.Dataset(path_, "r+" if apply else "r") as ds:
        for name, v in list(ds.variables.items()):
            if str(getattr(v, "calendar", "")).lower() not in NOLEAP:
                continue
            targets = [v]
            b = getattr(v, "bounds", None)
            if b in ds.variables:
                targets.append(ds.variables[b])
            if apply:
                for t in targets:
                    t[:] = reencode(t, v.units, v.calendar)
                v.calendar = "gregorian"
            done.append(name)
    return done

files = biol_files_2023()
print(f"{len(files)} biol files ({'APPLY' if apply else 'dry run'})")
n = 0
for f in files:
    if not os.path.exists(f):
        print("MISSING", f); continue
    try:
        m = fix(f)
        if m:
            n += 1
            print("FIXED   " if apply else "WOULD FIX", f, m)
    except Exception as e:
        print("ERROR", f, repr(e))
print(f"\n{n} files {'fixed' if apply else 'to fix'}")

print("\nCheck: biol files whose time_counter still decodes as non-datetime64")
bad = 0
for f in files:
    if os.path.exists(f):
        with xr.open_dataset(f) as ds:
            if "datetime64" not in str(ds["time_counter"].dtype):
                bad += 1
                print("STILL CFTIME", f)
print(f"{bad} problem files")