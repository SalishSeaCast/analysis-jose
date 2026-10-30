#!/usr/bin/env python
"""Convert noleap/365_day time axes to gregorian, in place.
Dry run by default. Usage:
    python fix_calendar.py '/path/to/2023/*_grid_T.nc' '/path/to/winds/*.nc'
    python fix_calendar.py --apply '/path/to/2023/*.nc'
"""
import sys, glob, argparse
import numpy as np
import netCDF4 as nc
import cftime

BAD = {"noleap", "365_day"}
NEW = "gregorian"

def convert_var(ds, v, apply):
    """Re-encode one time variable. Returns (changed_values, skipped_reason)."""
    vals = np.ma.filled(v[:], np.nan).astype("float64")
    ok = np.isfinite(vals)
    dates = cftime.num2date(vals[ok], v.units, v.calendar)
    # a noleap model year that is a leap year has no Feb 29, so dates drift: don't guess
    for d in (dates[0], dates[-1]):
        y = d.year
        if (y % 4 == 0 and (y % 100 != 0 or y % 400 == 0)) and d.month >= 3:
            return False, f"leap year {y} after Feb 28 (needs manual look)"
    new_dates = [cftime.datetime(d.year, d.month, d.day, d.hour, d.minute, d.second,
                                 d.microsecond, calendar=NEW) for d in dates]
    new_vals = cftime.date2num(new_dates, v.units, NEW)
    changed = not np.allclose(new_vals, vals[ok])
    if apply:
        out = vals.copy()
        out[ok] = new_vals
        v[:] = out.astype(v.dtype) if np.issubdtype(v.dtype, np.floating) else np.round(out).astype(v.dtype)
        v.calendar = NEW
    return changed, None

def process(path, apply):
    with nc.Dataset(path, "r+" if apply else "r") as ds:
        todo = [v for v in ds.variables.values() if getattr(v, "calendar", None) in BAD]
        if not todo:
            return "ok (no noleap)"
        msgs = []
        for v in todo:
            changed, skip = convert_var(ds, v, apply)
            if skip:
                return f"SKIPPED {v.name}: {skip}"
            msgs.append(f"{v.name}: {'values shifted' if changed else 'label only'}")
            b = getattr(v, "bounds", None)          # e.g. time_counter_bounds
            if b in ds.variables and "units" in ds[b].ncattrs():
                pass  # bounds usually inherit units/calendar from parent; see note below
        return ("FIXED " if apply else "WOULD FIX ") + ", ".join(msgs)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("patterns", nargs="+")
    ap.add_argument("--apply", action="store_true")
    a = ap.parse_args()
    files = sorted({f for p in a.patterns for f in glob.glob(p)})
    print(f"{len(files)} files ({'APPLY' if a.apply else 'dry run'})")
    for f in files:
        try:
            print(f, "->", process(f, a.apply))
        except Exception as e:
            print(f, "-> ERROR", repr(e))