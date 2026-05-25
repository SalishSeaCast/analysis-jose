from pathlib import Path
from multiprocessing import Pool
import xarray as xr
import numpy as np
import pandas as pd


MASK_PATH = "/home/jvalenti/MOAD/grid/mesh_mask202108.nc"
VAR_NAME = "heterotrophic_bacteria"  # change this


def load_volume():
    mask = xr.open_dataset(MASK_PATH)

    e1t = mask["e1t"][0, :, :]
    e2t = mask["e2t"][0, :, :]
    e3t = mask["e3t_0"][0, :, :, :]
    tmask = mask["tmask"][0, :, :, :]

    vol = (e1t * e2t * e3t).where(tmask == 1)

    return vol


def process_one_day(path):
    path = Path(path)

    vol = load_volume()
    vol_sum = vol.sum(("deptht", "y", "x"))

    ds = xr.open_dataset(path)
    hbac = ds[VAR_NAME]

    hbac_day = hbac.mean("time_counter")

    hbac_total = (hbac_day * vol).sum(("deptht", "y", "x"))
    hbac_vmean = hbac_total / vol_sum

    hbac_vstd = np.sqrt(
        (((hbac_day - hbac_vmean) ** 2) * vol).sum(("deptht", "y", "x"))
        / vol_sum
    )

    date = pd.to_datetime(ds.time_counter.values[0]).date()

    result = {
        "date": date,
        "hbac_total": float(hbac_total.values),
        "hbac_vmean": float(hbac_vmean.values),
        "hbac_vstd": float(hbac_vstd.values),
        "file": str(path),
    }

    ds.close()

    return result


if __name__ == "__main__":

    files = sorted(Path("/path/to/results").glob("*/SHEM_1h_*_biol_T.nc"))

    nproc = 24

    with Pool(processes=nproc) as pool:
        results = pool.map(process_one_day, files)

    df = pd.DataFrame(results).sort_values("date")
    df.to_csv("hbac_daily_diagnostics.csv", index=False)