import pyshtools as psh
from nastro import types as nt, graphics as ng, constants as nc
import numpy as np
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
from pathlib import Path


def colormap_data(field) -> tuple[nt.Vector, nt.Vector, nt.Vector]:
    lons, lats = np.meshgrid(field.lons(), field.lats())
    return lons, lats, field.data


LMAX = 300

if __name__ == "__main__":

    # Reference parameters from WSG84
    a = psh.constants.Earth.wgs84.a.value
    b = psh.constants.Earth.wgs84.b.value
    f = psh.constants.Earth.wgs84.f.value
    u0 = psh.constants.Earth.wgs84.u0.value
    rho_crust = 2670.0

    # Import gravity field data [GOCO06S - Only satellite data]
    gcoeffs = psh.datasets.Earth.GOCO06S()
    gfield = gcoeffs.expand(a=a, f=f).total * 1e5

    # Re-import gravity field data for full acceleration field
    gfcoeffs = psh.datasets.Earth.GOCO06S()
    gfcoeffs.set_coeffs(ls=0, ms=0, values=0)
    gfcoeffs.set_coeffs(ls=2, ms=0, values=0)
    gffield = gfcoeffs.expand(a=a, f=f, lmax=LMAX).total * (-1)

    # Import topography [Earth2012 - Rock equivalent]
    tcoeffs = psh.datasets.Earth.Earth2012.ret(lmax=LMAX)
    assert tcoeffs is not None
    tfield = tcoeffs.expand() / 1e3

    # Compute Bouguer anomalies
    scoeffs = psh.datasets.Earth.Earth2012.shape_ret(lmax=LMAX)
    bcoeffs = psh.SHGravCoeffs.from_shape(scoeffs, rho_crust, gm=gcoeffs.gm, lmax=LMAX)
    bcoeffs = bcoeffs.change_ref(r0=gcoeffs.r0)
    bcoeffs.set_coeffs(ls=0, ms=0, values=0)
    bcoeffs.set_coeffs(ls=2, ms=0, values=0)
    bcoeffs = gcoeffs.pad(lmax=LMAX) - bcoeffs
    bfield = bcoeffs.expand(a=a, f=f, lmax=LMAX).total * 1e5

    # Generate figures
    mediadir = Path(__file__).parents[2] / "report/media/q2"
    mediadir.mkdir(parents=True, exist_ok=True)

    figure, axes = plt.subplots(2, 2, figsize=(7, 5), layout="tight")

    gffield.plot(
        projection=ccrs.PlateCarree(central_longitude=0.1),
        cmap="RdBu_r",
        colorbar="bottom",
        cb_triangles="both",
        cb_label="$g$ [m/s$^2$]",
        tick_interval=[60, 30],
        ax=axes[0, 0],
    )
    gfield.plot(
        projection=ccrs.PlateCarree(central_longitude=0.1),
        cmap_limits=[-200, 200],
        cmap="RdBu_r",
        colorbar="bottom",
        cb_triangles="both",
        cb_label="Gravity anomaly [mGal]",
        tick_interval=[60, 30],
        ax=axes[0, 1],
    )
    tfield.plot(
        projection=ccrs.PlateCarree(central_longitude=0.1),
        cmap="RdBu_r",
        colorbar="bottom",
        cb_triangles="both",
        cb_label="Rock-equivalent topography [km]",
        tick_interval=[60, 30],
        ax=axes[1, 0],
    )
    bfield.plot(
        projection=ccrs.PlateCarree(central_longitude=0.1),
        cmap="RdBu_r",
        colorbar="bottom",
        cb_triangles="both",
        cb_label="Bouguer anomaly [mGal]",
        cmap_limits=[-600, 600],
        tick_interval=[60, 30],
        ax=axes[1, 1],
    )

    figure.savefig(mediadir / "fields.png")
