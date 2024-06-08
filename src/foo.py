import pyshtools as psh
from nastro import types as nt
import numpy as np
from matplotlib import pyplot as plt
import cartopy.crs as ccrs


def colormap_data(field) -> tuple[nt.Vector, nt.Vector, nt.Vector]:
    lons, lats = np.meshgrid(field.lons(), field.lats())
    return lons, lats, field.data


LMAX = 400

if __name__ == "__main__":

    # Import gravity field data
    gcoeffs = psh.datasets.Earth.EGM2008(lmax=LMAX)
    # gcoeffs.coeffs[0, 0, 0] = 0.0
    # gcoeffs.coeffs[0, 2, 0] = 0.0

    # Use WSG84 ellipsoid as reference to expand gravity field
    gfield = gcoeffs.expand(
        a=psh.constants.Earth.wgs84.a.value,
        f=psh.constants.Earth.wgs84.f.value,
        lmax=LMAX,
    ).total

    plt.show()
