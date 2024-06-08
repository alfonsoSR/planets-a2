from pathlib import Path
from pyshtools.datasets.Earth import Earth2014, EGM2008, Earth2012
from pyshtools.constants import Earth
from cartopy import crs as ccrs
from nastro import graphics as ng
import numpy as np

outdir = Path(__file__).parents[1] / "output"
sourcedir = outdir.parent / "input"
use_water = False


def colormap_data(field):

    lons, lats = np.meshgrid(field.lons(), field.lats())
    return lons, lats, field.data


if __name__ == "__main__":

    # CONFIGURATION
    LMAX = 400
    TH_GRAV = 200

    gcoeffs = EGM2008(lmax=LMAX)
    tcoeffs_surface = Earth2014.surface(lmax=LMAX)
    tcoeffs_bedrock = Earth2014.bedrock(lmax=LMAX)

    # Retrieve reference values
    a = Earth.wgs84.a.value
    f = Earth.wgs84.f.value

    print(f"Reference radius: {a} m")
    print(f"Reference flattening: {f}")
    print(gcoeffs.info())

    gcoeffs.to_file(f"./egm2008-{LMAX}.txt")

    gfield = gcoeffs.expand(a=a, f=f, lmax=LMAX)
    tfield = tcoeffs_surface.expand(lmax=LMAX)

    setup = ng.PlotSetup(aspect="equal", grid=False, xscilimits=(0, 3))

    with ng.Mosaic("ab") as figure:

        with figure.subplot(setup) as tmap:
            tmap.colormap(*colormap_data(tfield))

        with figure.subplot(setup) as gmap:
            gmap.colormap(*colormap_data(gfield))

    # with ng.SingleAxis() as fig:
    #     fig.imshow(pot.data)

    # PROJECTION = ccrs.PlateCarree(central_longitude=0.1)
    # fig, ax = gfield.expand(a=a, f=f, lmax=LMAX).plot(
    #     projection=PROJECTION, cmap="GnBu", show=False
    # )
    # plt.show()
    # plt.show()

    # gfield = GravityField(outdir / "Crust10_crust.mat")
    # gdata = General(outdir / "data_Crust10_crust_0_179_26-May-2024 11:50:36.mat")
