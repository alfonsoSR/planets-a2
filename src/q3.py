import pyshtools as psh
from nastro import types as nt, graphics as ng, constants as nc
import numpy as np
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
from pathlib import Path
from scipy.interpolate import griddata


def colormap_data(field) -> tuple[nt.Vector, nt.Vector, nt.Vector]:
    lons, lats = np.meshgrid(field.lons(), field.lats())
    return lons, lats, field.data


def save_to_gmt(field, filename) -> None:

    if Path(filename).exists():
        Path(filename).unlink()
    Path(filename).parent.mkdir(parents=True, exist_ok=True)

    # From corners to centers
    _lons = field.lons()
    _lats = field.lats()
    lons = 0.5 * (_lons[:-1] + _lons[1:])
    lats = 0.5 * (_lats[:-1] + _lats[1:])
    __data = field.data
    _data = 0.5 * (__data[:-1] + __data[1:])
    data = 0.5 * (_data[:, :-1] + _data[:, 1:])

    with open(filename, "w") as file:
        for idx, lat in enumerate(lats):
            for jdx, lon in enumerate(lons):
                file.write(f"{lon} {lat} {data[idx, jdx]}\n")

    # with open(filename, "w") as file:
    #     for idx, lat in enumerate(field.lats()):
    #         for jdx, lon in enumerate(field.lons()):
    #             file.write(f"{lon} {lat} {data[idx, jdx]}\n")

    return None


LMAX = 90
outdir = Path(__file__).parent / "GSH/input"

if __name__ == "__main__":

    # Reference parameters from WSG84
    a = psh.constants.Earth.wgs84.a.value
    b = psh.constants.Earth.wgs84.b.value
    f = psh.constants.Earth.wgs84.f.value
    u0 = psh.constants.Earth.wgs84.u0.value
    rho_crust = 2670.0

    # Topography
    retcoeffs = psh.datasets.Earth.Earth2012.ret(lmax=LMAX)
    assert retcoeffs is not None
    retfield = retcoeffs.expand()
    save_to_gmt(retfield, outdir / "custom/topography.gmt")

    # Bouguer anomalies
    gcoeffs = psh.datasets.Earth.GOCO06S(lmax=LMAX)
    # with open(outdir / "custom/gravity_coeffs.gmt", "w") as file:

    #     for n in gcoeffs.degrees():
    #         for m in range(n, gcoeffs.lmax + 1):
    #             file.write(
    #                 f"{m} {n} {gcoeffs.coeffs[0, m, n]} {gcoeffs.coeffs[1, m, n]}\n"
    #             )

    # gcoeffs.to_file(outdir / "custom/gravity_coeffs.txt")
    gfield = gcoeffs.expand(a=a, f=f, lmax=LMAX).total

    save_to_gmt(gfield, outdir / "custom/gravity_anomaly.gmt")
    # exit(0)

    scoeffs = psh.datasets.Earth.Earth2012.shape_ret(lmax=LMAX)
    bccoeffs = psh.SHGravCoeffs.from_shape(scoeffs, rho_crust, gm=gcoeffs.gm, lmax=LMAX)
    bccoeffs = bccoeffs.change_ref(r0=gcoeffs.r0)
    bccoeffs.set_coeffs(ls=0, ms=0, values=0)
    bccoeffs.set_coeffs(ls=2, ms=0, values=0)
    bcfield = bccoeffs.expand(a=a, f=f, lmax=LMAX).total  # Bouguer correction field
    bcoeffs = gcoeffs.pad(lmax=LMAX) - bccoeffs
    bfield = bcoeffs.expand(a=a, f=f, lmax=LMAX).total  # Bouguer anomalies field

    save_to_gmt(bcfield, outdir / "custom/bouguer_correction.gmt")
    save_to_gmt(bfield, outdir / "custom/bouguer_anomalies.gmt")
    exit(0)

    # Import gravity field data [GOCO06S - Only satellite data]
    gcoeffs = psh.datasets.Earth.GOCO06S()
    gfield = gcoeffs.expand(a=a, f=f).total

    # Compute Bouguer anomalies
    scoeffs = psh.datasets.Earth.Earth2012.shape_ret(lmax=LMAX)
    bccoeffs = psh.SHGravCoeffs.from_shape(scoeffs, rho_crust, gm=gcoeffs.gm, lmax=LMAX)
    bccoeffs = bccoeffs.change_ref(r0=gcoeffs.r0)
    bccoeffs.set_coeffs(ls=0, ms=0, values=0)
    bccoeffs.set_coeffs(ls=2, ms=0, values=0)
    bcfield = bccoeffs.expand(a=a, f=f, lmax=LMAX).total  # Bouguer correction field
    bcoeffs = gcoeffs.pad(lmax=LMAX) - bccoeffs
    bfield = bcoeffs.expand(a=a, f=f, lmax=LMAX).total  # Bouguer anomalies field

    # Save Bouguer correction field to file
    save_to_gmt(bcfield, outdir / "bcfield.gmt")
    save_to_gmt(bfield, outdir / "bfield.gmt")
    save_to_gmt(gfield, outdir / "gfield.gmt")

    # Generate topography to compare with Bart
    btcoeffs = psh.datasets.Earth.Earth2012.topo_bathy(lmax=LMAX)
    assert btcoeffs is not None
    btfield = btcoeffs.expand() * 1e-3
    save_to_gmt(btfield, outdir / "btfield.gmt")

    retcoeffs = psh.datasets.Earth.Earth2012.ret(lmax=LMAX)
    assert retcoeffs is not None
    retfield = retcoeffs.expand() * 1e-3
    save_to_gmt(retfield, outdir / "retfield.gmt")

    # Generate lower boundary for the model
    lgrid = np.loadtxt(outdir / "retfield.gmt").T
    lgrid[-1] = -200 * np.ones_like(lgrid[-1])
    np.savetxt(outdir / "lgrid.gmt", lgrid.T)

    exit(0)

    # scoeffs = psh.datasets.Earth.Earth2012.shape_ret(lmax=LMAX)
    # bcoeffs = psh.SHGravCoeffs.from_shape(scoeffs, rho_crust, gm=gcoeffs.gm, lmax=LMAX)
    # bcoeffs = bcoeffs.change_ref(r0=gcoeffs.r0)
    # bcoeffs.set_coeffs(ls=0, ms=0, values=0)
    # bcoeffs.set_coeffs(ls=2, ms=0, values=0)
    # bcoeffs = gcoeffs.pad(lmax=LMAX) - bcoeffs
    # bfield = bcoeffs.expand(a=a, f=f, lmax=LMAX).total

    print(bfield.data.shape)
    print(bfield.lats().shape)
    print(bfield.lons().shape)
    exit(0)
    assert tcoeffs is not None
    tfield = tcoeffs.expand()
    h = tfield.data
    print(h.shape)
    print(h)
    tfield.plot(colorbar="bottom", projection=ccrs.PlateCarree(central_longitude=0.1))
    plt.show()

    # lmax = 200
    # alpha = ncat.Earth.mu / (ncat.Earth.R**2)
    # coeffs_array = np.zeros((2, lmax + 1, lmax + 1))
    # coeffs_array[0, 0, 0] = 1.0
    # coeffs_array[0, 2, 0] = ncat.Earth.j2
    # coeffs = psh.SHCoeffs.from_array(coeffs_array)
    # cfield = coeffs.expand(grid="DH2") * alpha
    # cfield.plot(colorbar="bottom")
    # plt.show()
    exit(0)

    lats = np.linspace(-0.5 * np.pi, 0.5 * np.pi, 300)
    lons = np.linspace(0, 2.0 * np.pi, 600)
    LON, LAT = np.meshgrid(lons, lats)
    model = psh.datasets.Earth.GOCO06S()
    shape_data = model.r0 * np.ones_like(LON)
    shape = psh.SHGrid.from_array(shape_data)
    gcoeffs = psh.SHGravCoeffs.from_shape(shape, rho=2670.0, gm=model.gm, lmax=1)
    gfield = gcoeffs.expand()
    gfield.plot(colorbar="bottom")
    # grid = psh.SHGrid.from_array(shape)
    # grid.plot(colorbar="bottom")
    plt.show()
    exit(0)
    a = psh.constants.Earth.wgs84.a.value
    f = psh.constants.Earth.wgs84.f.value
    gcoeffs = psh.SHGravCoeffs.from_zeros(
        2, gm=model.gm, r0=model.r0, omega=model.omega
    )
    gcoeffs.set_coeffs(values=model.coeffs[0, 0, 0], ls=0, ms=0)
    # gcoeffs.set_coeffs(values=model.coeffs[0, 2, 0], ls=2, ms=0)
    # gcoeffs.set_coeffs(values=model.coeffs[0, 2, 1], ls=2, ms=1)
    # gcoeffs.set_coeffs(values=model.coeffs[0, 2, 2], ls=2, ms=2)
    zfield = gcoeffs.expand(r=model.r0)
    zfield.plot_total(colorbar="bottom")
    plt.show()
    exit()

    h = 6370e3 * (1 + np.cos(LAT))

    data = np.random.uniform(-1, 1, size=LON.shape)
    grid = psh.SHGrid.from_array(h)
    coeffs = psh.SHGravCoeffs.from_shape(
        grid,
        rho=2670.0,
        omega=model.omega,
        gm=model.gm,
        lmax=2,
        nmax=9,
    )
    gfield = coeffs.expand(a=a, f=f, lmax=2)
    gfield.plot(colorbar="bottom")
    print(coeffs)
    grid.plot(colorbar="bottom")

    # # Compute Bouguer correction
    # scoeffs = psh.datasets.Earth.Earth2012.shape_ret(lmax=400)
    # shape = scoeffs.expand()

    # __gcoeffs = psh.datasets.Earth.GOCO06S()

    # gcoeffs = psh.SHGravCoeffs.from_shape(shape, rho=2670.0, gm=__gcoeffs.gm, lmax=400)
    # gfield = gcoeffs.expand()

    # shape.plot(colorbar="bottom")
    # gfield.plot(colorbar="bottom")
    plt.show()
