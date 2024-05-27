from pathlib import Path
from processing import GravityField, General

outdir = Path(__file__).parents[1] / "output"
sourcedir = outdir.parent / "input"


if __name__ == "__main__":

    gfield = GravityField(outdir / "Crust10_crust.mat")
    gdata = General(outdir / "data_Crust10_crust_0_179_26-May-2024 11:50:36.mat")
