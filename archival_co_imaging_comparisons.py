
from glob import glob
from pathlib import Path
import os

data_path = Path("/reduction11/erickoch/M31/sma/per_mosaic/")

# M31_154

this_target = 'M31_154'

os.chdir(data_path / this_target)


ms_files = list(glob("*.ms"))

if len(ms_files) == 0:
    raise ValueError(f"Found no MSs in {this_target}")

tclean(vis=ms_files,
    imagename=str(f"{this_target}_co21_test_5kms"),
    imsize=512,
    cell='0.5arcsec',
    spw='',
    width='5km/s',
    nchan=20,
    start='-100km/s',
    specmode='cube',
    gridder='standard',
    deconvolver='multiscale',
    scales=[0, 3, 5, 10],
    niter=5000,
    cyclefactor=5.0,
    nsigma=3.0,
    perchanweightdensity=False,
    mosweight=False,
    pbcor=False,
    pblimit=0.15,
    restfreq='230.538GHz',
    phasecenter='',
    outframe='LSRK',
    )

impbcor(str(f"{this_target}_co21_test_5kms.image"),
        pbimage=str(f"{this_target}_co21_test_5kms.pb"),
        outfile=str(f"{this_target}_co21_test_5kms.image.pbcor"))



# M31_157

this_target = 'M31_157'

os.chdir(data_path / this_target)


ms_files = list(glob("*.ms"))

if len(ms_files) == 0:
    raise ValueError(f"Found no MSs in {this_target}")

tclean(vis=ms_files,
    imagename=str(f"{this_target}_co21_test_5kms"),
    imsize=512,
    cell='0.5arcsec',
    spw='',
    width='5km/s',
    nchan=20,
    start='-100km/s',
    specmode='cube',
    gridder='standard',
    deconvolver='multiscale',
    scales=[0, 3, 5, 10],
    niter=5000,
    cyclefactor=5.0,
    nsigma=3.0,
    perchanweightdensity=False,
    mosweight=False,
    pbcor=False,
    pblimit=0.15,
    restfreq='230.538GHz',
    phasecenter='',
    outframe='LSRK',
    )

impbcor(str(f"{this_target}_co21_test_5kms.image"),
        pbimage=str(f"{this_target}_co21_test_5kms.pb"),
        outfile=str(f"{this_target}_co21_test_5kms.image.pbcor"))
