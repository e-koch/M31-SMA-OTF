
from casatasks import tclean, split, concat, rmtables

import os
from glob import glob

os.chdir('/reduction/erickoch/M31/sma/per_mosaic/NGC7293')

final_ms = 'NGC7293_spw6.ms'

if not os.path.exists(final_ms):
    # Split out SPW 6 then concat

    ms_files = glob('*_NGC7293.ms')

    print(ms_files)

    split_ms_files = []


    for ms_file in ms_files:
        print(ms_file)
        split(vis=ms_file,
            outputvis=ms_file.replace('.ms', '_spw6.ms'),
            spw='6', datacolumn='data', keepflags=False)

        split_ms_files.append(ms_file.replace('.ms', '_spw6.ms'))

    concat(vis=split_ms_files,
        concatvis='NGC7293_spw6.ms')

rmtables('ngc7293_co21_test_1p5kms.*')

tclean(vis='NGC7293_spw6.ms',
    imagename=str("ngc7293_co21_test_1p5kms"),
    imsize=1024,
    cell='1.3arcsec',
    spw='0',
    width='1.5km/s',
    nchan=50,
    start='-60km/s',
    specmode='cube',
    gridder='mosaic',
    niter=0,
    perchanweightdensity=False,
    mosweight=False,
    pbcor=False,
    pblimit=0.15,
    restfreq='230.538GHz',
    phasecenter='ICRS 22:29:38.545 -20.50.13.75',
    outframe='LSRK',
    )

