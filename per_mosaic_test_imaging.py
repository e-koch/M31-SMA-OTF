
import numpy as np

from pathlib import Path
import sys

import os

from casatasks import tclean, rmtables, imstat


data_path = Path('/reduction/erickoch/M31/sma/per_mosaic')

this_target = sys.argv[-1]

if this_target == "all":
    this_target = None

skip_existing = False

os.chdir(str(data_path))

rms_values = []

# Iterate through data for all existing individual mosaics.
# These are indicated by the folder names in the data_path folder
# As created when running `apply_compass_solutions.py`

all_target_paths = list(data_path.glob("M31*"))
all_targets = [target.name for target in all_target_paths]

if this_target is not None:

    this_target_path = data_path / this_target

    # Check this_target is in all_targets
    if this_target not in all_targets:
        raise ValueError(f"Target {this_target} not found in {all_targets}")
    else:
        all_target_paths = [this_target_path]

print("All targets to image:")
for this_target_path in all_target_paths:
    print(this_target_path.name)

for this_target_path in all_target_paths:

    this_target = this_target_path.name

    print(f"Test imaging on {this_target}.")

    this_target_ms = [str(filename) for filename in this_target_path.glob("*.ms")]

    print(f"Found {len(this_target_ms)} MSs to image.")

    # Check each MS is complete:
    bad_ms = []
    for this_ms in this_target_ms:

        try:
            if ms.open(this_ms, check=True):
                ms.close()
        except Exception:
            print(f"Issue with {this_ms}. Skipping.")
            bad_ms.append(this_ms)
    # Remove any bad MSs
    this_target_ms = list(set(this_target_ms) - set(bad_ms))
    print(f"Imaging with {len(this_target_ms)} MSs.")


    if len(this_target_ms) == 0:
        print(f"Found no MSs for {this_target}. Skipping")
        continue



    if (this_target_path / f"{this_target}_co21_test_cube.image").exists():
        if skip_existing:
            print(f"Existing image found. Skipping")
            continue
        else:
            rmtables(str(this_target_path / f"{this_target}_co21_test_cube.*"))


    tclean(vis=this_target_ms,
           imagename=str(this_target_path / f"{this_target}_co21_test_cube"),
           imsize=1024 if "Brick-E" not in this_target else 1600,
           cell='1.3arcsec',
           spw='6',
           width='5.0km/s',
           nchan=120,
           start='-600km/s',
           specmode='cube',
           gridder='mosaic',
           niter=0,
           perchanweightdensity=False,
           mosweight=False,
           pbcor=False,
           pblimit=0.2,
           restfreq='230.538GHz',
           )

    impbcor(str(this_target_path / f"{this_target}_co21_test_cube.image"),
            pbimage=str(this_target_path / f"{this_target}_co21_test_cube.pb"),
            outfile=str(this_target_path / f"{this_target}_co21_test_cube.image.pbcor"))

    immoments(str(this_target_path / f"{this_target}_co21_test_cube.image.pbcor"),
              moments=[8],
              outfile=str(this_target_path / f"{this_target}_co21_test_cube.image.pbcor.tpeak"))

