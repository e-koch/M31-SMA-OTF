
import numpy as np

from pathlib import Path
import sys

import os

from casatasks import tclean, rmtables, imstat


skip_existing = False

data_path = Path('/reduction/erickoch/M31/sma/per_mosaic')
os.chdir(str(data_path))


this_target = sys.argv[-1]

# Iterate through data for all existing individual mosaics.
# These are indicated by the folder names in the data_path folder
# As created when running `apply_compass_solutions.py`

all_target_paths = list(data_path.glob("M31*"))
all_targets = [target.name for target in all_target_paths]

print(f"Searching on input target: {this_target}")

if this_target == "all":
    # Loop over all.
    save_final = True

    pass

# Or check for imaging of a maps in a single Brick
elif "Brick" in this_target and not "Row" in this_target:
    all_target_paths = [this_target_path for this_target_path in all_target_paths
                        if this_target in this_target_path.name]
    save_final = True

else:
    this_target_path = data_path / this_target

    # Check this_target is in all_targets
    if this_target not in all_targets:
        raise ValueError(f"Target {this_target} not found in {all_targets}")
    else:
        all_target_paths = [this_target_path]

    save_final = False

rms_values = []

print("All targets to image:")
for this_target_path in all_target_paths:
    print(this_target_path.name)


for this_target_path in all_target_paths:

    this_target = this_target_path.name

    print(f"Test imaging on {this_target}.")

    this_target_ms = [str(filename) for filename in this_target_path.glob("*.ms")]

    print(f"Found {len(this_target_ms)}")

    if len(this_target_ms) == 0:
        print(f"Found no MSs for {this_target}. Skipping")
        continue

    # if (this_target_path / f"{this_target}_co21_test_channel").exists():
    #     if skip_existing:
    #         print(f"Existing image found. Skipping")
    #         continue
    #     else:
    rmtables(str(this_target_path / f"{this_target}_co21_test_channel.*"))


    try:
        tclean(vis=this_target_ms,
            imagename=str(this_target_path / f"{this_target}_co21_test_channel"),
            imsize=1024,
            cell='1.3arcsec',
            spw='6',
            width='5.0km/s',
            nchan=1,
            start='-550km/s',
            specmode='cube',
            gridder='mosaic',
            niter=0,
            perchanweightdensity=False,
            mosweight=False,
            pbcor=True,
            pblimit=0.2,
            restfreq='230.538GHz',
            )
    except RuntimeError:
        continue

    # Estimate the rms:
    # Most of the field
    this_pb_name = str(this_target_path / f"{this_target}_co21_test_channel.pb")
    output = imstat(imagename=str(this_target_path / f"{this_target}_co21_test_channel.image"),
                    mask=f'"{this_pb_name}">0.6')

    this_rms_value = [this_target, output['rms'][0]] # , this_target_ms]

    print(f"Found {this_rms_value}")

    # Save the noise per map
    output_name_target = this_target_path / f"{this_target}_co21_test_channel_rms.npy"
    if output_name_target.exists():
        output_name_target.unlink()

    np.save(output_name_target, np.array(this_rms_value))

    # Append to aggregated list over targets:
    rms_values.append(this_rms_value)

# Save rms over all targets run:
if save_final:
    output_name_all = data_path / "all_rms_values.npy"
    np.save(output_name_all, rms_values)
