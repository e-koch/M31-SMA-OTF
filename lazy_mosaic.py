from pathlib import Path
from astropy.io import fits
from reproject import reproject_interp
from reproject.mosaicking import find_optimal_celestial_wcs, reproject_and_coadd
import numpy as np
from spectral_cube import SpectralCube

from astropy.stats import mad_std


from tqdm import tqdm

# Directory containing subdirectories with FITS files
base_dir = Path("/reduction/erickoch/M31/sma/per_mosaic")

# Find all FITS files in subdirectories (e.g. subdir/file.fits)
# input_files = list(base_dir.glob("*/*.pbcor.tpeak"))
all_filenames = list(base_dir.glob("M31*/*cube*.pbcor"))

if not all_filenames:
    raise FileNotFoundError(f"No tpeak files found under {base_dir}")

print(f"Mosaicking these files: {[this_file.name for this_file in all_filenames]}")

# Collect headers and data shapes for all input FITS files
# input_projs = [SpectralCube.read(str(p), format='casa')[0] for p in all_filenames]
# headers = [h.header for h in input_projs]
# shapes = [h.data.shape for h in input_projs]

# # Determine an optimal output WCS/header that covers all inputs
# # The function returns (output_header, footprint) but we only need header here
# out_wcs, out_shape = find_optimal_celestial_wcs([(shape, header) for header, shape in zip(headers, shapes)])

# # Reproject and coadd. reproject_and_coadd will handle reprojection and weighting.
# # It expects an iterable of (data, header) tuples. We convert hdus to arrays.
# array_header_iter = [(proj.unitless_filled_data[:], proj.header) for proj in input_projs]

# # reproject_and_coadd returns (array, footprint) where array is the mosaicked data
# mosaic, out_footprint = reproject_and_coadd(array_header_iter, out_wcs, shape_out=out_shape,
#                                             reproject_function=reproject_interp)

# # Save the mosaic to a new FITS file
# mosaic_hdu = fits.PrimaryHDU(mosaic, out_wcs.to_header())
# mosaic_hdu.writeto(base_dir / "tpeak_mosaic.fits", overwrite=True)

ref_cube = SpectralCube.read(base_dir / "M31-Brick-D-Row-1-Col-1/M31-Brick-D-Row-1-Col-1_co21_test_cube.image.pbcor",
                             format='casa').with_spectral_unit(u.km/u.s, 'radio')


# Save out as FITS.
orig_cubes = []

all_fits_filenames = []
for filename in tqdm(all_filenames):
    this_cube = SpectralCube.read(filename, format='casa').with_spectral_unit(u.km/u.s, 'radio')
    fits_filename = str(filename) + ".fits"
    this_cube.write(fits_filename, overwrite=True)
    all_fits_filenames.append(fits_filename)

all_filenames = all_fits_filenames

orig_cubes = [SpectralCube.read(fits.open(filename)[0]) for filename in all_filenames]

cubes = orig_cubes


# Skipping this. All cubes should have the same spectral axis by construction

# cubes = []
# for this_cube in orig_cubes:
#     cubes.append(this_cube.spectral_interpolate(ref_cube.spectral_axis))

# for this_filename, this_cube in zip(all_filenames, cubes):
#     this_cube.write(this_filename, overwrite=True)


wcs_out, shape_out = find_optimal_celestial_wcs([cube[0].hdu for cube in cubes])

output_cube = np.zeros((ref_cube.shape[0],) + shape_out)

# Assume same weighting per channel
output_weight = np.zeros(shape_out)

# Only grab the first channel.
# weights = [SpectralCube.read(filename.replace('.image.pbcor', '.weight'), format='casa')[0]
#            for filename in all_filenames]

# weights = [SpectralCube.read(filename.replace('.image.pbcor.fits', '.weight'), format='casa')[0]
#            for filename in all_filenames]

pbs = [SpectralCube.read(filename.replace('.image.pbcor.fits', '.pb'), format='casa')[0]
           for filename in all_filenames]

# No CO emission expected at ~ -25 km/s. Use this for rms estimate:
this_chan = ref_cube.closest_spectral_channel(-25*u.km/u.s)


weights = []
for this_pb, this_cube in tqdm(zip(pbs, cubes)):

    # These are already pb corrected. Multiply by the pb first.
    this_rms = mad_std(this_cube[this_chan] * this_pb, ignore_nan=True)

    this_weight = this_pb / this_rms**2
    weights.append(this_weight)

for ii in tqdm(range(ref_cube.shape[0])):
    channel_reproj, channel_weight = reproject_and_coadd([cube[ii].hdu for cube in cubes], wcs_out,
                                                         shape_out=shape_out,
                                                         input_weights=[weight.hdu for weight in weights],
                                                         reproject_function=reproject_interp)
    output_cube[ii] = channel_reproj

    # Just keep a single channel of weights since we assume they're constant.
    if ii == 0:
        output_weight = channel_weight

output_hdr = wcs_out.to_header()
for hdr_key in ref_cube.header:
    if "3" in hdr_key:
        output_hdr[hdr_key] = ref_cube.header[hdr_key]

cube_hdu = fits.PrimaryHDU(output_cube, output_hdr)
cube_hdu.verify('fix')
cube_hdu.writeto('test_full_mosaic_cube.fits', overwrite=True)


cube_hdu = fits.PrimaryHDU(output_weight, output_hdr)
cube_hdu.verify('fix')
cube_hdu.writeto('test_full_mosaic_weights.fits', overwrite=True)
