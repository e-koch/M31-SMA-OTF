import numpy as np
import matplotlib.pyplot as plt

from spectral_cube import SpectralCube
from radio_beam import Beam
from astropy import units as u

from astropy.convolution import Gaussian1DKernel

from uvcombine import feather_compare
from uvcombine.scale_factor import find_scale_factor

from cube_analysis.feather_cubes import feather_compare_cube
from cube_analysis.register_cubes import cube_registration

from uvcombine import feather_simple_cube


import seaborn as sns
from pathlib import Path
from copy import deepcopy

smt_data_path = Path("//Volumes/Expansion/storage/M31/SMT/")
sma_data_path = Path("/Volumes/Expansion/storage/M31/SMA/")

def taper_weights(pb_plane,
                erosion_interations=5):

    # mask = np.logical_or(np.isfinite(pb_plane), pb_plane > 0.)
    mask = pb_plane > 0.

    if erosion_interations > 0:
        mask = nd.binary_erosion(mask, iterations=erosion_interations)

    smoothed_weights = nd.gaussian_filter(mask.astype(float), sigma=erosion_interations/2.)

    weight_arr = pb_plane * smoothed_weights

    return weight_arr


smt_cube = SpectralCube.read(smt_data_path / "pilot2+3_map_cube.fits")
smt_cube = smt_cube.with_spectral_unit(u.km/u.s, 'radio')
smt_cube._unit = u.K

# SMT cube is gridded per pixel with a top-hat.
orig_beam = Beam(31 * u.arcsec)
pix_scale = (smt_cube.wcs.wcs.cdelt[1] * u.deg).to(u.arcsec)
eff_beam = Beam((pix_scale**2 + orig_beam.major**2)**0.5)

print("Effective beam FWHM:", eff_beam.major)
#

smt_cube._beam = Beam(0*u.arcsec)
smt_cube = smt_cube.convolve_to(orig_beam,
                                preserve_nan=True)


# Two SMA cubes that overlap the SMT pilot.
sma_cube_row1_col8 = SpectralCube.read(sma_data_path / "M31-Brick-A-Row-1-Col-8_co21_deconv_cube.image.pbcor",
                                       format='casa')

sma_cube_row2_col8 = SpectralCube.read(sma_data_path / "M31-Brick-A-Row-2-Col-8_co21_deconv_cube.image.pbcor",
                                       format='casa')

sma_cube_row1_col8 = sma_cube_row1_col8.with_spectral_unit(u.km/u.s, 'radio')
sma_cube_row2_col8 = sma_cube_row2_col8.with_spectral_unit(u.km/u.s, 'radio')

# Convert to K
sma_cube_row1_col8 = sma_cube_row1_col8.to(u.K)
sma_cube_row2_col8 = sma_cube_row2_col8.to(u.K)

# Spectral smooth and downsample the SMT cube.

# SMA spectral axes should match. Interpolate the SMT cube to it.
fwhm_factor = np.sqrt(8*np.log(2))
smt_velchan = np.abs(np.diff(smt_cube.spectral_axis)[0])
sma_velchan = np.abs(np.diff(sma_cube_row1_col8.spectral_axis)[0])
gaussian_width = ((sma_velchan**2 - smt_velchan**2)**0.5 / smt_velchan / fwhm_factor)
kernel = Gaussian1DKernel(gaussian_width.value)

smt_cube_specsmooth = smt_cube.spectral_smooth(kernel)

smt_cube_specinterp = smt_cube_specsmooth.spectral_interpolate(sma_cube_row1_col8.spectral_axis)

# smt_cube_specinterp.write(smt_data_path / "pilot2+3_map_cube_sma_specinterp.fits", overwrite=True)
smt_cube_specinterp.write(smt_data_path / "pilot2+3_map_cube_sma_specinterp_spatialsmooth.fits", overwrite=True)


##########
# Plot mesn spectra:
##########

# sma_cube_row1_col8.mean(axis=(1, 2)).quicklook()
# smt_cube_specinterp.mean(axis=(1, 2)).quicklook()



##########
# Registration checks:
##########

# spec_slicer = slice(None)


# spatial_offsets = cube_registration(smt_cube_specinterp[spec_slicer],
#                                     sma_cube_row1_col8[spec_slicer],
#                                     verbose=True, num_cores=1,
#                                     restfreq=230.538 * u.GHz,
#                                     check_specaxis=True,)

# cmap = sns.color_palette("icefire", as_cmap=True)

# fig, ax = plt.subplots()
# points = ax.scatter(spatial_offsets[:, 1],
#                     spatial_offsets[:, 0],
#                     c=smt_cube_specinterp[spec_slicer].spectral_axis.to(u.km / u.s).value,
#                     s=50,
#                     cmap=cmap)
# fig.colorbar(points, label="Velocity (km/s)")

# ax.set_ylabel("Dec pixel offset")
# ax.set_xlabel("RA pixel offset")

# plt.close()

# fig.savefig(smt_data_path / f"{this_gal}_spatialoffsets.pdf")
# plt.close(fig)

# Match spatial areas:

sma_reproj_plane = sma_cube_row1_col8[0].reproject(smt_cube_specinterp[0].header)
smt_cube_specinterp_masked_1_8 = smt_cube_specinterp.with_mask(np.isfinite(sma_reproj_plane)).minimal_subcube(spatial_only=True)

# Hign S/N velocity range:
chan_min = sma_cube_row1_col8.closest_spectral_channel(-90 * u.km/u.s)
chan_max = sma_cube_row1_col8.closest_spectral_channel(-40 * u.km/u.s)

if chan_min > chan_max:
    chan_min, chan_max = chan_max, chan_min

print(f"Channels to slice: {chan_min}, {chan_max}")

spec_slicer = slice(chan_min, chan_max)


co_freq = 230.538 * u.GHz
min_baseline = 7 * u.m
las = (co_freq.to(u.cm, u.spectral()) / min_baseline).to(u.arcsec, u.dimensionless_angles())

# 5 iterations is sufficient to downweight noisy edges with low PB values
weights = taper_weights(np.isfinite(sma_cube_row1_col8[0]),
                        erosion_interations=5)

radii, ratios, high_pts, low_pts, chan_out = \
            feather_compare_cube(sma_cube_row1_col8,
                                smt_cube_specinterp_masked_1_8,
                                las,
                                lowresfwhm=None,
                                num_cores=1,
                                chunk=250,
                                verbose=True,
                                weights=weights,
                                relax_spectral_check=False,
                                spec_check_kwargs={'rtol': 0.01})


sc_factor, sc_err = find_scale_factor(np.hstack(low_pts[:-1]),
                                    np.hstack(high_pts[:-1]),
                                    method='distrib',
                                    verbose=True)
print(f"Scaling factor (all channels): {sc_factor:.2f}+/-{sc_err:.2f}")
plt.title(f"Scaling factor (all channels): {sc_factor:.2f}+/-{sc_err:.2f}")

plt.savefig(smt_data_path / f"smt_sma_ro1_col8_sdfactor_allchans.pdf")
plt.close()


sc_factor, sc_err = find_scale_factor(np.hstack(low_pts[spec_slicer]),
                                    np.hstack(high_pts[spec_slicer]),
                                    method='distrib',
                                    verbose=True)
print(f"Scaling factor (gal channels): {sc_factor:.2f}+/-{sc_err:.2f}")
plt.title(f"Scaling factor (gal channels): {sc_factor:.2f}+/-{sc_err:.2f}")
# plt.xlim([-1, 1])

plt.savefig(smt_data_path / f"smt_sma_ro1_col8_sdfactor_galchans.pdf")
plt.close()


# Now for Row 2 Col 8


sma_reproj_plane = sma_cube_row2_col8[0].reproject(smt_cube_specinterp[0].header)
smt_cube_specinterp_masked_2_8 = smt_cube_specinterp.with_mask(np.isfinite(sma_reproj_plane)).minimal_subcube(spatial_only=True)


weights = taper_weights(np.isfinite(sma_cube_row2_col8[0]),
                        erosion_interations=5)

radii, ratios, high_pts, low_pts, chan_out = \
            feather_compare_cube(sma_cube_row2_col8,
                                smt_cube_specinterp_masked_2_8,
                                las,
                                lowresfwhm=None,
                                num_cores=1,
                                chunk=250,
                                verbose=True,
                                weights=weights,
                                relax_spectral_check=False,
                                spec_check_kwargs={'rtol': 0.01})


sc_factor, sc_err = find_scale_factor(np.hstack(low_pts[:-1]),
                                    np.hstack(high_pts[:-1]),
                                    method='distrib',
                                    verbose=True)
print(f"Scaling factor (all channels): {sc_factor:.2f}+/-{sc_err:.2f}")
plt.title(f"Scaling factor (all channels): {sc_factor:.2f}+/-{sc_err:.2f}")

plt.savefig(smt_data_path / f"smt_sma_row2_col8_sdfactor_allchans.pdf")
plt.close()


sc_factor, sc_err = find_scale_factor(np.hstack(low_pts[spec_slicer]),
                                    np.hstack(high_pts[spec_slicer]),
                                    method='distrib',
                                    verbose=True)
print(f"Scaling factor (gal channels): {sc_factor:.2f}+/-{sc_err:.2f}")
plt.title(f"Scaling factor (gal channels): {sc_factor:.2f}+/-{sc_err:.2f}")
# plt.xlim([-1, 1])

plt.savefig(smt_data_path / f"smt_sma_row2_col8_sdfactor_galchans.pdf")
plt.close()


# Noise only ratio comparison


sc_factor, sc_err = find_scale_factor(np.hstack(low_pts[:5]),
                                    np.hstack(high_pts[:5]),
                                    method='distrib',
                                    verbose=True)
print(f"(Noise only) Scaling factor (gal channels): {sc_factor:.2f}+/-{sc_err:.2f}")
plt.title(f"(Noise only) Scaling factor (gal channels): {sc_factor:.2f}+/-{sc_err:.2f}")
# plt.xlim([-1, 1])

plt.savefig(smt_data_path / f"smt_sma_row2_col8_sdfactor_noisechans.pdf")
plt.close()



smt_channel = smt_cube_specinterp_masked_1_8[smt_cube_specinterp_masked_1_8.closest_spectral_channel(-60 * u.km/u.s)]
sma_channel = sma_cube_row1_col8[sma_cube_row1_col8.closest_spectral_channel(-60 * u.km/u.s)]

# Test smooth:
smt_channel._beam = Beam(0*u.arcsec)
smt_channel_smooth = smt_channel.convolve_to(smt_cube.beam,
                                             preserve_nan=True)

smt_channel_reproj = smt_channel_smooth.reproject(sma_channel.header, order='bilinear')

# smt_cube_specinterp_masked_1_8_reproj = smt_cube_specinterp_masked_1_8.reproject(sma_cube_row1_col8.header)

plt.figure()
plt.subplot(121)
smt_channel_reproj.quicklook()
# smt_channel_smooth.quicklook()
plt.subplot(122)
sma_channel.quicklook()

plt.close()

# Is the SMT convolved with the correct beam? The noise looks like it has no spatial structure
from turbustat.statistics import PowerSpectrum

plt.close()

pspec_smt = PowerSpectrum(smt_channel_reproj.hdu)
pspec_smt.run(verbose=True, apodize_kernel='tukey', fit_2D=False)
plt.close()



# Smooth and reproject the SMA cubes to SMT:

# sma_cube_row2_col8_conv = sma_cube_row2_col8.convolve_to(smt_cube.beam)
# sma_cube_row2_col8_reproj = sma_cube_row2_col8_conv.reproject(smt_cube.header)


##########
# Naive feathering with no flux correction:
##########


feathered_cube = feather_simple_cube(sma_cube_row2_col8,
                                smt_cube_specinterp,
                                allow_lo_reproj=True,
                                allow_spectral_resample=False,
                                lowresscalefactor=1.0)

# NaN out blank areas post-FFT.
feathered_cube = feathered_cube.with_mask(sma_cube_row2_col8.mask)

# this_feathered_filename = vla_data_path / f"{this_vla_filename[:-5]}_feathered.fits"
this_feathered_filename = sma_data_path / f"M31-Brick-A-Row-2-Col-8_co21_deconv_cube.image.pbcor_feathered.K.fits"

feathered_cube.write(this_feathered_filename, overwrite=True)


feathered_cube = feather_simple_cube(sma_cube_row1_col8,
                                smt_cube_specinterp,
                                allow_lo_reproj=True,
                                allow_spectral_resample=False,
                                lowresscalefactor=1.0)

# NaN out blank areas post-FFT.
feathered_cube = feathered_cube.with_mask(sma_cube_row1_col8.mask)

# this_feathered_filename = vla_data_path / f"{this_vla_filename[:-5]}_feathered.fits"
this_feathered_filename = sma_data_path / f"M31-Brick-A-Row-1-Col-8_co21_deconv_cube.image.pbcor_feathered.K.fits"

feathered_cube.write(this_feathered_filename, overwrite=True)
