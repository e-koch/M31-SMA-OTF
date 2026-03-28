
import numpy as np
import json

from pyuvdata.uvdata.mir_parser import MirParser

from astropy.coordinates import SkyCoord
from astropy.table import Table, Column, vstack
import astropy.units as u

import matplotlib.pyplot as plt


def split_by_map_and_row(mjd_s, tint=0.6, trow_diff=3.0):


    time_diffs = np.unique(np.round(np.diff(mjd_s), 1))
    # print(f"Time diffs (s): {time_diffs}")

    if time_diffs.size == 0:
        raise ValueError("No integrations returned.")


    map_idx = [0]
    row_idx = [0]

    row_id = 0
    map_id = 0

    for ii in range(len(mjd_s)-1):

        diff_s = mjd_s[ii+1] - mjd_s[ii]

        if diff_s > tint and diff_s <= trow_diff * 2:
            row_id += 1

        if diff_s > 2 * trow_diff:
            row_id = 0
            map_id += 1

        map_idx.append(map_id)
        row_idx.append(row_id)


    map_idx = np.array(map_idx)
    row_idx = np.array(row_idx)

    # Add start times by map
    start_times = np.array([mjd_s[map_idx == ii][0] for ii in np.unique(map_idx)])

    return map_idx, row_idx, start_times



def get_observation_segments(souids, mjd_s, source_names):
    """
    Given numpy arrays souids and mjd_s, return a list of (start, stop, souid) tuples
    for each contiguous block of the same souid.
    """
    if len(souids) == 0 or len(mjd_s) == 0:
        return []

    # Find indices where souid changes
    change_indices = np.where(souids[:-1] != souids[1:])[0] + 1
    # Add start and end indices
    segment_starts = np.concatenate(([0], change_indices))
    segment_ends = np.concatenate((change_indices, [len(souids)]))

    segments = []
    for start, end in zip(segment_starts, segment_ends):
        duration = mjd_s[end-1] - mjd_s[start]
        segments.append((source_names[souids[start]], souids[start], mjd_s[start], mjd_s[end-1], duration))

    return segments

if __name__ == "__main__":

    import argparse
    from pathlib import Path

    scripts_path = '/mnt/COMPASS9/sma/scripts/%s.json'


    parser = argparse.ArgumentParser(description="Extract spatial coverage from MIR file or JSON config.")
    parser.add_argument('--mir', type=str, help='Path to MIR file')
    parser.add_argument('--obsid', type=int, help='ObsID to read in JSON config')
    parser.add_argument('--plot', action='store_true',
                        help='Enable plot creation (overrides JSON)')
    args = parser.parse_args()

    # Determine input source and plot flag
    if args.obsid:
        json_file = scripts_path % args.obsid
        with open(json_file, 'r') as f:
            config = json.load(f)
        all_data_names= config["dataDirList"]
        mir_filename_list = [Path(filename) for filename in all_data_names]
    elif args.mir:
        mir_filename_list = [Path(args.mir)]
    else:
        raise ValueError('Either --mir or --json must be provided.')

    do_plots = args.plot

    for mir_filename in mir_filename_list:

        print(f"Reading {mir_filename}")
        mir_data = MirParser(mir_filename)

        print("Selecting data")
        mir_data.select([("rinteg","lt", 2), ("flags", "eq", 0)])
        offx, offy, rar, decr, souids, mjd = mir_data.in_data[["offx","offy","rar","decr","souid","mjd"]]

        mjd_s = mjd * 24 * 3600

        # Dictionary mapping source names to souids
        source_names = mir_data.codes_data['source']

        # Group by common rar values, then split and calculate with the offsets.
        unique_rar, unique_indices_rar = np.unique(rar, return_index=True)
        unique_decr, unique_indices_decr = np.unique(decr, return_index=True)

        # Sort indices to maintain order
        sorted_indices_rar = np.sort(unique_indices_rar)
        sorted_indices_decr = np.sort(unique_indices_decr)

        # Use sorted indices to get the unique values in the original order
        unique_rar_ordered = rar[sorted_indices_rar]
        unique_decr_ordered = decr[sorted_indices_decr]

        unique_souids = souids[sorted_indices_rar]
        unique_source_names = []
        for this_souid in unique_souids:
            unique_source_names.append(source_names[this_souid])

        assert unique_rar_ordered.size == unique_souids.size

        coord_target = SkyCoord(unique_rar_ordered * u.rad,
                                unique_decr_ordered * u.rad)

        if do_plots:
            plt.close('all')
            plt.figure()

        coords_target_otf_tabs = []
        for this_rar, this_decr, this_coord, this_souid in zip(unique_rar_ordered, unique_decr_ordered, coord_target, unique_souids):

            this_source = source_names[this_souid]

            print(f"Processing {this_coord.to_string('hmsdms')}")
            print(f"This is source {this_souid} with name {this_source}")

            # Identify the offsets for this source.
            this_mask = (rar == this_rar) & (decr == this_decr)
            print(f"Found {np.sum(this_mask)} offsets")

            this_offx = offx[this_mask]
            this_offy = offy[this_mask]

            coords_target_otf = this_coord.spherical_offsets_by(this_offx * u.arcsec, this_offy * u.arcsec)


            coords_target_otf_tabs.append(Table([Column(coords_target_otf.ra, name='ra'),
                                                Column(coords_target_otf.dec, name='dec'),
                                                Column([this_source] * coords_target_otf.ra.size, name='source')]))

            if do_plots:
                if "M31" not in this_source:
                    print(f"Only plotting M31 targets. This is {this_source}. Skipping")
                    continue
                plt.scatter(coords_target_otf.ra.value, coords_target_otf.dec.value, label=this_source)

        if do_plots:
            plt.legend()
            plt.savefig(f"{mir_filename.name}_target_otf_coords.png", bbox_inches='tight')
            plt.close()


        # Combine into a single table
        combined_table = vstack(coords_target_otf_tabs)
        combined_table.write(f"{mir_filename.name}_target_otf_coords.fits", overwrite=True)

        # Make per target hexbin heatmaps to assess uniformity of the coverage.
        if do_plots:

            for this_source in np.unique(combined_table['source']):
                this_mask = combined_table['source'] == this_source
                this_coords = combined_table[this_mask]

                plt.close('all')
                plt.figure()
                # Gridsize is number of hexagons along x axis. Set to 15 so each will be ~0.5 HPBW
                # for the 14 x 9' SMA M31 mosaics
                plt.hexbin(this_coords['ra'], this_coords['dec'], gridsize=15, cmap='inferno')
                plt.gca().invert_xaxis()
                plt.axis('equal')

                plt.savefig(f"{mir_filename.name}_{this_source}_hexbin.png", bbox_inches='tight')
                plt.close()

            # And then a version with them all!
            plt.hexbin(combined_table['ra'], combined_table['dec'], gridsize=15*3, cmap='inferno')
            plt.gca().invert_xaxis()
            plt.axis('equal')

            plt.savefig(f"{mir_filename.name}_all_hexbin.png", bbox_inches='tight')
            plt.close()


        # Save the unique rar and decr values. These should be the inputs for the observe cmd
        field_centers = Table([Column(coord_target.ra, name='ra'),
                            Column(coord_target.dec, name='dec'),
                            Column(unique_source_names, name='source')])
        field_centers.write(f"{mir_filename.name}_field_centers.fits", overwrite=True)



