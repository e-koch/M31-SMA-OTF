

import numpy as np

from pathlib import Path

from pyuvdata.uvdata.mir_parser import MirParser

from astropy.coordinates import SkyCoord
from astropy.table import Table, Column, vstack
import astropy.units as u

import matplotlib.pyplot as plt



def get_observation_segments(souids, mjd_s, source_names, expected_duration=1877 * u.s):
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
        duration = (mjd_s[end-1] - mjd_s[start]) * u.s
        # Execution fraction. Round down to nearest 0.5
        exec_frac = np.floor((duration / expected_duration).to(u.one).value * 2) / 2
        segments.append((source_names[souids[start]], souids[start], mjd_s[start], mjd_s[end-1], duration.value, exec_frac))

    return segments


executions = {'A': [], 'B': [], 'C': [], 'D': []}

# track ID, [mir filenames], nants

# executions['A'] = \
#     [
#     #  ['13757', 'MIRFILE'],  # swarm issues.
#      ['13763', ['250905_07:15:37'], 6],
#      ['13768', ['250908_07:09:38'], 6],
#      ['13773', ['250910_06:55:13'], 6],
#      ['13777', ['250912_06:45:28'], 5],
#      ['13783', ['250915_06:31:44'], 5],
#      ['13801', ['250923_06:00:47'], 6],
#      ['13806', ['250925_05:45:34'], 6],
#      ['13813', ['250929_05:39:24'], 5],
#      ]

# executions['B'] = \
#     [
#     #  ['13785', ['250916_06:50:33', '250916_14:30:05'], 4],
#      ['13815', ['250930_05:36:19'], 5],
#      ['13819', ['251002_05:22:58'], 6],
#      ['13821', ['251003_05:26:05'], 6],
#      ['13824', ['251004_05:11:22'], 6],
#      ['13826', ['251005_05:08:39'], 6],
#      ['13828', ['251006_05:06:47'], 6],
#      ['13831', ['251007_05:11:59'], 6],
#      ['13834', ['251008_05:05:55'], 6],
#     ]

# executions['C'] = \
#     [
#      ['13837', ['251009_04:51:37'], 6],
#      ['13839', ['251010_05:00:27'], 6],
#      ['13862', ['251022_06:50:06'], 6],
#      ['13865', ['251023_06:42:04'], 6],
#      ['13887', ['251102_03:28:23'], 6],
#      ['13895', ['251106_03:52:38'], 6],
#      ['13898', ['251107_04:11:52'], 6],
#      ['13900', ['251108_04:21:21'], 5],
#      ['13926', ['251122_03:07:34'], 6],
#     ]

executions['D'] = \
    [
     ['13901', ['251109_03:41:14'], 5],
     ['13902', ['251110_03:17:01'], 5],
     ['13904', ['251111_03:46:44'], 5],
     ['13907', ['251112_03:30:34'], 6],
     ['13914', ['251116_03:18:30'], 6],
     ['13921', ['251121_03:12:08'], 6],
     ['13927', ['251123_02:59:59'], 6],
     ['13931', ['251125_02:55:12'], 6],
     ['13935', ['251127_03:23:52'], 6],
     ['13936', ['251128_03:11:50'], 6],
     ['13938', ['251129_04:07:00'], 6],
     ['13939', ['251130_04:33:41'], 6],
     ['13947', ['251204_05:37:44'], 6],
     ['13949', ['251205_02:40:00'], 6],
     ['13951', ['251206_03:51:55'], 6],
    ]

# executions['E'] = \
#     [
#      ['13886', ['251101_03:50:59'], 6],
#      ['13928', ['251124_03:23:42'], 6],
#      ['13933', ['251126_03:15:54'], 6],
#     ]

execution_summary = {'A': [], 'B': [], 'C': [], 'D': [], 'E': []}


# Time to do both interleaved OTF maps for a single target
# Used to assess completed map.
expected_duration = 1800 * u.s

min_ants = 6

# For SUB
maps_per_track = 15
nmaps_expected = 22

# For COM
# maps_per_track = 19
# nmaps_expected = 29


total_map_anthr = min_ants * nmaps_expected * expected_duration


data_path = Path('/sma/data/science/mir_data/')

for this_brick in executions.keys():

    print(f"On brick {this_brick}.")

    brick_executions = executions[this_brick]

    for obsid, mir_filenames, num_ants in brick_executions:

        print(f"On obsid {obsid}.")

        for ii, mir_filename in enumerate(mir_filenames):

            mir_filename = data_path / mir_filename

            mir_data = MirParser(mir_filename)
            mir_data.select([("rinteg","lt", 2), ("flags", "eq", 0)])

            parsed_data = mir_data.in_data[["offx","offy","rar","decr","souid","mjd"]]
            this_offx, this_offy, this_rar, this_decr, this_souids, this_mjd = parsed_data

            this_mjd_s = this_mjd * 24 * 3600

            # Dictionary mapping source names to souids
            this_source_names = mir_data.codes_data['source']

            if ii == 0:
                offx = this_offx
                offy = this_offy
                rar = this_rar
                decr = this_decr
                souids = this_souids
                mjd_s = this_mjd_s
                source_names = this_source_names
            else:
                offx = np.append(offx, this_offx)
                offy = np.append(offy, this_offy)
                rar = np.append(rar, this_rar)
                decr = np.append(decr, this_decr)
                souids = np.append(souids, this_souids)
                mjd_s = np.append(mjd_s, this_mjd_s)

                # Dictionary union
                source_names = source_names | this_source_names

        segments = get_observation_segments(souids, mjd_s, source_names,
                                            expected_duration)


        execution_summary[this_brick].append([obsid] + segments)



# Per execution summary: number of maps per target, number of maps per track, and total number of maps per track
from collections import defaultdict

maps_per_target = defaultdict(lambda: defaultdict(int))  # {brick: {target: count}}
maps_per_track = defaultdict(int)
maps_per_brick = defaultdict(int)  # {brick: total count}

for brick, executions in execution_summary.items():
    for execution in executions:
        obsid = execution[0]
        segments = execution[1:]
        for seg in segments:
            target_name = seg[0]
            maps_per_target[brick][target_name] += seg[-1]
            maps_per_brick[brick] += seg[-1]
            maps_per_track[obsid] += seg[-1]

# Print summary
print("\nSummary of maps per target (per brick):")
for brick in maps_per_target:
    print(f"Brick {brick}:")
    for target, count in maps_per_target[brick].items():
        print(f"  Target {target}: {count} maps out of {nmaps_expected}. Completion fraction: {count / nmaps_expected}")
    print(f"  Total maps in brick {brick}: {maps_per_brick[brick]}")

for obsid in maps_per_track:
    print(f"Execution {obsid}: {maps_per_track[obsid]} maps per track")
