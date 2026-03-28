
import os
import sys
import json
from pathlib import Path
import numpy as np

from astropy.time import Time
from astropy.table import Table

import matplotlib.pyplot as plt

osjoin = os.path.join


scripts_path = '/mnt/COMPASS9/sma/scripts/%s.json'
flag_path = '/home/erickoch/M31_SMA/manual_flags/%s_flags.txt'

# output_path = '/reduction11/erickoch/M31/sma/per_mosaic'
output_path = '/reduction11/erickoch/M31/sma_12CO/per_mosaic'
# output_path = '/reduction11/erickoch/IC342/sma/per_mosaic'

print("All args:", sys.argv)


script_arg_pos = np.where(["apply_compass_solutions.py" in val for val in np.array(sys.argv)])[0][0]
print("Script arg pos:", script_arg_pos)
script_args = sys.argv[script_arg_pos+1:]

print("Script args:", script_args)

obs_id = script_args[0]

# Otherwise specify the line to extract. This is to create
# the correct SPW mapping from the MS (containing a subset of SPW) and
# the full list written out in the calibration tables
if len(script_args) > 1:
    line_selec = script_args[1]

    if line_selec == '12CO':
        # 1 SPW expected: upper SB, chunk 1.
        # Maps to SPW 6 in the full list.
        spw_map = [6]
    elif line_selec == '13CO':
        # 2 SPWs expected: lower SB chunks 1 and 2.
        spw_map = [4, 5]
    else:
        raise ValueError("Current line selection must be 12CO or 13CO.")

else:
    line_selec = None
    spw_map = []


# Find the MS:
this_path = Path(".")

search_str = "_bin*.ms"

ms_name = None
for file in this_path.glob(f"*{search_str}"):
    ms_name = file
    print(f"Found {ms_name}.")
    break

if ms_name is None:
    print("No MS found.")
    sys.exit(1)

# ms_name = Path(sys.argv[-2])

data_path = ms_name.parent

json_file = scripts_path % obs_id

with open(json_file, 'r') as file_handle:
    json_dict = json.load(file_handle)


tb.open(f"{ms_name}/SPECTRAL_WINDOW")
nspws = tb.getcol('NUM_CHAN').size
tb.close()

# If there's 14 spw not 12, we need to drop the pseudo-cont outputs.
if nspws == 14:
    split(vis=str(ms_name),
          outputvis=f"{ms_name}.temp",
          datacolumn='data',
          spw='1~12',)

    rmtables(str(ms_name))
    os.system(f"mv {ms_name}.temp {ms_name}")


# Otherwise check for SPW sub-sets when extracting a single line
if line_selec == '12CO':
    if nspws != 1:
        raise ValueError("Only 1 SPW expected for 12CO.")
elif line_selec == '13CO':
    if nspws != 2:
        raise ValueError("Only 2 SPWs expected for 13CO.")
else:
    pass


# Grab the target and cal names:
sciTargs = [str(val) for val in json_dict["sciTargs"]]
gainCals = [str(val) for val in json_dict["gainCals"]]
fluxCals = [str(val) for val in json_dict["fluxCals"]]
bpCals = [str(val) for val in json_dict["bpCals"]]

# Check the names in the ms:
tb.open(f"{ms_name}/FIELD")
names = [str(val) for val in np.unique(tb.getcol("NAME"))]
tb.close()

name_dict = dict.fromkeys(sciTargs+gainCals+fluxCals+bpCals)
for name in name_dict.keys():

    if name in names:
        name_dict[name] = name
        continue

    # Check case.
    name_cape = name.capitalize()
    if name_cape in names:
        name_dict[name] = name_cape
        continue

    # Check lower case.
    name_lower = name.lower()
    if name_lower in names:
        name_dict[name] = name_lower
        continue

    # Check upper case.
    name_upper = name.upper()
    if name_upper in names:
        name_dict[name] = name_upper
        continue

    # raise ValueError(f"Name {name} not found in {ms_name}")


# Pre-apply manual flags if the file exists

if os.path.exists(flag_path % obs_id):

    flagmanager(vis=str(ms_name),
                mode='save',
                versionname='initial_flags')

    flagdata(vis=str(ms_name),
             mode='list',
             inpfile=flag_path % obs_id)

    flagmanager(vis=str(ms_name),
                mode='save',
                versionname='manual_flagging_file')

# Apply the self-cal solutions for the bp and flux cals:
for this_cal in bpCals:  #  + fluxCals:

    selfcal_table = data_path / f"{obs_id}_{this_cal}_selfcal_solns.ms"

    this_name = name_dict[this_cal]

    applycal(vis=str(ms_name),
            field=str(this_name),
            spw="",
            gaintable=[str(selfcal_table)],
            interp=['nearest,linear'],
            spwmap=[spw_map],
            gainfield=[str(this_name)],
            flagbackup=False,
            calwt=False)

# Apply the general amp/phase solutions from the gain cals:
amp_table = data_path / f"{obs_id}_amp_solns.ms"
pha_table = data_path / f"{obs_id}_pha_solns.ms"

gain_names = []
for this_cal in gainCals:

    this_name = name_dict[this_cal]

    applycal(vis=str(ms_name),
            field=this_name,
            spw="",
            gaintable=[str(amp_table), str(pha_table)],
            interp=['nearest,linear', 'nearest,linear'],
            spwmap=[spw_map, spw_map],
            gainfield=[this_name, this_name],
            flagbackup=False,
            calwt=False)

    gain_names.append(this_name)

# Apply to all science targets:

gainCal_names = [str(name_dict[name]) for name in gainCals]

sci_names = [str(name_dict[name]) for name in sciTargs]

applycal(vis=str(ms_name),
        field=",".join(sci_names),
        spw="",
        gaintable=[str(amp_table),
                   str(pha_table)],
        interp=['linearPD,linear',
                'linear,linear'],
        spwmap=[spw_map,
                spw_map],
        gainfield=[",".join(gain_names),
                   ",".join(gain_names)],
        flagbackup=False,
        calwt=False,
        applymode='calflagstrict')

flagmanager(vis=str(ms_name),
            mode='save',
            versionname='applycal',
            comment='applycal to science targets'
            )



# Early OTF tracks ran into issues with datacatcher and otf not failing gracefully
# this led to some cases where the data source is mislabeled (e.g., science target labeling
# is actually a gain cal)
# Quick fix for this is to flag average amplitudes that are much higher than expected

### visstat is incredibly slow with our data for some reason.

# Export time series with visstat
# Returns times in MJD seconds
# test = visstat(vis=str(ms_name),
#                field=",".join(sci_names),
#                datacolumn='corrected',
#                reportingaxes='integration',
#                useflags=True,
#                timeaverage=True,
#                timebin='10s',)

# amp = np.array([test[key]['median'] for key in test.keys()])
# amp_rms = np.array([test[key]['medabsdevmed'] for key in test.keys()])

# times_mjd = [float(key.split("TIME=")[-1]) for key in test.keys()]


# Using plotms exporting a txt file instead
output_name_amp = data_path / f"{ms_name.name}_amp_vs_time.txt"
if output_name_amp.exists():
    output_name_amp.unlink()

plotms(vis=str(ms_name),
       xaxis='time',
       yaxis='amp',
    ydatacolumn='corrected',
    selectdata=True,
    field=",".join(sci_names),
    avgchannel=str(1e8),
    avgtime='',
    correlation="",
    averagedata=True,
    avgbaseline=True,
    avgspw=True,
    transform=False,
    extendflag=False,
    plotfile=str(output_name_amp),
    overwrite=True,
    showgui=False)

output_name_wts = data_path / f"{ms_name.name}_wt_vs_time.txt"
if output_name_wts.exists():
    output_name_wts.unlink()

msmd.open(str(ms_name))
nchan = msmd.nchan(0) # Grab 0th SPW by default
msmd.close()

plotms(vis=str(ms_name),
       xaxis='time',
       yaxis='wt',
    ydatacolumn='corrected',
    selectdata=True,
    field=",".join(sci_names),
    spw=f'0:{nchan//4}',
    avgchannel=str(1e8),
    avgtime='',
    correlation="",
    averagedata=True,
    avgbaseline=False,
    avgspw=False,
    transform=False,
    extendflag=False,
    plotfile=str(output_name_wts),
    overwrite=True,
    showgui=False)


def read_casa_txt(filename):

    # Grab the meta-data from the header
    meta_lines = skim_header_metadata(filename)

    # Grab the meta-data from the header
    meta_dict = make_meta_dict(meta_lines)

    # After the plot 0 line
    header_start = len(meta_lines) + 1
    # One for column names, another for units.
    data_start = len(meta_lines) + 3

    try:
        tab = Table.read(filename,
                        format='ascii.commented_header',
                        header_start=header_start,
                        data_start=data_start)
    except Exception as e:
        print(f"Failured reading {filename} with exception {e}.")
        tab = Table()

    return tab, meta_dict


def skim_header_metadata(filename):
    '''
    Search for "From plot 0"
    '''
    search_str = "# From plot 0"

    # Should be close to ~10 or below, I think
    # This just stops reading too far if something
    # goes wrong.
    max_line = 50

    meta_lines = []

    with open(filename, 'r') as f:

        for i, line in enumerate(f):
            if search_str in line:
                break

            meta_lines.append(line)

            if i > max_line:
                raise ValueError(f"Could not find header in {filename}")

    return meta_lines


def make_meta_dict(meta_lines):
    '''
    Convert the meta lines into something nice.
    '''

    data_dict = {}

    for line in meta_lines:

        # Skip "# "
        line = line[2:]

        # Some plotms output will have multiple name:value pairs
        num_names = len(line.split(": ")) // 2

        for ii in range(num_names):

            name, value = line.split(": ")[2*ii:2*(ii)+2]

            name = name.strip(" ")
            value = value.strip(" ")
            value = value.strip("\n")

            data_dict[name] = value

    # CASA 6.6 uses 'file' instead of 'vis'.
    if 'file' in data_dict:
        data_dict['vis'] = data_dict['file']
        del data_dict['file']

    return data_dict



def find_consecutive_ranges(mask,
                            min_length=1,
                            max_gap=10,
                            pad_size=2):

    # Find the start and end points of True regions
    diff = np.diff(np.concatenate(([False], mask, [False])).astype(int))
    start_indices = np.where(diff == 1)[0]
    end_indices = np.where(diff == -1)[0] - 1

    if len(start_indices) == 0:
        return []

    # Join regions within max_gap
    if max_gap > 0:
        gaps = start_indices[1:] - end_indices[:-1] - 1
        merge_mask = gaps <= max_gap
        merge_indices = np.where(merge_mask)[0]

        for idx in merge_indices[::-1]:
            start_indices = np.delete(start_indices, idx + 1)
            end_indices = np.delete(end_indices, idx)

    # Apply padding
    if pad_size > 0:
        start_indices = np.maximum(start_indices - pad_size, 0)
        end_indices = np.minimum(end_indices + pad_size, len(mask) - 1)

    # Calculate lengths of regions
    lengths = end_indices - start_indices + 1

    # Filter based on min_length
    valid = lengths >= min_length
    start_indices = start_indices[valid]
    end_indices = end_indices[valid]

    # Return as list of tuples
    return list(zip(start_indices, end_indices))


# MAD approach to identify outliers
def find_outliers(data, sigma=5.):
    d = np.abs(data - np.nanmedian(data))
    mdev = np.nanmedian(d)
    s = d/mdev
    return s >= sigma


# Preserve initial flag version
flagmanager(vis=str(ms_name),
            mode='save',
            versionname='initial_flags_preoutlier',
            comment='Initial flags')

apply_flagging = True


# Loop through outlier checks in amplitude and weight:
for output_name, outlier_type in zip([output_name_wts, output_name_amp],
                       ['amp', 'wt']):

    print(f"Running outlier checks on {outlier_type}...")

    if not os.path.exists(output_name):
        print(f"Could not find file {output_name}. Skipping.")
        continue

    tab, meta_dict = read_casa_txt(output_name)

    # Iterate separately through polarizations
    these_pols = np.unique(np.array(tab['corr']))

    print(f"Found polarizations: {these_pols}")

    for this_pol in these_pols:

        mask = tab['corr'] == this_pol

        times_mjd = tab['time'][mask]
        amp = tab['y'][mask]

        # At this point we should have 1 amplitude per integration.
        # assert np.unique(times_mjd).size == len(times_mjd)

        # Only reject significant outliers by default (sigma > 5)
        outliers = find_outliers(amp, sigma=10.)

        plt.scatter(times_mjd, amp)

        if np.sum(outliers) == 0:
            print(f"No outliers found for pol {this_pol}. Skipping additional time flagging.")
        else:
            print(f"Flagging outliers for {this_pol}...")

            outlier_ranges = find_consecutive_ranges(outliers,
                                                    min_length=2,
                                                    max_gap=10,
                                                    pad_size=1)

            if len(outlier_ranges) == 0:
                print(f"No outliers found for pol {this_pol}. Skipping additional time flagging.")
                continue

            for start, end in outlier_ranges:
                plt.axvspan(times_mjd[start], times_mjd[end],
                            alpha=0.5, color='gray')

            for start, end in outlier_ranges:
                start_time = Time((times_mjd[start] - 0.5) / (24 * 3600.), format='mjd').iso.replace("-", "/").replace(" ", "/")
                end_time = Time((times_mjd[end] + 0.5) / (24 * 3600.), format='mjd').iso.replace("-", "/").replace(" ", "/")
                print(f"Flagging time range {start_time} - {end_time} (ISO)")

                if not apply_flagging:
                    # Apply the flag
                    continue

                flagdata(vis=str(ms_name),
                        field=",".join(sci_names),
                        correlation=this_pol,
                        mode='manual',
                        timerange=f"{start_time}~{end_time}",
                        flagbackup=False,)

                # Append to flagging file.
                with open(f"{data_path}/{ms_name.name}_flagging_{this_pol}.txt", 'a') as f:
                    f.write(f"mode='manual' timerange='{start_time}~{end_time}' correlation='{this_pol}'\n")

            if apply_flagging:
                flagmanager(vis=str(ms_name),
                            mode='save',
                            versionname=f'outlier_flagging_{outlier_type}_{this_pol}',
                            comment=f'Flagged outliers in OTF track for {outlier_type} with pol {this_pol}')

        plt.savefig(f"{data_path}/{ms_name.name}_{outlier_type}_vs_time_{this_pol}.png")
        plt.close()


# Export per sci target MSs to the output directory

for this_target in sci_names:

    print(f"Splitting {this_target}...")

    this_target_output_dir = osjoin(output_path, this_target)

    if not os.path.exists(this_target_output_dir):
        os.makedirs(this_target_output_dir)

    if line_selec is not None:
        line_str = f"_{line_selec}"
    else:
        line_str = ""

    this_outputvis = osjoin(this_target_output_dir, f"{obs_id}_{ms_name.name}_{this_target}{line_str}.ms")

    # Clean up any existing versions.
    rmtables(this_outputvis)

    split(vis=str(ms_name),
          outputvis=this_outputvis,
          datacolumn='corrected',
          field=this_target)
