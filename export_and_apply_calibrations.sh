#!/bin/bash


# export DATADIR="/reduction/erickoch/M31/sma/"
export DATADIR="/reduction11/erickoch/M31/sma_12CO/"
# export DATADIR="/reduction11/erickoch/IC342/sma/"
export COMPASSDIR="/mnt/COMPASS9/sma/scripts/"

export THIS_PYTHON="/home/gkeating/.conda/envs/sma-compass/bin/python"
# export THIS_SCRIPT="/home/gkeating/repo/SMA-COMPASS/sma_compass/ms_output.py"
export THIS_SCRIPT="/home/erickoch/SMA-COMPASS/sma_compass/ms_output.py"
# export THIS_SCRIPT="/home/erickoch/M81_SMA/ms_output_skipints.py"

export THIS_COORDS_SCRIPT="/home/erickoch/M31_SMA/extract_spatial_coverage.py"


# export THIS_CASA="/py3opt/casa-6.6.3-22-py3.8.el8/bin/casa --nogui --nologger"
export THIS_CASA="xvfb-run -a /py3opt/casa-6.6.1-17-pipeline-2024.1.0.8/bin/casa --nogui --nologger"
export THIS_CASA_SCRIPT="/home/erickoch/M31_SMA/apply_compass_solutions.py"


# All SPWs
# export RECHUNK=16
# export RECHUNK=4


# 12CO only
export RECHUNK=8

export LINENAME="12CO"
export SBAND='u'
export CHUNK='1'



# 13CO only
# export RECHUNK=8

# export LINENAME="13CO"
# export SBAND='l'
# export CHUNK='1 2'


export NUMJOBS=10


mkdir -p "$DATADIR"


# Processing function:
# Process a single track
#
# Parameters:
#   TRACKNUM - the track number to process (required)
#
# Description:
#   This function processes a single track by:
#   1. Generating the compass solutions
#   2. Applying the compass solutions in CASA
#   3. Outputting the processed data in a subdirectory
process_track() {
    local TRACKNUM=$1

    if [ -z "$TRACKNUM" ]; then
        echo "Error: TRACKNUM is required"
        return 1
    fi

    local LOGFILE="${DATADIR}/${TRACKNUM}.log"

    {
        echo "Processing started for TRACKNUM: $TRACKNUM"

        mkdir -p "$DATADIR"
        cd "$DATADIR" || return 1

        # Export to MS
        "$THIS_PYTHON" "$THIS_SCRIPT" -j "$COMPASSDIR/$TRACKNUM.json" -r "$RECHUNK" -sb "$SBAND" -c "$CHUNK"

        cd $TRACKNUM

        # Create spatial coverage table and plots
        # "$THIS_PYTHON" "$THIS_COORDS_SCRIPT" --obsid "$TRACKNUM" --plot

        $THIS_CASA -c "$THIS_CASA_SCRIPT" "$TRACKNUM" "$LINENAME"

        echo "Processing completed for TRACKNUM: $TRACKNUM"

        cd $DATADIR

    } 2>&1 | tee "$LOGFILE"

}


# Process multiple tracks in parallel using GNU parallel
#
# Parameters:
#   -j JOBS     - number of jobs to run in parallel (optional, default: unlimited)
#   TRACKNUMS   - array of track numbers to process (required)
#
# Description:
#   This function processes multiple tracks in parallel using GNU parallel
#
# Usage:
#   process_tracks_parallel "13757" "13763" "13768"
#   process_tracks_parallel -j 4 "13757" "13763" "13768"
#   or
#   tracks=("13757" "13763" "13768")
#   process_tracks_parallel -j 4 "${tracks[@]}"
process_tracks_parallel() {
    local jobs=""
    local tracknums=()

    # Parse options
    while [[ $# -gt 0 ]]; do
        case $1 in
            -j)
                jobs="$2"
                shift 2
                ;;
            *)
                tracknums+=("$1")
                shift
                ;;
        esac
    done

    if [ ${#tracknums[@]} -eq 0 ]; then
        echo "Error: At least one track number is required"
        return 1
    fi

    echo "Starting parallel processing of ${#tracknums[@]} tracks: ${tracknums[*]}"
    if [ -n "$jobs" ]; then
        echo "Running with max $jobs jobs in parallel"
    fi

    # Export the function so parallel can use it
    export -f process_track
    export DATADIR THIS_PYTHON THIS_SCRIPT THIS_COORDS_SCRIPT THIS_CASA THIS_CASA_SCRIPT COMPASSDIR RECHUNK

    # Use parallel to process all tracks
    if [ -n "$jobs" ]; then
        printf "%s\n" "${tracknums[@]}" | parallel -j "$jobs" process_track {}
    else
        printf "%s\n" "${tracknums[@]}" | parallel process_track {}
    fi

    echo "All tracks completed"
}


# M31

###############
# Brick A
###############

brick_a_tracks=(
    "13757"  # Track 1
    "13763"  # Track 1
    "13768"  # Track 2
    "13773"  # Track 3
    "13777"  # Track 4
    "13783"  # Track 4
    "13813"  # Track 4
    "13801"  # Track 5
    "13806"  # Track 6
)

process_tracks_parallel -j $NUMJOBS "${brick_a_tracks[@]}"



###############
# Brick B
###############

brick_b_tracks=(
    "13785"  # Track 8
    "13815"  # Track 8
    "13819"  # Track 8
    "13821"  # Track 9
    "13824"  # Track 10
    "13826"  # Track 11
    "13828"  # Track 12
    "13831"  # Track 13
    "13834"  # Track 14
)

# process_tracks_parallel -j $NUMJOBS "${brick_b_tracks[@]}"

###############
# Brick C
###############

brick_c_tracks=(
    "13837"  # Track 16
    "13839"  # Track 17
    "13862"  # Track 18
    "13865"  # Track 19
    "13887"  # Track 20
    "13895"  # Track 21
    "13898"  # Track 22
    "13900"  # Track 23
    "13926"  # Track 23
)

# process_tracks_parallel -j $NUMJOBS "${brick_c_tracks[@]}"


###############
# Brick D
###############

brick_d_tracks=(
    # "13901"  # Track 24
    # "13921"  # Track 24
    # "13936"  # Track 24
    "13951"  # Track 24
    # "13902"  # Track 25
    "13935"  # Track 25
    # "13904"  # Track 26
    "13938"  # Track 26
    "13949"  # Track 26
    # "13907"  # Track 27
    "13914"  # Track 28
    "13939"  # Track 28
    "13947"  # Track 28
    # "13927"  # Track 29
    "13931"  # Track 30
)

# process_tracks_parallel -j $NUMJOBS "${brick_d_tracks[@]}"


###############
# Brick E
###############

brick_e_tracks=(
    # "13886"  # Track 7
    "13928"  # Track 15
    "13933"  # Track 15
)

# process_tracks_parallel -j $NUMJOBS "${brick_e_tracks[@]}"

# process_track "13933"


# IC342
ic342_tracks=(
    # "13981"  # Track 1
    "13985"  # Track 2
    "13990"  # Track 3
    "13999"  # Track 4
)

# process_tracks_parallel -j $NUMJOBS "${ic342_tracks[@]}"

# process_track "13981"

# Tests on 2025A-S014
# NGC7293
# process_track "13790"
# process_track "13808"
# process_track "13867"
