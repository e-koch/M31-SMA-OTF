#!/bin/bash
#SBATCH --time=100:00:00
#SBATCH --cpus-per-task=64
#SBATCH --mem=64G
#SBATCH --job-name=m31-sma-co21-imaging-%J
#SBATCH --output=m31-sma-co21-imaging-%J.out
#SBATCH --mail-user=ekoch@ualberta.ca
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


# For now, just pass the target name from the cmd line with sbatch.
# Will want the ability to pass an array of target names later.

export this_target=$1
# export this_line_product='hilores'

echo $SLURM_CPUS_PER_TASK

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

module load StdEnv
module load qt

source /home/ekoch/.bashrc

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/usr/

export CASALD_LIBRARY_PATH=$LD_LIBRARY_PATH


# Ensure no time overlap in job start times
python3 -c "import time, random; time.sleep(random.randint(2, 120))"

export data_path="/home/ekoch/scratch/VLAXL_imaging/MeasurementSets/"

export casa_executable="/home/ekoch/casa-6.6.1-17-pipeline-2024.1.0.8/bin/casa"
export casa_script="/home/ekoch/M31-SMA-OTF/fir_imaging/fir_per_mosaic_test_deconv_imaging.py"


export script_args="$this_target"
echo "Args passed to script: $script_args"
xvfb-run -a $casa_executable --rcdir ~/.casa --nologger --nogui --log2term -c $casa_script $script_args
