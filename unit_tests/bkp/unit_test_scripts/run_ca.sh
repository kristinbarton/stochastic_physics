#!/bin/sh
#SBATCH -e err
#SBATCH -o out
#SBATCH --account=ufs-artic
#SBATCH --qos=debug
#SBATCH --nodes=10
#SBATCH --ntasks-per-node=40
#SBATCH --time=20
#SBATCH --job-name="ca_tests"


RESDIR="RESTART"
if [ ! -d "$RESDIR" ]; then
    mkdir -p "$RESDIR"
fi

INDIR="INPUT"
if [ ! -d "$INDIR" ]; then
    mkdir -p "$INDIR"
fi

INFILE="ca_data.tile1.nc"
if [ ! -f "$INDIR/$INFILE" ]; then
    conda_env="/scratch4/BMC/ufs-artic/Kristin.Barton/envs/ufs-arctic"
    module load rdhpcs-conda
    conda activate ${conda_env}
    python gen_restart_files.py
    conda deactivate

    if [ ! -f "$INDIR/$INFILE" ]; then
        echo "ERROR GENERATING $INDIR/$INFILE"
        exit 1
    fi
fi

GRIDSPEC="C768_grid_spec.nc"
MOSAICDIR="/scratch4/BMC/ufs-artic/Kristin.Barton/repos/ufs-community/UFS_UTILS/build/fix/orog/C768/"
if [ ! -d "$MOSAICDIR" ]; then
    echo "DOES NOT EXIST: $MOSAIC"
    exit 1
fi
if [ ! -L "$INDIR/$GRIDSPEC" ]; then
    ln -s "$MOSAICDIR"/* "$INDIR"/.
    ln -s "$MOSAICDIR"/C768_mosaic.nc "$INDIR/$GRIDSPEC"
fi

EXEC=standalone_ca.x

if [ ! -f "$EXEC" ]; then
    sh compile_standalone_ca.ursa.intel
    if [ ! -f "$EXEC" ]; then
        echo "ERROR COMPILING $EXEC"
        exit 1
    fi
fi

# Spack stack conflicts with python environmnet, so this must go after the previous setup
module purge
module use /scratch4/BMC/ufs-artic/Kristin.Barton/repos/kristinbarton/ufs-arctic-workflow/main/ufs-weather-model/modulefiles/
module load ufs_ursa.intel.lua

ulimit -s unlimited
export OMP_STACKSIZE=512M
export KMP_AFFINITY=scatter
export OMP_NUM_THREADS=1

cp input.nml.noise input.nml
#sed -i -e "s/NOISE/0/g" input.nml
echo "option 0 run 1"
sleep 5
time srun --label -n 384 $EXEC  >& stdout_option_0
mkdir option_0
mv ca_out* option_0

#cp input.nml.noise input.nml
#sed -i -e "s/NOISE/1/g" input.nml
#echo "option 1 run 1"
#sleep 5
#time srun --label -n 384 $EXEC  >& stdout_option_1
#mkdir option_1
#mv ca_out* option_1
#
#cp input.nml.noise input.nml
#sed -i -e "s/NOISE/2/g" input.nml
#echo "option 2 run 1"
#sleep 5
#time srun --label -n 384 $EXEC  >& stdout_option_2
#mkdir option_2
#mv ca_out* option_2
#exit
#
#cp input.nml.noise input.nml
#sed -i -e "s/NOISE/2/g" input.nml
#echo "option 2 run 2"
#sleep 5
#time srun --label -n 384 $EXEC  >& stdout_option_2b
#mkdir option_2b
#mv ca_out* option_2b
#
#cp input.nml.noise input.nml
#sed -i -e "s/NOISE/1/g" input.nml
#echo "option 1 run 2"
#sleep 5
#time srun --label -n 384 $EXEC  >& stdout_option_1b
#mkdir option_1b
#mv ca_out* option_1b
#
#cp input.nml.noise input.nml
#sed -i -e "s/NOISE/0/g" input.nml
#echo "option 0 run 2"
#sleep 5
#time srun --label -n 384 $EXEC  >& stdout_option_0b
#mkdir option_0b
#mv ca_out* option_0b
