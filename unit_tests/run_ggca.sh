#!/bin/bash
#SBATCH -e err
#SBATCH -o out
#SBATCH --account=ufs-artic
#SBATCH --qos=debug
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --time=10
#SBATCH --job-name="run_ggca"

# This is for compiling and testing the standalone version
# of the Gaussian Grid Cellular Automata code

# Compile the standalone_ggca.x code
EXEC=standalone_ggca.x
if [ ! -f "build/$EXEC" ]; then
    echo "Compiling $EXEC"
    source ./envs/env_ursa_intelllvm.sh
    make ggca
    if [ ! -f "build/$EXEC" ]; then
        echo "ERROR COMPILING $EXEC"
        exit 1
    fi
fi

# Use low resolution grid for testing...
RES=96
NPX=`expr $RES + 1`
NPY=`expr $RES + 1`
TEST="run_ggca"

# Restart & INput directories must exist already or executable will fail
mkdir -p $TEST && cd $TEST
mkdir -p "INPUT"
mkdir -p "RESTART"

# Gather cubed sphere grid files
GRIDSPEC="C${RES}_grid_spec.nc"
MOSAICDIR="/scratch4/BMC/ufs-artic/Kristin.Barton/repos/ufs-community/UFS_UTILS/build/fix/orog/C${RES}/"
if [ ! -d "$MOSAICDIR" ]; then
    echo "DOES NOT EXIST: $MOSAIC"
    exit 1
fi
if [ ! -L "INPUT/$GRIDSPEC" ]; then
    ln -s "$MOSAICDIR"/* INPUT/.
    ln -s "$MOSAICDIR"/C${RES}_mosaic.nc "INPUT/$GRIDSPEC"
fi

# Populate namelist template
cp ../templates/input.nml.gg_template input.nml
sed -i -e "s/LOX/1/g" input.nml
sed -i -e "s/LOY/1/g" input.nml
sed -i -e "s/NPX/$NPX/g" input.nml
sed -i -e "s/NPY/$NPY/g" input.nml
sed -i -e "s/RES/$RES/g" input.nml

# Run executable
ln -s ../build/$EXEC .
export OMP_NUM_THREADS=1
echo "Running executable"
time srun --label -n 6 $EXEC >& ggca.stdout
