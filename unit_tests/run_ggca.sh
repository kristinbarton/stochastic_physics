#!/bin/bash
#SBATCH -e err
#SBATCH -o out
#SBATCH --account=ufs-artic
#SBATCH --qos=debug
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=40
#SBATCH --time=20
#SBATCH --job-name="run_ggca"

# This is for compiling and testing the standalone version
# of the Gaussian Grid Cellular Automata code

# Compile the standalone_ca.x code
EXEC=standalone_ca.x
if [ ! -f "build/$EXEC" ]; then
    source ./env_ursa_intelllvm.sh
    make ca
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
cp ../templates/input.nml.ca_template input.nml
sed -i -e "s/LOX/1/g" input.nml
sed -i -e "s/LOY/1/g" input.nml
sed -i -e "s/NPX/$NPX/g" input.nml
sed -i -e "s/NPY/$NPY/g" input.nml
sed -i -e "s/RES/$RES/g" input.nml
sed -i -e "s/CA_SGS/.false./g" input.nml
sed -i -e "s/CA_GLOBAL/.true./g" input.nml
sed -i -e "s/WARM_START/.false./g" input.nml

# Run executable
ln -s ../build/$EXEC .
export OMP_NUM_THREADS=1
time srun --label -n 6 $EXEC >& ggca.stdout
