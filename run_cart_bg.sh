#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=1:00:00
#SBATCH --partition=normal

currentdir=$PWD
echo $currentdir

# create tmp dir and copy run file
tmp_dir="$(mktemp -d -p /scratch-shared)"
echo $tmp_dir
cp covidga $tmp_dir
cd $tmp_dir
# init conda and execute script
module load 2020
module load GCC/9.3.0
module load GSL/2.6-GCC-9.3.0
bash covidga

# copy the contents in the tmp folder to the output folder in this directory
mkdir $currentdir/output
cp -r $tmp_dir $currentdir/output
cd $currentdir
rm -r $tmp_dir
