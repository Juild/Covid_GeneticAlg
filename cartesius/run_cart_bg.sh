#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=10:00:00
#SBATCH --partition=normal

currentdir=$PWD
echo $currentdir

# create tmp dir and copy run file
tmp_dir="$(mktemp -d -p /scratch-shared)"
echo $tmp_dir
cp full_sci.py $tmp_dir
cd $tmp_dir
# init conda and execute script
module load 2019
module load Miniconda3
source /sw/arch/RedHatEnterpriseServer7/EB_production/2020/software/Miniconda3/4.7.12.1/etc/profile.d/conda.sh
conda activate "netsquid"
python3 full_sci.py
conda deactivate

# copy the contents in the tmp folder to the output folder in this directory
mkdir $currentdir/output
cp -r $tmp_dir $currentdir/output
cd $currentdir
rm -r $tmp_dir
