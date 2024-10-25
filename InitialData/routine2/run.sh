#!/bin/bash

exec 3>&1 1>>"../Logs/LogFile_$(date +"%B_%d_%Y_%H-%M").txt" 2>&1

echo "======BLACK HOLE DISK: INITIAL DATA GENERATION ROUTINE======"

# Check if a file name is provided as an argument
if [[ $# -eq 0 ]]
then
  echo "USAGE: $0 sho100.f95" 
  exit 1 
fi

echo "$(date +"%B_%d_%Y_%H-%M")"
mkdir "$(date +"%B_%d_%Y_%H-%M")"

file_name=$1

cp $1 "$(date +"%B_%d_%Y_%H-%M")"

cd "$(date +"%B_%d_%Y_%H-%M")"

# Create necessary folders
mkdir 2D_data
mkdir 3D_data
mkdir processed_grids
mkdir plots

# Run the gfortran command with the provided file name
#file_name=$1
echo "------------------------"
echo "generating 2D initial data..."
gfortran -fopenmp -O2 "$file_name" -o sho100 -L/usr/lib/x86_64-linux-gnu/lapack -llapack

# Create 2D data folder and execute sho100 in order to generate x-z data
mv sho100 2D_data/
cd 2D_data/
./sho100
cd ../../
 

# Generate 3D data and perform necessary variable calculations
#echo "------------------------"
#echo "generating 3D initial data..."
#python calc_vars_uj.py "$(date +"%B_%d_%Y_%H-%M")" -v
#echo "generating interpolated data"
#python interpolate_grid.py "$(date +"%B_%d_%Y_%H-%M")" 

echo "====END===="
