#!/bin/bash

# Check if a file name is provided as an argument
if [ $# -eq 0 ]; then
    echo "Usage: $0 <FILE_NAME>"
    exit 1
fi

# Create necessary folders
mkdir 2D_data
mkdir 3D_data
mkdir processed_grids

# Run the gfortran command with the provided file name
file_name=$1
gfortran -fopenmp -O2 "$file_name" -o sho100 -L/usr/lib/x86_64-linux-gnu/lapack -llapack

# Create 2D data folder and execute sho100 in order to generate x-z data
mv sho100 2D_data/
cd 2D_data/
./sho100
cd ..

# Generate 3D data and perform necessary variable calculations
python calc_vars.py
python calc_horizon.py
python combine.py
