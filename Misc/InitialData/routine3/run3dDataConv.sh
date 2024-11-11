#!/bin/bash 

exec 3>&1 1>>"cowboyoutput.txt" 2>&1


echo "======RUN THIS SHIT COWBOY====="
# Check if a file name is provided as an argument
python calc_vars_uj.py October_15_2024_22-05/ 2D_data/   
python calc_vars_uj.py October_15_2024_22-05/ 2D_data_100000/   
python calc_vars_uj.py October_15_2024_22-05/ 2D_data_1k/

echo "DONEEEEEE"
