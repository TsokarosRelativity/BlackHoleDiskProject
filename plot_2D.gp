# Set common plot settings
set pm3d map
set size ratio 1.0

# Create the subdirectory if it doesn't exist
system("mkdir -p plots")

# Loop through each .dat file in the directory
# filelist = system("ls *.dat")
filelist = "all_vars.dat"
do for [file in filelist] {
    # Extract the filename without extension
    # filename = strcol(1)
    filename_noext = file[:strlen(file)-4]

    # Set output filename with subdirectory path
    set terminal pngcairo enhanced
    # set output 'plots/'.filename_noext.'.png'
    set output 'plots/b2.png'

    # Add title
    # set title sprintf("%s: Columns 1, 2, and 15", filename_noext)
    set title "b2 Plot in Positive x-z Plane" 
    set xlabel "x"
    set ylabel "z" 
    
    # Plot the data
    splot [0:35][0:35] file u 1:2:15
    

    # Reset terminal
    set output
}



