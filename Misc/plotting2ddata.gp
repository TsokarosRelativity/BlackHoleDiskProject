set pm3d map
set size ratio 1.0

# Get current timestamp
current_time = system("date +%Y%m%d_%H%M%S")

# Create a directory with the timestamp
directory = sprintf("plots_%s", current_time)

system(sprintf("mkdir -p %s", directory))


file = "September_20_2024_16-53/2D_data/all_vars.dat"
set terminal pngcairo size 1000,800 enhanced font "Arial, 24"

array colnames[16] = ["x", "z", "r", "theta", "alpha", "psi", "q", "betak", "betat", "dbetatdr", "dbetatdt", "He", "Hf", "Omega", "b2", "rho"] 

do for [i=5:16] {
  set output sprintf("%s/cartesian_2dinitialdata_%s.png", directory, colnames[i])
  set title sprintf("Initial 2d Data: %s", colnames[i])
  set xlabel colnames[1]
  set ylabel colnames[2]
  set zlabel colnames[i]
  splot [0:35][0:35] file u 2:3:i with pm3d title ""
}

system(sprintf("mv %s %s/plots/", directory, "September_20_2024_16-53"))
