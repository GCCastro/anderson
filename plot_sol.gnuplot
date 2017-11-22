set pm3d
set palette defined (0 "blue", 0.001 "yellow", 0.002 "red")

splot "sol.dat" w pm3d using 1:2:(1/($3))
