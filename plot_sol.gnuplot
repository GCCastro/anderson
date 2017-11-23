set pm3d
set palette defined (0 "blue", 0.001 "yellow", 0.002 "red")

set zrange[0.0001:]

splot "sol.dat" using 1:2:(($3)) w pm3d
