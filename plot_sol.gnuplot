set pm3d
#set palette defined (0 "blue", 0.001 "yellow", 0.002 "red")
set dgrid3d 300,300,2
set palette model HSV
set palette rgb 3,2,2

set zrange[0.:0.0015]

splot "sol.dat" using 1:2:(($3)) w pm3d
