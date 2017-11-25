set contour base
set zrange[0:0.0005]
set cntrparam level incremental 0, 0.00002, 0.0005
unset surface
set table 'cont.dat'
splot "sol.dat" u 1:2:3
unset table

reset
unset key
#set palette model HSV
#set palette rgb 3,2,2
set palette rgb 21,22,23
set zrange[0:0.0005]
p 'sol.dat' u 1:2:3 with image, 'cont.dat' w l lt -1 lw 1.5
pause -1
