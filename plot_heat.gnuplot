set xrange[0:1.6]
set yrange[0:1.6]
set logscale cb

#set view map
#set dgrid3d
set pm3d map
set palette defined (0 0 0 0.5, 1 0 0 1, 2 0 0.5 1, 3 0 1 1, 4 0.5 1 0.5, 5 1 1 0, 6 1 0.5 0, 7 1 0 0, 8 0.5 0 0)
splot "eigenvectors0.dat" u 1:2:(abs($3)), "vale.dat" u 1:2:(($3)<1/1006.07 ? (1.0) : 1/0) w points ps 0.01 lc "white"
#splot "vale.dat" u 1:2:(0.0) w points
#splot "eigenvectors0.dat" u 1:2:(abs($3))

pause -1
