set zrange[0.0001:]

splot "eigenvectors0.dat" using 1:2:(abs($3)) ps 0.1,"eigenvectors1.dat" using 1:2:(abs($3)) ps 0.1,"eigenvectors2.dat" using 1:2:(abs($3)) ps 0.1,"eigenvectors3.dat" using 1:2:(abs($3)) ps 0.1,"eigenvectors4.dat" using 1:2:(abs($3)) ps 0.1, "vale3.dat" ps 0.1
