
set zrange[0.0001:]
splot "vale.dat" using 1:2:(8e-6*($3)) ps 0.000001,"eigenvectors0.dat" using 1:2:(abs($3)) ps 0.01,"eigenvectors1.dat" using 1:2:(abs($3)) ps 0.01,"eigenvectors2.dat" using 1:2:(abs($3)) ps 0.01,"eigenvectors3.dat" using 1:2:(abs($3)) ps 0.01,"eigenvectors4.dat" using 1:2:(abs($3)) ps 0.01
