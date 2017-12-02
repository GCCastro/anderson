set term png

E4=4305.15
E3=4180.24
E2=4120.68
E1=4038.77
E0=3984.51


set output "1d_300_80_0.png"
plot "eigenvectors0_1d.dat" u 1:(abs($2))ps 0.0001, "sol1D.dat" u 1:(($2)*E0)ps 0.00001

set output "1d_300_80_1.png"
plot "eigenvectors1_1d.dat" u 1:(abs($2))ps 0.0001, "sol1D.dat" u 1:(($2)*E1)ps 0.00001

set output "1d_300_80_2.png"
plot "eigenvectors2_1d.dat" u 1:(abs($2))ps 0.0001, "sol1D.dat" u 1:(($2)*E2)ps 0.00001

set output "1d_300_80_3.png"
plot "eigenvectors3_1d.dat" u 1:(abs($2))ps 0.0001, "sol1D.dat" u 1:(($2)*E3)ps 0.00001

set output "1d_300_80_4.png"
plot "eigenvectors4_1d.dat" u 1:(abs($2))ps 0.0001, "sol1D.dat" u 1:(($2)*E4)ps 0.00001

