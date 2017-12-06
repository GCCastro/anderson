set term png

set logscale y

set output "agmon_1d_0.png"
plot "eigenvectors0_1d.dat" u 1:(abs($2)) w l, "agmon0_1d.dat" w l
set output "agmon_1d_1.png"
plot "eigenvectors1_1d.dat" u 1:(abs($2)) w l, "agmon1_1d.dat" w l
set output "agmon_1d_2.png"
plot "eigenvectors2_1d.dat" u 1:(abs($2)) w l, "agmon2_1d.dat" w l
set output "agmon_1d_3.png"
plot "eigenvectors3_1d.dat" u 1:(abs($2)) w l, "agmon3_1d.dat" w l
set output "agmon_1d_4.png"
plot "eigenvectors4_1d.dat" u 1:(abs($2)) w l, "agmon4_1d.dat" w l

