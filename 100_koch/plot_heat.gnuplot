set term png
set samples 10000000

file(n)  = sprintf("eigenvectors%d.dat",n)
graph(n) = sprintf("graf%d.png",n)

print file(1)

set xrange[0:1.6]
set yrange[0:1.6]
set logscale cb
set format cb "%T"

array max[100]

max[100]=4370.28
max[99]=4347.35
max[98]=4337.62
max[97]=4320.13
max[96]=4302.23
max[95]=4273.64
max[94]=4270.19
max[93]=4253.48
max[92]=4210.42
max[91]=4198.69
max[90]=4172.46
max[89]=4157.7
max[88]=4126.35
max[87]=4100.8
max[86]=4099.16
max[85]=4026.85
max[84]=4011.31
max[83]=3941.82
max[82]=3935.05
max[81]=3931.53
max[80]=3923.63
max[79]=3916.11
max[78]=3887.84
max[77]=3871.29
max[76]=3870.11
max[75]=3859.66
max[74]=3858.43
max[73]=3832.85
max[72]=3810.26
max[71]=3807.45
max[70]=3779.93
max[69]=3770.35
max[68]=3728.12
max[67]=3718.78
max[66]=3703.77
max[65]=3695.4
max[64]=3662.75
max[63]=3662.44
max[62]=3646.49
max[61]=3632.78
max[60]=3609.37
max[59]=3585.51
max[58]=3526.54
max[57]=3493.14
max[56]=3473.88
max[55]=3457.71
max[54]=3412.01
max[53]=3352.04
max[52]=3311.75
max[51]=3307.16
max[50]=3268.47
max[49]=3254.1
max[48]=3250.93
max[47]=3239.35
max[46]=3223.69
max[45]=3142.02
max[44]=3135.71
max[43]=3111.07
max[42]=3080.55
max[41]=3067.95
max[40]=3059.17
max[39]=3045.79
max[38]=2954.88
max[37]=2942.22
max[36]=2884.39
max[35]=2853.16
max[34]=2807.21
max[33]=2745.8
max[32]=2675.43
max[31]=2674.63
max[30]=2642.42
max[29]=2615.69
max[28]=2556.72
max[27]=2545.9
max[26]=2540.21
max[25]=2502.11
max[24]=2458.48
max[23]=2434.73
max[22]=2355.05
max[21]=2319.87
max[20]=2306.02
max[19]=2273.8
max[18]=2272.88
max[17]=2117.24
max[16]=2093.1
max[15]=2011.16
max[14]=1976.7
max[13]=1921.
max[12]=1831.97
max[11]=1822.4
max[10] =1767.13
max[9] =1739.89
max[8] =1581.14
max[7] =1551.58
max[6] =1504.89
max[5] =1468.53
max[4] =1451.27
max[3] =1436.7
max[2] =1252.47
max[1] =1185.89

#set view map
#set dgrid3d
set pm3d map
set palette defined (0 0 0 0.5, 1 0 0 1, 2 0 0.5 1, 3 0 1 1, 4 0.5 1 0.5, 5 1 1 0, 6 1 0.5 0, 7 1 0 0, 8 0.5 0 0)
do for [i=0:99]{
set output graph(i)
splot file(i) u 1:2:(abs($3)), "vale.dat" u 1:2:(($3)<1/max[i+1] ? (1.0) : 1/0) w points ps 0.01 lc "white"
}
#splot "vale.dat" u 1:2:(0.0) w points
#splot "eigenvectors0.dat" u 1:2:(abs($3))

pause -1
