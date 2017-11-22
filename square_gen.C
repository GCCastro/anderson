#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

using namespace std;

int main()
{
  double L = 1.7;
  double h = 0.0049;
  int N = int(L/h);
  int leng = N*N;
  int bleng = 4*N-4;

  cout << N << endl;

  double *boundary = new double[bleng];

  for(int i=0;i<N;i++)
    boundary[i] = i;
  for(int i=0;i<N-2;i++)
  {
    boundary[2*i+N] = (i+1)*N;
    boundary[2*i+N+1] = (i+2)*N-1;
  }
  for(int i=0;i<N;i++)
    boundary[3*N-4+i] = N*(N-1)+i;

  int *grid = new int[leng];
  bool bflag = false;
  for(int i=0; i<leng; i++)
  {
    for(int j=0; j<bleng; j++)
    {
      if(i==boundary[j])
      {
        grid[i] = 1;
	bflag = true;
	break;
      }
    }
    if(bflag)
    {
      bflag=false;
      continue;
    }
    grid[i] = 2;
  }


  ofstream outfile2;
  outfile2.open("square.dat");

  for(int i=0; i<leng;i++)
  {
    div_t coord;
    coord = div(i,N);
    double x = coord.rem*h;
    double y = coord.quot*h;

    outfile2 << x << " " << y << " " << grid[i] << endl;
  }

  outfile2.close();

  delete[] boundary;

  return 0;
}
