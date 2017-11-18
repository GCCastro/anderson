#include <iostream>
#include <fstream>
#include <stdio.h>

#include <vector>
#include <math.h>

using namespace std;

int main()
{

  double L = 1.6;
  double l = 1.;

//Verificar se cada ponto esta dentro ou fora do quadrado de Koch

  double h = 0.0039;
  int N = (int)(L/h);
  int leng = N*N;

  int *grid = new int[leng];
  for(int i=0; i<leng; i++)
    grid[i] = 0;

  for(int i=0; i<leng; i++)
  {
    div_t coord;
    coord = div(i,N);
    double x = coord.rem*h;
    double y = coord.quot*h;
    if(x>2./5.*L && y<2./5.*L && y>1./3.*L)
      grid[i] = 0;
    else if(x>2./5.*L && x< 1./2.*L && y<1./3.*L && y>1./5.*L)
      grid[i] = 0;
    else if(x<h || x>L-2*h || y>L-2*h || y<h)
      grid[i] = 1;
    else if(x>2./5.*L-5*h && x<2./5.*L+5*h && y>3./4.*L-5*h && y<3./4.*L+5*h)
      grid[i] = 1;
    else if(x>1./5.*L-5*h && x<1./5.*L+5*h && y>3./5.*L-5*h && y<3./5.*L+5*h)
      grid[i] = 1;
    else if(x>4./5.*L && y>2./3.*L-5*h && y<2./3.*L+5*h)
      grid[i] = 1;
    else if(x>2./5.*L && y>2./5.*L-5*h && y<2./5.*L+5*h)
      grid[i] = 1;
    else if(x>1./2.*L && y>1./3.*L-5*h && y<1./3.*L+5*h)
      grid[i] = 1;
    else if(x>2./5.*L && x<1./2.*L && y>1./5.*L-5*h && y<1./5.*L+5*h)
      grid[i] = 1;
    else if(x>2./5.*L-5*h && x<2./5.*L+5*h && y>1./5.*L-5*h && y<2./5.*L+5*h)
      grid[i] = 1;
    else if(x>1./2.*L-5*h && x<1./2.*L+5*h && y>1./5.*L-5*h && y<1./3.*L+5*h)
      grid[i] = 1;
    else
      grid[i] = 2;

  }


  ofstream outfile2;
  outfile2.open("estranho.dat");

  for(int i=0; i<leng;i++)
  {
    div_t coord;
    coord = div(i,N);
    double x = coord.rem*h;
    double y = coord.quot*h;

    outfile2 << x << " " << y << " " << grid[i] << endl;
  }

  outfile2.close();



  return 0;
}
