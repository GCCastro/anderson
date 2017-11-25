#include <iostream>
#include <fstream>
#include <stdio.h>

#include <vector>
#include <math.h>

using namespace std;

int main()
{

  double L = 1.7;
  double l = 1.;
  int iters = 2;
  double i_square = L/2-l/2;
  double f_square = L/2+l/2;

  vector<pair<double,double>> verts (2);

//Definicao dos pontos inicial e final da linha de Koch

  verts[0].first = i_square;
  verts[0].second = f_square;
  verts[1].first = f_square;
  verts[1].second = f_square;

//Criacao dos vertices da iteracao iters da linha de Koch 

  for (int j=0; j<iters; j++)
  {
    double Lit = L/pow(3,j);
    //int leng = verts.size();
    for (int i=0; i<verts.size()-1; i++)
    {
      std::vector<pair<double,double>>::iterator it = verts.begin() + i;
      double sinth = (verts[i+1].second-verts[i].second)/Lit;
      double costh = (verts[i+1].first-verts[i].first)/Lit; 

      pair<double,double> x1(verts[i].first+Lit/3*costh,verts[i].second+Lit/3*sinth);
      pair<double,double> x2(verts[i].first+Lit/2*costh-Lit/sqrt(12)*sinth,verts[i].second+Lit/2*sinth+Lit/sqrt(12)*costh);
      pair<double,double> x3(verts[i].first+2*Lit/3*costh,verts[i].second+2*Lit/3*sinth);

      verts.insert(it+1,x1);
      it = verts.begin()+i;
      verts.insert(it+2,x2);
      it = verts.begin()+i;
      verts.insert(it+3,x3);
      i+=3;
    }
  }

//Verificar se cada ponto esta dentro ou fora do quadrado de Koch

  double h = 0.0049;
  int N = (int)(L/h);
  int leng = N*N;
  int npoints = verts.size();

  int *grid = new int[leng];
  for(int i=0; i<leng; i++)
    grid[i] = 0;

  for(int i=0; i<leng; i++)
  {
    div_t coord;
    coord = div(i,N);
    double x = coord.rem*h;
    double y = coord.quot*h;
    //Dentro do quadrado menor
    if(x>i_square && x<f_square)
    {
      if(y>i_square && y<f_square)
      {
	grid[i]=2;      
	continue;
      }
      else
      {
        if(y>f_square)  //acima do quadrado menor
	{
          int n_cross = 0;
          for(int j=0; j<npoints-1; j++)
	  {
            double m = (verts[j+1].second-verts[j].second)/(verts[j+1].first-verts[j].first);
            double b = verts[j].second - m*verts[j].first;
            double mx = y - b;
            bool inv = verts[j+1].first>verts[j].first;
            if(abs(mx)<abs(m*max(verts[j+1].first,verts[j].first)) && abs(mx)>abs(m*min(verts[j].first,verts[j+1].first)))
              if((abs(mx) <= abs(m*x) && !inv) || (abs(mx)>= abs(m*x) && inv))
                n_cross++;
          }

          if(n_cross%2)
          {
            grid[i] = 2;
	    continue;
	  }
	}

	else if(y<i_square)  //abaixo do quadrado menor
	{
          int n_cross = 0;
          for(int j=0; j<npoints-1; j++)
	  {
            double m = -(verts[j+1].second-verts[j].second)/(verts[j+1].first-verts[j].first);
            double b = i_square-abs(verts[j].second-f_square) - m*verts[j].first;
            double mx = y - b;
            bool inv = verts[j+1].first>verts[j].first;
            if(abs(mx)<abs(m*max(verts[j+1].first,verts[j].first)) && abs(mx)>abs(m*min(verts[j].first,verts[j+1].first)))
              if((abs(mx) <= abs(m*x) && !inv) || (abs(mx)>= abs(m*x) && inv))
                n_cross++;

          }

          if(n_cross%2)
          {
            grid[i] = 2;
	    continue;
	  }
	}


      }
    }
    else
    {
      if(x>f_square)  //a direita do quadrado menor
      {
        int n_cross = 0;
        for(int j=0; j<npoints-1; j++)
        {
          double m = (verts[j+1].second-verts[j].second)/(verts[j+1].first-verts[j].first);
          double b = verts[j].second - m*verts[j].first;
          double my = x - b;
          bool inv = verts[j+1].first>verts[j].first;
          if(abs(my)<abs(m*max(verts[j+1].first,verts[j].first)) && abs(my)>abs(m*min(verts[j].first,verts[j+1].first)))
            if((abs(my) <= abs(m*y) && !inv) || (abs(my)>= abs(m*y) && inv))
              n_cross++;
        }

        if(n_cross%2)
        {
          grid[i] = 2;
          continue;
        }
      }

      if(x<i_square)  //a direita do quadrado menor
      {
        int n_cross = 0;
        for(int j=0; j<npoints-1; j++)
        {
          double m = -(verts[j+1].second-verts[j].second)/(verts[j+1].first-verts[j].first);
          double b = i_square-abs(verts[j].second-f_square) - m*verts[j].first;
          double my = x - b;
          bool inv = verts[j+1].first>verts[j].first;
          if(abs(my)<abs(m*max(verts[j+1].first,verts[j].first)) && abs(my)>abs(m*min(verts[j].first,verts[j+1].first)))
            if((abs(my) <= abs(m*y) && !inv) || (abs(my)>= abs(m*y) && inv))
              n_cross++;
        }

        if(n_cross%2)
        {
          grid[i] = 2;
          continue;
        }
      }


    }
  }


  for(int i=0; i<leng; i++)
  {
    if(grid[i]==2)
    {
      //if(grid[i+1]==0 || grid[i-1]==0 || grid[i+N]==0 || grid[i-N]==0 || grid[i+N+1]==0 || grid[i-N-1]==0 || grid[i-N+1]==0 || grid[i+N-1]==0)
      if(grid[i+1]==0 || grid[i-1]==0 || grid[i+N]==0 || grid[i-N]==0)
        grid[i]=1;
    }
  }



  ofstream outfile;
  outfile.open("koch.dat");

  for(int i=0; i<verts.size();i++)
    outfile << verts[i].first << " " << verts[i].second << endl;

  outfile.close();

  ofstream outfile2;
  outfile2.open("kochsquare.dat");

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
