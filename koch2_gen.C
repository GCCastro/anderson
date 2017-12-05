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
  int iters = 3;
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
    double Lit = l/pow(4,j+1);
    //int leng = verts.size();
    for (int i=0; i<verts.size()-1; i++)
    {
      std::vector<pair<double,double>>::iterator it = verts.begin() + i;
      int sinth = round((verts[i+1].second-verts[i].second)/Lit/4);
      int costh = round((verts[i+1].first-verts[i].first)/Lit/4);

      pair<double,double> x1;
      pair<double,double> x2;
      pair<double,double> x3;
      pair<double,double> x4;
      pair<double,double> x5;
      pair<double,double> x6;
      pair<double,double> x7;

      if(sinth==0 && costh==1)
      {
        x1 = make_pair(verts[i].first+Lit  ,verts[i].second    );
        x2 = make_pair(verts[i].first+Lit  ,verts[i].second+Lit);
        x3 = make_pair(verts[i].first+2*Lit,verts[i].second+Lit);
        x4 = make_pair(verts[i].first+2*Lit,verts[i].second    );
        x5 = make_pair(verts[i].first+2*Lit,verts[i].second-Lit);
        x6 = make_pair(verts[i].first+3*Lit,verts[i].second-Lit);
        x7 = make_pair(verts[i].first+3*Lit,verts[i].second    );
      }
      if (sinth==1 && costh==0)
      {
        x1 = make_pair(verts[i].first    ,verts[i].second+Lit  );
        x2 = make_pair(verts[i].first+Lit,verts[i].second+Lit  );
        x3 = make_pair(verts[i].first+Lit,verts[i].second+2*Lit);
        x4 = make_pair(verts[i].first    ,verts[i].second+2*Lit);
        x5 = make_pair(verts[i].first-Lit,verts[i].second+2*Lit);
        x6 = make_pair(verts[i].first-Lit,verts[i].second+3*Lit);
        x7 = make_pair(verts[i].first    ,verts[i].second+3*Lit);
      }

      if (sinth==-1 && costh==0)
      {
        x1 = make_pair(verts[i].first    ,verts[i].second-Lit  );
        x2 = make_pair(verts[i].first-Lit,verts[i].second-Lit  );
        x3 = make_pair(verts[i].first-Lit,verts[i].second-2*Lit);
        x4 = make_pair(verts[i].first    ,verts[i].second-2*Lit);
        x5 = make_pair(verts[i].first+Lit,verts[i].second-2*Lit);
        x6 = make_pair(verts[i].first+Lit,verts[i].second-3*Lit);
        x7 = make_pair(verts[i].first    ,verts[i].second-3*Lit);
      }

      if (sinth==0 && costh==-1)
      {
        x1 = make_pair(verts[i].first-Lit  ,verts[i].second    );
        x2 = make_pair(verts[i].first-Lit  ,verts[i].second-Lit);
        x3 = make_pair(verts[i].first-2*Lit,verts[i].second-Lit);
        x4 = make_pair(verts[i].first-2*Lit,verts[i].second    );
        x5 = make_pair(verts[i].first-2*Lit,verts[i].second+Lit);
        x6 = make_pair(verts[i].first-3*Lit,verts[i].second+Lit);
        x7 = make_pair(verts[i].first-3*Lit,verts[i].second    );
      }

      verts.insert(it+1,x1);
      it = verts.begin()+i;
      verts.insert(it+2,x2);
      it = verts.begin()+i;
      verts.insert(it+3,x3);
      it = verts.begin()+i;
      verts.insert(it+4,x4);
      it = verts.begin()+i;
      verts.insert(it+5,x5);
      it = verts.begin()+i;
      verts.insert(it+6,x6);
      it = verts.begin()+i;
      verts.insert(it+7,x7);
      i+=7;
    }
  }

  ofstream outfile;
  outfile.open("koch2.dat");

  for(int i=0; i<verts.size();i++)
    outfile << verts[i].first << " " << verts[i].second << endl;

  outfile.close();

  cout << "passei por aqui" << endl;

//Verificar se cada ponto esta dentro ou fora do quadrado de Koch

  double h = 0.009;
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
    if(y-L/2 > abs(x-L/2))  //acima do quadrado menor
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
      cout << n_cross << endl;
      if(n_cross%2)
      {
        grid[i] = 2;
        continue;
      }
    }

    else if(y-L/2 < abs(x-L/2))  //abaixo do quadrado menor
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

    else if(x-L/2>abs(y-L/2))  //a direita do quadrado menor
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

    if(x-L/2<abs(y-L/2))  //a direita do quadrado menor
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


  for(int i=0; i<leng; i++)
  {
    if(grid[i]==2)
    {
      //if(grid[i+1]==0 || grid[i-1]==0 || grid[i+N]==0 || grid[i-N]==0 || grid[i+N+1]==0 || grid[i-N-1]==0 || grid[i-N+1]==0 || grid[i+N-1]==0)
      if(grid[i+1]==0 || grid[i-1]==0 || grid[i+N]==0 || grid[i-N]==0)
        grid[i]=1;
    }
  }




  ofstream outfile2;
  outfile2.open("koch2square.dat");

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
