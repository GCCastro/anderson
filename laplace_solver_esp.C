#include <eigen3/Eigen/Eigenvalues>
#include <eigen3/Eigen/Eigen>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Core>

#include <GenEigsSolver.h>
#include <MatOp/SparseGenMatProd.h> 


#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <vector>
#include <tuple>
#include <set>

#include <gsl/gsl_rng.h>  //random generator
#include "rtnorm.hpp"     //truncated normal distribution

typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> T;

using namespace std;
using namespace Spectra;

double f(double h, int xi, int yj, int N, double **a, int D)
{
  int x = (int)((double)xi/((double)N)*D);
  int y = (int)((double)yj/((double)N)*D);

  return a[x][y];
  //return 0.;
}


bool sort_function(pair<int,int> pair1,pair<int,int> pair2)
{
  return (pair1.first<pair2.first);
}


int main()
{
  cout << "Hail Eris" << endl;

  cout << "sigma da distribuicao?" << endl;
  double sig;
  cin >> sig;


  double Vini = 0.;
  double Vmed = 4000.;
  double Vend = 8000.;
  int D = 20;
  double **a = new double*[D];
  for(int i=0; i<D; i++)
    a[i] = new double[D];

  //--- GSL random init ---
  gsl_rng_env_setup();                          // Read variable environnement
  const gsl_rng_type* type = gsl_rng_default;   // Default algorithm 'twister'
  gsl_rng *gen = gsl_rng_alloc (type);          // Rand generator allocation
  gsl_rng_set(gen,time(NULL));

  for(int i=0; i<D; i++)
    for(int j=0; j<D; j++)
    {
      a[i][j] = rtnorm(gen,Vini,Vend,Vmed,sig).first;
    }


  //isto devia vir directamente do ficheiro
  double L = 1.7;
  double h = 0.009;
  int Ntot = int(L/h);
  cout << Ntot << endl;

//  cout << "condicoes fronteira?" << endl
//  string bound;
//  cin >> bound;

  cout << "ficheiro da regiao?" << endl;
  string filename;
  cin >> filename;

//  if((bound == "esp") && (filename != "square.dat"))
//  {
//    cout << "isto assim nao vai funcionar" << endl;
//    return 0;
//  }

  ifstream infile;
  infile.open(filename.c_str());
  vector<vector<pair<int,int>>> grid;
  vector<pair<int,int>> line;
  double xf,yf;
  int c;
  int leng=0;  //o leng neste ficheiro contem a fronteira
  double xold=0.;


  while(infile >> xf >> yf >> c)
  {
    if(xf<xold)
    {
      grid.push_back(line);
      line.clear();
    }
    if(c>=1)
    {
      pair<int,int> entry(leng,c);
      line.push_back(entry);
      leng++;
    }
    else if(c==0)
    {
      pair<int,int> entry(0,c);
      line.push_back(entry);
    }
    else
    {
      cout << "Estas a fazer alguma coisa mal, utilizador deste programa" << endl;
    }
    xold=xf;
  }
  grid.push_back(line);

  cout << "BATATA" << endl;

  int Nmodos = 3;

  vector<T> coefficients;

  Eigen::VectorXd b(leng);

  for(int i=0;i<Ntot;i++) //loop ao longo das colunas
  {
    for(int j=0;j<Ntot;j++) //loop dentro da linha
    {
      if(grid[i][j].second>=1)
      {
        T entry1(grid[i][j].first,grid[i][j].first,4./h/h+f(h,j,i,Ntot,a,D));
        coefficients.push_back(entry1);

        int tmp1 = j+1;
        if(tmp1==Ntot)
          tmp1 = 0;
        if(grid[i][tmp1].second>=1)
	{
          T entry2(grid[i][j].first,grid[i][tmp1].first,-1./h/h);
          coefficients.push_back(entry2);
	}

        int tmp2 = j-1;
        if(tmp2==-1)
          tmp2 = Ntot-1;
        if(grid[i][tmp2].second>=1)
	{
          T entry3(grid[i][j].first,grid[i][tmp2].first,-1./h/h);
          coefficients.push_back(entry3);
	}

        int tmp3 = i+1;
        if(tmp3==Ntot)
          tmp3 = 0;
        if(grid[tmp3][j].second>=1)
	{
          T entry4(grid[i][j].first,grid[tmp3][j].first,-1./h/h);
          coefficients.push_back(entry4);
	}

        int tmp4 = i-1;
        if(tmp4==-1)
          tmp4 = Ntot-1;
        if(grid[tmp4][j].second>=1)
	{
          T entry5(grid[i][j].first,grid[tmp4][j].first,-1./h/h);
          coefficients.push_back(entry5);
	}
  
        b[grid[i][j].first] = 1.;
      }
    }
  }

  cout << "BATATA" << endl;

  SpMat A(leng,leng);
  A.setFromTriplets(coefficients.begin(),coefficients.end());

//  for(int i=0; i<leng; i++)
//  {
//    for(int j=0; j<leng; j++)
//      cout << A.coeff(i,j) << " ";
//    cout << endl;
//  }
//
//  for(int i=0; i<leng; i++)
//    cout << b(i) << endl;

  Eigen::SparseLU<SpMat> solver(A);

  Eigen::VectorXd x = solver.solve(b);

  cout << "sobrevivi a calcular o u" << endl;


//Algoritmo para descobrir minimos da solucao
  bool vizflag = false;
  vector<pair<int,int>> minpoints;
  for(int i=0; i<Ntot; i++)
  {
    for(int j=0; j<Ntot; j++)
    {
      if(grid[j][i].second >= 1)
      {
        int tmp1 = j+1;
        if(tmp1==Ntot)
          tmp1 = 0;
        int tmp2 = j-1;
        if(tmp2==-1)
          tmp2 = Ntot-1;
        int tmp3 = i+1;
        if(tmp3==Ntot)
          tmp3 = 0;
        int tmp4 = i-1;
        if(tmp4==-1)
          tmp4 = Ntot-1;

        pair<int,int> viz[8] = {make_pair(i,tmp1),make_pair(i,tmp2),make_pair(tmp4,tmp1),make_pair(tmp4,tmp2),make_pair(tmp3,tmp1),make_pair(tmp3,tmp2),make_pair(tmp4,j),make_pair(tmp3,j)};
        for(int k=0; k<8; k++)
          if(x(grid[viz[k].second][viz[k].first].first)<x(grid[j][i].first)) //se estiver < estou a descobrir minimos; isto nao vai nada correr mal...
          {
            vizflag = true;
            break;
          }
        if(vizflag)
        {
          vizflag=false;
          continue;
        }
        // se quiser os maximos/minimos uso este
        minpoints.push_back(make_pair(i,j));
        //se quiser os pontos a volta dos maximos/minimos uso este
        //for(int v=0; v<8; v++)
        //  minpoints.push_back(viz[v]);
      }
    }
  }

//  int imin=0;
//  pair<int,int> minabs=minpoints[0];
//  for(int i=1; i<minpoints.size(); i++)
//  {
//    if(x(grid[minabs.second][minabs.first].first)>x(grid[minpoints[i].second][minpoints[i].first].first))
//    {
//      imin = i;
//      minabs = minpoints[i];
//    }
//  }
//  cout << x(grid[minpoints[imin].second][minpoints[imin].first].first) << endl;
//  minpoints.erase(minpoints.begin()+imin);
//  for(int i=1; i<minpoints.size(); i++)
//  {
//    if(x(grid[minabs.second][minabs.first].first)>x(grid[minpoints[i].second][minpoints[i].first].first))
//    {
//      imin = i;
//      minabs = minpoints[i];
//    }
//  }

//Imprimir maximos/minimos para ficheiro
  ofstream outfilemin;
  outfilemin.open("min.dat");
  for(int i=0; i<minpoints.size(); i++)
  {
    outfilemin << minpoints[i].first*h << "   " << minpoints[i].second*h << endl;
  }
  outfilemin.close();


//Algoritmo para descobrir maximos da solucao
  vizflag = false;
  vector<pair<int,int>> maxpoints;
  for(int i=0; i<Ntot; i++)
  {
    for(int j=0; j<Ntot; j++)
    {
      if(grid[j][i].second > 1)
      {
        int tmp1 = j+1;
        if(tmp1==Ntot)
          tmp1 = 0;
        int tmp2 = j-1;
        if(tmp2==-1)
          tmp2 = Ntot-1;
        int tmp3 = i+1;
        if(tmp3==Ntot)
          tmp3 = 0;
        int tmp4 = i-1;
        if(tmp4==-1)
          tmp4 = Ntot-1;

        pair<int,int> viz[8] = {make_pair(i,tmp1),make_pair(i,tmp2),make_pair(tmp4,tmp1),make_pair(tmp4,tmp2),make_pair(tmp3,tmp1),make_pair(tmp3,tmp2),make_pair(tmp4,j),make_pair(tmp3,j)};
        for(int k=0; k<8; k++)
          if(x(grid[viz[k].second][viz[k].first].first)>x(grid[j][i].first)) //se estiver < estou a descobrir minimos; isto nao vai nada correr mal...
          {
            vizflag = true;
            break;
          }
        if(vizflag)
        {
          vizflag=false;
          continue;
        }
        // se quiser os maximos/minimos uso este
        maxpoints.push_back(make_pair(i,j));
        //se quiser os pontos a volta dos maximos/minimos uso este
        //for(int v=0; v<8; v++)
        //  maxpoints.push_back(viz[v]);
      }
    }
  }


  double *W = new double[leng];
  for(int i=0; i<leng; i++)
    W[i]=1/x(i);


//Algoritmo de watershed + anterior
  map<int,pair<int,int>> idtocoord;
  for(int i=0; i<Ntot; i++)
  {
    for(int j=0; j<Ntot; j++)
    {
      if(grid[j][i].second>=1)
        idtocoord[grid[j][i].first] = make_pair(i,j);
    }
  }
  vector<pair<int,int>> vals;
  for(int i=0; i<leng; i++)
  {
    vals.push_back(make_pair(i,0));
  }

  for(int i=0; i<maxpoints.size(); i++)
  {
    vals[grid[maxpoints[i].second][maxpoints[i].first].first].second = i+1;
  }

  for(int i=0; i<vals.size(); i++)
  {
    bool notmin = true;
    pair<int,int> pos = idtocoord[vals[i].first];
    int posxi = pos.first;
    int posyi = pos.second;

    if(grid[pos.second][pos.first].second>=1)
    {
      while(notmin)
      {
        int posx = pos.first;
        int posy = pos.second;
        if(grid[posy][posx].second==0)
        {
          break;
        }
        pair<int,int> newpos = pos;

        int tmp1 = posy+1;
        if(tmp1==Ntot)
          tmp1 = 0;
        int tmp2 = posy-1;
        if(tmp2==-1)
          tmp2 = Ntot-1;
        int tmp3 = posx+1;
        if(tmp3==Ntot)
          tmp3 = 0;
        int tmp4 = posx-1;
        if(tmp4==-1)
          tmp4 = Ntot-1;

        pair<int,int> viz[8] = {make_pair(posx,tmp1),make_pair(posx,tmp2),make_pair(tmp4,tmp1),make_pair(tmp4,tmp2),make_pair(tmp3,tmp1),make_pair(tmp3,tmp2),make_pair(tmp4,posy),make_pair(tmp3,posy)};

        for(int j=0; j<8; j++)
        {
          if(W[grid[newpos.second][newpos.first].first]>W[grid[viz[j].second][viz[j].first].first])
            newpos = viz[j];
        }
        if(pos == newpos)
        {
          notmin = false;
        }

        if(!notmin)
	{
          for(int k=0; k<maxpoints.size(); k++)
	  {
            if(maxpoints[k]==pos)
	    {
              vals[i].second = vals[grid[posy][posx].first].second;
              break;
            }
          }
          //vals[i].second = 1000;
	  //cout << posx << "  " << posy << endl;
          break;
	}
        pos = newpos;
      }
    }
  }


  vector<tuple<int,int,double>> vale;
  for(int i=0; i<vals.size(); i++)
  {
    pair<int,int> pos = idtocoord[vals[i].first];
    int posx = pos.first;
    int posy = pos.second;

    int tmp1 = posy+1;
    if(tmp1==Ntot)
      tmp1 = 0;
    int tmp2 = posy-1;
    if(tmp2==-1)
      tmp2 = Ntot-1;
    int tmp3 = posx+1;
    if(tmp3==Ntot)
      tmp3 = 0;
    int tmp4 = posx-1;
    if(tmp4==-1)
      tmp4 = Ntot-1;

    pair<int,int> viz[8] = {make_pair(posx,tmp1),make_pair(posx,tmp2),make_pair(tmp4,tmp1),make_pair(tmp4,tmp2),make_pair(tmp3,tmp1),make_pair(tmp3,tmp2),make_pair(tmp4,posy),make_pair(tmp3,posy)};

    set<int> labs;
    labs.insert(vals[i].second);
    for(int j=0; j<8; j++)
      labs.insert(vals[grid[viz[j].second][viz[j].first].first].second);
    if(labs.size()>=2)
      vale.push_back(make_tuple(posx,posy,x(grid[posy][posx].first)));
  }



//escrever coisas em ficheiros

//Escrever solucao para ficheiro
  ofstream outfile;
  outfile.open("sol.dat");
  int contador=0;
  for(int i=0; i<Ntot; i++)
  {
    for(int j=0; j<Ntot; j++)
    {
      if(grid[j][i].second>=1)
        outfile << h*i << "   " << h*j << "   " << x[grid[j][i].first] << endl;
      else if(grid[j][i].second<=1)
        outfile << h*i << "   " << h*j << "   " << 0 << endl;
      contador++;
    }
    outfile << endl;
  }
  outfile.close();
  cout << "escrevi a solucao para um ficheiro" << endl;
  cout << contador << endl;

//Imprimir maximos/minimos para ficheiro
  ofstream outfilemax;
  outfilemax.open("max.dat");
  for(int i=0; i<maxpoints.size(); i++)
  {
    outfilemax << maxpoints[i].first*h << "   " << maxpoints[i].second*h << endl;
  }
  outfilemax.close();




//  Eigen::MatrixXd Adense;
//  Adense = Eigen::MatrixXd(A);
//
//  for(int i=0; i<leng; i++)
//  {
//    for(int j=0; j<leng; j++)
//    {
//      cout << Adense(i,j) << " ";
//    }
//    cout << "     " << b[i] << "     " << x[i] << endl;
//  }
//
//  for(int i=0; i<bleng; i++)
//  {
//    cout << x[boundary[i]] << endl;
//  }


//Escrever potencial para ficheiro
  ofstream outfilef;
  outfilef.open("pot.dat");

  for(int i=0; i<Ntot; i++)
    for(int j=0; j<Ntot; j++)
    {
      outfilef << h*i << "   " << h*j << "   " << f(h,i,j,Ntot,a,D) << endl;
    }
  outfilef.close();
  cout << "escrevi o potencial para um ficheiro" << endl;

//Escrever vale para ficheiro - ultima versao
  ofstream outfilevale;
  outfilevale.open("vale.dat");

  for(int i=0; i<vale.size(); i++)
  {
    outfilevale << get<0>(vale[i])*h << "   " << get<1>(vale[i])*h << "   " << get<2>(vale[i]) << endl;
  }
  outfilevale.close();


//Calculo de valores e vectores proprios

  SparseGenMatProd<double> op(A);
  cout << "BATATA1" << endl;
  GenEigsSolver< double, SMALLEST_MAGN, SparseGenMatProd<double> > eigs(&op, Nmodos, 2*Nmodos+100);
  cout << "BATATA2" << endl;
  eigs.init();
  cout << "BATATA3" << endl;
  int nconv = eigs.compute();
  cout << "BATATA4" << endl;
  Eigen::VectorXcd evalues;
  cout << "BATATA5" << endl;
  if(eigs.info() == SUCCESSFUL)
  {
    evalues = eigs.eigenvalues();
    cout << "fui bem sucedido" << endl;
  }
  else
    cout << "falhei na vida" << endl;


  cout << "sobrevivi a calcular os valores proprios" << endl;


  ofstream outfile_ev;
  outfile_ev.open("eigenvalues.dat");

  outfile_ev << evalues << endl;

  outfile_ev.close();

  cout << "valor proprio fundamental: " << evalues(Nmodos-1) << endl;


  Eigen::MatrixXcd evectors = eigs.eigenvectors();


  cout << "sobrevivi a calcular os vectores proprios" << endl;

  //cout << evectors(0,0).real() << endl;

  double max0 = 0.;
  double max1 = 0.;
  double max2 = 0.;
//  double max3 = 0.;
//  double max4 = 0.;
//  double max5 = 0.;
//  double max6 = 0.;
//  double max7 = 0.;
//  double max8 = 0.;
//  double max9 = 0.;

  for(int i=0; i<leng; i++)
  {
    if(abs(evectors(i,Nmodos-1).real()) > max0)
      max0 = abs(evectors(i,Nmodos-1).real());
    if(abs(evectors(i,Nmodos-2).real()) > max1)
      max1 = abs(evectors(i,Nmodos-2).real());
    if(abs(evectors(i,Nmodos-3).real()) > max2)
      max2 = abs(evectors(i,Nmodos-3).real());
//    if(abs(evectors(i,Nmodos-4).real()) > max3)
//      max3 = abs(evectors(i,Nmodos-4).real());
//    if(abs(evectors(i,Nmodos-5).real()) > max4)
//      max4 = abs(evectors(i,Nmodos-5).real());
//    if(abs(evectors(i,Nmodos-6).real()) > max5)
//      max5 = abs(evectors(i,Nmodos-6).real());
//    if(abs(evectors(i,Nmodos-7).real()) > max6)
//      max6 = abs(evectors(i,Nmodos-7).real());
//    if(abs(evectors(i,Nmodos-8).real()) > max7)
//      max7 = abs(evectors(i,Nmodos-8).real());
//    if(abs(evectors(i,Nmodos-9).real()) > max8)
//      max8 = abs(evectors(i,Nmodos-9).real());
//    if(abs(evectors(i,Nmodos-10).real()) > max9)
//      max9 = abs(evectors(i,Nmodos-10).real());

  }


  ofstream outfile_evec0;
  outfile_evec0.open("eigenvectors0.dat");
  for(int i=0; i<Ntot; i++)
  {
    for(int j=0; j<Ntot; j++)
    {
      if(grid[i][j].second>=1)
        outfile_evec0 << h*j << "   " << h*i << "   " << evectors(grid[i][j].first,Nmodos-1).real()/max0 << endl;
    }
  }

  outfile_evec0.close();

  ofstream outfile_evec1;
  outfile_evec1.open("eigenvectors1.dat");

  for(int i=0; i<Ntot; i++)
  {
    for(int j=0; j<Ntot; j++)
    {
      if(grid[i][j].second>=1)
        outfile_evec1 << h*j << "   " << h*i << "   " << evectors(grid[i][j].first,Nmodos-2).real()/max1 << endl;
    }
  }

  outfile_evec1.close();

  cout << "escrevi o segundo vector proprio" << endl;

  ofstream outfile_evec2;
  outfile_evec2.open("eigenvectors2.dat");

  for(int i=0; i<Ntot; i++)
  {
    for(int j=0; j<Ntot; j++)
    {
      if(grid[i][j].second>=1)
        outfile_evec2 << h*j << "   " << h*i << "   " << evectors(grid[i][j].first,Nmodos-3).real()/max2 << endl;
    }
  }

  outfile_evec2.close();

  cout << "escrevi o terceiro vector proprio" << endl;

/*  ofstream outfile_evec3;
  outfile_evec3.open("eigenvectors3.dat");

  for(int i=0; i<Ntot; i++)
  {
    for(int j=0; j<Ntot; j++)
    {
      if(grid[i][j].second>=1)
        outfile_evec3 << h*j << "   " << h*i << "   " << evectors(grid[i][j].first,Nmodos-4).real()/max3 << endl;
    }
  }

  outfile_evec3.close();

  cout << "escrevi o quarto vector proprio" << endl;

  ofstream outfile_evec4;
  outfile_evec4.open("eigenvectors4.dat");

  for(int i=0; i<Ntot; i++)
  {
    for(int j=0; j<Ntot; j++)
    {
      if(grid[i][j].second>=1)
        outfile_evec4 << h*j << "   " << h*i << "   " << evectors(grid[i][j].first,Nmodos-5).real()/max4 << endl;
    }
  }

  outfile_evec4.close();

  cout << "escrevi o quinto vector proprio" << endl;


  ofstream outfile_evec5;
  outfile_evec5.open("eigenvectors5.dat");
  for(int i=0; i<Ntot; i++)
  {
    for(int j=0; j<Ntot; j++)
    {
      if(grid[i][j].second>=1)
        outfile_evec5 << h*j << "   " << h*i << "   " << evectors(grid[i][j].first,Nmodos-6).real()/max5 << endl;
    }
  }

  outfile_evec5.close();


  ofstream outfile_evec6;
  outfile_evec6.open("eigenvectors6.dat");
  for(int i=0; i<Ntot; i++)
  {
    for(int j=0; j<Ntot; j++)
    {
      if(grid[i][j].second>=1)
        outfile_evec6 << h*j << "   " << h*i << "   " << evectors(grid[i][j].first,Nmodos-7).real()/max6 << endl;
    }
  }

  outfile_evec6.close();


  ofstream outfile_evec7;
  outfile_evec7.open("eigenvectors7.dat");
  for(int i=0; i<Ntot; i++)
  {
    for(int j=0; j<Ntot; j++)
    {
      if(grid[i][j].second>=1)
        outfile_evec7 << h*j << "   " << h*i << "   " << evectors(grid[i][j].first,Nmodos-8).real()/max7 << endl;
    }
  }

  outfile_evec7.close();

  ofstream outfile_evec8;
  outfile_evec8.open("eigenvectors8.dat");
  for(int i=0; i<Ntot; i++)
  {
    for(int j=0; j<Ntot; j++)
    {
      if(grid[i][j].second>=1)
        outfile_evec8 << h*j << "   " << h*i << "   " << evectors(grid[i][j].first,Nmodos-9).real()/max8 << endl;
    }
  }

  outfile_evec8.close();


  ofstream outfile_evec9;
  outfile_evec9.open("eigenvectors9.dat");
  for(int i=0; i<Ntot; i++)
  {
    for(int j=0; j<Ntot; j++)
    {
      if(grid[i][j].second>=1)
        outfile_evec9 << h*j << "   " << h*i << "   " << evectors(grid[i][j].first,Nmodos-10).real()/max9 << endl;
    }
  }

  outfile_evec9.close();*/

  for(int i=0; i<D; i++)
    delete[] a[i];
  delete[] a;

  return 0;
}
