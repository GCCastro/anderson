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
typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> T;

using namespace std;
using namespace Spectra;

double f(double h, int xi, int yj, int N, int **a, int D)
{
  int x = (int)((double)xi/((double)N)*D);
  int y = (int)((double)yj/((double)N)*D);

  return a[x][y];
}

int main()
{


  cout << "Hail Eris" << endl;

  srand (time(NULL));
  int VMax = 8000;
  int D = 20;
  int **a = new int*[D];
  for(int i=0; i<D; i++)
    a[i] = new int[D];

  for(int i=0; i<D; i++)
    for(int j=0; j<D; j++)
    {
      a[i][j] = double(rand())/RAND_MAX * VMax;
    }


  //isto devia vir directamente do ficheiro
  double L = 1.7;
  double h = 0.0049;
  int Ntot = int(L/h);
  cout << Ntot << endl;

  cout << "ficheiro de regiao?" << endl;
  string filename;
  cin >> filename;

  ifstream infile;
  infile.open(filename.c_str());
  vector<vector<pair<int,int>>> grid;
  vector<pair<int,int>> line;
  double xf,yf;
  int c;
  int leng=0;
  double xold=0.;


  while(infile >> xf >> yf >> c)
  {
    if(xf<xold)
    {
      grid.push_back(line);
      line.clear();
    }
    if(c==2)
    {
      pair<int,int> entry(leng,c);
      line.push_back(entry);
      leng++;
    }
    else if(c<=1)
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


  int Nmodos = 5;

  vector<T> coefficients;

  //double *boundary = new double[bleng];

  Eigen::VectorXd b(leng);

  /*for(int i=0;i<Ntot;i++) //loop ao longo das colunas
  {
    for(int j=0;j<Ntot;j++) //loop dentro da linha
    {
      if(grid[i][j].second==1)
      {
        T entry(grid[i][j].first,grid[i][j].first,1);
        coefficients.push_back(entry);
        b[grid[i][j].first] = 0;
      }
      else if(grid[i][j].second==2)
      {
        T entry1(grid[i][j].first,grid[i][j].first,4+h*h*f(h,j,i,Ntot,a,D));
        T entry2(grid[i][j].first,grid[i][j+1].first,-1);
        T entry3(grid[i][j].first,grid[i][j-1].first,-1);
        T entry4(grid[i][j].first,grid[i+1][j].first,-1);
        T entry5(grid[i][j].first,grid[i-1][j].first,-1);
  
        coefficients.push_back(entry1);
        coefficients.push_back(entry2);
        coefficients.push_back(entry3);
        coefficients.push_back(entry4);
        coefficients.push_back(entry5);
  
        b[grid[i][j].first] = h*h;

      }
    }
  }*/

  for(int i=0;i<Ntot;i++) //loop ao longo das colunas
  {
    for(int j=0;j<Ntot;j++) //loop dentro da linha
    {
      if(grid[i][j].second==2)
      {
        T entry1(grid[i][j].first,grid[i][j].first,4+h*h*f(h,j,i,Ntot,a,D));
        coefficients.push_back(entry1);
        if(grid[i][j+1].second==2)
	{
          T entry2(grid[i][j].first,grid[i][j+1].first,-1);
          coefficients.push_back(entry2);
	}
        if(grid[i][j-1].second==2)
	{
          T entry3(grid[i][j].first,grid[i][j-1].first,-1);
          coefficients.push_back(entry3);
	}
        if(grid[i+1][j].second==2)
	{
          T entry4(grid[i][j].first,grid[i+1][j].first,-1);
          coefficients.push_back(entry4);
	}
        if(grid[i-1][j].second==2)
	{
          T entry5(grid[i][j].first,grid[i-1][j].first,-1);
          coefficients.push_back(entry5);
	}
  
        b[grid[i][j].first] = h*h;

      }
    }
  }

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

//  Eigen::SimplicialCholesky<SpMat> chol(A);
//  Eigen::VectorXd x = chol.solve(b);

  Eigen::SparseLU<SpMat> solver(A);

  Eigen::VectorXd x = solver.solve(b);

  cout << "sobrevivi a calcular o u" << endl;


//Algoritmo para descobrir maximos/minimos da solucao
  bool vizflag = false;
  vector<pair<int,int>> maxpoints;
  for(int i=0; i<Ntot; i++)
  {
    for(int j=0; j<Ntot; j++)
    {
      if(grid[j][i].second > 1)
      {
        pair<int,int> viz[8] = {make_pair(i,j+1),make_pair(i,j-1),make_pair(i-1,j+1),make_pair(i-1,j-1),make_pair(i+1,j+1),make_pair(i+1,j-1),make_pair(i-1,j),make_pair(i+1,j)};
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
        //maxpoints.push_back(i);
        //se quiser os pontos a volta dos maximos/minimos uso este
        for(int v=0; v<8; v++)
          maxpoints.push_back(viz[v]);
      }
    }
  }


//Algoritmo para descobrir o sistema de vales da solucao; o algoritmo percorre os percursos de maior declive desde a vizinhanca dos maximos ate a minimos locais
/*  int vale[Ntot][Ntot];

  for(int i=0; i<Ntot; i++)
    for(int j=0; j<Ntot; j++)
      vale[i][j] = 0;

  for(int i=0; i<maxpoints.size(); i++)
  {
    bool notmin = true;
    pair<int,int> pos = maxpoints[i];
    if(grid[pos.second][pos.first].second>=1)
    {
      while(notmin)
      {
        int posx = pos.first;
        int posy = pos.second;
	if(grid[posy][posx].second!=2)
	{
          break;
	}
        pair<int,int> newpos = pos;

        pair<int,int> viz[8] = {make_pair(posx,posy+1),make_pair(posx,posy-1),make_pair(posx-1,posy+1),make_pair(posx-1,posy-1),make_pair(posx+1,posy+1),make_pair(posx+1,posy-1),make_pair(posx-1,posy),make_pair(posx+1,posy)};

        for(int j=0; j<8; j++)
        {
          if(x(grid[newpos.second][newpos.first].first)>x(grid[viz[j].second][viz[j].first].first))
            newpos = viz[j];
        }
        if(pos == newpos)
        {
          notmin = false;
          break;
        }

        if(!notmin)
          break;
        pos = newpos;
        vale[pos.first][pos.second] += 1;
      }
    }
  }*/


//Algoritmo para vale versao 2: fazer o mesmo mas para todos os pontos com valores nao nulos
  int vale[Ntot][Ntot];
  double *W = new double[leng];
  for(int i=0; i<leng; i++)
    W[i]=1/x(i);


  for(int i=0; i<Ntot; i++)
    for(int j=0; j<Ntot; j++)
      vale[i][j] = 0;

  for(int i=0; i<Ntot; i++)
  {
    for(int k=0; k<Ntot; k++)
    {
      if(grid[k][i].second==2)
      {
        bool notmin = true;
        pair<int,int> pos = make_pair(i,k);
        if(grid[pos.second][pos.first].second>=1)
        {
          while(notmin)
          {
            int posx = pos.first;
            int posy = pos.second;
          if(grid[posy][posx].second!=2)
          {
              break;
          }
            pair<int,int> newpos = pos;
  
            pair<int,int> viz[8] = {make_pair(posx,posy+1),make_pair(posx,posy-1),make_pair(posx-1,posy+1),make_pair(posx-1,posy-1),make_pair(posx+1,posy+1),make_pair(posx+1,posy-1),make_pair(posx-1,posy),make_pair(posx+1,posy)};
  
            for(int j=0; j<8; j++)
            {
              if(W[grid[newpos.second][newpos.first].first]>W[grid[viz[j].second][viz[j].first].first])
                newpos = viz[j];
            }
            if(pos == newpos)
            {
              notmin = false;
              break;
            }
  
            if(!notmin)
              break;
            pos = newpos;
            vale[pos.first][pos.second] += 1;
          }
        }
      }
    }
  }




//Algoritmo de \"watershed\"


//  double *W = new double(leng);
//  for(int i=0; i<leng; i++)
//    W[i]=x(i);

  


//escrever coisas em ficheiros

//Escrever solucao para ficheiro
  ofstream outfile;
  outfile.open("sol.dat");
  int contador=0;
  for(int i=0; i<Ntot; i++)
  {
    for(int j=0; j<Ntot; j++)
    {
      if(grid[j][i].second==2)
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

//Escrever vale para ficheiro
  ofstream outfilevale;
  outfilevale.open("vale.dat");

  for(int i=0; i<Ntot; i++)
    for(int j=0; j<Ntot; j++)
    {
      if(vale[i][j] > 1)
        outfilevale << double(i)*h << "   " << double(j)*h << "   " << vale[i][j] << endl;
    }
  outfilevale.close();



//Calculo de valores e vectores proprios

  SparseGenMatProd<double> op(A);
  cout << "BATATA1" << endl;
  GenEigsSolver< double, LARGEST_MAGN, SparseGenMatProd<double> > eigs(&op, Nmodos, 2*Nmodos+3);
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

  cout << "valor proprio 1 estado excitado: " << evalues(0) << endl;
  cout << "valor proprio fundamental: " << evalues(1) << endl;

  Eigen::MatrixXcd evectors = eigs.eigenvectors();


  cout << "sobrevivi a calcular os vectores proprios" << endl;

  //cout << evectors(0,0).real() << endl;

  ofstream outfile_evec0;
  outfile_evec0.open("eigenvectors0.dat");


  for(int i=0; i<Ntot; i++)
  {
    for(int j=0; j<Ntot; j++)
    {
      if(grid[i][j].second>=1)
        outfile_evec0 << h*j << "   " << h*i << "   " << evectors(grid[i][j].first,Nmodos-1).real() << endl;
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
        outfile_evec1 << h*j << "   " << h*i << "   " << evectors(grid[i][j].first,Nmodos-2).real() << endl;
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
        outfile_evec2 << h*j << "   " << h*i << "   " << evectors(grid[i][j].first,Nmodos-3).real() << endl;
    }
  }

  outfile_evec2.close();

  cout << "escrevi o terceiro vector proprio" << endl;


  ofstream outfile_evec3;
  outfile_evec3.open("eigenvectors3.dat");

  for(int i=0; i<Ntot; i++)
  {
    for(int j=0; j<Ntot; j++)
    {
      if(grid[i][j].second>=1)
        outfile_evec3 << h*j << "   " << h*i << "   " << evectors(grid[i][j].first,Nmodos-4).real() << endl;
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
        outfile_evec4 << h*j << "   " << h*i << "   " << evectors(grid[i][j].first,Nmodos-5).real() << endl;
    }
  }

  outfile_evec4.close();

  cout << "escrevi o quinto vector proprio" << endl;



  for(int i=0; i<D; i++)
    delete[] a[i];
  delete[] a;

  return 0;
}
