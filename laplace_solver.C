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
  double L = 1.6;
  double h = 0.0039;
  int Ntot = int(L/h);

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


  int Nmodos = 2;

  vector<T> coefficients;

  //double *boundary = new double[bleng];

  Eigen::VectorXd b(leng);

  for(int i=0;i<Ntot;i++) //loop ao longo das colunas
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

//  Eigen::SimplicialCholesky<SpMat> chol(A);
//  Eigen::VectorXd x = chol.solve(b);


  Eigen::SparseLU<SpMat> solver(A);
  Eigen::VectorXd x = solver.solve(b);

  cout << "sobrevivi a calcular o u" << endl;


//Algoritmo para descobrir maximos/minimos da solucao (nao adaptado ao metodo geral)
  /*bflag = false;
  bool vizflag = false;
  vector<double> maxpoints;
  for(int i=0; i<leng; i++)
  {
    for(int j=0; j<bleng; j++)
      if(i==boundary[j])
      {
        bflag = true;
      }
    if(bflag==true)
    {
      bflag = false;
      continue;
    }

    int viz[8] = {i-N-1,i-N,i-N+1,i-1,i+N-1,i+N,i+N+1,i+1};
    for(int j=0; j<8; j++)
      if(x(viz[j])<x(i)) //se estiver < estou a descobrir minimos; isto nao vai nada correr mal...
      {
        vizflag = true;
        break;
      }
    if(vizflag)
    {
      vizflag=false;
      continue;
    }
    // se quiser os maximos uso este
    //maxpoints.push_back(i);
    //se quiser os ponts a volta do maximo uso este
    for(int v=0; v<8; v++)
      maxpoints.push_back(viz[v]);
  }*/


//Algoritmo para descobrir o sistema de vales da solucao (nao adaptado ao metodo geral)
  /*int *vale = new int[leng];
  for(int i=0; i<leng; i++)
    vale[i] = 0;

  for(int i=0; i<maxpoints.size(); i++)
  {
    bool notmin = true;
    int pos = maxpoints[i];
    while(notmin)
    {
      cout << pos << endl;
      int newpos = pos;
      int viz[8] = {pos-1,pos+1,pos-N,pos+N,pos-1-N,pos+1+N,pos-N+1,pos+N-1};
      for(int j=0; j<8; j++)
      {
        if(x(newpos)>x(viz[j]))
          newpos = viz[j];
      }
      if(pos == newpos)
      {
        notmin = false;
        break;
      }

      for(int k=0; k<bleng; k++)
        if(newpos==boundary[k])
        {
          notmin = false;
          break;
        }
      if(!notmin)
        break;
      pos = newpos;
      vale[pos] += 1;
    }
  }*/


//escrever coisas em ficheiros

//Escrever solucao para ficheiro
  ofstream outfile;
  outfile.open("sol.dat");
  for(int i=0; i<Ntot; i++)
  {
    for(int j=0; j<Ntot; j++)
    {
      if(grid[i][j].second>=1)
        outfile << h*j << "   " << h*i << "   " << x[grid[i][j].first] << endl;
    }
  }
  outfile.close();
  cout << "escrevi a solucao para um ficheiro" << endl;

//Imprimir maximos/minimos para ficheiro
  /*ofstream outfilemax;
  outfilemax.open("max.dat");
  for(int i=0; i<maxpoints.size(); i++)
  {

    div_t coord;
    coord = div(maxpoints[i],N);

    outfilemax << h*coord.rem << "   " << h*coord.quot << endl;
  }
  outfilemax.close();*/

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
  /*ofstream outfilevale;
  outfilevale.open("vale.dat");

  for(int i=0; i<leng; i++)
  {
    div_t coord;
    coord = div(i,N);
    if(vale[i] > 1)
      outfilevale << h*coord.rem << "   " << h*coord.quot << "   " << vale[i] << endl;
  }
  outfilevale.close();*/



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

  ofstream outfile_evec1;
  outfile_evec1.open("eigenvectors0.dat");


  for(int i=0; i<Ntot; i++)
  {
    for(int j=0; j<Ntot; j++)
    {
      if(grid[i][j].second>=1)
        outfile_evec1 << h*j << "   " << h*i << "   " << evectors(grid[i][j].first,1).real() << endl;
    }
  }

  outfile_evec1.close();


  ofstream outfile_evec2;
  outfile_evec2.open("eigenvectors1.dat");

  for(int i=0; i<Ntot; i++)
  {
    for(int j=0; j<Ntot; j++)
    {
      if(grid[i][j].second>=1)
        outfile_evec1 << h*j << "   " << h*i << "   " << evectors(grid[i][j].first,0).real() << endl;
    }
  }

  outfile_evec2.close();




  for(int i=0; i<D; i++)
    delete[] a[i];
  delete[] a;

  return 0;
}
