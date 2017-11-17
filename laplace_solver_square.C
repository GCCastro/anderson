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
      a[i][j] = rand() % VMax;

  double L = 1.;
  double h = 0.09;
  int N = (int)L/h;
  int leng = N*N;
  int bleng = 4*N-4;
  int Nmodos = 2;

  vector<T> coefficients;

  //double *grid = new double[leng];
  double *boundary = new double[bleng];

  Eigen::VectorXd b(leng);

  for(int i=0;i<N;i++)
    boundary[i] = i;
  for(int i=0;i<N-2;i++)
  {
    boundary[2*i+N] = (i+1)*N;
    boundary[2*i+N+1] = (i+2)*N-1;
  }
  for(int i=0;i<N;i++)
    boundary[3*N-4+i] = N*(N-1)+i;

  for(int i=0;i<4*N-4;i++)
    cout << boundary[i] << endl;



  int n = 0;
  bool bflag = false;
  for(int i=0;i<leng;i++)
  {
    for(int j=0;j<bleng;j++)
    {
      if(i==boundary[j])
      {
        bflag = true;
        T entry(i,i,1);
        coefficients.push_back(entry);
        b[i] = 0;
      }
    }
    if(bflag==true)
    {
      bflag = false;
      continue;
    }

    div_t coord;
    coord = div(i,N);

    T entry1(i,i,4+h*h*f(h,coord.rem,coord.quot,N,a,D));
    T entry2(i,i+1,-1);
    T entry3(i,i-1,-1);
    T entry4(i,i+N,-1);
    T entry5(i,i-N,-1);

    coefficients.push_back(entry1);
    coefficients.push_back(entry2);
    coefficients.push_back(entry3);
    coefficients.push_back(entry4);
    coefficients.push_back(entry5);


    b[i] = h*h;
//    if (i%(N/10)==0)
//      cout << i << endl;
  }

  SpMat A(leng,leng);
  A.setFromTriplets(coefficients.begin(),coefficients.end());

//  Eigen::SimplicialCholesky<SpMat> chol(A);
//  Eigen::VectorXd x = chol.solve(b);

  Eigen::SparseLU<SpMat> solver(A);
  cout << solver.determinant() << endl;
  Eigen::VectorXd x = solver.solve(b);

  cout << "sobrevivi a calcular o u" << endl;

  bflag = false;
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
  }


  int *vale = new int[leng];
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
  }


//escrever coisas em ficheiros
  ofstream outfile;
  ofstream outfilemax;
  outfile.open("sol.dat");
  outfilemax.open("max.dat");
  for(int i=0; i<leng; i++)
  {

    div_t coord;
    coord = div(i,N);

    outfile << h*coord.rem << "   " << h*coord.quot << "   " << x[i] << endl;
  }
  for(int i=0; i<maxpoints.size(); i++)
  {

    div_t coord;
    coord = div(maxpoints[i],N);

    outfilemax << h*coord.rem << "   " << h*coord.quot << endl;
  }

  outfile.close();
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

  ofstream outfilef;
  outfilef.open("pot.dat");

  for(int i=0; i<D; i++)
    for(int j=0; j<D; j++)
    {
      outfilef << h*i << "   " << h*j << "   " << a[i][j] << endl;
    }
  outfilef.close();


  ofstream outfilevale;
  outfilevale.open("vale.dat");

  for(int i=0; i<leng; i++)
  {
    div_t coord;
    coord = div(i,N);
    if(vale[i] > 1)
      outfilevale << h*coord.rem << "   " << h*coord.quot << "   " << vale[i] << endl;
  }
  outfilevale.close();



//Calculo de valores e vectores proprios

  /*SparseGenMatProd<double> op(A);
  GenEigsSolver< double, LARGEST_MAGN, SparseGenMatProd<double> > eigs(&op, Nmodos, 2*Nmodos+3);
  eigs.init();
  int nconv = eigs.compute();

  Eigen::VectorXcd evalues;
  if(eigs.info() == SUCCESSFUL)
    evalues = eigs.eigenvalues();

  cout << "sobrevivi a calcular os valores proprios" << endl;


  Eigen::MatrixXcd evectors = eigs.eigenvectors();


  cout << "sobrevivi a calcular os vectores proprios" << endl;

  ofstream outfile_ev;
  outfile_ev.open("eigenvalues.dat");

  outfile_ev << evalues << endl;

  outfile_ev.close();

  cout << "valor proprio fundamental: " << evalues(1) << endl;
  cout << "valor proprio 1 estado excitado: " << evalues(0) << endl;


  //cout << evectors(0,0).real() << endl;

  ofstream outfile_evec1;
  outfile_evec1.open("eigenvectors0.dat");

  for(int i=0; i<leng; i++)
  {
    div_t coord;
    coord = div(i,N);

    outfile_evec1 << h*coord.rem << "   " << h*coord.quot << "   " << evectors(i,1).real() << endl;
  }

  outfile_evec1.close();


  ofstream outfile_evec2;
  outfile_evec2.open("eigenvectors1.dat");

  for(int i=0; i<leng; i++)
  {
    div_t coord;
    coord = div(i,N);

    outfile_evec2 << h*coord.rem << "   " << h*coord.quot << "   " << evectors(i,0).real() << endl;
  }

  outfile_evec2.close();
*/



  delete[] boundary;
  for(int i=0; i<D; i++)
    delete[] a[i];
  delete[] a;

  return 0;
}
