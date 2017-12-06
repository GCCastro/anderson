#include <eigen3/Eigen/Eigenvalues>
#include <eigen3/Eigen/Eigen>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Core>

#include <SymEigsSolver.h>
#include <MatOp/SparseGenMatProd.h> 


#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <vector>
#include <tuple>
#include <random>

#include <gsl/gsl_rng.h>  //random generator
#include "rtnorm.hpp"     //truncated normal distribution

typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> T;

using namespace std;
using namespace Spectra;

double f(double h, int xi, int N, double *a, int D)
{
  int x = (int)((double)xi/((double)N)*D);
  return a[x];
}


int main()
{

  cout << "Hail Eris" << endl;

  cout << "sigma da distribuicao?" << endl;
  double sig;
  cin >> sig;

  double Vmed = 4000.;
  double Vini = 0.;
  double Vend = 8000.;
  int D = 40;
  double *a = new double[D];

  //--- GSL random init ---
  gsl_rng_env_setup();                          // Read variable environnement
  const gsl_rng_type* type = gsl_rng_default;   // Default algorithm 'twister'
  gsl_rng *gen = gsl_rng_alloc (type);          // Rand generator allocation
  gsl_rng_set(gen,time(NULL));


  for(int i=0; i<D; i++)
    a[i] = rtnorm(gen,Vini,Vend,Vmed,sig).first;
   

  //isto devia vir directamente do ficheiro
  double L = 1.;
  double h = 0.0001;
  int N = int(L/h);
  cout << N << endl;
  int Nmodos = 5;

  vector<T> coefficients;

  Eigen::VectorXd b(N);

  for(int i=1; i<N-1; i++)
  {
    T entry1(i,i,2./h/h+f(h,i,N,a,D));
    coefficients.push_back(entry1);
    T entry2(i+1,i,-1./h/h);
    coefficients.push_back(entry2);
    T entry3(i-1,i,-1./h/h);
    coefficients.push_back(entry3);
    b[i] = 1.;
  }
  T entry00(0,0,2./h/h+f(h,0,N,a,D));
  coefficients.push_back(entry00);
  T entry01(1,0,-1./h/h);
  coefficients.push_back(entry01);
  T entryN0(N-1,N-1,2./h/h+f(h,N-1,N,a,D));
  coefficients.push_back(entryN0);
  T entryN1(N-2,N-1,-1./h/h);
  coefficients.push_back(entryN1);
  b[0] = 0;
  b[N-1] = 0;

  SpMat A(N,N);
  A.setFromTriplets(coefficients.begin(),coefficients.end());
//  for(int i=0; i<N; i++)
//  {
//    for(int j=0; j<N; j++)
//    {
//      cout << A.coeff(i,j) << " ";
//    }
//    cout << endl;
//  }
  Eigen::SparseLU<SpMat> solver(A);
  Eigen::VectorXd x = solver.solve(b);
  cout << "sobrevivi a calcular o u" << endl;



  SparseGenMatProd<double> op(A);
  SymEigsSolver< double, SMALLEST_MAGN, SparseGenMatProd<double> > eigs(&op, Nmodos, 100);
  eigs.init();
  int nconv = eigs.compute();
  Eigen::VectorXcd evalues;
  if(eigs.info() == SUCCESSFUL)
  {
    evalues = eigs.eigenvalues();
    cout << "fui bem sucedido" << endl;
  }
  else
    cout << "falhei na vida" << endl;

  cout << "sobrevivi a calcular os valores proprios" << endl;


  ofstream outfile_ev;
  outfile_ev.open("eigenvalues_1d.dat");
  outfile_ev << evalues << endl;
  outfile_ev.close();

  Eigen::MatrixXcd evectors = eigs.eigenvectors();
  cout << "sobrevivi a calcular os vectores proprios" << endl;

  double agmon[Nmodos][N];

  for(int i=0; i<Nmodos; i++)
  {
    double max=0;
    int imax=0;
    for(int j=0; j<N; j++)
    {
      if(abs(evectors(j,Nmodos-(i+1)).real()) > max)
      {
        max = abs(evectors(j,Nmodos-(i+1)).real());
        imax = j;
      }
    }
    double ximax = 0;
    for(int k=1; k<N-1; k++)
    {
      ximax += 1./2.*(k*(evectors(k,Nmodos-(i+1)).real())*(evectors(k,Nmodos-(i+1)).real())+(k+1)*(evectors(k+1,Nmodos-(i+1)).real())*(evectors(k+1,Nmodos-(i+1)).real()));
    }
    cout << "maximo antes" << endl;
    cout << imax << endl;
    imax = (int)round(ximax);
    cout << "maximo depois" << endl;
    cout << imax << endl;
    for(int k=1; k<N-1; k++)
    {
      double rho = 0.;
      if(k==imax)
      {
        agmon[i][k]=1.;
        continue;
      }
      int step = (k-imax)/abs(k-imax);
      for(int l=imax; l!=k-1; l+=step)
      {
        double arg1 = 1./x(l)-evalues(Nmodos-(i+1)).real();
        double arg2 = 1./x(l+step)-evalues(Nmodos-(i+1)).real();

        if(arg1<0.)
          arg1=0.;
        if(arg2<0.)
          arg2=0.;
        rho += h/2.*(sqrt(arg1)+sqrt(arg2));
      }
      agmon[i][k] = exp(-rho);
    }

    ofstream outfile_evec;
    string file_vec = string("eigenvectors") + to_string(i);
    file_vec += string("_1d.dat");
    outfile_evec.open(file_vec.c_str());
    for(int jj=0; jj<N; jj++)
    {
      outfile_evec << h*jj << "   " << evectors(jj,Nmodos-(i+1)).real()/max << endl;
    }
    outfile_evec.close();
    cout << "escrevi o " << i << "º vector proprio" << endl;
    ofstream outfile_agmon;
    string file_agmon = string("agmon") + to_string(i);
    file_agmon += string("_1d.dat");
    outfile_agmon.open(file_agmon.c_str());
    for(int jj=0; jj<N; jj++)
    {
      outfile_agmon << h*jj << "   " << agmon[i][jj] << endl;
    }
    outfile_agmon.close();
    cout << "escrevi a " << i << "ª estimativa de Agmon" << endl;
  }


//Escrever solucao para ficheiro
  ofstream outfile;
  outfile.open("sol1D.dat");
  for(int i=0; i<N; i++)
  {
    outfile << h*i << "   " << x(i) << endl;
  }
  outfile.close();
  cout << "escrevi a solucao para um ficheiro" << endl;


  ofstream outfilepot;
  outfilepot.open("pot1D.dat");
  for(int i=0; i<N; i++)
  {
    outfilepot << h*i << "  " << f(h,i,N,a,D) << endl;
  }
  outfilepot.close();
  cout << "escrevi o potencial para um ficheiro" << endl;



  gsl_rng_free(gen);                            // GSL rand generator deallocation

  return 0;
}
