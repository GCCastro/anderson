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

  srand (time(NULL));
  double VMax = 8000.;
  int D = 20;
  double *a = new double[D];

  for(int i=0; i<D; i++)
    a[i] = double(rand())/RAND_MAX * VMax;


  //isto devia vir directamente do ficheiro
  double L = 1.;
  double h = 0.0009;
  int N = int(L/h);
  cout << N << endl;
  int Nmodos = 5;

  vector<T> coefficients;

  Eigen::VectorXd b(N);

  for(int i=1; i<N-1; i++)
  {
    T entry1(i,i,2+h*h*f(h,i,N,a,D));
    coefficients.push_back(entry1);
    T entry2(i+1,i,-1);
    coefficients.push_back(entry2);
    T entry3(i-1,i,-1);
    coefficients.push_back(entry3);
    b[i] = h*h;
  }
  T entry00(0,0,2+h*h*f(h,0,N,a,D));
  coefficients.push_back(entry00);
  T entry01(1,0,-1);
  coefficients.push_back(entry01);
  T entryN0(N-1,N-1,2+h*h*f(h,N-1,N,a,D));
  coefficients.push_back(entryN0);
  T entryN1(N-2,N-1,-1);
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

  double integral=1.;
//  integral = 0.;
//  for(int i=1; i<N; i++)
//  {
//    integral+=(evectors(i,1).real()+evectors(i-1,1).real())*h/2;
//  }
//  cout << "integral=" << integral << endl;


  ofstream outfile_evec0;
  outfile_evec0.open("eigenvectors0_1d.dat");
  for(int i=0; i<N; i++)
  {
    outfile_evec0 << h*i << "   " << evectors(i,Nmodos-1).real()/integral << endl;
  }
  outfile_evec0.close();
  ofstream outfile_evec1;
  outfile_evec1.open("eigenvectors1_1d.dat");
  for(int i=0; i<N; i++)
  {
    outfile_evec1 << h*i << "   " << evectors(i,Nmodos-2).real() << endl;
  }
  outfile_evec1.close();
  ofstream outfile_evec2;
  outfile_evec2.open("eigenvectors2_1d.dat");
  for(int i=0; i<N; i++)
  {
    outfile_evec2 << h*i << "   " << evectors(i,Nmodos-3).real() << endl;
  }
  outfile_evec2.close();

  ofstream outfile_evec3;
  outfile_evec3.open("eigenvectors3_1d.dat");
  for(int i=0; i<N; i++)
  {
    outfile_evec3 << h*i << "   " << evectors(i,Nmodos-4).real() << endl;
  }
  outfile_evec3.close();

  ofstream outfile_evec4;
  outfile_evec4.open("eigenvectors4_1d.dat");
  for(int i=0; i<N; i++)
  {
    outfile_evec4 << h*i << "   " << evectors(i,Nmodos-5).real() << endl;
  }
  outfile_evec4.close();


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

}