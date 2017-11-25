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
typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> T;

using namespace std;
using namespace Spectra;

double f(double h, int xi, int yj, int N, int **a, int D)
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

  cout << "ficheiro da regiao?" << endl;
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
        T entry1(grid[i][j].first,grid[i][j].first,4./h/h+f(h,j,i,Ntot,a,D));
        coefficients.push_back(entry1);
        if(grid[i][j+1].second==2)
	{
          T entry2(grid[i][j].first,grid[i][j+1].first,-1./h/h);
          coefficients.push_back(entry2);
	}
        if(grid[i][j-1].second==2)
	{
          T entry3(grid[i][j].first,grid[i][j-1].first,-1./h/h);
          coefficients.push_back(entry3);
	}
        if(grid[i+1][j].second==2)
	{
          T entry4(grid[i][j].first,grid[i+1][j].first,-1./h/h);
          coefficients.push_back(entry4);
	}
        if(grid[i-1][j].second==2)
	{
          T entry5(grid[i][j].first,grid[i-1][j].first,-1./h/h);
          coefficients.push_back(entry5);
	}
  
        b[grid[i][j].first] = 1.;


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


//Algoritmo para descobrir minimos da solucao
  bool vizflag = false;
  vector<pair<int,int>> minpoints;
  for(int i=0; i<Ntot; i++)
  {
    for(int j=0; j<Ntot; j++)
    {
      if(grid[j][i].second > 1)
      {
        pair<int,int> viz[8] = {make_pair(i,j+1),make_pair(i,j-1),make_pair(i-1,j+1),make_pair(i-1,j-1),make_pair(i+1,j+1),make_pair(i+1,j-1),make_pair(i-1,j),make_pair(i+1,j)};
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

//  for(int i=0; i<minpoints.size(); i++)
//  {
//    cout << x(grid[minpoints[i].second][minpoints[i].first].first) << endl;
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
        maxpoints.push_back(make_pair(i,j));
        //se quiser os pontos a volta dos maximos/minimos uso este
        //for(int v=0; v<8; v++)
        //  maxpoints.push_back(viz[v]);
      }
    }
  }


//Algoritmo para descobrir o sistema de vales da solucao; o algoritmo percorre os percursos de maior declive desde a vizinhanca dos maximos ate a minimos locais
  /*int vale[Ntot][Ntot];

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
  /*cout << "o algoritmo de watershed comeca aqui" << endl;

//  double *W = new double(leng);  //Feito mais acima
//  for(int i=0; i<leng; i++)
//    W[i]=x(i);

  map<int,pair<int,int>> idtocoord;
  for(int i=0; i<Ntot; i++)
  {
    for(int j=0; j<Ntot; j++)
    {
      if(grid[j][i].second==2)
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
  sort(vals.begin(),vals.end(),[&W](pair<int,int> pair1,pair<int,int> pair2) -> bool {return W[pair2.first] > W[pair1.first];});

  bool minflag = false;

  for(int i=0; i<vals.size(); i++)
  {
    int id = vals[i].first;
    int coordx = idtocoord[id].first;
    int coordy = idtocoord[id].second;
    for(int j=0; j<maxpoints.size(); j++)
    {
      if(id == grid[maxpoints[j].second][maxpoints[j].first].first)
      {
        minflag = true;
        break;
      }
    }
    if(minflag)
    {
      cout << "BATATA MINIMA" << endl;
      minflag=false;
      continue;
    }
    pair<int,int> vizs[8] = {make_pair(coordx+1,coordy+1),make_pair(coordx,coordy+1),make_pair(coordx-1,coordy+1),make_pair(coordx+1,coordy),make_pair(coordx-1,coordy),make_pair(coordx+1,coordy-1),make_pair(coordx,coordy-1),make_pair(coordx-1,coordy-1)};

    vector<int> labels;
    bool lab_flag = false;
    for(int j=0; j<8; j++)
    {
      if(grid[vizs[j].second][vizs[j].first].second==2)
      {
        int lab = vals[grid[vizs[j].second][vizs[j].first].first].second;
        if(lab==0 || lab==1000){continue;}
        for(int k=0; k<labels.size(); k++)
        {
          if(labels[k]==lab)
          {
            lab_flag=true;
            break;
          }
        }
        if(lab_flag)
        {
          lab_flag = false;
          continue;
        }
        else
          labels.push_back(lab);
      }
    }
    if(labels.size()>=2)
    {
      vals[i].second=0;
    }
    else if(labels.empty())
    {
      vals[i].second=1000;
    }
    else
    {
      //cout << "BATATA  " << labels[0] << endl;
      vals[i].second=labels[0];
    }
    //if(i<50)
    //{
    //  cout << coordx << "  " << coordy << "  " << W[id] << "  " << vals[i].second << endl;
    //}

  }
  cout << "o algoritmo de watershed acaba aqui" << endl;*/









//Algoritmo de watershed + anterior
  map<int,pair<int,int>> idtocoord;
  for(int i=0; i<Ntot; i++)
  {
    for(int j=0; j<Ntot; j++)
    {
      if(grid[j][i].second==2)
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

    if(grid[pos.second][pos.first].second==2)
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


  vector<tuple<int,int,double>> valley;
  for(int i=0; i<vals.size(); i++)
  {
    pair<int,int> pos = idtocoord[vals[i].first];
    int posx = pos.first;
    int posy = pos.second;
    pair<int,int> viz[8] = {make_pair(posx,posy+1),make_pair(posx,posy-1),make_pair(posx-1,posy+1),make_pair(posx-1,posy-1),make_pair(posx+1,posy+1),make_pair(posx+1,posy-1),make_pair(posx-1,posy),make_pair(posx+1,posy)};
    set<int> labs;
    labs.insert(vals[i].second);
    for(int j=0; j<8; j++)
      labs.insert(vals[grid[viz[j].second][viz[j].first].first].second);
    if(labs.size()>=2)
      valley.push_back(make_tuple(posx,posy,x(grid[posy][posx].first)));
  }



//o outro algoritmo para o vale
  vector<tuple<int,int,double>> vale3;
  for(int i=2; i<Ntot-2; i++)
  {
    for(int j=2; j<Ntot-2; j++)
    {
      if(x(grid[j][i].first)<x(grid[j-1][i].first) && x(grid[j][i].first)<x(grid[j+1][i].first))
        vale3.push_back(make_tuple(i,j,x(grid[j][i].first)));
      if(x(grid[j][i].first)<x(grid[j][i-1].first) && x(grid[j][i].first)<x(grid[j][i+1].first))
        vale3.push_back(make_tuple(i,j,x(grid[j][i].first)));
      if(x(grid[j][i].first)<x(grid[j-1][i-1].first) && x(grid[j][i].first)<x(grid[j+1][i+1].first))
        vale3.push_back(make_tuple(i,j,x(grid[j][i].first)));
      if(x(grid[j][i].first)<x(grid[j-1][i+1].first) && x(grid[j][i].first)<x(grid[j+1][i-1].first))
        vale3.push_back(make_tuple(i,j,x(grid[j][i].first)));

      if(x(grid[j][i].first)<x(grid[j-2][i].first) && x(grid[j][i].first)<x(grid[j+2][i].first))
        vale3.push_back(make_tuple(i,j,x(grid[j][i].first)));
      if(x(grid[j][i].first)<x(grid[j][i-2].first) && x(grid[j][i].first)<x(grid[j][i+2].first))
        vale3.push_back(make_tuple(i,j,x(grid[j][i].first)));

      if(x(grid[j][i].first)<x(grid[j-2][i-1].first) && x(grid[j][i].first)<x(grid[j+2][i+1].first))
        vale3.push_back(make_tuple(i,j,x(grid[j][i].first)));
      if(x(grid[j][i].first)<x(grid[j-1][i-2].first) && x(grid[j][i].first)<x(grid[j+1][i+2].first))
        vale3.push_back(make_tuple(i,j,x(grid[j][i].first)));
      if(x(grid[j][i].first)<x(grid[j+2][i-1].first) && x(grid[j][i].first)<x(grid[j-2][i+1].first))
        vale3.push_back(make_tuple(i,j,x(grid[j][i].first)));
      if(x(grid[j][i].first)<x(grid[j-1][i+2].first) && x(grid[j][i].first)<x(grid[j+1][i-2].first))
        vale3.push_back(make_tuple(i,j,x(grid[j][i].first)));


    }
  }

//Ainda outro algoritmo par ao vale -> atraves de derivadas
  vector<tuple<int,int,double>> vale4;
  for(int i=1; i<Ntot-1; i++)
  {
    for(int j=1; j<Ntot-1; j++)
    {
      if((abs((x(grid[j+1][i].first)-x(grid[j-1][i].first))/x(grid[j][i].first))<0.05 || abs((x(grid[j][i+1].first)-x(grid[j][i-1].first))/x(grid[j][i].first))<0.05) && !((x(grid[j+1][i].first)<x(grid[j][i].first) && x(grid[j-1][i].first)<x(grid[j][i].first)) || (x(grid[j][i-1].first)<x(grid[j][i].first) && x(grid[j][i-1].first)<x(grid[j][i].first))))
        vale4.push_back(make_tuple(i,j,x(grid[j][i].first)));
    }
  }


//  for(int j=0; j<Ntot; j++)
//  {
//    for(int i=1; i<Ntot-1; i++)
//    {
//      if(x(grid[j][i].first)<x(grid[j][i-1].first) && x(grid[j][i].first)<x(grid[j][i+1].first))
//        vale3.push_back(make_tuple(i,j,x(grid[j][i].first)));
//    }
//  }



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



//Escrever vale para ficheiro
  ofstream outfilevale2;
  outfilevale2.open("vale2.dat");

  for(int i=0; i<Ntot; i++)
    for(int j=0; j<Ntot; j++)
    {
      if(grid[j][i].second == 2)
        outfilevale2 << double(i)*h << "   " << double(j)*h << "   " << vals[grid[j][i].first].second << endl;
      else
        outfilevale2 << double(i)*h << "   " << double(j)*h << "   " << 0 << endl;
    }
  outfilevale2.close();


//Escrever vale para ficheiro
  ofstream outfilevale3;
  outfilevale3.open("vale3.dat");

  for(int i=0; i<vale3.size(); i++)
  {
    outfilevale3 << get<0>(vale3[i])*h << "   " << get<1>(vale3[i])*h << "   " << get<2>(vale3[i]) << endl;
  }
  outfilevale3.close();


//Escrever vale para ficheiro
  ofstream outfilevale4;
  outfilevale4.open("vale4.dat");

  for(int i=0; i<vale4.size(); i++)
  {
    outfilevale4 << get<0>(vale4[i])*h << "   " << get<1>(vale4[i])*h << "   " << get<2>(vale4[i]) << endl;
  }
  outfilevale4.close();

//Escrever vale para ficheiro - ultima versao
  ofstream outfilevale5;
  outfilevale5.open("valley.dat");

  for(int i=0; i<valley.size(); i++)
  {
    outfilevale5 << get<0>(valley[i])*h << "   " << get<1>(valley[i])*h << "   " << get<2>(valley[i]) << endl;
  }
  outfilevale5.close();


//Calculo de valores e vectores proprios

  SparseGenMatProd<double> op(A);
  cout << "BATATA1" << endl;
  GenEigsSolver< double, SMALLEST_MAGN, SparseGenMatProd<double> > eigs(&op, Nmodos, 2*Nmodos+30);
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

//  cout << "valor proprio 1 estado excitado: " << evalues(0) << endl;
//  cout << "valor proprio fundamental: " << evalues(1) << endl;
  cout << "valor proprio fundamental: " << evalues(0) << endl;


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
