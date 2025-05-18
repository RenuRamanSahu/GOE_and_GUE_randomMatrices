#include<iostream>
#include<random>
#include<ctime>
#include<cstdlib>
#include<armadillo>
#include<fstream>
#include<iomanip>

using namespace std;
using namespace arma;


void run_simulation()
{
  time_t current_time=time(NULL);
  int seed=(int)current_time;
  seed=1573230545;
  cout<<seed<<"\n";

  mt19937 gen(seed);
  double mu=0, var=1.0;
  normal_distribution<double> d(mu, var);
  int n=8;
  mat A(n, n, fill::zeros);
  for(int i=0; i<n; i++)
    {
      for(int j=0; j<n; j++)
	{
	  A(i,j)=d(gen);
	}

    }
//A=(A+A.t())/2.0;
  cout<<A<<"\n";


  ofstream file;
  file.open("gaussianRM.txt");
  for(int i=0; i<n; i++)
    {
      for(int j=0; j<n; j++)
	{
	  if(j!=(n-1))
	    {file<<setprecision(4)<<A(i, j)<<"  &   ";}
	  else
	    {file<<setprecision(4)<<A(i, j)<<" \\\\ \n";}
	}
    }

  file<<"\n"<<real(eig_gen(A))<<"\n";
  file.close();

  
  
}

int main()
{
  run_simulation();
  return 0;
}
