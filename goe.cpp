#include<iostream>
#include<ctime>
#include<cstdlib>
#include<random>
#include<fstream>
#include<armadillo>


using namespace arma;
using namespace std;

void run_simulation()
{
  cout<<"Hi\n";
  // Seed is generated for the random number generator
  time_t current_time=time(NULL); //number of seconds since 00:00 1 Jan 1970
  int seed = (int)current_time;
  cout<<"The seed used for this simulation is "<<seed<<"\n";
  //Mersenne twister engine is called
  mt19937 gen(seed);
  double mu=0.0, sigma_sq=1.0;
  normal_distribution<double> distribution(mu, sigma_sq);

  int N;//dimension of the matrix
  int T;//number of matrices
  int m;
  ofstream file;

  N=200;
  m=N*N;
  T=500;
  mat A(N, N);
  cx_vec eig_val;


  file.open("eig_vals_goe.dat");
  file<<"";
  file.close();
  
  file.open("eig_vals_goe.dat", ios::app);
    for(int i=0; i<T; i++)
      {
        for(int j=0; j<m; j++)
	  {
	    A[j]=distribution(gen);
	  }
//	cout<<A<<"\n";

	cout<<"After symmetrising...Matrix "<<i<<"\n";
	A=(A+A.t())/2.0;
	eig_val = eig_gen(A);
	file<<real(eig_val);
//	cout<<A<<"\n";
	cout<<"The eigen values are \n"<<real(eig_val)<<"\n";
       
	
      }
    file.close();
  
}


int main()
{
  run_simulation();
  return 0;
}
